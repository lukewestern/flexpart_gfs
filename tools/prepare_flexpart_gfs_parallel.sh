#!/usr/bin/env bash
#
# Parallel version of prepare_flexpart_gfs.sh
#
# Each timestep is processed independently in a subprocess.
# Uses GNU parallel when available, otherwise falls back to xargs -P.
#
# Additional env var:
#   NPROC   number of parallel workers (default: 8)

set -euo pipefail

INPUT_DIR=${1:-.}
OUTPUT_DIR=${2:-./flexpart_inputs}
AVAILABLE_FILE=${3:-$OUTPUT_DIR/AVAILABLE}
REPACK_TO_SIMPLE=${REPACK_TO_SIMPLE:-1}
OVERWRITE_EXISTING=${OVERWRITE_EXISTING:-0}
STEP_FILTER_MODE=${STEP_FILTER_MODE:-nonzero_max}
NPROC=${NPROC:-8}
REQUIRED_FIELDS=${REQUIRED_FIELDS:-"2t@heightAboveGround@2 2r@heightAboveGround@2 10u@heightAboveGround@10 10v@heightAboveGround@10 sp@surface@0"}

mkdir -p "$OUTPUT_DIR"

echo "Input dir:  $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "AVAILABLE:  $AVAILABLE_FILE"
echo "Repack to grid_simple: $REPACK_TO_SIMPLE"
echo "Overwrite existing GF files: $OVERWRITE_EXISTING"
echo "Step filter mode: $STEP_FILTER_MODE"
echo "Workers: $NPROC"
echo "Required fields: $REQUIRED_FIELDS"

if [[ "$REPACK_TO_SIMPLE" == "1" ]] && ! command -v grib_set >/dev/null 2>&1; then
    echo "ERROR: REPACK_TO_SIMPLE=1 but grib_set is not in PATH"
    exit 2
fi

if [[ "$STEP_FILTER_MODE" != "off" ]] && ! command -v grib_copy >/dev/null 2>&1; then
    echo "ERROR: STEP_FILTER_MODE=$STEP_FILTER_MODE but grib_copy is not in PATH"
    exit 2
fi

shopt -s nullglob

# Build a sorted list of unique valid times from pgrbh files.
pfiles=("$INPUT_DIR"/*.cdas1.*.pgrbh.grb2)
if (( ${#pfiles[@]} == 0 )); then
    echo "ERROR: no pgrbh files found in $INPUT_DIR"
    exit 1
fi

mapfile -t timesteps < <(
    printf '%s\n' "${pfiles[@]}" |
    sed -E 's#.*/([0-9]{12})\.cdas1\.[0-9]{8}\.pgrbh\.grb2#\1#' |
    sort -u
)

if (( ${#timesteps[@]} == 0 )); then
    echo "ERROR: no valid timestamps could be parsed"
    exit 1
fi

echo "Found ${#timesteps[@]} timesteps to process"

# Temp directory for per-timestep output: AVAILABLE entries and status markers.
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

# Validate one required field triplet in a GRIB file.
_has_required_triplet() {
    local file="$1"
    local short_name="$2"
    local level_type="$3"
    local level="$4"

    grib_ls -p shortName,typeOfLevel,level "$file" 2>/dev/null | \
        awk -v s="$short_name" -v t="$level_type" -v l="$level" '
            NR > 2 && $1 == s && $2 == t && $3 == l {
                found = 1
                exit
            }
            END { exit(found ? 0 : 1) }
        '
}

# Ensure all required fields are present in the prepared output file.
_validate_required_fields() {
    local file="$1"
    local spec short_name level_type level
    local missing=()

    for spec in $REQUIRED_FIELDS; do
        IFS='@' read -r short_name level_type level <<< "$spec"
        if ! _has_required_triplet "$file" "$short_name" "$level_type" "$level"; then
            missing+=("$spec")
        fi
    done

    if (( ${#missing[@]} > 0 )); then
        echo "  WARNING: missing required fields in $(basename "$file"): ${missing[*]}"
        return 1
    fi

    return 0
}

# ---------------------------------------------------------------------------
# Per-timestep worker function.
# Writes:
#   $work_dir/<ts>.avail   — AVAILABLE file line (only on success)
#   $work_dir/<ts>.status  — one of: written | reused | skipped | filtered | filter_skipped
# ---------------------------------------------------------------------------
_process_one() {
    local ts="$1"

    echo "Processing timestep: $ts"

    # Find available pgrbh cycles for this valid time.
    local ts_pfiles
    mapfile -t ts_pfiles < <(printf '%s\n' "$INPUT_DIR"/${ts}.cdas1.*.pgrbh.grb2 2>/dev/null || true)
    # nullglob may not be set in subshell; filter missing
    local real_pfiles=()
    for f in "${ts_pfiles[@]}"; do [[ -f "$f" ]] && real_pfiles+=("$f"); done

    if (( ${#real_pfiles[@]} == 0 )); then
        echo "  WARNING: no pgrbh files for $ts"
        echo "skipped" > "$work_dir/${ts}.status"
        return 0
    fi

    mapfile -t cycles < <(
        printf '%s\n' "${real_pfiles[@]}" |
        sed -E 's#.*\.cdas1\.([0-9]{8})\.pgrbh\.grb2#\1#' |
        sort -ru
    )

    local yymmddhh=${ts:2:8}
    local outfile="$OUTPUT_DIR/GF${yymmddhh}"

    echo "  Writing $outfile"

    local status="written"
    local selected_cycle=""

    if [[ -f "$outfile" && "$OVERWRITE_EXISTING" != "1" ]]; then
        if _validate_required_fields "$outfile"; then
            echo "  Reusing existing $outfile"
            status="reused"
        else
            echo "  Existing $outfile failed validation; rebuilding"
        fi
    fi

    if [[ "$status" != "reused" ]]; then
        local cycle
        local candidate_ok=0

        for cycle in "${cycles[@]}"; do
            local pfile="$INPUT_DIR/${ts}.cdas1.${cycle}.pgrbh.grb2"
            local ipfile="$INPUT_DIR/${ts}.cdas1.${cycle}.ipvgrbh.grb2"

            if [[ ! -f "$pfile" ]]; then
                continue
            fi

            echo "  Trying cycle: $cycle"

            local tmpfile
            tmpfile=$(mktemp)
            local workfile="$tmpfile"
            local candidate_status="written"

            cp "$pfile" "$tmpfile"
            if [[ -f "$ipfile" ]]; then
                cat "$ipfile" >> "$tmpfile"
            else
                echo "  WARNING: Missing ipvgrbh file for $ts cycle $cycle (continuing with pgrbh only)"
            fi

            if [[ "$STEP_FILTER_MODE" == "nonzero_max" ]]; then
                local step_info
                mapfile -t step_info < <(
                    grib_ls -p stepRange "$tmpfile" 2>/dev/null |
                    awk '
                        NR > 2 && NF == 1 && $1 ~ /^[0-9]+$/ {
                            if (!seen[$1]++) {
                                unique_count++
                            }
                            if (!have_any || $1 > max_any) {
                                max_any = $1
                                have_any = 1
                            }
                            if ($1 > 0 && (!have_nonzero || $1 > max_nonzero)) {
                                max_nonzero = $1
                                have_nonzero = 1
                            }
                        }
                        END {
                            print unique_count + 0
                            if (have_nonzero) {
                                print max_nonzero
                            } else if (have_any) {
                                print max_any
                            } else {
                                print ""
                            }
                        }
                    '
                )
                local unique_steps=${step_info[0]:-0}
                local target_step=${step_info[1]:-}
                if [[ -n "$target_step" && "$unique_steps" -gt 1 ]]; then
                    local filtered_tmp
                    filtered_tmp=$(mktemp)
                    grib_copy -w stepRange="$target_step" "$tmpfile" "$filtered_tmp"
                    if [[ -s "$filtered_tmp" ]]; then
                        echo "  Keeping stepRange=$target_step"
                        workfile="$filtered_tmp"
                        candidate_status="filtered"
                    else
                        echo "  WARNING: step filter produced empty output; keeping all steps"
                        rm -f "$filtered_tmp"
                    fi
                elif [[ -n "$target_step" ]]; then
                    echo "  Step filtering not needed (single numeric stepRange=$target_step)"
                    candidate_status="filter_skipped"
                else
                    echo "  WARNING: no numeric stepRange values found; keeping all steps"
                fi
            elif [[ "$STEP_FILTER_MODE" != "off" ]]; then
                echo "ERROR: unsupported STEP_FILTER_MODE=$STEP_FILTER_MODE"
                rm -f "$tmpfile"
                return 2
            fi

            local candidate_out
            candidate_out=$(mktemp)
            if [[ "$REPACK_TO_SIMPLE" == "1" ]]; then
                grib_set -r -s packingType=grid_simple "$workfile" "$candidate_out"
            else
                cp "$workfile" "$candidate_out"
            fi

            if _validate_required_fields "$candidate_out"; then
                mv "$candidate_out" "$outfile"
                status="$candidate_status"
                selected_cycle="$cycle"
                candidate_ok=1
                if [[ "$workfile" != "$tmpfile" ]]; then
                    rm -f "$workfile"
                fi
                rm -f "$tmpfile"
                break
            fi

            rm -f "$candidate_out"
            if [[ "$workfile" != "$tmpfile" ]]; then
                rm -f "$workfile"
            fi
            rm -f "$tmpfile"
            echo "  Cycle $cycle rejected due to missing required fields"
        done

        if (( candidate_ok == 0 )); then
            echo "  WARNING: no valid cycle found for $ts; skipping timestep"
            echo "skipped" > "$work_dir/${ts}.status"
            return 0
        fi

        echo "  Selected cycle: $selected_cycle"
    fi

    echo "$status" > "$work_dir/${ts}.status"

    local yyyy=${ts:0:4}
    local mm=${ts:4:2}
    local dd=${ts:6:2}
    local hh=${ts:8:2}
    printf '%s%s%s %s0000      %-12s ON DISK\n' "$yyyy" "$mm" "$dd" "$hh" "GF${yymmddhh}" \
        > "$work_dir/${ts}.avail"
}

export -f _has_required_triplet
export -f _validate_required_fields
export -f _process_one
export INPUT_DIR OUTPUT_DIR REPACK_TO_SIMPLE OVERWRITE_EXISTING STEP_FILTER_MODE REQUIRED_FIELDS work_dir

# ---------------------------------------------------------------------------
# Dispatch in parallel.
# ---------------------------------------------------------------------------
if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${timesteps[@]}" | parallel --will-cite -j "$NPROC" _process_one {}
else
    printf '%s\n' "${timesteps[@]}" | \
        xargs -P "$NPROC" -I{} bash -c '_process_one "$@"' _ {}
fi

# ---------------------------------------------------------------------------
# Assemble AVAILABLE file in timestep order.
# ---------------------------------------------------------------------------
available_tmp=$(mktemp)

cat > "$available_tmp" <<'EOF'
XXXXXX EMPTY LINES XXXXXXXXX
XXXXXX EMPTY LINES XXXXXXXX
YYYYMMDD HHMMSS   name of the file(up to 80 characters)
EOF

count_written=0
count_skipped=0
count_reused=0
count_filtered=0
count_filter_skipped=0

for ts in "${timesteps[@]}"; do
    status_file="$work_dir/${ts}.status"
    avail_file="$work_dir/${ts}.avail"

    if [[ -f "$status_file" ]]; then
        status=$(< "$status_file")
        case "$status" in
            written)        ((count_written+=1)) ;;
            reused)         ((count_reused+=1)) ;;
            skipped)        ((count_skipped+=1)) ;;
            filtered)       ((count_written+=1)); ((count_filtered+=1)) ;;
            filter_skipped) ((count_written+=1)); ((count_filter_skipped+=1)) ;;
        esac
    fi

    if [[ -f "$avail_file" ]]; then
        cat "$avail_file" >> "$available_tmp"
    fi
done

mv "$available_tmp" "$AVAILABLE_FILE"

echo "Done."
echo "Wrote $count_written meteorological files"
echo "Reused $count_reused existing meteorological files"
echo "Applied step filtering to $count_filtered files"
echo "Skipped step filtering for $count_filter_skipped single-step files"
echo "Skipped $count_skipped timesteps"
echo "AVAILABLE written to $AVAILABLE_FILE"
