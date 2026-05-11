#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR=${1:-.}
OUTPUT_DIR=${2:-./flexpart_inputs}
AVAILABLE_FILE=${3:-$OUTPUT_DIR/AVAILABLE}
REPACK_TO_SIMPLE=${REPACK_TO_SIMPLE:-1}
OVERWRITE_EXISTING=${OVERWRITE_EXISTING:-0}
STEP_FILTER_MODE=${STEP_FILTER_MODE:-nonzero_max}

mkdir -p "$OUTPUT_DIR"

echo "Input dir:  $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "AVAILABLE:  $AVAILABLE_FILE"
echo "Repack to grid_simple: $REPACK_TO_SIMPLE"
echo "Overwrite existing GF files: $OVERWRITE_EXISTING"
echo "Step filter mode: $STEP_FILTER_MODE"

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
    echo "Processing timestep: $ts"

    # Find available pgrbh cycles for this valid time.
    ts_pfiles=("$INPUT_DIR"/${ts}.cdas1.*.pgrbh.grb2)
    if (( ${#ts_pfiles[@]} == 0 )); then
        echo "  WARNING: no pgrbh files for $ts"
        ((count_skipped+=1))
        continue
    fi

    mapfile -t cycles < <(
        printf '%s\n' "${ts_pfiles[@]}" |
        sed -E 's#.*\.cdas1\.([0-9]{8})\.pgrbh\.grb2#\1#' |
        sort -u
    )

    # Pick the latest cycle (by YYYYMMDD lexical order).
    latest_cycle=${cycles[-1]}

    echo "  Using cycle: $latest_cycle"

    # Locate selected cycle files.
    pfile="$INPUT_DIR/${ts}.cdas1.${latest_cycle}.pgrbh.grb2"
    ipfile="$INPUT_DIR/${ts}.cdas1.${latest_cycle}.ipvgrbh.grb2"

    if [[ ! -f "$pfile" ]]; then
        echo "  WARNING: Missing pgrbh file for $ts"
        ((count_skipped+=1))
        continue
    fi

    # FLEXPART naming: GFyymmddhh
    yymmddhh=${ts:2:8}
    outfile="$OUTPUT_DIR/GF${yymmddhh}"

    echo "  Writing $outfile"

    if [[ -f "$outfile" && "$OVERWRITE_EXISTING" != "1" ]]; then
        echo "  Reusing existing $outfile"
        ((count_reused+=1))
    else
        # Build a temporary merged GRIB stream.
        tmpfile=$(mktemp)
        workfile="$tmpfile"
        cp "$pfile" "$tmpfile"
        if [[ -f "$ipfile" ]]; then
            cat "$ipfile" >> "$tmpfile"
        else
            echo "  WARNING: Missing ipvgrbh file for $ts (continuing with pgrbh only)"
        fi

        if [[ "$STEP_FILTER_MODE" == "nonzero_max" ]]; then
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
            unique_steps=${step_info[0]:-0}
            target_step=${step_info[1]:-}
            if [[ -n "$target_step" && "$unique_steps" -gt 1 ]]; then
                filtered_tmp=$(mktemp)
                grib_copy -w stepRange="$target_step" "$tmpfile" "$filtered_tmp"
                if [[ -s "$filtered_tmp" ]]; then
                    echo "  Keeping stepRange=$target_step"
                    workfile="$filtered_tmp"
                    ((count_filtered+=1))
                else
                    echo "  WARNING: step filter produced empty output; keeping all steps"
                    rm -f "$filtered_tmp"
                fi
            elif [[ -n "$target_step" ]]; then
                echo "  Step filtering not needed (single numeric stepRange=$target_step)"
                ((count_filter_skipped+=1))
            else
                echo "  WARNING: no numeric stepRange values found; keeping all steps"
            fi
        elif [[ "$STEP_FILTER_MODE" != "off" ]]; then
            echo "ERROR: unsupported STEP_FILTER_MODE=$STEP_FILTER_MODE"
            rm -f "$tmpfile"
            exit 2
        fi

        if [[ "$REPACK_TO_SIMPLE" == "1" ]]; then
            # Repack to simple packing so FLEXPART ecCodes does not require JPEG support.
            grib_set -r -s packingType=grid_simple "$workfile" "$outfile"
        else
            cp "$workfile" "$outfile"
        fi
        if [[ "$workfile" != "$tmpfile" ]]; then
            rm -f "$workfile"
        fi
        rm -f "$tmpfile"
        ((count_written+=1))
    fi

    yyyy=${ts:0:4}
    mm=${ts:4:2}
    dd=${ts:6:2}
    hh=${ts:8:2}
    printf '%s%s%s %s0000      %-12s ON DISK\n' "$yyyy" "$mm" "$dd" "$hh" "GF${yymmddhh}" >> "$available_tmp"

done

mv "$available_tmp" "$AVAILABLE_FILE"

echo "Done."
echo "Wrote $count_written meteorological files"
echo "Reused $count_reused existing meteorological files"
echo "Applied step filtering to $count_filtered files"
echo "Skipped step filtering for $count_filter_skipped single-step files"
echo "Skipped $count_skipped timesteps"
echo "AVAILABLE written to $AVAILABLE_FILE"

