#!/bin/bash
set -euo pipefail

# One-shot setup script for Svante workflows:
# - Build ecbuild and ecCodes (compiler/toolchain configurable)
# - Build FLEXPART binary with portable flags
# - Optionally create/update a stable postprocess conda env
#
# Run on an interactive EDR/FDR compute node.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ECCODES_VERSION="2.23.0"
ECCODES_PREFIX="${HOME}/local/eccodes-gcc11"
ECBUILD_PREFIX="${HOME}/local/ecbuild"
GCC11_PREFIX="/home/software/rhel/8/gcc/11.3.0"
FORTRAN_COMPILER=""
C_COMPILER=""
CXX_COMPILER=""
POSTPROCESS_ENV_NAME="flexpart-post"
BUILD_JOBS="8"
DRY_RUN="0"
SKIP_ECCODES="0"
SKIP_FLEXPART="0"
SETUP_POSTPROCESS_ENV="1"
BUILD_MODULE_INIT="module load cmake/3.26.4 || true"

usage() {
  cat <<'EOF'
Usage: run_scripts/setup_svante_environment.sh [options]

Options:
  --eccodes-version VER        ecCodes version (default: 2.23.0)
  --eccodes-prefix DIR         ecCodes install prefix (default: $HOME/local/eccodes-gcc11)
  --ecbuild-prefix DIR         ecbuild install prefix (default: $HOME/local/ecbuild)
  --postprocess-env NAME       Conda env name for postprocess stack (default: flexpart-post)
  --gcc-prefix DIR             GCC toolchain prefix containing bin/{gcc,g++,gfortran}
                                (default: /home/software/rhel/8/gcc/11.3.0)
  --fortran-compiler PATH      Explicit gfortran path for ecCodes/FLEXPART build
  --c-compiler PATH            Explicit C compiler path for ecCodes build
  --cxx-compiler PATH          Explicit C++ compiler path for ecCodes build
  --jobs N                     Parallel build jobs (default: 8)
  --skip-eccodes               Skip ecbuild/ecCodes build step
  --skip-flexpart              Skip FLEXPART build step
  --no-postprocess-env         Do not create/update postprocess conda env
  --dry-run                    Print commands without executing them
  -h, --help                   Show this help message

Notes:
- Run on Svante compute nodes (EDR/FDR), not login/fs nodes.
- This script follows the README build flow and keeps binary flags portable.
EOF
}

log() {
  echo "[setup] $*"
}

run_cmd() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "+ $*"
    return 0
  fi
  echo "+ $*"
  "$@"
}

run_shell() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "+ $*"
    return 0
  fi
  echo "+ bash -lc $*"
  bash -lc "$*"
}

run_with_build_modules() {
  local cmd="$*"
  run_shell "source ~/.bashrc >/dev/null 2>&1 || true; ${BUILD_MODULE_INIT}; ${cmd}"
}

run_clean_conda_shell() {
  local cmd="$*"
  run_shell "env -u CONDA_PREFIX -u CONDA_DEFAULT_ENV -u CONDA_EXE -u CONDA_PYTHON_EXE -u CONDA_SHLVL -u _CE_M -u _CE_CONDA CONDA_AUTO_ACTIVATE_BASE=false CONDA_NO_PLUGINS=true bash -lc 'source ~/.bashrc >/dev/null 2>&1 || true; ${cmd}'"
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: missing ${label}: ${path}" >&2
    exit 2
  fi
}

require_executable() {
  local path="$1"
  local label="$2"
  if [[ ! -x "${path}" ]]; then
    echo "ERROR: missing ${label} executable: ${path}" >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --eccodes-version)
      ECCODES_VERSION="$2"
      shift 2
      ;;
    --eccodes-prefix)
      ECCODES_PREFIX="$2"
      shift 2
      ;;
    --ecbuild-prefix)
      ECBUILD_PREFIX="$2"
      shift 2
      ;;
    --postprocess-env)
      POSTPROCESS_ENV_NAME="$2"
      shift 2
      ;;
    --gcc-prefix)
      GCC11_PREFIX="$2"
      shift 2
      ;;
    --fortran-compiler)
      FORTRAN_COMPILER="$2"
      shift 2
      ;;
    --c-compiler)
      C_COMPILER="$2"
      shift 2
      ;;
    --cxx-compiler)
      CXX_COMPILER="$2"
      shift 2
      ;;
    --jobs)
      BUILD_JOBS="$2"
      shift 2
      ;;
    --skip-eccodes)
      SKIP_ECCODES="1"
      shift
      ;;
    --skip-flexpart)
      SKIP_FLEXPART="1"
      shift
      ;;
    --no-postprocess-env)
      SETUP_POSTPROCESS_ENV="0"
      shift
      ;;
    --dry-run)
      DRY_RUN="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

compiler_can_link() {
  local lang="$1"
  local compiler="$2"
  local src bin logf

  if [[ "${DRY_RUN}" == "1" ]]; then
    return 0
  fi

  case "${lang}" in
    fortran)
      src="$(mktemp /tmp/fp_link_XXXXXX.f90)"
      cat > "${src}" <<'EOF'
program fp_link_test
  print *, 1
end program fp_link_test
EOF
      ;;
    c)
      src="$(mktemp /tmp/fp_link_XXXXXX.c)"
      cat > "${src}" <<'EOF'
int main(void) { return 0; }
EOF
      ;;
    cxx)
      src="$(mktemp /tmp/fp_link_XXXXXX.cpp)"
      cat > "${src}" <<'EOF'
int main() { return 0; }
EOF
      ;;
    *)
      echo "ERROR: unsupported language for compiler_can_link: ${lang}" >&2
      exit 2
      ;;
  esac

  bin="${src%.*}"
  logf="${bin}.log"
  if "${compiler}" "${src}" -o "${bin}" >"${logf}" 2>&1; then
    rm -f "${src}" "${bin}" "${logf}"
    return 0
  fi

  log "Compiler link test failed for ${lang} compiler: ${compiler}"
  head -n 20 "${logf}" || true
  rm -f "${src}" "${bin}" "${logf}"
  return 1
}

compiler_supports_pthread() {
  local lang="$1"
  local compiler="$2"
  local src bin logf

  if [[ "${DRY_RUN}" == "1" ]]; then
    return 0
  fi

  case "${lang}" in
    c)
      src="$(mktemp /tmp/fp_pth_XXXXXX.c)"
      cat > "${src}" <<'EOF'
#include <pthread.h>
void* worker(void* x) { return x; }
int main(void) {
  pthread_t t;
  return pthread_create(&t, 0, worker, 0);
}
EOF
      ;;
    cxx)
      src="$(mktemp /tmp/fp_pth_XXXXXX.cpp)"
      cat > "${src}" <<'EOF'
#include <pthread.h>
void* worker(void* x) { return x; }
int main() {
  pthread_t t;
  return pthread_create(&t, 0, worker, 0);
}
EOF
      ;;
    *)
      echo "ERROR: unsupported language for compiler_supports_pthread: ${lang}" >&2
      exit 2
      ;;
  esac

  bin="${src%.*}"
  logf="${bin}.log"
  if "${compiler}" -pthread "${src}" -o "${bin}" >"${logf}" 2>&1; then
    rm -f "${src}" "${bin}" "${logf}"
    return 0
  fi

  log "Compiler pthread test failed for ${lang} compiler: ${compiler}"
  head -n 20 "${logf}" || true
  rm -f "${src}" "${bin}" "${logf}"
  return 1
}

toolchain_can_link() {
  local f="$1"
  local c="$2"
  local cxx="$3"

  compiler_can_link fortran "${f}" || return 1
  compiler_can_link c "${c}" || return 1
  compiler_can_link cxx "${cxx}" || return 1
  compiler_supports_pthread c "${c}" || return 1
  compiler_supports_pthread cxx "${cxx}" || return 1
  return 0
}

resolve_compilers() {
  if [[ -n "${FORTRAN_COMPILER}" ]]; then
    GF="${FORTRAN_COMPILER}"
  elif [[ -x "${GCC11_PREFIX}/bin/gfortran" ]]; then
    GF="${GCC11_PREFIX}/bin/gfortran"
  else
    GF="$(command -v gfortran || true)"
  fi

  if [[ -z "${GF}" || ! -x "${GF}" ]]; then
    echo "ERROR: no usable Fortran compiler found. Set --fortran-compiler PATH or --gcc-prefix DIR." >&2
    exit 2
  fi

  if [[ -n "${C_COMPILER}" ]]; then
    CC="${C_COMPILER}"
  elif [[ -x "${GCC11_PREFIX}/bin/gcc" ]]; then
    CC="${GCC11_PREFIX}/bin/gcc"
  else
    CC="$(command -v gcc || true)"
  fi

  if [[ -n "${CXX_COMPILER}" ]]; then
    CXX="${CXX_COMPILER}"
  elif [[ -x "${GCC11_PREFIX}/bin/g++" ]]; then
    CXX="${GCC11_PREFIX}/bin/g++"
  else
    CXX="$(command -v g++ || true)"
  fi

  if [[ -z "${CC}" || ! -x "${CC}" ]]; then
    echo "ERROR: no usable C compiler found. Set --c-compiler PATH or --gcc-prefix DIR." >&2
    exit 2
  fi
  if [[ -z "${CXX}" || ! -x "${CXX}" ]]; then
    echo "ERROR: no usable C++ compiler found. Set --cxx-compiler PATH or --gcc-prefix DIR." >&2
    exit 2
  fi

  if ! toolchain_can_link "${GF}" "${CC}" "${CXX}"; then
    explicit_compiler_overrides=0
    if [[ -n "${FORTRAN_COMPILER}" || -n "${C_COMPILER}" || -n "${CXX_COMPILER}" ]]; then
      explicit_compiler_overrides=1
    fi

    if [[ "${explicit_compiler_overrides}" == "1" ]]; then
      echo "ERROR: selected compilers cannot link test programs on this node." >&2
      echo "This is often a GLIBC/toolchain mismatch (e.g., libgfortran built against newer GLIBC than node provides)." >&2
      echo "Try different compiler paths for --fortran-compiler/--c-compiler/--cxx-compiler." >&2
      exit 2
    fi

    fallback_gf="$(command -v gfortran || true)"
    fallback_cc="$(command -v gcc || true)"
    fallback_cxx="$(command -v g++ || true)"

    if [[ -n "${fallback_gf}" && -n "${fallback_cc}" && -n "${fallback_cxx}" ]] && \
       [[ "${fallback_gf}" != "${GF}" || "${fallback_cc}" != "${CC}" || "${fallback_cxx}" != "${CXX}" ]]; then
      log "Primary toolchain failed link test; trying system compilers as fallback."
      if toolchain_can_link "${fallback_gf}" "${fallback_cc}" "${fallback_cxx}"; then
        GF="${fallback_gf}"
        CC="${fallback_cc}"
        CXX="${fallback_cxx}"
        log "Using fallback compiler toolchain from PATH."
      else
        echo "ERROR: both preferred and fallback compilers failed link tests." >&2
        echo "Likely GLIBC/toolchain mismatch on this node. Try running on a compatible build node or pass explicit compiler paths." >&2
        exit 2
      fi
    else
      echo "ERROR: selected compilers cannot link test programs and no alternate toolchain was found in PATH." >&2
      echo "Likely GLIBC/toolchain mismatch on this node. Try a different --gcc-prefix or explicit compiler paths." >&2
      exit 2
    fi
  fi

  export GF CC CXX
  log "Compiler selection:"
  log "  Fortran: ${GF}"
  log "  C:       ${CC}"
  log "  C++:     ${CXX}"
}

check_fortran_abi_compat() {
  local ecc_f90_lib="$1"

  if [[ "${DRY_RUN}" == "1" ]]; then
    log "Dry-run: skipping Fortran ABI compatibility check."
    return 0
  fi

  if [[ ! -f "${ecc_f90_lib}" ]]; then
    echo "ERROR: expected ecCodes Fortran library not found: ${ecc_f90_lib}" >&2
    exit 2
  fi

  local ecc_soname compiler_soname
  ecc_soname="$(ldd "${ecc_f90_lib}" | awk '/libgfortran\.so\./ {print $1; exit}')"
  if [[ -z "${ecc_soname}" ]]; then
    echo "ERROR: could not determine libgfortran soname required by ${ecc_f90_lib}" >&2
    exit 2
  fi

  local test_src test_bin
  test_src="$(mktemp /tmp/fp_abi_XXXXXX.f90)"
  test_bin="${test_src%.f90}"
  cat > "${test_src}" <<'EOF'
program fp_abi_test
  print *, 1
end program fp_abi_test
EOF

  "${GF}" "${test_src}" -o "${test_bin}" >/dev/null 2>&1
  compiler_soname="$(ldd "${test_bin}" | awk '/libgfortran\.so\./ {print $1; exit}')"
  rm -f "${test_src}" "${test_bin}"

  if [[ -z "${compiler_soname}" ]]; then
    echo "ERROR: could not determine libgfortran soname produced by compiler ${GF}" >&2
    exit 2
  fi

  if [[ "${ecc_soname}" != "${compiler_soname}" ]]; then
    echo "ERROR: Fortran ABI mismatch detected." >&2
    echo "  ecCodes Fortran library needs: ${ecc_soname}" >&2
    echo "  compiler ${GF} produces:      ${compiler_soname}" >&2
    echo "Select a different compiler/toolchain or rebuild ecCodes with this compiler." >&2
    exit 2
  fi

  log "Fortran ABI check passed: ${compiler_soname}"
}

resolve_compilers

if [[ "${GF}" == "${GCC11_PREFIX}"/* || "${CC}" == "${GCC11_PREFIX}"/* || "${CXX}" == "${GCC11_PREFIX}"/* ]]; then
  BUILD_MODULE_INIT="module load cmake/3.26.4 || true; module load gcc/11.3.0 || true"
  log "Build modules: cmake/3.26.4 + gcc/11.3.0"
else
  log "Build modules: cmake/3.26.4 (gcc module not loaded; using selected compilers directly)"
fi

HOST_SHORT="$(hostname -s || true)"
if [[ "${HOST_SHORT}" != c* && "${HOST_SHORT}" != c* ]]; then
  log "WARNING: hostname '${HOST_SHORT}' does not look like an EDR/FDR compute node."
  log "         Continue only if this node has the expected build/runtime libraries."
fi

if [[ "${SKIP_ECCODES}" != "1" ]]; then
  log "Building ecbuild + ecCodes ${ECCODES_VERSION}"

  run_with_build_modules "cmake --version | head -n 1"

  if [[ "${DRY_RUN}" != "1" ]]; then
    # Avoid recursive filesystem scans on shared storage (can appear hung for minutes).
    # Use known toolchain/system lib paths and include only directories that exist.
    # If we fell back to non-GCC11 compilers, avoid forcing GCC11 runtime libs.
    extra_lib_dirs=()
    if [[ "${GF}" == "${GCC11_PREFIX}"/* || "${CC}" == "${GCC11_PREFIX}"/* || "${CXX}" == "${GCC11_PREFIX}"/* ]]; then
      [[ -d "${GCC11_PREFIX}/lib64" ]] && extra_lib_dirs+=("${GCC11_PREFIX}/lib64")
      [[ -d "${GCC11_PREFIX}/lib" ]] && extra_lib_dirs+=("${GCC11_PREFIX}/lib")
    fi
    [[ -d "/usr/lib64" ]] && extra_lib_dirs+=("/usr/lib64")
    [[ -d "/usr/lib" ]] && extra_lib_dirs+=("/usr/lib")
    EXTRA_LIBS="$(printf '%s\n' "${extra_lib_dirs[@]}" | awk 'NF && !seen[$0]++' | paste -sd: -)"
    if [[ -n "${EXTRA_LIBS}" ]]; then
      export LD_LIBRARY_PATH="${EXTRA_LIBS}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
    fi
    log "Using LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
  fi

  SRC_BASE="${HOME}/local/src"
  BUILD_BASE="${HOME}/local/build"
  run_cmd mkdir -p "${SRC_BASE}" "${BUILD_BASE}"

  if [[ ! -d "${SRC_BASE}/ecbuild" ]]; then
    run_cmd git clone --depth 1 --branch 3.8.2 https://github.com/ecmwf/ecbuild.git "${SRC_BASE}/ecbuild"
  else
    log "ecbuild source already exists: ${SRC_BASE}/ecbuild"
  fi

  run_cmd mkdir -p "${BUILD_BASE}/ecbuild-build"
  run_with_build_modules "cmake -S '${SRC_BASE}/ecbuild' -B '${BUILD_BASE}/ecbuild-build' -DCMAKE_INSTALL_PREFIX='${ECBUILD_PREFIX}'"
  run_with_build_modules "cmake --build '${BUILD_BASE}/ecbuild-build' -j '${BUILD_JOBS}'"
  run_with_build_modules "cmake --install '${BUILD_BASE}/ecbuild-build'"

  ECCODES_SRC_DIR="${SRC_BASE}/eccodes-${ECCODES_VERSION}"
  ECCODES_TARBALL="${SRC_BASE}/eccodes-${ECCODES_VERSION}-Source.tar.gz"
  if [[ ! -d "${ECCODES_SRC_DIR}" ]]; then
    if [[ ! -f "${ECCODES_TARBALL}" ]]; then
      log "Downloading ecCodes source tarball"
      if command -v curl >/dev/null 2>&1; then
        run_cmd curl -fL -o "${ECCODES_TARBALL}" "https://confluencehpc.ecmwf.int/eccodes-${ECCODES_VERSION}-Source.tar.gz"
      elif command -v wget >/dev/null 2>&1; then
        run_cmd wget -O "${ECCODES_TARBALL}" "https://confluencehpc.ecmwf.int/eccodes-${ECCODES_VERSION}-Source.tar.gz"
      else
        echo "ERROR: neither curl nor wget found to download ecCodes." >&2
        exit 2
      fi
    fi
    run_cmd tar xzf "${ECCODES_TARBALL}" -C "${SRC_BASE}"
    if [[ -d "${SRC_BASE}/eccodes-${ECCODES_VERSION}-Source" ]]; then
      run_cmd mv "${SRC_BASE}/eccodes-${ECCODES_VERSION}-Source" "${ECCODES_SRC_DIR}"
    fi
  else
    log "ecCodes source already exists: ${ECCODES_SRC_DIR}"
  fi

  if [[ "${GF}" == "${GCC11_PREFIX}"/* || "${CC}" == "${GCC11_PREFIX}"/* || "${CXX}" == "${GCC11_PREFIX}"/* ]]; then
    ECCODES_BUILD_DIR="${BUILD_BASE}/eccodes-gcc11-build"
  else
    ECCODES_BUILD_DIR="${BUILD_BASE}/eccodes-system-build"
  fi
  run_cmd mkdir -p "${ECCODES_BUILD_DIR}"

  # Always clear CMake cache before configure to avoid stale compiler/toolchain
  # state when switching between GCC11 and fallback system compilers.
  if [[ -f "${ECCODES_BUILD_DIR}/CMakeCache.txt" ]]; then
    run_cmd rm -f "${ECCODES_BUILD_DIR}/CMakeCache.txt"
  fi
  if [[ -d "${ECCODES_BUILD_DIR}/CMakeFiles" ]]; then
    run_cmd rm -rf "${ECCODES_BUILD_DIR}/CMakeFiles"
  fi

  run_with_build_modules "cmake -S '${ECCODES_SRC_DIR}' -B '${ECCODES_BUILD_DIR}' -DCMAKE_INSTALL_PREFIX='${ECCODES_PREFIX}' -DCMAKE_C_COMPILER='${CC}' -DCMAKE_CXX_COMPILER='${CXX}' -DCMAKE_Fortran_COMPILER='${GF}' -DCMAKE_PREFIX_PATH='${ECBUILD_PREFIX}' -DENABLE_FORTRAN=ON -DBUILD_SHARED_LIBS=ON -DENABLE_AEC=OFF -DENABLE_NETCDF=OFF -DENABLE_JPG=ON -DENABLE_PNG=OFF"

  run_with_build_modules "cmake --build '${ECCODES_BUILD_DIR}' -j '${BUILD_JOBS}'"
  run_with_build_modules "cmake --install '${ECCODES_BUILD_DIR}'"

  if [[ "${DRY_RUN}" != "1" ]]; then
    require_file "${ECCODES_PREFIX}/include/eccodes.h" "ecCodes header"
    if [[ -f "${ECCODES_PREFIX}/lib64/libeccodes.so" ]]; then
      require_file "${ECCODES_PREFIX}/lib64/libeccodes_f90.so" "ecCodes Fortran library"
      log "Verified ecCodes artifacts under ${ECCODES_PREFIX}/lib64"
    else
      require_file "${ECCODES_PREFIX}/lib/libeccodes.so" "ecCodes shared library"
      require_file "${ECCODES_PREFIX}/lib/libeccodes_f90.so" "ecCodes Fortran library"
      log "Verified ecCodes artifacts under ${ECCODES_PREFIX}/lib"
    fi
  fi
fi

if [[ "${SKIP_FLEXPART}" != "1" ]]; then
  log "Building FLEXPART in ${REPO_ROOT}/src"

  export ECCODES_PREFIX
  export CPATH="/usr/include:${ECCODES_PREFIX}/include:/usr/lib64/gfortran/modules"
  export LIBRARY_PATH="/usr/lib64:${ECCODES_PREFIX}/lib64:${ECCODES_PREFIX}/lib"
  unset CONDA_PREFIX
  unset LD_LIBRARY_PATH

  ECCODES_F90_LIB=""
  if [[ -f "${ECCODES_PREFIX}/lib64/libeccodes_f90.so" ]]; then
    ECCODES_F90_LIB="${ECCODES_PREFIX}/lib64/libeccodes_f90.so"
  elif [[ -f "${ECCODES_PREFIX}/lib/libeccodes_f90.so" ]]; then
    ECCODES_F90_LIB="${ECCODES_PREFIX}/lib/libeccodes_f90.so"
  else
    ECCODES_F90_LIB="${ECCODES_PREFIX}/lib64/libeccodes_f90.so"
  fi

  check_fortran_abi_compat "${ECCODES_F90_LIB}"

  run_cmd make -C "${REPO_ROOT}/src" -f makefile_svante cleanall
  run_cmd make -C "${REPO_ROOT}/src" -f makefile_svante eta=no ncf=yes SERIAL=yes arch=x86-64 GFLAG= -j "${BUILD_JOBS}" FC="${GF}" F90="${GF}"

  if [[ "${DRY_RUN}" != "1" ]]; then
    require_executable "${REPO_ROOT}/src/FLEXPART" "FLEXPART"
    log "Verified FLEXPART executable: ${REPO_ROOT}/src/FLEXPART"
  fi
fi

if [[ "${SETUP_POSTPROCESS_ENV}" == "1" ]]; then
  log "Creating/updating conda postprocess env: ${POSTPROCESS_ENV_NAME}"

  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH; cannot create postprocess env." >&2
    exit 2
  fi

  run_clean_conda_shell "conda create -y -n '${POSTPROCESS_ENV_NAME}' python=3.12"
  run_clean_conda_shell "conda run -n '${POSTPROCESS_ENV_NAME}' python -m pip install --upgrade pip && conda run -n '${POSTPROCESS_ENV_NAME}' python -m pip install numpy pandas xarray netCDF4"
fi

log "Setup complete."
log "Recommended config values in run_scripts/slurm_array_config.sh:"
log "  DISABLE_AUTO_POSTPROCESS=\"1\""
log "  PRUNE_TO_GRID_FILES=\"1\""
log "  POSTPROCESS_PYTHON_CMD=\"/home/${USER}/.conda/envs/${POSTPROCESS_ENV_NAME}/bin/python\""
log "Then run:"
log "  ./run_scripts/run_slurm_array_backward.sh"
log "  ./run_scripts/run_slurm_postprocess_all.sh"
