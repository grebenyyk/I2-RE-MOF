#!/usr/bin/env bash
set -euo pipefail

QE_ROOT=/workspace/qe-7.5
QE_BIN="$QE_ROOT/bin/pw.x"
PSEUDO_DIR="$QE_ROOT/pseudo"
MPI_ROOT=/workspace/openmpi
NVHPC_ROOT=/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11
CUDA_VER=13.0

export PATH="$MPI_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_ROOT/lib:$NVHPC_ROOT/compilers/lib:$NVHPC_ROOT/math_libs/$CUDA_VER/lib64:$NVHPC_ROOT/cuda/$CUDA_VER/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export TERM=xterm-256color

NP="${NP:-4}"
PW_FLAGS="${PW_FLAGS:--nk 1 -nd 1 -nb 1 -nt ${OMP_NUM_THREADS}}"
DETACH=0

if [[ "${1:-}" == "--detach" ]]; then
  DETACH=1
  shift
fi

usage() {
  cat <<'EOF'
Usage:
  ./run.sh [--detach] input.in [output.out]

Environment variables:
  NP                MPI ranks to launch (default: 4)
  OMP_NUM_THREADS   OpenMP threads per rank (default: 1)
  PW_FLAGS          Extra pw.x parallel flags (default: -nk 1 -nd 1 -nb 1 -nt 1)
  LAUNCH_LOG        File to capture detached launcher output

Examples:
  ./run.sh cofe2o4_co3_tetra_relax.in
  ./run.sh --detach cofe2o4_co3_octa_relax.in
  NP=8 ./run.sh cofe2o4_co3_tetra_relax.in tetra.out
EOF
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="${2:-${INPUT_FILE%.in}.out}"

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file not found: $INPUT_FILE" >&2
  exit 1
fi

if [[ ! -x "$QE_BIN" ]]; then
  echo "pw.x not found or not executable: $QE_BIN" >&2
  exit 1
fi

if ! command -v mpirun >/dev/null 2>&1; then
  echo "mpirun not found in PATH after adding $MPI_ROOT/bin" >&2
  exit 1
fi

if [[ "$DETACH" -eq 1 && "${RUN_SH_DETACHED:-0}" != "1" ]]; then
  DEFAULT_LAUNCH_LOG="${OUTPUT_FILE%.out}.launch.log"
  LAUNCH_LOG="${LAUNCH_LOG:-$DEFAULT_LAUNCH_LOG}"

  nohup env \
    RUN_SH_DETACHED=1 \
    NP="$NP" \
    OMP_NUM_THREADS="$OMP_NUM_THREADS" \
    PW_FLAGS="$PW_FLAGS" \
    LAUNCH_LOG="$LAUNCH_LOG" \
    bash "$0" "$INPUT_FILE" "$OUTPUT_FILE" \
    > "$LAUNCH_LOG" 2>&1 < /dev/null &

  echo "Detached background job started"
  echo "PID           : $!"
  echo "Input         : $INPUT_FILE"
  echo "Output        : $OUTPUT_FILE"
  echo "Launch log    : $LAUNCH_LOG"
  exit 0
fi

echo "QE executable : $QE_BIN"
echo "MPI launcher  : $(command -v mpirun)"
echo "Pseudo dir    : $PSEUDO_DIR"
echo "Input         : $INPUT_FILE"
echo "Output        : $OUTPUT_FILE"
echo "MPI ranks     : $NP"
echo "OMP threads   : $OMP_NUM_THREADS"
echo "pw.x flags    : $PW_FLAGS"

mpirun --allow-run-as-root -np "$NP" "$QE_BIN" $PW_FLAGS -inp "$INPUT_FILE" | tee "$OUTPUT_FILE"