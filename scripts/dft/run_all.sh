#!/bin/bash
# DFT pipeline: scf -> critical k-points (bands) -> projwfc  (Quantum ESPRESSO 7.x)
# Run from this folder after setup.sh. Outputs are copied to ../../results/dft/.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
WORK=~/kolB_dft
mkdir -p "$WORK" && cd "$WORK"
cp "$SRC/scf.in" "$SRC/kpts.in" "$SRC/projwfc.in" .
BIN=${QE_BIN:-$(dirname "$(command -v pw.x)")}
MPIRUN="${MPIRUN:-mpirun}"
export OMP_NUM_THREADS=1

echo "=== SCF: $(date) ==="
$MPIRUN -np 8 $BIN/pw.x -nk 4 -in scf.in > scf.out
grep -E "!|highest occupied|convergence has" scf.out | tail -5

echo "=== critical k-points (bands) ==="
$MPIRUN -np 8 $BIN/pw.x -nk 1 -in kpts.in > kpts.out
echo "=== projwfc ==="
$MPIRUN -np 1 $BIN/projwfc.x -in projwfc.in > proj.out

cp scf.out kpts.out proj.out "$SRC/../../results/dft/"
echo "=== done: $(date) ==="
