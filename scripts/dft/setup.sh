#!/bin/bash
# Kol B DFT: ortam kurulumu + pseudo indirme
set -e
nproc
mkdir -p ~/kolB_dft/pseudo
cd ~/kolB_dft/pseudo
for f in Mo.pbe-spn-kjpaw_psl.1.0.0.UPF S.pbe-n-kjpaw_psl.1.0.0.UPF; do
  if [ ! -f "$f" ]; then
    curl -sfLo "$f" "https://pseudopotentials.quantum-espresso.org/upf_files/$f"
  fi
done
ls -la
