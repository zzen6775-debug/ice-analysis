#!/bin/bash
#SBATCH --job-name=FRAMES
#SBATCH --partition=pub
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/xujianning/zsmx/_frames/extract_%j.out
#SBATCH --error=/scratch/xujianning/zsmx/_frames/extract_%j.err
set -uo pipefail
R=/scratch/xujianning/zsmx
OUT=$R/_frames
mkdir -p "$OUT"
cd "$OUT" || exit 1
module load gromacs/2025-cpu 2>/dev/null || module load gromacs 2>/dev/null || true
if command -v gmx >/dev/null 2>&1; then GMX=gmx
elif command -v gmx_mpi >/dev/null 2>&1; then GMX=gmx_mpi
else echo "no gmx found"; exit 1; fi
echo "using $GMX"

echo "=== $(date)  extracting frames ==="
for pair in "HEMA_on_ice_slab:HEMA" "MMA_on_ice_slab:MMA" "PMMA_on_ice_slab:PMMA10" "PMMA20_on_ice_slab:PMMA20"; do
  sys="${pair%%:*}"; tag="${pair##*:}"
  d="$R/$sys/npt"
  for t in 0 5000 10000 20000; do
    out="$OUT/${tag}_${t}ps.gro"
    [ -s "$out" ] && { echo "  $out exists, skip"; continue; }
    printf 'System\n' | $GMX trjconv -f "$d/npt_260K.xtc" -s "$d/npt_260K.tpr" \
      -pbc whole -dump "$t" -o "$out" > /dev/null 2>&1
    if [ -s "$out" ]; then
      echo "  $tag @ ${t} ps -> $(basename "$out")  ($(head -2 "$out" | tail -1 | tr -d ' ') atoms)"
    else
      echo "  $tag @ ${t} ps FAILED"
    fi
  done
done

echo "=== $(date)  packing ==="
tar czf "$OUT/frames.tar.gz" -C "$OUT" $(cd "$OUT" && ls *.gro 2>/dev/null) 2>/dev/null
ls -lh "$OUT/frames.tar.gz" 2>/dev/null
echo "=== $(date)  done ==="
