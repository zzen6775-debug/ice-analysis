#!/bin/bash
#SBATCH --job-name=FRAMES2
#SBATCH --partition=pub
#SBATCH --account=liukang
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/xujianning/zsmx/_frames2/extract2_%j.out
#SBATCH --error=/scratch/xujianning/zsmx/_frames2/extract2_%j.err
set -uo pipefail
R=/scratch/xujianning/zsmx
OUT=$R/_frames2
mkdir -p "$OUT"
cd "$OUT" || exit 1
module load gromacs/2025-cpu 2>/dev/null || module load gromacs 2>/dev/null || true
if command -v gmx >/dev/null 2>&1; then GMX=gmx
elif command -v gmx_mpi >/dev/null 2>&1; then GMX=gmx_mpi
else echo "no gmx"; exit 1; fi
echo "using $GMX"
echo "=== $(date) ==="

for pair in "HEMA_on_ice_slab:HEMA" "MMA_on_ice_slab:MMA" "PMMA_on_ice_slab:PMMA10" "PMMA20_on_ice_slab:PMMA20"; do
  sys="${pair%%:*}"; tag="${pair##*:}"
  d="$R/$sys/npt"
  ndx="$OUT/ow_${tag}.ndx"
  out="$OUT/frames_${tag}.gro"
  [ -s "$ndx" ] || { echo "  $tag: index file missing"; continue; }
  [ -s "$out" ] && { echo "  $tag: already done, skip"; continue; }
  echo "  -- $tag --"
  printf 'OW\n' | $GMX trjconv -f "$d/npt_260K.xtc" -s "$d/npt_260K.tpr" -n "$ndx" \
    -dt 2000 -o "$out" > /dev/null 2>&1
  if [ -s "$out" ]; then
    nfr=$(grep -c '^ *OW' "$out" 2>/dev/null || true)
    echo "     $(basename "$out"): $(du -h "$out" | awk '{print $1}'), $(grep -c 't=' "$out" 2>/dev/null) frames"
  else
    echo "     FAILED"
  fi
done

echo "=== $(date)  packing ==="
tar czf "$OUT/frames2.tar.gz" -C "$OUT" $(cd "$OUT" && ls frames_*.gro 2>/dev/null) 2>/dev/null
ls -lh "$OUT/frames2.tar.gz" 2>/dev/null
echo "=== $(date)  done ==="
