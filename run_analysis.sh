#!/bin/bash
# One-shot: fetch the frame tarball from the cluster, extract, analyse, fit.
# Run inside WSL:  bash /mnt/g/MD_useful/ice_analysis/run_analysis.sh [REMOTE_TARBALL]
set -uo pipefail
KEY_G=/mnt/g/MD_useful/_hpc_tools/.ssh/whu_hpc
REMOTE="${1:-/scratch/xujianning/zsmx/_frames/frames.tar.gz}"
HOST=xujianning@202.114.96.180
WORK=/tmp/icework

echo "=== $(date) ==="
cp "$KEY_G" /tmp/k 2>/dev/null && chmod 600 /tmp/k
echo "-- fetching $REMOTE"
ssh -i /tmp/k -o BatchMode=yes -o StrictHostKeyChecking=accept-new "$HOST" "cat $REMOTE" > /tmp/f.tar.gz 2>/dev/null
rc=$?
sz=$(stat -c %s /tmp/f.tar.gz 2>/dev/null || echo 0)
echo "   ssh_rc=$rc  size=$((sz/1024/1024)) MB"
[ "$sz" -gt 1000 ] || { echo "fetch failed"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK"
tar xzf /tmp/f.tar.gz -C "$WORK" || exit 1
echo "-- frames: $(ls "$WORK"/*.gro 2>/dev/null | wc -l)"

echo "-- computing q(z)"
python3 /mnt/g/MD_useful/ice_analysis/ice_front.py --frames "$WORK" --out /tmp/fronts.csv
cp /tmp/fronts.csv /mnt/g/MD_useful/ice_analysis/fronts.csv 2>/dev/null && echo "   CSV copied to the workspace"
cp /tmp/fronts.csv.qprofile /mnt/g/MD_useful/ice_analysis/fronts.csv.qprofile 2>/dev/null

echo
echo "-- fitting growth rates"
python3 /mnt/g/MD_useful/ice_analysis/growth_fit.py /tmp/fronts.csv --plot /tmp/growth.png
cp /tmp/growth.png /mnt/g/MD_useful/ice_analysis/growth.png 2>/dev/null && echo "   plot copied to the workspace"
echo "=== done ==="
