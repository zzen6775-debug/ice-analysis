#!/bin/bash
# Full pipeline, from the frame tarball to the Origin tables.
#
#   bash run_pipeline.sh <frames.tar.gz> [cache-dir]
#
# Steps
#   1. unpack the tarball
#   2. split the multi-frame .gro into one file per frame
#   3. q(z) -> sub-bin front positions with anchored tracking + ice amount
#   4. linear / window fits
#   5. Origin tables (thickness, ice amount, both fronts, centred fronts)
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAR="${1:?usage: bash run_pipeline.sh <frames.tar.gz> [cache-dir]}"
CACHE="${2:-$HOME/ice_cache}"

mkdir -p "$CACHE"
echo "=== $(date) ==="

var ext
cp "$TAR" "$CACHE/frames.tar.gz" || exit 1
rm -rf "$CACHE/raw" "$CACHE/split"
mkdir -p "$CACHE/raw"
tar xzf "$CACHE/frames.tar.gz" -C "$CACHE/raw" || exit 1
echo "-- unpacked: $(ls "$CACHE"/raw/*.gro 2>/dev/null | wc -l) multi-frame file(s)"

echo "-- splitting"
python3 "$HERE/split_frames.py" "$CACHE/raw" "$CACHE/split"

echo
echo "-- fronts (0.2 nm grid, anchored tracking)"
python3 "$HERE/fronts_track.py" --frames "$CACHE/split" --bin 0.2 --window 1.2 \
    --out "$CACHE/fronts.csv" 2>&1 | tail -6

echo
echo "-- fits"
python3 "$HERE/growth_fit.py" "$CACHE/fronts.csv" --tube 4 --windows 0 20 40 60

echo
echo "-- Origin tables"
OUT="$CACHE/origin_data"
python3 "$HERE/export_for_origin.py" "$CACHE/fronts.csv" --outdir "$OUT" 2>&1 | head -8
python3 "$HERE/centre_fronts.py" "$OUT/origin_fronts_both.csv" \
    --out "$OUT/origin_fronts_centered.csv"

echo
echo "results in $OUT"
echo "  origin_thickness.csv          thickness vs time"
echo "  origin_ice_amount.csv         continuous ice amount vs time"
echo "  origin_fronts_both.csv        8 fronts, absolute z"
echo "  origin_fronts_centered.csv    8 fronts, centred on the slab"
echo "=== done ==="
