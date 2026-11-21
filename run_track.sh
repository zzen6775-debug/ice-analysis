#!/bin/bash
# Anchored-tracking front analysis + the 8-curve (both fronts) Origin table.
set -uo pipefail
ANA=/mnt/g/MD_useful/ice_analysis
CACHE=${CACHE:-$HOME/ice_cache}

[ -d "$CACHE/split" ] || { echo "frame cache missing at $CACHE/split"; exit 1; }

echo "=== $(date) ==="
echo "-- anchored tracking (0.2 nm grid, +/-1.2 nm window)"
python3 "$ANA/fronts_track.py" --frames "$CACHE/split" --bin 0.2 --window 1.2 \
    --out "$CACHE/fronts_track.csv" 2>&1 | tail -6

echo
echo "-- linear fits"
python3 "$ANA/growth_fit.py" "$CACHE/fronts_track.csv" --tube 4 --windows 0 20 40 60 2>&1

echo
echo "-- export for Origin"
python3 "$ANA/export_for_origin.py" "$CACHE/fronts_track.csv" --outdir "$ANA/origin_data_track" 2>&1 | head -8

cp "$CACHE/fronts_track.csv" "$ANA/fronts_track.csv" 2>/dev/null
echo
echo "-- sanity: how many distinct front positions per system"
python3 - "$CACHE/fronts_track.csv" <<'PY'
import csv, sys, collections
d = collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[1])):
    try:
        d[r['system']].append(float(r['ice_lo']))
    except (TypeError, ValueError):
        pass
for s in sorted(d):
    v = d[s]
    jumps = [abs(v[i+1]-v[i]) for i in range(len(v)-1)]
    print(f"  {s:<8} n={len(v):3d}  distinct={len(set(round(x,2) for x in v)):3d}  "
          f"max step={max(jumps):.3f} nm  total move={v[0]-v[-1]:.2f} nm")
PY
