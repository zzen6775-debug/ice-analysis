#!/bin/bash
# One-shot pipeline: fetch a frames tarball from the cluster, extract the
# frames, compute q(z) front positions, then fit the growth rates.
#
# Usage:
#   bash run_analysis.sh <frames.tar.gz> [workdir]
#
#   <frames.tar.gz>  can be a local file OR an scp-style remote path
#   workdir          scratch directory, default /tmp/icework
#
# Cluster settings are read from the environment, so nothing is hard-coded:
#   HPC_HOST   ssh target, e.g. user@host      (required only for remote fetch)
#   HPC_KEY    path to the private key         (optional)
#
# Example:
#   HPC_HOST=xujianning@202.114.96.180 \
#   HPC_KEY=/path/to/key \
#   bash run_analysis.sh /scratch/user/zsmx/_frames/frames.tar.gz
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARBALL_IN="${1:-}"
WORK="${2:-/tmp/icework}"

if [ -z "$TARBALL_IN" ]; then
  echo "usage: bash run_analysis.sh <frames.tar.gz> [workdir]" >&2
  exit 2
fi

echo "=== $(date) ==="
mkdir -p "$WORK"

# ---- obtain the tarball -----------------------------------------------------
if [ -f "$TARBALL_IN" ]; then
  echo "-- using local tarball $TARBALL_IN"
  cp "$TARBALL_IN" "$WORK/frames.tar.gz"
else
  : "${HPC_HOST:?set HPC_HOST=user@host to fetch $TARBALL_IN}"
  echo "-- fetching $TARBALL_IN from $HPC_HOST"
  if [ -n "${HPC_KEY:-}" ]; then
    cp "$HPC_KEY" "$WORK/k" && chmod 600 "$WORK/k"
    KEYOPT="-i $WORK/k"
  else
    KEYOPT=""
  fi
  # shellcheck disable=SC2086
  ssh $KEYOPT -o BatchMode=yes -o StrictHostKeyChecking=accept-new "$HPC_HOST" \
      "cat $TARBALL_IN" > "$WORK/frames.tar.gz" 2>/dev/null
fi

SZ=$(stat -c %s "$WORK/frames.tar.gz" 2>/dev/null || echo 0)
echo "   size: $((SZ/1024)) KB"
[ "$SZ" -gt 1000 ] || { echo "tarball is empty or missing" >&2; exit 1; }

# ---- extract ----------------------------------------------------------------
rm -rf "$WORK/frames"; mkdir -p "$WORK/frames"
tar xzf "$WORK/frames.tar.gz" -C "$WORK/frames" || exit 1
N=$(find "$WORK/frames" -name '*.gro' | wc -l)
echo "-- frames: $N"
[ "$N" -gt 0 ] || { echo "no .gro files inside" >&2; exit 1; }

# ---- analyse ----------------------------------------------------------------
echo "-- computing q(z)"
python3 "$HERE/ice_front.py" --frames "$WORK/frames" --out "$WORK/fronts.csv" || exit 1

# ---- fit --------------------------------------------------------------------
echo
echo "-- fitting growth rates"
python3 "$HERE/growth_fit.py" "$WORK/fronts.csv"

echo
echo "results: $WORK/fronts.csv and $WORK/fronts.csv.qprofile"
echo "next:    python3 $HERE/export_for_origin.py $WORK/fronts.csv --outdir $WORK/origin_data"
echo "=== done ==="
