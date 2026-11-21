#!/usr/bin/env python3
"""Sub-bin front positions with anchored tracking.

The naive "longest contiguous ice band" criterion jumps between candidate bands
when the profile is noisy on a fine grid. Here each system's frames are walked in
time order and every front is searched only inside a small window around its
previous position, which makes the positions continuous without quantising them.

Also reports `ice_amount` (integral of the ice fraction), which needs no front at
all and is the most robust growth measure.

Usage:
  python3 fronts_track.py --frames DIR --out fronts_track.csv [--bin 0.2] [--window 1.2]
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ice_front import load_gro, tetra_q            # noqa: E402


def profile(O, q, box_z, binw, qcut, mincount=8):
    nbin = int(box_z // binw)
    idx = np.clip((O[:, 2] // binw).astype(int), 0, nbin - 1)
    prof = np.full(nbin, np.nan)
    frac = np.full(nbin, np.nan)
    for b in range(nbin):
        sel = idx == b
        if sel.sum() >= mincount:
            prof[b] = np.nanmean(q[sel])
            frac[b] = np.mean(q[sel] > qcut)
    return prof, frac


def edge_near(prof, qcut, binw, guess, window, lower):
    """Find the ice/water boundary nearest to `guess`, inside +/- window."""
    lo_b = max(0, int((guess - window) / binw))
    hi_b = min(len(prof) - 1, int((guess + window) / binw))
    seg = prof[lo_b:hi_b + 1]
    ok = np.isfinite(seg)
    if ok.sum() < 3:
        return None
    inside = np.where(ok & (seg > qcut))[0]
    if inside.size == 0:
        return None

    # bin inside the ice nearest to the guess (convert back to an absolute bin)
    rel = inside[np.argmin(np.abs((lo_b + inside) * binw + binw / 2 - guess))]
    b_in = int(lo_b + rel)
    if lower:                      # walk down while still ice
        b_edge = b_in
        while b_edge - 1 >= 0 and np.isfinite(prof[b_edge - 1]) and prof[b_edge - 1] > qcut:
            b_edge -= 1
        b_out = b_edge - 1
        if b_out < 0 or not np.isfinite(prof[b_out]):
            return float(b_edge * binw)
    else:                          # walk up while still ice
        b_edge = b_in
        while b_edge + 1 < len(prof) and np.isfinite(prof[b_edge + 1]) and prof[b_edge + 1] > qcut:
            b_edge += 1
        b_out = b_edge + 1
        if b_out >= len(prof) or not np.isfinite(prof[b_out]):
            return float((b_edge + 1) * binw)

    q_in, q_out = prof[b_edge], prof[b_out]
    if q_out == q_in:
        return float((b_edge + 0.5) * binw)
    frac = (qcut - q_out) / (q_in - q_out)     # 0 at b_out, 1 at b_edge
    z_out = (b_out + 0.5) * binw
    z_in = (b_edge + 0.5) * binw
    return float(z_out + frac * (z_in - z_out))


def analyse_system(files, binw, cut, qcut, window, mincount=8):
    files = sorted(files, key=lambda p: int(re.search(r"_(\d+)ps", p).group(1)))
    out, prev_lo, prev_hi = [], None, None
    for path in files:
        xyz, names, res, box = load_gro(path)
        tag = os.path.basename(path)[:-4]
        m = re.search(r"_(\d+)ps$", tag)
        t_ns = int(m.group(1)) / 1000.0
        system = tag[:m.start()]

        O = xyz[names == "OW"] % box
        q = tetra_q(O, box, cut)
        prof, frac = profile(O, q, box[2], binw, qcut, mincount)
        amount = float(np.nansum(frac[np.isfinite(frac)]) * binw)

        if prev_lo is None:
            # bootstrap from the longest contiguous run of ice-like bins,
            # not from the first/last one (stray high-q bins exist near the box edges)
            good = np.where(np.isfinite(prof) & (prof > qcut))[0]
            if good.size == 0:
                lo = hi = None
            else:
                runs, start = [], good[0]
                for a, b in zip(good[:-1], good[1:]):
                    if b != a + 1:
                        runs.append((start, a))
                        start = b
                runs.append((start, good[-1]))
                lo_r, hi_r = max(runs, key=lambda r: (r[1] - r[0], r[1]))
                lo = float((lo_r + 0.5) * binw)
                hi = float((hi_r + 0.5) * binw)
        else:
            lo = edge_near(prof, qcut, binw, prev_lo, window, lower=True)
            hi = edge_near(prof, qcut, binw, prev_hi, window, lower=False)
            if lo is None:
                lo = prev_lo
            if hi is None:
                hi = prev_hi

        if lo is not None:
            prev_lo, prev_hi = lo, hi

        out.append(dict(system=system, t_ns=t_ns, n_waters=int((names == "OW").sum()),
                        box_z=round(float(box[2]), 4),
                        ice_lo=round(lo, 3) if lo is not None else None,
                        ice_hi=round(hi, 3) if hi is not None else None,
                        thickness=round(hi - lo, 3) if (lo is not None and hi is not None) else None,
                        ice_amount=round(amount, 3),
                        chainA_z=None, chainB_z=None, gapA=None, gapB=None))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", required=True)
    ap.add_argument("--out", default="fronts_track.csv")
    ap.add_argument("--bin", type=float, default=0.2)
    ap.add_argument("--cut", type=float, default=0.35)
    ap.add_argument("--qcut", type=float, default=0.85)
    ap.add_argument("--window", type=float, default=1.2, help="search half-window in nm")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.frames, "*.gro")))
    if not files:
        sys.exit(f"no .gro files in {a.frames}")

    by_system = {}
    for f in files:
        m = re.search(r"_(\d+)ps\.gro$", os.path.basename(f))
        if m:
            by_system.setdefault(os.path.basename(f)[:m.start()], []).append(f)

    rows = []
    for system, fl in sorted(by_system.items()):
        rows.extend(analyse_system(fl, a.bin, a.cut, a.qcut, a.window))
        print(f"  {system}: {len(fl)} frames")

    cols = ["system", "t_ns", "n_waters", "box_z", "ice_lo", "ice_hi",
            "thickness", "ice_amount", "chainA_z", "chainB_z", "gapA", "gapB"]
    with open(a.out, "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
