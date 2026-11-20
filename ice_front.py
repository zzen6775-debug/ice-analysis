#!/usr/bin/env python3
"""Tetrahedral-order q(z) analysis of an ice/water slab, per frame.

Computes, for every frame in a directory of .gro files:
  * mean tetrahedral order parameter q in 1 nm z bins
  * the ice band (bins with q > qcut), its thickness and both front positions
  * the polymer position (two chains) and the gap to each front

Writes two CSVs:
  <out>            one row per frame: fronts, thickness, polymer z, gaps
  <out>.qprofile   one row per (frame, z bin): q and ice fraction

Frame naming convention: <tag>_<time>ps.gro   (e.g. HEMA_20000ps.gro)

Usage:
  python3 ice_front.py --frames DIR --out fronts.csv [--bin 1.0] [--cut 0.35] [--qcut 0.85]
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree


def load_gro(path):
    """Read a .gro file using fixed-width columns (atom names can fuse otherwise)."""
    raw = open(path, errors="replace").read().splitlines()
    n = int(raw[1].split()[0])
    box = np.array([float(v) for v in raw[-1].split()[:3]])
    names, res, xyz = [], [], []
    for ln in raw[2:2 + n]:
        names.append(ln[10:15].strip())
        res.append(ln[5:10].strip())
        xyz.append((float(ln[20:28]), float(ln[28:36]), float(ln[36:44])))
    return np.array(xyz), np.array(names), np.array(res), box


def tetra_q(O, box, cut):
    """Tetrahedral order parameter q = 1 - 3/8 sum_{j<k} (cos(theta_jik) + 1/3)^2."""
    O = O % box                       # -pbc whole can leave atoms outside the box
    tree = cKDTree(O, boxsize=box)
    pairs = tree.query_ball_point(O, cut)
    q = np.full(len(O), np.nan)
    for i, nb in enumerate(pairs):
        nb = [j for j in nb if j != i]
        if len(nb) < 3:
            continue
        v = O[nb] - O[i]
        v -= box * np.round(v / box)
        nrm = np.linalg.norm(v, axis=1)
        nrm[nrm == 0] = 1.0
        v /= nrm[:, None]
        c = v @ v.T
        iu = np.triu_indices(len(nb), 1)
        q[i] = 1.0 - 3.0 / 8.0 * np.sum((c[iu] + 1.0 / 3.0) ** 2)
    return q


def ice_band(prof_q, qcut, binw):
    """Contiguous run of bins above qcut with the largest extent."""
    good = np.where(np.isfinite(prof_q) & (prof_q > qcut))[0]
    if good.size == 0:
        return None, None
    runs, start = [], good[0]
    for a, b in zip(good[:-1], good[1:]):
        if b != a + 1:
            runs.append((start, a))
            start = b
    runs.append((start, good[-1]))
    lo, hi = max(runs, key=lambda r: (r[1] - r[0], r[1]))
    return lo * binw, (hi + 1) * binw


def chain_centres(poly):
    """Split the polymer atoms into two chains by a z gap and return their medians."""
    if poly.size == 0:
        return None, None
    z = np.sort(poly[:, 2])
    gaps = np.diff(z)
    if gaps.size and gaps.max() > 3.0:
        cut = z[np.argmax(gaps)]
        a, b = poly[poly[:, 2] < cut], poly[poly[:, 2] >= cut]
        if a.size and b.size:
            return float(np.median(a[:, 2])), float(np.median(b[:, 2]))
    half = len(poly) // 2
    return (float(np.median(poly[:half, 2])), float(np.median(poly[half:, 2])))


def analyse(path, binw, cut, qcut):
    xyz, names, res, box = load_gro(path)
    tag = os.path.basename(path)[:-4]
    m = re.search(r"_(\d+)ps$", tag)
    t_ns = int(m.group(1)) / 1000.0 if m else float("nan")
    system = tag[:m.start()] if m else tag

    isow = names == "OW"
    O = xyz[isow] % box
    poly = xyz[res != "SOL"] % box
    q = tetra_q(O, box, cut)

    nbin = int(box[2] // binw)
    idx = np.clip((O[:, 2] // binw).astype(int), 0, nbin - 1)
    prof_q = np.full(nbin, np.nan)
    prof_i = np.full(nbin, np.nan)
    for b in range(nbin):
        m_b = idx == b
        if m_b.sum() >= 20:
            prof_q[b] = np.nanmean(q[m_b])
            prof_i[b] = np.mean(q[m_b] > qcut)

    lo, hi = ice_band(prof_q, qcut, binw)
    cA, cB = chain_centres(poly)
    row = dict(system=system, t_ns=t_ns, n_waters=int(isow.sum()),
               n_polymer=int((res != "SOL").sum()),
               box_z=round(float(box[2]), 4),
               ice_lo=lo, ice_hi=hi,
               thickness=(round(hi - lo, 3) if lo is not None else None),
               chainA_z=cA, chainB_z=cB,
               gapA=(round(cA - lo, 3) if (cA is not None and lo is not None) else None),
               gapB=(round(hi - cB, 3) if (cB is not None and hi is not None) else None))
    prof = [(system, t_ns, round(b * binw, 2), prof_q[b], prof_i[b]) for b in range(nbin)]
    return row, prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", required=True, help="directory of <tag>_<time>ps.gro")
    ap.add_argument("--out", default="fronts.csv")
    ap.add_argument("--bin", type=float, default=1.0, help="z bin width in nm")
    ap.add_argument("--cut", type=float, default=0.35, help="first-shell cutoff in nm")
    ap.add_argument("--qcut", type=float, default=0.85, help="ice threshold on q")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.frames, "*.gro")))
    if not files:
        sys.exit(f"no .gro files in {a.frames}")

    rows, profs = [], []
    for f in files:
        r, p = analyse(f, a.bin, a.cut, a.qcut)
        rows.append(r)
        profs.extend(p)
        print(f"  {r['system']:<10} t={r['t_ns']:6.2f} ns  ice {r['ice_lo']}–{r['ice_hi']} nm "
              f"({r['thickness']} nm)  chains {r['chainA_z']}/{r['chainB_z']}")

    cols = ["system", "t_ns", "n_waters", "n_polymer", "box_z", "ice_lo", "ice_hi",
            "thickness", "chainA_z", "chainB_z", "gapA", "gapB"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    with open(a.out + ".qprofile", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "t_ns", "z_nm", "q", "ice_fraction"])
        w.writerows(profs)
    print(f"\nwrote {a.out} and {a.out}.qprofile")


if __name__ == "__main__":
    main()
