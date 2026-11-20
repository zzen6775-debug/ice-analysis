#!/usr/bin/env python3
"""Fit ice-growth rates from fronts.csv produced by ice_front.py.

Two models per system:

  linear     thickness(t) = 2*v*t + t0          -> single growth velocity v per front
  piecewise  thickness(t) = 2*v1*t + t0 for t<t*,  2*v2*(t-t*) + ... for t>=t*
             (use --piecewise when the front has reached the polymer and the
              rate is expected to change)

Also fits each front separately (lower goes down, upper goes up) so you can see
whether growth is symmetric.

Usage:
  python3 growth_fit.py fronts.csv                 # linear + per-front
  python3 growth_fit.py fronts.csv --piecewise     # also try one breakpoint
  python3 growth_fit.py fronts.csv --plot growth.png
"""
import argparse
import csv
from collections import defaultdict

import numpy as np


def read_csv(path):
    data = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            def num(k):
                v = r.get(k, "")
                try:
                    return float(v)
                except (TypeError, ValueError):
                    return np.nan
            data[r["system"]].append((num("t_ns"), num("thickness"), num("ice_lo"),
                                      num("ice_hi"), num("chainA_z"), num("chainB_z")))
    for k in data:
        data[k].sort()
    return data


def linfit(x, y):
    x, y = np.asarray(x), np.asarray(y)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return np.nan, np.nan, np.nan
    p, cov = np.polyfit(x[ok], y[ok], 1, cov=True)
    resid = y[ok] - np.polyval(p, x[ok])
    ss_res = (resid ** 2).sum()
    ss_tot = ((y[ok] - y[ok].mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return p[0], float(np.sqrt(cov[0, 0])), r2


def piecewise_fit(t, thick):
    """Fit thickness = 2*v1*t + c for t<t*, and 2*v2*(t-t*) + thick(t*) after."""
    best = None
    ts = np.asarray(t)
    for tstar in np.unique(ts)[1:-1]:
        left = ts < tstar
        right = ~left
        if left.sum() < 2 or right.sum() < 2:
            continue
        s1, _, _ = linfit(ts[left], thick[left])
        s2, _, _ = linfit(ts[right], thick[right])
        if not np.isfinite(s1) or not np.isfinite(s2):
            continue
        # residual of the continuous broken line
        c1 = np.nanmean(thick[left] - s1 * ts[left])
        pred = np.where(left, s1 * ts + c1, s2 * (ts - tstar) + (s1 * tstar + c1))
        err = np.nansum((thick - pred) ** 2)
        if best is None or err < best[0]:
            best = (err, tstar, s1, s2)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--piecewise", action="store_true")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--tube", type=float, default=0.5,
                    help="only fit t >= this value (skip the initial transient)")
    a = ap.parse_args()

    data = read_csv(a.csv)
    print(f"{'system':<10}{'v_total (nm/ns)':>18}{'v_per_front':>14}{'R^2':>8}"
          f"{'v_lower':>10}{'v_upper':>10}{'n':>4}")
    print("-" * 76)
    results = {}
    for sys_, rows in data.items():
        t = np.array([r[0] for r in rows])
        th = np.array([r[1] for r in rows])
        lo = np.array([r[2] for r in rows])
        hi = np.array([r[3] for r in rows])
        sel = t >= a.tube
        slope, err, r2 = linfit(t[sel], th[sel])
        v_lo, _, _ = linfit(t[sel], lo[sel])      # should be negative
        v_hi, _, _ = linfit(t[sel], hi[sel])
        print(f"{sys_:<10}{slope:>14.4f}±{err:<5.3f}{slope/2:>12.4f}{r2:>8.3f}"
              f"{-v_lo:>10.4f}{v_hi:>10.4f}{sel.sum():>4}")
        results[sys_] = (t, th, slope)
        if a.piecewise:
            best = piecewise_fit(t, th)
            if best:
                _, tstar, s1, s2 = best
                print(f"{'':10}  piecewise: v_total {s1:.4f} nm/ns before t*={tstar:.1f} ns, "
                      f"{s2:.4f} nm/ns after")

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available, skipping plot")
            return
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
        for sys_, (t, th, slope) in results.items():
            ax.plot(t, th, "o-", ms=4, label=f"{sys_} ({slope*1000:.0f} pm/ns)")
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("ice thickness (nm)")
        ax.set_title("Ice slab growth at 260 K")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"wrote {a.plot}")


if __name__ == "__main__":
    main()
