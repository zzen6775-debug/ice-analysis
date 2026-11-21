#!/usr/bin/env python3
"""Centre the two fronts on each system's own slab centre.

Absolute z is meaningless across systems: the box length differs by several nm
after the cooling stage (67.8 - 73.6 nm), so the same slab sits at a different
absolute z in each. Subtracting the t=0 slab centre puts all four systems on a
common reference, where the lower front is negative and the upper one positive.

Usage:
  python3 centre_fronts.py origin_fronts_both.csv --out origin_fronts_centered.csv
"""
import argparse
import csv
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", default="origin_fronts_centered.csv")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.csv)))
    systems = sorted({c[:-len("_lower")] for c in rows[0] if c.endswith("_lower")})

    # slab centre at the earliest time, per system
    first = min(rows, key=lambda r: float(r["t_ns"]))
    centre = {}
    for s in systems:
        lo = float(first[f"{s}_lower"])
        hi = float(first[f"{s}_upper"])
        centre[s] = 0.5 * (lo + hi)
    print("slab centre at t=0 (nm):")
    for s in systems:
        print(f"  {s:<8} {centre[s]:.3f}")

    head = ["t_ns"]
    for s in systems:
        head += [f"{s}_lower", f"{s}_upper"]
    with open(a.out, "w", newline="\n") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerow(["# z relative to each system's t=0 slab centre (nm)"])
        for r in sorted(rows, key=lambda r: float(r["t_ns"])):
            line = [r["t_ns"]]
            for s in systems:
                for key in (f"{s}_lower", f"{s}_upper"):
                    v = r.get(key, "")
                    line.append(f"{float(v) - centre[s]:.3f}" if v not in ("", None) else "")
            w.writerow(line)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
