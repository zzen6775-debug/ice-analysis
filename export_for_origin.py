#!/usr/bin/env python3
"""Turn fronts.csv into Origin-friendly wide tables.

Origin imports these directly (File > Import > CSV); the first column is X
(time in ns), every other column is one system.

Usage:
  python3 export_for_origin.py fronts.csv [--outdir .]
"""
import argparse
import csv
import os
from collections import defaultdict


def load(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            def num(k):
                try:
                    return float(r[k])
                except (TypeError, ValueError):
                    return float("nan")
            rows.append(dict(system=r["system"], t=num("t_ns"),
                             thick=num("thickness"), lo=num("ice_lo"), hi=num("ice_hi"),
                             amount=num("ice_amount"),
                             cA=num("chainA_z"), cB=num("chainB_z"),
                             gapA=num("gapA"), gapB=num("gapB")))
    rows.sort(key=lambda r: (r["system"], r["t"]))
    return rows


def pivot(rows, key):
    systems = sorted({r["system"] for r in rows})
    times = sorted({r["t"] for r in rows})
    table = {s: {r["t"]: r[key] for r in rows if r["system"] == s} for s in systems}
    return systems, times, table


def write_wide(rows, key, path, header_note):
    systems, times, table = pivot(rows, key)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([header_note] + systems)
        for t in times:
            w.writerow([f"{t:g}"] + [("" if table[s].get(t) is None else f"{table[s][t]:.4f}")
                                     for s in systems])
    print(f"  {path}   ({len(times)} rows x {len(systems)} systems)")


def load_profiles(path):
    prof = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                prof[(r["system"], float(r["t_ns"]))].append(
                    (float(r["z_nm"]), float(r["q"]), float(r["ice_fraction"])))
            except (TypeError, ValueError):
                continue
    return prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rows = load(a.csv)
    print("writing Origin tables:")

    write_wide(rows, "thick", os.path.join(a.outdir, "origin_thickness.csv"), "t_ns")
    write_wide(rows, "amount", os.path.join(a.outdir, "origin_ice_amount.csv"), "t_ns")
    write_wide(rows, "lo", os.path.join(a.outdir, "origin_front_lower.csv"), "t_ns")
    write_wide(rows, "hi", os.path.join(a.outdir, "origin_front_upper.csv"), "t_ns")
    write_wide(rows, "gapA", os.path.join(a.outdir, "origin_gap_chainA.csv"), "t_ns")
    write_wide(rows, "gapB", os.path.join(a.outdir, "origin_gap_chainB.csv"), "t_ns")

    # both fronts of every system in one table: 8 series against time
    systems, times, _ = pivot(rows, "lo")
    table = {}
    for r in rows:
        table.setdefault(r["system"], {})[r["t"]] = (r["lo"], r["hi"])
    with open(os.path.join(a.outdir, "origin_fronts_both.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        head = ["t_ns"]
        for s in systems:
            head += [f"{s}_lower", f"{s}_upper"]
        w.writerow(head)
        for t in times:
            line = [f"{t:g}"]
            for s in systems:
                pair = table[s].get(t)
                if pair is None:
                    line += ["", ""]
                else:
                    line += [f"{pair[0]:.3f}", f"{pair[1]:.3f}"]
            w.writerow(line)
    print(f"  origin_fronts_both.csv   ({len(times)} rows x {2*len(systems)} series)")

    prof_path = a.csv + ".qprofile"
    if os.path.exists(prof_path):
        prof = load_profiles(prof_path)
        for t in sorted({k[1] for k in prof}):
            out = os.path.join(a.outdir, f"origin_qprofile_{t:g}ns.csv")
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh)
                systems = sorted({k[0] for k in prof if k[1] == t})
                w.writerow(["z_nm"] + systems)
                zs = sorted({z for k in prof if k[1] == t for z, _, _ in prof[k]})
                table = {s: {z: q for z, q, _ in prof.get((s, t), [])} for s in systems}
                for z in zs:
                    w.writerow([f"{z:g}"] + [("" if table[s].get(z) is None else f"{table[s][z]:.4f}")
                                             for s in systems])
            print(f"  {out}")

    print("\nOrigin 用法：File > Import > CSV 导入；第一列设为 X，"
          "其余列 Plot > Line+Symbol；拟合用 Analysis > Fitting > Linear Fit。")


if __name__ == "__main__":
    main()
