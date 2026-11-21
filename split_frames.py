#!/usr/bin/env python3
"""Split multi-frame .gro files (as written by trjconv -dt) into one file per frame.

Output names follow the convention ice_front.py expects: <tag>_<time>ps.gro

Usage:
  python3 split_frames.py <indir> <outdir>
"""
import glob
import os
import re
import sys


def split(path, outdir):
    lines = open(path, errors="replace").read().splitlines()
    tag = re.sub(r"^frames_", "", os.path.splitext(os.path.basename(path))[0])
    i, n = 0, 0
    while i < len(lines):
        title = lines[i]
        try:
            nat = int(lines[i + 1].split()[0])
        except (IndexError, ValueError):
            break
        block = lines[i:i + nat + 3]
        if len(block) < nat + 3:
            break
        m = re.search(r"t=\s*([0-9.]+)", title)
        t_ps = int(round(float(m.group(1)))) if m else n * 2000
        out = os.path.join(outdir, f"{tag}_{t_ps}ps.gro")
        with open(out, "w", newline="\n") as fh:
            fh.write("\n".join(block) + "\n")
        n += 1
        i += nat + 3
    print(f"  {os.path.basename(path)}: {n} frames -> {outdir}")
    return n


if __name__ == "__main__":
    indir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    total = 0
    for f in sorted(glob.glob(os.path.join(indir, "frames_*.gro"))):
        total += split(f, outdir)
    print(f"total frames: {total}")
