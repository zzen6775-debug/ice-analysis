"""Generate GROMACS index files selecting only the water oxygens.

The topology order is SOL first (OW, HW1, HW2, MW repeating), polymer last,
so the OW atoms are 1, 5, 9, ... up to 4*N_waters - 3.
"""
import os

SYSTEMS = {
    "HEMA":   (23070, "HEMA_on_ice_slab"),
    "MMA":    (23077, "MMA_on_ice_slab"),
    "PMMA10": (22869, "PMMA_on_ice_slab"),
    "PMMA20": (22655, "PMMA20_on_ice_slab"),
}

OUT = r"G:\MD_useful\_frames2"
os.makedirs(OUT, exist_ok=True)

for tag, (nwat, _dir) in SYSTEMS.items():
    path = os.path.join(OUT, f"ow_{tag}.ndx")
    atoms = list(range(1, 4 * nwat - 3 + 1, 4))
    with open(path, "w", newline="\n") as fh:
        fh.write("[ OW ]\n")
        for i in range(0, len(atoms), 15):
            fh.write(" ".join(f"{a:5d}" for a in atoms[i:i + 15]) + "\n")
    print(f"{path}: {len(atoms)} atoms, last={atoms[-1]}")
