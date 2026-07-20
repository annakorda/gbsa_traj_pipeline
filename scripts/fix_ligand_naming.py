#!/usr/bin/env python3
"""
Fix a generic SDF (e.g. downloaded from PubChem) so it can be used to
build a ligand topology: assigns unique, indexed atom names in the same
convention as typical MD/PDB ligand naming (C1, C2, ..., H1, H2, ..., O1, ...),
in the order atoms appear in the file, and writes a proper mol2 with the
SDF's own bond table (not re-derived from geometry).

Pure python, no openbabel/rdkit dependency - a plain SDF (MDL V2000) atom
block has no name field at all, which is why a raw obabel conversion gives
you bare, non-unique element symbols instead of indexed names.

IMPORTANT CAVEAT (see NOTES.md): this only makes names unique and gives
them the right *style*. It does NOT prove atom N in this SDF is the same
physical atom as the identically-styled name in your MD frames - that
requires the file's atom order to already correspond to your frame's
numbering (e.g. because both trace back to the same PDB Chemical Component
Dictionary / CGenFF definition). Check the match-count printout from
build_trajectory.py before trusting any downstream number.
"""
import argparse

BOND_TYPE = {"1": "1", "2": "2", "3": "3", "4": "ar"}


def parse_sdf_v2000(path):
    with open(path) as f:
        lines = f.readlines()

    counts_line = lines[3]
    n_atoms = int(counts_line[0:3])
    n_bonds = int(counts_line[3:6])

    atoms = []
    for line in lines[4:4 + n_atoms]:
        x = float(line[0:10])
        y = float(line[10:20])
        z = float(line[20:30])
        element = line[31:34].strip()
        atoms.append({"x": x, "y": y, "z": z, "element": element})

    bonds = []
    for line in lines[4 + n_atoms:4 + n_atoms + n_bonds]:
        a1 = int(line[0:3])
        a2 = int(line[3:6])
        btype = line[6:9].strip()
        bonds.append((a1, a2, BOND_TYPE.get(btype, "1")))

    return atoms, bonds


def assign_names(atoms):
    counters = {}
    for a in atoms:
        el = a["element"]
        counters[el] = counters.get(el, 0) + 1
        a["name"] = f"{el}{counters[el]}"
    return atoms


def write_mol2(atoms, bonds, outfile, molname="MOL"):
    with open(outfile, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{molname}\n")
        f.write(f" {len(atoms)} {len(bonds)} 1 0 0\n")
        f.write("SMALL\nUSER_CHARGES\n\n")

        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i:>7} {a['name']:<8} {a['x']:>10.4f} {a['y']:>10.4f} {a['z']:>10.4f} "
                    f"{a['element']:<5} 1 {molname:<8} 0.0000\n")

        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, btype) in enumerate(bonds, 1):
            f.write(f"{i:>6} {a1:>5} {a2:>5} {btype:>4}\n")

        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"1 {molname} 1\n")


def main():
    parser = argparse.ArgumentParser(description="Assign unique, indexed atom names to an SDF and write a named mol2.")
    parser.add_argument("-i", dest="infile", required=True, help="Input SDF (MDL V2000).")
    parser.add_argument("-o", dest="outfile", required=True, help="Output mol2 with unique atom names.")
    parser.add_argument("-name", dest="molname", default="MOL", help="Residue/molecule name (default: MOL, unigbsa's own default).")
    args = parser.parse_args()

    atoms, bonds = parse_sdf_v2000(args.infile)
    assign_names(atoms)
    write_mol2(atoms, bonds, args.outfile, molname=args.molname)

    counts = {}
    for a in atoms:
        counts[a["element"]] = counts.get(a["element"], 0) + 1
    print(f"Wrote {args.outfile}: {len(atoms)} atoms {counts}, {len(bonds)} bonds")
    print("Names assigned in file order: " + ", ".join(a["name"] for a in atoms))


if __name__ == "__main__":
    main()
