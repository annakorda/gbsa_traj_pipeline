#!/usr/bin/env python3
"""
Remap N raw per-frame complex PDBs (protein+ligand, from whatever MD
engine) onto the fixed atom order of the reference topology built by
build_topology_once.py, and concatenate them into a multi-MODEL
trajectory PDB for unigbsa-traj.

Matching rule:
  - protein atoms: matched by (chain, resnum, atom name)
  - ligand atoms:   matched by atom name only (there's one ligand copy)

Any reference atom with no match in a given frame (typically the
receptor hydrogens, which pdb2gmx rebuilt fresh and which never
existed in the raw MD frame) keeps its reference-build coordinate.
Check the printed match counts per frame before trusting the output -
a low match count on the protein or ligand means the naming/numbering
assumption below doesn't hold for your files.
"""
import argparse


def parse_pdb_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            rec = line[0:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            atoms.append({
                "line": line.rstrip("\n"),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resseq": line[22:26].strip(),
                "x": line[30:38], "y": line[38:46], "z": line[46:54],
            })
    return atoms


def remap_frame(reference_atoms, frame_atoms, ligand_resname, ref_ligand_resname):
    protein_lookup = {}
    ligand_lookup = {}
    for a in frame_atoms:
        if a["resname"].upper() == ligand_resname.upper():
            ligand_lookup[a["name"].upper()] = a
        else:
            protein_lookup[(a["chain"], a["resseq"], a["name"].upper())] = a

    out_lines = []
    n_matched = 0
    n_kept = 0
    for ref in reference_atoms:
        if ref["resname"].upper() == ref_ligand_resname.upper():
            src = ligand_lookup.get(ref["name"].upper())
        else:
            src = protein_lookup.get((ref["chain"], ref["resseq"], ref["name"].upper()))
        line = ref["line"]
        if src:
            line = line[:30] + src["x"] + src["y"] + src["z"] + line[54:]
            n_matched += 1
        else:
            n_kept += 1
        out_lines.append(line)
    return out_lines, n_matched, n_kept


def main():
    parser = argparse.ArgumentParser(description="Build a multi-frame trajectory PDB matched to a reference topology.")
    parser.add_argument("-ref", required=True, help="Reference complex.pdb from build_topology_once.py")
    parser.add_argument("-frames", nargs="+", required=True, help="Raw per-frame complex PDBs, in order")
    parser.add_argument("-ligand-resname", required=True, dest="ligand_resname",
                         help="Ligand residue name as it appears in the RAW frame files")
    parser.add_argument("-ref-ligand-resname", default="MOL", dest="ref_ligand_resname",
                         help="Ligand residue name in the built reference topology (unigbsa default: MOL)")
    parser.add_argument("-o", dest="outfile", default="trajectory.pdb")
    args = parser.parse_args()

    reference_atoms = parse_pdb_atoms(args.ref)
    print(f"Reference has {len(reference_atoms)} atoms (this defines the fixed atom order).")

    with open(args.outfile, "w") as out:
        for i, framefile in enumerate(args.frames, 1):
            frame_atoms = parse_pdb_atoms(framefile)
            lines, n_matched, n_kept = remap_frame(
                reference_atoms, frame_atoms, args.ligand_resname, args.ref_ligand_resname
            )
            print(f"Frame {i} ({framefile}): {n_matched} matched, {n_kept} kept-from-reference "
                  f"(expect kept ~= receptor H count only)")
            out.write(f"MODEL{i:>9d}\n")
            for line in lines:
                out.write(line + "\n")
            out.write("ENDMDL\n")
        out.write("END\n")

    print(f"Wrote {len(args.frames)}-frame trajectory: {args.outfile}")


if __name__ == "__main__":
    main()
