#!/usr/bin/env python3
"""
Fix a generic SDF (e.g. downloaded from GPCRmd/PubChem) so it can be used
to build a ligand topology: writes a mol2 with the SDF's own bond table
(not re-derived from geometry) and generic sequential atom names (O1, C1,
C2, ...). Also writes two sidecars next to the output mol2:
  <outfile>.net_charge  - net formal charge, from the SDF's own M CHG
                          records (see NOTES.md).
  <outfile>.atom_map    - generic name -> real raw-frame atom name, only
                          written if -frame-pdb is given and direct name
                          matching against it isn't already 1:1 (see
                          NOTES.md). Consumed by build_trajectory_v2.py's
                          -ligand-atom-map.
"""
import argparse
import re

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

    # M CHG property lines: one or more "(atom_index, charge)" pair lists.
    net_charge = 0
    for line in lines[4 + n_atoms + n_bonds:]:
        tokens = line.split()
        if len(tokens) >= 2 and tokens[0] == "M" and tokens[1] == "CHG":
            pairs = tokens[3:]
            net_charge += sum(int(c) for c in pairs[1::2])
        elif tokens[:2] == ["M", "END"]:
            break

    return atoms, bonds, net_charge


def assign_names(atoms):
    counters = {}
    for a in atoms:
        el = a["element"]
        counters[el] = counters.get(el, 0) + 1
        a["name"] = f"{el}{counters[el]}"
    return atoms


def parse_ligand_from_pdb(pdb_path, resname):
    """Atom name/element/coords for one residue name in a raw MD frame PDB."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line[0:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip() != resname:
                continue
            name = line[12:16].strip()
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                m = re.match(r"[A-Za-z]+", name)
                element = m.group(0)[0] if m else "C"
            atoms.append({
                "name": name, "element": element.strip().capitalize(),
                "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
            })
    return atoms


def direct_name_map(atoms, frame_atoms):
    """Cheap path: if the generic names already exist verbatim in the frame,
    no graph matching is needed at all. Returns None if not a full match."""
    frame_names = {a["name"].upper() for a in frame_atoms}
    if all(a["name"].upper() in frame_names for a in atoms):
        return {a["name"]: a["name"] for a in atoms}
    return None


def _graph_only_mol(elements, bond_pairs):
    """Element identity + connectivity only - no valence/charge chemistry,
    which is all a topology match needs."""
    from rdkit import Chem
    rw = Chem.RWMol()
    for el in elements:
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    for a1, a2 in bond_pairs:
        rw.AddBond(a1, a2, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^
                      Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_PROPERTIES)
    return mol


def _perceive_bonds_from_geometry(frame_atoms):
    """Connectivity from real 3D coordinates via RDKit's distance-based
    bond perception (rdDetermineBonds) - no OpenBabel."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(frame_atoms))
    for i, a in enumerate(frame_atoms):
        rw.AddAtom(Chem.Atom(a["element"]))
        conf.SetAtomPosition(i, Point3D(a["x"], a["y"], a["z"]))
    mol = rw.GetMol()
    mol.AddConformer(conf)
    rdDetermineBonds.DetermineConnectivity(mol)
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]


def match_atoms_to_frame(sdf_elements, sdf_bond_pairs, frame_atoms):
    """Fallback for when direct_name_map() fails: match the SDF's authoritative
    bond graph against the frame's geometry-perceived one (element +
    connectivity only) via RDKit substructure matching. Returns, for each SDF
    atom index, the matching index into frame_atoms. Raises rather than
    returning a partial match - see NOTES.md."""
    template = _graph_only_mol(sdf_elements, sdf_bond_pairs)

    frame_elements = [a["element"] for a in frame_atoms]
    frame_bond_pairs = _perceive_bonds_from_geometry(frame_atoms)
    frame_mol = _graph_only_mol(frame_elements, frame_bond_pairs)

    if template.GetNumAtoms() != frame_mol.GetNumAtoms():
        raise ValueError(
            f"SDF has {template.GetNumAtoms()} atoms but the frame ligand has "
            f"{frame_mol.GetNumAtoms()} - refusing to guess a name mapping."
        )

    match = frame_mol.GetSubstructMatch(template, useChirality=False)
    if len(match) != template.GetNumAtoms():
        raise ValueError(
            f"Could not match the SDF's bond graph onto the frame ligand's "
            f"geometry-perceived bond graph ({len(match)}/{template.GetNumAtoms()} "
            f"atoms matched) - refusing to guess. Investigate before trusting "
            f"this ligand's results."
        )
    return list(match)


def resolve_atom_map(atoms, bonds, frame_atoms):
    """Direct name match first; only fall back to graph isomorphism if that
    doesn't fully resolve. Returns (atom_map, used_fallback)."""
    direct = direct_name_map(atoms, frame_atoms)
    if direct is not None:
        return direct, False

    sdf_elements = [a["element"] for a in atoms]
    sdf_bond_pairs = [(a1 - 1, a2 - 1) for a1, a2, _bt in bonds]
    match = match_atoms_to_frame(sdf_elements, sdf_bond_pairs, frame_atoms)
    atom_map = {atoms[i]["name"]: frame_atoms[match[i]]["name"] for i in range(len(atoms))}
    return atom_map, True


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
    parser = argparse.ArgumentParser(description="Assign atom names to an SDF ligand and write a named mol2.")
    parser.add_argument("-i", dest="infile", required=True, help="Input SDF (MDL V2000).")
    parser.add_argument("-o", dest="outfile", required=True, help="Output mol2.")
    parser.add_argument("-name", dest="molname", default="MOL", help="Residue/molecule name (default: MOL).")
    parser.add_argument("-frame-pdb", dest="frame_pdb", default=None,
                         help="A raw MD frame PDB containing this ligand, used to verify/resolve real atom names.")
    parser.add_argument("-frame-resname", dest="frame_resname", default=None,
                         help="Ligand residue name in -frame-pdb (default: same as -name).")
    args = parser.parse_args()

    atoms, bonds, net_charge = parse_sdf_v2000(args.infile)
    assign_names(atoms)
    write_mol2(atoms, bonds, args.outfile, molname=args.molname)

    if args.frame_pdb:
        frame_resname = args.frame_resname or args.molname
        frame_atoms = parse_ligand_from_pdb(args.frame_pdb, frame_resname)
        if not frame_atoms:
            raise SystemExit(f"No atoms with resname '{frame_resname}' found in {args.frame_pdb}.")
        atom_map, used_fallback = resolve_atom_map(atoms, bonds, frame_atoms)
        map_file = args.outfile + ".atom_map"
        with open(map_file, "w") as f:
            for generic_name, real_name in atom_map.items():
                f.write(f"{generic_name} {real_name}\n")
        method = "graph-isomorphism fallback" if used_fallback else "direct name match"
        print(f"Wrote {map_file}: {len(atom_map)}/{len(atoms)} atoms resolved via {method}.")
    else:
        print("WARNING: no -frame-pdb given, no atom_map written - see NOTES.md before trusting "
              "downstream ligand matching.")

    charge_file = args.outfile + ".net_charge"
    with open(charge_file, "w") as f:
        f.write(f"{net_charge}\n")

    counts = {}
    for a in atoms:
        counts[a["element"]] = counts.get(a["element"], 0) + 1
    print(f"Wrote {args.outfile}: {len(atoms)} atoms {counts}, {len(bonds)} bonds")
    print(f"Wrote {charge_file}: net formal charge {net_charge:+d}")


if __name__ == "__main__":
    main()
