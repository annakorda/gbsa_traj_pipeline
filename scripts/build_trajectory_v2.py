#!/usr/bin/env python3
"""
Remap N raw per-frame complex PDBs (protein+ligand) onto the fixed atom
order of the reference topology built by build_topology_once.py, and
concatenate them into a multi-MODEL trajectory PDB for unigbsa-traj.

Matching: protein atoms by (residue ordinal position, atom name) - not
literal chain/resnum, since pdb2gmx drops/renumbers those. Ligand atoms
by name, translated through -ligand-atom-map first if given (see
NOTES.md - SDF and raw-frame ligand naming rarely coincide directly).

Requires a reference built by clean_receptor_v2.py (first residue ACE,
last NME) - refuses to run otherwise. See NOTES.md for why.

Any reference atom with no match in a given frame (typically receptor
hydrogens pdb2gmx rebuilt fresh) keeps its reference-build coordinate -
check the printed match counts per frame before trusting the output.
"""
import argparse

# AMBER (pdb2gmx -ff amber99sb-ildn) atom name -> CHARMM equivalents, for the
# handful of protein atoms whose names differ between force fields. Found by
# actually diagnosing a clash: 304/337 unmatched atoms in a real test were
# this exact AMBER "H" vs CHARMM "HN" backbone amide hydrogen. A frozen,
# unmatched atom on a common name like this WILL clash with something once
# the surrounding sidechains move even slightly (VDW repulsion ~1/r^12).
AMBER_TO_CHARMM_ALIASES = {
    "H": ["HN"],                 # backbone amide H
    "H1": ["HT1"], "H2": ["HT2"], "H3": ["HT3"],   # N-terminal NH3+
    "O": ["OT1"], "OXT": ["OT2"], "OC1": ["OT1"], "OC2": ["OT2"],  # C-terminal COO-
    "HG": ["HG1"],                # SER/CYS single sidechain H
    "HE2": ["HE21"],               # fallback, residue-dependent
}

# Raw CHARMM cap atom name -> AMBER equivalent (amber99sb-ildn's ACE/NME rtp
# entries), renamed directly during the frame-side split below so no alias
# lookup is needed for these - see split_frame_terminal_caps().
NTERM_CAP_RENAME = {"CAY": "CH3", "CY": "C", "OY": "O",
                     "HY1": "HH31", "HY2": "HH32", "HY3": "HH33"}
CTERM_CAP_RENAME = {"NT": "N", "CAT": "CH3", "HNT": "H",
                     "HT1": "HH31", "HT2": "HH32", "HT3": "HH33"}


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
                "icode": line[26:27].strip(),
                "x": line[30:38], "y": line[38:46], "z": line[46:54],
            })
    return atoms


def assign_residue_ordinals(atoms, ligand_resname):
    """Ordinal index of each non-ligand atom's residue, in file order (None for ligand atoms)."""
    ordinals = []
    current_key = None
    idx = -1
    for a in atoms:
        if a["resname"].upper() == ligand_resname.upper():
            ordinals.append(None)
            continue
        key = (a["chain"], a["resseq"], a["icode"])
        if key != current_key:
            idx += 1
            current_key = key
        ordinals.append(idx)
    return ordinals


def _residue_blocks(atoms):
    """Group atoms into contiguous (chain, resseq, icode) residue blocks, in file order."""
    blocks, keys = [], []
    i = 0
    while i < len(atoms):
        key = (atoms[i]["chain"], atoms[i]["resseq"], atoms[i]["icode"])
        block = []
        while i < len(atoms) and (atoms[i]["chain"], atoms[i]["resseq"], atoms[i]["icode"]) == key:
            block.append(atoms[i])
            i += 1
        blocks.append(block)
        keys.append(key)
    return blocks, keys


def split_frame_terminal_caps(atoms, ligand_resname):
    """
    Mirror clean_receptor_v2.py's ACE/NME cap split on a raw MD frame, so
    frame residue ordinals line up with a reference built by that script.
    Renames the CHARMM cap atom names to their AMBER equivalents in place.
    Ligand atoms are skipped - they never count as a protein residue here.

    Handles any number of independently-capped fragments per chain, matching
    clean_receptor_v2.py's split_terminal_caps() - see NOTES.md.
    """
    blocks, _ = _residue_blocks(atoms)
    protein_idx = {i for i, b in enumerate(blocks) if b[0]["resname"].upper() != ligand_resname.upper()}
    if not protein_idx:
        return atoms

    nterm_atoms_per_block = {i: [a for a in blocks[i] if a["name"].strip().upper() in NTERM_CAP_RENAME]
                              for i in protein_idx}
    cterm_atoms_per_block = {i: [a for a in blocks[i] if a["name"].strip().upper() in CTERM_CAP_RENAME]
                              for i in protein_idx}

    out = []
    expect_fragment_start = {}  # chain -> bool
    for i, block in enumerate(blocks):
        if i not in protein_idx:
            out.extend(block)
            continue

        chain = block[0]["chain"]
        nterm = nterm_atoms_per_block[i]
        cterm = cterm_atoms_per_block[i]
        is_fragment_start = expect_fragment_start.get(chain, True)  # unseen chain -> starts a fragment

        if nterm and not is_fragment_start:
            raise ValueError(
                f"Frame has ACE-cap-style atom names on chain '{chain}' resseq {block[0]['resseq']}, "
                f"but the preceding residue of this chain didn't close its own fragment with an "
                f"NME cap first - unexpected, investigate."
            )

        if not nterm and not cterm:
            out.extend(block)
            expect_fragment_start[chain] = False
            continue

        orig_resseq = block[0]["resseq"]  # cache before any mutation - block[0] IS one of the cap
                                           # atoms being renamed below, so reading it fresh inside
                                           # the loop after the first iteration would read its own
                                           # already-mutated value and compound on every atom after it
        rest = [a for a in block if a not in nterm and a not in cterm]
        if nterm:
            for a in nterm:
                a["name"] = NTERM_CAP_RENAME[a["name"].strip().upper()]
                a["resname"] = "ACE"
                a["resseq"] = f"{orig_resseq}~ACE"
            out.extend(nterm)
        out.extend(rest)
        if cterm:
            for a in cterm:
                a["name"] = CTERM_CAP_RENAME[a["name"].strip().upper()]
                a["resname"] = "NME"
                a["resseq"] = f"{orig_resseq}~NME"
            out.extend(cterm)
            expect_fragment_start[chain] = True  # this fragment just closed; next block starts a new one
        else:
            expect_fragment_start[chain] = False

    return out


def remap_frame(reference_atoms, ref_ordinals, frame_atoms, ligand_resname, ref_ligand_resname,
                 ligand_atom_map=None):
    frame_ordinals = assign_residue_ordinals(frame_atoms, ligand_resname)

    protein_lookup = {}
    ligand_lookup = {}
    for a, ordinal in zip(frame_atoms, frame_ordinals):
        if ordinal is None:
            ligand_lookup[a["name"].upper()] = a
        else:
            protein_lookup[(ordinal, a["name"].upper())] = a

    out_lines = []
    counts = {"protein_matched": 0, "protein_kept": 0, "ligand_matched": 0, "ligand_kept": 0}
    for ref, ordinal in zip(reference_atoms, ref_ordinals):
        is_ligand = ref["resname"].upper() == ref_ligand_resname.upper()
        if is_ligand:
            # ligand_atom_map translates the reference's generic acpype name
            # to the raw frame's real atom name - see NOTES.md.
            name = ref["name"].upper()
            if ligand_atom_map:
                name = ligand_atom_map.get(name, name)
            src = ligand_lookup.get(name)
        else:
            name = ref["name"].upper()
            src = protein_lookup.get((ordinal, name))
            if src is None:
                for alias in AMBER_TO_CHARMM_ALIASES.get(name, []):
                    src = protein_lookup.get((ordinal, alias))
                    if src:
                        break
        line = ref["line"]
        bucket = "ligand" if is_ligand else "protein"
        if src:
            line = line[:30] + src["x"] + src["y"] + src["z"] + line[54:]
            counts[f"{bucket}_matched"] += 1
        else:
            counts[f"{bucket}_kept"] += 1
        out_lines.append(line)
    return out_lines, counts


def main():
    parser = argparse.ArgumentParser(description="Build a multi-frame trajectory PDB matched to a reference topology.")
    parser.add_argument("-ref", required=True, help="Reference complex.pdb from build_topology_once.py")
    parser.add_argument("-frames", nargs="+", required=True, help="Raw per-frame complex PDBs, in order")
    parser.add_argument("-ligand-resname", required=True, dest="ligand_resname",
                         help="Ligand residue name as it appears in the RAW frame files")
    parser.add_argument("-ref-ligand-resname", default="MOL", dest="ref_ligand_resname",
                         help="Ligand residue name in the built reference topology "
                              "(unigbsa hardcodes this to MOL in build_lignad(), regardless of your input file's resname)")
    parser.add_argument("-ligand-atom-map", dest="ligand_atom_map", default=None,
                         help="Path to the <ligand>.atom_map sidecar written by fix_ligand_naming.py.")
    parser.add_argument("-o", dest="outfile", default="trajectory.pdb")
    args = parser.parse_args()

    ligand_atom_map = None
    if args.ligand_atom_map:
        ligand_atom_map = {}
        with open(args.ligand_atom_map) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    ligand_atom_map[parts[0].upper()] = parts[1].upper()
        print(f"Loaded {len(ligand_atom_map)} ligand atom name mappings from {args.ligand_atom_map}")

    reference_atoms = parse_pdb_atoms(args.ref)
    ref_ordinals = assign_residue_ordinals(reference_atoms, args.ref_ligand_resname)
    n_protein = len([o for o in ref_ordinals if o is not None])
    n_residues = (max([o for o in ref_ordinals if o is not None], default=-1) + 1)
    n_ligand_ref = len(reference_atoms) - n_protein
    print(f"Reference has {len(reference_atoms)} atoms (this defines the fixed atom order), "
          f"{n_protein} in {n_residues} protein residues, {n_ligand_ref} ligand atoms.")

    ref_blocks, _ = _residue_blocks(reference_atoms)
    ref_protein_blocks = [b for b in ref_blocks if b[0]["resname"].upper() != args.ref_ligand_resname.upper()]
    first_resname = ref_protein_blocks[0][0]["resname"].upper()
    last_resname = ref_protein_blocks[-1][0]["resname"].upper()
    if first_resname != "ACE" or last_resname != "NME":
        raise SystemExit(
            f"This script requires a reference topology built by clean_receptor_v2.py, whose first "
            f"residue must be ACE and last must be NME (got '{first_resname}' / '{last_resname}'). "
            f"An unsplit reference will silently misalign every residue - refusing to run."
        )

    any_ligand_mismatch = False
    with open(args.outfile, "w") as out:
        for i, framefile in enumerate(args.frames, 1):
            frame_atoms = parse_pdb_atoms(framefile)
            frame_atoms = split_frame_terminal_caps(frame_atoms, args.ligand_resname)
            lines, counts = remap_frame(
                reference_atoms, ref_ordinals, frame_atoms, args.ligand_resname, args.ref_ligand_resname,
                ligand_atom_map=ligand_atom_map,
            )
            print(f"Frame {i} ({framefile}): "
                  f"protein {counts['protein_matched']} matched / {counts['protein_kept']} kept-from-reference "
                  f"(expect kept ~= receptor H count only); "
                  f"ligand {counts['ligand_matched']} matched / {counts['ligand_kept']} kept-from-reference "
                  f"(expect 0 kept - there's exactly one copy, matched by name only)")
            if counts["ligand_kept"] > 0:
                any_ligand_mismatch = True
                print(f"  WARNING: {counts['ligand_kept']} ligand atom(s) did not match by name in this frame - "
                      f"they are frozen at the reference build's coordinate. Pass -ligand-atom-map "
                      f"(see NOTES.md) if you haven't already.")
            out.write(f"MODEL{i:>9d}\n")
            for line in lines:
                out.write(line + "\n")
            out.write("ENDMDL\n")
        out.write("END\n")

    print(f"Wrote {len(args.frames)}-frame trajectory: {args.outfile}")
    if any_ligand_mismatch:
        print("WARNING: at least one frame had unmatched ligand atoms - do not trust the ligand's "
              "contribution to TOTAL until this is resolved.")


if __name__ == "__main__":
    main()
