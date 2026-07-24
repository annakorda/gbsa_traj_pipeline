#!/usr/bin/env python3
"""
Remap N raw per-frame complex PDBs (protein+ligand, from whatever MD
engine) onto the fixed atom order of the reference topology built by
build_topology_once.py, and concatenate them into a multi-MODEL
trajectory PDB for unigbsa-traj.

v2 - changes vs build_trajectory.py:

1. Reporting/safety-check only, no matching-logic change: protein and
   ligand match/kept counts are now reported separately per frame
   instead of combined (v1 combined them, which could mask a
   ligand-specific matching problem inside an overall count that still
   "looks fine" because the much larger protein atom count dominates
   it), and a warning prints if any ligand atom fails to match (should
   always be 0 or all - there's exactly one copy, matched by name
   only). Ligand naming (fix_ligand_naming.py) itself is untouched.

2. Terminal-cap handling, paired with clean_receptor_v2.py: that script
   now splits the CHARMM N-terminal ACE cap and C-terminal NME cap into
   standalone residues in the REFERENCE topology (2 attempts at doing
   this on the topology side only were tried and reverted - see that
   file's docstring for why: it changes the reference's residue count
   relative to the raw frame, which breaks pure ordinal matching,
   either by shifting every downstream ordinal, or by creating a
   reference ordinal with nothing to match). The fix that actually
   works: split the SAME way on the FRAME side too, right here, before
   computing ordinals - so reference and frame residue counts match
   throughout, at every ordinal, not just count-for-count overall.
   split_frame_terminal_caps() below does that, renaming the raw
   CHARMM cap atom names directly to their AMBER equivalents so no new
   alias-table entries are needed.
   Verified end-to-end: built the topology, ran this against all 3
   test frames, and directly re-checked that reference ordinal 1 still
   matches frame ordinal 1's real MET atoms (not a neighboring
   residue's) - the exact check that caught the two reverted attempts.

   This version of the script REQUIRES a reference built by
   clean_receptor_v2.py (first residue literally named ACE, last
   literally named NME) - it will refuse to run against an
   unsplit/v1-style reference rather than silently misalign everything.

Matching rule (unchanged from v1):
  - protein atoms: matched by (residue ORDINAL POSITION, atom name).
    Not by literal chain/resnum - pdb2gmx drops the chain ID and
    renumbers residues (confirmed: the ligand occupies residue slot 1
    in the built topology, shifting every protein residue by +1), so
    literal (chain, resnum) never lines up. Ordinal position (1st
    residue encountered, 2nd, ...) survives that renumbering as long
    as residue order/count is unchanged, which it is here.
  - ligand atoms: matched by atom name only (there's one copy).
    Also note: unigbsa's build_lignad() hardcodes the ligand residue
    name to "MOL" regardless of what your input mol2 called it - that
    is why -ref-ligand-resname defaults to MOL, not your real ligand's
    resname.

Any reference atom with no match in a given frame (typically the
receptor hydrogens, which pdb2gmx rebuilt fresh and which never
existed in the raw MD frame) keeps its reference-build coordinate.
Check the printed match counts per frame before trusting the output -
a low match count on the protein or ligand means the matching
assumption above doesn't hold for your files.
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
    Renames the CHARMM cap atom names to their AMBER equivalents in place
    (no alias-table entries needed for them). Ligand atoms are skipped
    entirely - they never count as a protein residue here.
    """
    blocks, _ = _residue_blocks(atoms)
    protein_idx = [i for i, b in enumerate(blocks) if b[0]["resname"].upper() != ligand_resname.upper()]
    if not protein_idx:
        return atoms

    first_of_chain, last_of_chain = {}, {}
    for i in protein_idx:
        chain = blocks[i][0]["chain"]
        first_of_chain.setdefault(chain, i)
        last_of_chain[chain] = i

    out = []
    for i, block in enumerate(blocks):
        if i not in first_of_chain.values() and i not in last_of_chain.values():
            out.extend(block)
            continue
        chain = block[0]["chain"]
        nterm = [a for a in block if a["name"].strip().upper() in NTERM_CAP_RENAME]
        cterm = [a for a in block if a["name"].strip().upper() in CTERM_CAP_RENAME]
        if nterm and i != first_of_chain.get(chain):
            raise ValueError(f"Frame has ACE-cap-style atom names on a non-first residue of "
                              f"chain '{chain}' (resseq {block[0]['resseq']}) - unexpected, investigate.")
        if cterm and i != last_of_chain.get(chain):
            raise ValueError(f"Frame has NME-cap-style atom names on a non-last residue of "
                              f"chain '{chain}' (resseq {block[0]['resseq']}) - unexpected, investigate.")
        if not nterm and not cterm:
            out.extend(block)
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

    return out


def remap_frame(reference_atoms, ref_ordinals, frame_atoms, ligand_resname, ref_ligand_resname):
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
            src = ligand_lookup.get(ref["name"].upper())
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
    parser.add_argument("-o", dest="outfile", default="trajectory.pdb")
    args = parser.parse_args()

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
            f"Using this script against a v1-style (unsplit) reference will silently misalign every "
            f"residue - refusing to run. Use build_trajectory.py (v1) for a v1-style reference instead."
        )

    any_ligand_mismatch = False
    with open(args.outfile, "w") as out:
        for i, framefile in enumerate(args.frames, 1):
            frame_atoms = parse_pdb_atoms(framefile)
            frame_atoms = split_frame_terminal_caps(frame_atoms, args.ligand_resname)
            lines, counts = remap_frame(
                reference_atoms, ref_ordinals, frame_atoms, args.ligand_resname, args.ref_ligand_resname
            )
            print(f"Frame {i} ({framefile}): "
                  f"protein {counts['protein_matched']} matched / {counts['protein_kept']} kept-from-reference "
                  f"(expect kept ~= receptor H count only); "
                  f"ligand {counts['ligand_matched']} matched / {counts['ligand_kept']} kept-from-reference "
                  f"(expect 0 kept - there's exactly one copy, matched by name only)")
            if counts["ligand_kept"] > 0:
                any_ligand_mismatch = True
                print(f"  WARNING: {counts['ligand_kept']} ligand atom(s) did not match by name in this frame - "
                      f"they are frozen at the reference build's coordinate. Check that "
                      f"fix_ligand_naming.py's assigned names actually correspond to this frame's ligand atom "
                      f"names/order (see the caveat in that script's docstring).")
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
