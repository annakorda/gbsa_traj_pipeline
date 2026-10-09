#!/usr/bin/env python3
"""
Standalone receptor-only PDB cleaner for GBSA/GROMACS pdb2gmx input.
Takes a raw protein PDB and produces a GROMACS/amber-ff-compatible
protein-only PDB: translates CHARMM histidine tautomers to their AMBER
names, renames disulfide-bonded CYS to CYX, and splits N-/C-terminal
ACE/NME caps into standalone residues (any number of fragments per chain).

Must always be paired with build_trajectory_v2.py, which applies the
identical residue/cap handling to each raw frame - see NOTES.md for why.
"""
import re
import sys
import argparse

HIS_VARIANTS = {"HSD", "HSE", "HSP", "HID", "HIE", "HIP"}
HIS_TO_AMBER = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}  # already-AMBER names pass through unchanged
NTERM_CAP_ATOMS = {"CAY": "CH3", "CY": "C", "OY": "O"}   # split into a standalone ACE residue
CTERM_CAP_ATOMS = {"NT": "N", "CAT": "CH3"}              # split into a standalone NME residue

AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "CYX", "CYM", "ASH", "GLH", "LYN", "HSD", "HSE", "HSP", "HID", "HIE", "HIP",
}

DEFAULT_DROP_RESNAMES = {
    # water
    "HOH", "WAT", "H2O", "TIP3", "TIP4", "SPC",
    # ions
    "NA", "CL", "K", "CS", "LI", "MG", "CA", "ZN", "MN", "FE", "FE2", "CO",
    "CU", "NI", "CD", "BR", "IOD", "SO4", "PO4",
    # common crystallization additives
    "GOL", "EDO", "PEG", "PG4", "BME", "DMS", "ACT", "FMT", "MPD",
}

RESNAME_ELEMENT = {
    "NA": "NA", "CL": "CL", "K": "K", "CS": "CS", "LI": "LI", "MG": "MG",
    "CA": "CA", "ZN": "ZN", "MN": "MN", "FE": "FE", "FE2": "FE", "CO": "CO",
    "CU": "CU", "NI": "NI", "CD": "CD", "BR": "BR", "IOD": "I",
}


def parse_atom(line):
    rec = line[0:6].strip()
    if rec not in ("ATOM", "HETATM"):
        return None
    d = {
        "record": line[0:6],
        "serial": line[6:11],
        "name": line[12:16],
        "altloc": line[16:17],
        "resname": line[17:20],
        "chain": line[21:22],
        "resseq": line[22:26],
        "icode": line[26:27],
        "x": line[30:38],
        "y": line[38:46],
        "z": line[46:54],
        "occ": line[54:60],
        "temp": line[60:66],
        "element": line[76:78] if len(line) >= 78 else "  ",
    }
    try:
        d["_xf"] = float(d["x"])
        d["_yf"] = float(d["y"])
        d["_zf"] = float(d["z"])
    except ValueError:
        d["_xf"] = d["_yf"] = d["_zf"] = None
    return d


def is_hydrogen(d):
    el = d["element"].strip().upper()
    nm = d["name"].strip().upper()
    return el == "H" or bool(re.match(r"^\d*H", nm))


def guess_element(atom_name, resname):
    rn = resname.strip().upper()
    if rn in RESNAME_ELEMENT:
        return RESNAME_ELEMENT[rn]
    nm = atom_name.strip().upper()
    if nm == "CA" and rn in AA_NAMES:
        return "C"
    if nm.startswith("O"):
        return "O"
    if nm.startswith("N"):
        return "N"
    if nm.startswith("S"):
        return "S"
    for two in ("CL", "BR", "ZN", "FE", "MG", "MN", "CO", "CU", "NI", "SE", "NA", "K"):
        if nm.startswith(two):
            return two
    return "C"


def fmt_atom(d, new_serial):
    an = d["name"].strip()
    if len(an) == 4:
        atom_field = an
    elif len(an) == 3:
        atom_field = f" {an}"
    elif len(an) == 2:
        atom_field = f" {an} "
    else:
        atom_field = f" {an}  "
    occ = d["occ"] if d["occ"].strip() else "  1.00"
    temp = d["temp"] if d["temp"].strip() else " 20.00"
    element = d["element"].rjust(2)[:2]
    return (f"{d['record'][:6]:<6}{new_serial:5d} {atom_field} "
            f"{d['resname'].strip():>3s} {d['chain']}{d['resseq']}{d['icode']}   "
            f"{d['x']:>8}{d['y']:>8}{d['z']:>8}{occ:>6}{temp:>6}          {element}\n")


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


def _unique_adjacent_key(chain, base_resseq, used_keys):
    """A (chain, resseq, icode) key that doesn't collide with any existing
    residue, by reusing the adjacent real residue's resSeq with a fresh
    insertion code (see NOTES.md for why plain resseq+-1 isn't safe)."""
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        key = (chain, base_resseq, letter)
        if key not in used_keys:
            return key
    raise ValueError(f"Could not find a free insertion code near resseq {base_resseq!r} on chain '{chain}'")


def split_terminal_caps(atoms):
    """
    Split CHARMM N-/C-terminal cap atoms fused onto a residue into standalone
    ACE/NME residues (amber naming), matching amber99sb-ildn's own ACE/NME
    rtp entries. Must be paired with build_trajectory_v2.py, which applies
    the identical split to each raw frame.

    Handles any number of independently-capped fragments per chain (e.g. an
    ICL3-truncated construct), not just one ACE/NME pair at the chain's
    first/last residue - see NOTES.md.
    """
    blocks, keys = _residue_blocks(atoms)
    used_keys = set(keys)

    nterm_atoms_per_block = [[a for a in b if a["name"].strip().upper() in NTERM_CAP_ATOMS] for b in blocks]
    cterm_atoms_per_block = [[a for a in b if a["name"].strip().upper() in CTERM_CAP_ATOMS] for b in blocks]

    out = []
    n_ace, n_nme, n_fragments = 0, 0, 0
    expect_fragment_start = {}  # chain -> bool
    frag_id_of_chain = {}       # chain -> current fragment counter, for TER placement
    for idx, block in enumerate(blocks):
        chain = block[0]["chain"]
        nterm_atoms = nterm_atoms_per_block[idx]
        cterm_atoms = cterm_atoms_per_block[idx]
        is_fragment_start = expect_fragment_start.get(chain, True)  # unseen chain -> starts a fragment

        if nterm_atoms and not is_fragment_start:
            raise ValueError(
                f"ACE-cap atoms found on chain '{chain}' resseq {block[0]['resseq']}, but the "
                f"preceding residue of this chain didn't close its own fragment with an NME cap "
                f"first - refusing to split, investigate first (two open N-termini back to back "
                f"isn't a real construct)."
            )

        if is_fragment_start:
            frag_id_of_chain[chain] = frag_id_of_chain.get(chain, 0) + 1
            n_fragments += 1
        frag_id = frag_id_of_chain[chain]

        rest = [a for a in block if a not in nterm_atoms and a not in cterm_atoms]

        if nterm_atoms:
            ace_key = _unique_adjacent_key(chain, block[0]["resseq"], used_keys)
            used_keys.add(ace_key)
            for a in nterm_atoms:
                a["name"] = NTERM_CAP_ATOMS[a["name"].strip().upper()]
                a["resname"] = "ACE"
                a["resseq"], a["icode"] = ace_key[1], ace_key[2]
                a["_frag"] = (chain, frag_id)
            out.extend(nterm_atoms)
            n_ace += 1

        for a in rest:
            a["_frag"] = (chain, frag_id)
        out.extend(rest)

        if cterm_atoms:
            nme_key = _unique_adjacent_key(chain, block[0]["resseq"], used_keys)
            used_keys.add(nme_key)
            for a in cterm_atoms:
                a["name"] = CTERM_CAP_ATOMS[a["name"].strip().upper()]
                a["resname"] = "NME"
                a["resseq"], a["icode"] = nme_key[1], nme_key[2]
                a["_frag"] = (chain, frag_id)
            out.extend(cterm_atoms)
            n_nme += 1
            expect_fragment_start[chain] = True  # this fragment just closed; next block starts a new one
        else:
            expect_fragment_start[chain] = False

    return out, n_ace, n_nme, n_fragments


def clean_receptor_pdb(input_pdb, output_pdb, chains=None, extra_drop_resnames=None, verbose=True):
    """
    Clean a receptor-only PDB for GROMACS/amber pdb2gmx.

    chains: optional set of chain IDs to keep (default: all chains present).
    extra_drop_resnames: additional HETATM resnames to drop beyond DEFAULT_DROP_RESNAMES.
    """
    drop_resnames = set(DEFAULT_DROP_RESNAMES)
    if extra_drop_resnames:
        drop_resnames |= {r.upper() for r in extra_drop_resnames}

    atoms = []
    seen_res_altloc = {}
    n_dropped_hetatm = 0
    n_dropped_chain = 0
    n_dropped_altloc = 0
    n_dropped_h = 0
    his_translated_residues = set()
    dropped_resnames = set()
    lys_hz3 = set()  # LYS residues with a 3rd NZ hydrogen -> protonated, keep as LYS

    with open(input_pdb) as f:
        for line in f:
            d = parse_atom(line)
            if d is None:
                continue

            chain = d["chain"].strip() or "_"
            if chains and chain not in chains:
                n_dropped_chain += 1
                continue

            rn = d["resname"].strip().upper()
            if rn in HIS_VARIANTS:
                new_rn = HIS_TO_AMBER.get(rn, rn)  # HID/HIE/HIP already-AMBER names pass through
                if new_rn != rn:
                    his_translated_residues.add((chain, d["resseq"], d["icode"]))
                d["resname"] = new_rn
                rn = new_rn

            if rn in drop_resnames:
                n_dropped_hetatm += 1
                dropped_resnames.add(rn)
                continue

            nm = d["name"].strip().upper()
            if rn == "LYS" and nm == "HZ3":
                lys_hz3.add((chain, d["resseq"], d["icode"]))

            if is_hydrogen(d):
                n_dropped_h += 1
                continue

            altloc = d["altloc"].strip()
            if altloc and altloc != "A":
                n_dropped_altloc += 1
                continue

            res_key = (chain, d["resseq"], d["icode"], nm)
            if res_key in seen_res_altloc:
                continue
            seen_res_altloc[res_key] = True

            d["altloc"] = " "
            d["element"] = guess_element(d["name"], d["resname"]).rjust(2)
            atoms.append(d)

    if not atoms:
        raise ValueError(f"No atoms survived cleaning for {input_pdb} (check chain/resname filters)")

    atoms, n_ace, n_nme, n_fragments = split_terminal_caps(atoms)

    # LYS with no 3rd NZ hydrogen in the input -> neutral, rename to LYN
    # (amber99sb-ildn defines LYN directly; see NOTES.md).
    n_lyn = 0
    lys_keys = {(a["chain"], a["resseq"], a["icode"]) for a in atoms
                if a["resname"].strip().upper() == "LYS"}
    for key in lys_keys - lys_hz3:
        for a in atoms:
            if (a["chain"], a["resseq"], a["icode"]) == key and a["resname"].strip().upper() == "LYS":
                a["resname"] = "LYN"
        n_lyn += 1

    sg = []
    for i, d in enumerate(atoms):
        if d["resname"].strip().upper() == "CYS" and d["name"].strip().upper() == "SG" \
                and d["_xf"] is not None:
            try:
                resno = int(d["resseq"])
            except ValueError:
                resno = int(re.sub(r"\D", "", d["resseq"]) or 0)
            sg.append((i, d["chain"].strip() or "_", resno, d["_xf"], d["_yf"], d["_zf"]))

    conect_idx = set()
    disulfide_residues = set()
    thr2 = 2.2 * 2.2
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            a, b = sg[i], sg[j]
            d2 = (a[3] - b[3]) ** 2 + (a[4] - b[4]) ** 2 + (a[5] - b[5]) ** 2
            if d2 <= thr2:
                conect_idx.add(tuple(sorted((a[0], b[0]))))
                disulfide_residues.add(a[0])
                disulfide_residues.add(b[0])

    # v2: actually rename the bonded residues CYS -> CYX so pdb2gmx builds
    # the disulfide deterministically, instead of relying on its own
    # distance-based auto-detect + interactive prompt.
    n_cyx = 0
    for sg_idx in disulfide_residues:
        chain_, resseq_, icode_ = atoms[sg_idx]["chain"], atoms[sg_idx]["resseq"], atoms[sg_idx]["icode"]
        for a in atoms:
            if (a["chain"], a["resseq"], a["icode"]) == (chain_, resseq_, icode_) and \
                    a["resname"].strip().upper() == "CYS":
                a["resname"] = "CYX"
        n_cyx += 1

    chains_kept = sorted({d["chain"].strip() or "_" for d in atoms})
    serial_of = {}
    with open(output_pdb, "w") as out:
        serial = 1
        prev_frag = None
        for i, d in enumerate(atoms):
            # TER on chain change and on fragment change within a chain (see NOTES.md).
            cur_frag = d.get("_frag", (d["chain"].strip() or "_", 0))
            if prev_frag is not None and cur_frag != prev_frag:
                out.write("TER\n")
            serial_of[i] = serial
            out.write(fmt_atom(d, serial))
            serial += 1
            prev_frag = cur_frag
        out.write("TER\n")
        for i, j in sorted(conect_idx):
            a, b = serial_of[i], serial_of[j]
            out.write(f"CONECT{a:5d}{b:5d}\n")
            out.write(f"CONECT{b:5d}{a:5d}\n")
        out.write("END\n")

    if verbose:
        print(f"Cleaned PDB: {len(atoms)} atoms, {len(chains_kept)} chain(s) kept: {','.join(chains_kept)}")
        print(f"  Dropped: {n_dropped_h} H, {n_dropped_hetatm} HETATM atoms "
              f"({','.join(sorted(dropped_resnames)) or 'none'}), "
              f"{n_dropped_altloc} alt-conformers, {n_dropped_chain} wrong-chain atoms")
        print(f"  Histidine tautomers translated to AMBER naming (HSD/HSE/HSP -> HID/HIE/HIP): "
              f"{len(his_translated_residues)} residue(s)")
        print(f"  Neutral lysines (no HZ3) renamed LYS->LYN: {n_lyn} residue(s)")
        print(f"  Disulfide bonds detected (SG-SG <= 2.2 A) and renamed CYS->CYX: {n_cyx} residue(s), "
              f"{len(conect_idx)} bond(s)")
        print(f"  Terminal caps split into standalone residues: {n_ace} ACE, {n_nme} NME "
              f"(must be matched by build_trajectory_v2.py)")
        print(f"  Polypeptide fragments detected: {n_fragments} "
              f"({'single continuous chain' if n_fragments == len(chains_kept) else 'chain fused from multiple segments'})")


def main():
    parser = argparse.ArgumentParser(description="Clean a receptor-only PDB for GROMACS/amber pdb2gmx (v2).")
    parser.add_argument("-i", dest="input_pdb", required=True, help="Input receptor PDB.")
    parser.add_argument("-o", dest="output_pdb", required=True, help="Output cleaned PDB.")
    parser.add_argument("--chains", default=None,
                         help="Comma-separated chain IDs to keep (default: keep all chains present).")
    parser.add_argument("--keep-resname", dest="keep_resnames", default=None,
                         help="Comma-separated HETATM resnames to keep even though they're in the default drop list "
                              "(e.g. a bound Na+ or Zn2+ you want to retain).")
    parser.add_argument("--drop-resname", dest="extra_drop", default=None,
                         help="Comma-separated extra HETATM resnames to drop beyond the default list.")
    args = parser.parse_args()

    global DEFAULT_DROP_RESNAMES
    if args.keep_resnames:
        keep = {r.strip().upper() for r in args.keep_resnames.split(",")}
        DEFAULT_DROP_RESNAMES = DEFAULT_DROP_RESNAMES - keep

    chains = {c.strip() for c in args.chains.split(",")} if args.chains else None
    extra_drop = args.extra_drop.split(",") if args.extra_drop else None

    clean_receptor_pdb(args.input_pdb, args.output_pdb, chains=chains, extra_drop_resnames=extra_drop)


if __name__ == "__main__":
    main()
