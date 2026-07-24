#!/usr/bin/env python3
"""
Standalone receptor-only PDB cleaner for GBSA / GROMACS pdb2gmx input.
Takes a raw protein PDB (any source: crystal structure, MD frame, etc.)
and produces a GROMACS/amber-ff-compatible protein-only PDB.

v2 - changes vs clean_receptor.py, all confirmed against the actual
amber99sb-ildn.ff files installed in this env (aminoacids.rtp,
specbond.dat) before writing, not assumed:

1. Histidine tautomer is now PRESERVED, not discarded. v1 collapsed
   HSD/HSE/HSP (and HID/HIE/HIP) all to generic "HIS", after which
   pdb2gmx (called downstream with no -his flag) silently applies one
   uniform default tautomer to every histidine - wrong for whichever
   ones weren't actually that state, causing their sidechain HD1/HE2
   to permanently mismatch (see DEBUG_FROZEN_ATOM_CLASHES.md, cause 1).
   v2 translates CHARMM names to their literal AMBER equivalents
   instead: HSD->HID, HSE->HIE, HSP->HIP. Confirmed amber99sb-ildn.ff's
   aminoacids.rtp defines HID/HIE/HIP directly as residues, so pdb2gmx
   builds the correct sidechain per-residue automatically, no -his
   flag or interactive prompt needed.

2. Disulfide bonds are now actually communicated to pdb2gmx. v1
   detected SG-SG distances <=2.2 A and wrote CONECT records, but
   pdb2gmx does not read CONECT for bonding, and the residues were
   still left named plain CYS - relying entirely on pdb2gmx's own
   specbond.dat distance-based auto-detection, which (confirmed by
   reading specbond.dat) triggers an interactive per-pair prompt that
   the topology-build step's single piped "1" cannot reliably answer
   for every pair. v2 renames both residues of each detected disulfide
   pair from CYS to CYX directly. Confirmed amber99sb-ildn.ff's
   aminoacids.rtp defines CYX directly, so pdb2gmx builds the bonded
   topology deterministically with no auto-detection or prompt at all.

3. N-terminal ACE cap and C-terminal NME cap are both split into their
   own standalone residues (CAY/CY/OY -> CH3/C/O as ACE; NT/CAT -> CH3/N
   as NME - matching amber99sb-ildn's own ACE/NME entries). Two earlier
   attempts at this (splitting only here, on the topology-building side)
   were tried and reverted after actually testing them: splitting here
   alone changes the reference's residue COUNT relative to the raw
   frame (which still has the caps fused into their neighbor residue),
   which breaks build_trajectory.py's pure ordinal-based matching -
   either by shifting every downstream residue's ordinal (the ACE case:
   confirmed reference ordinal 1 = MET's real atoms silently getting
   matched against frame ordinal 1 = PRO's real atoms - undetectable,
   far worse than not fixing it), or, appended at the end (the NME
   case), by creating a reference ordinal with nothing to match at all
   (the entire NME residue froze instead of just 2 atoms).
   The fix that actually works, verified end-to-end (topology build +
   real matching against all 3 test frames, checked the same ordinal
   alignment that caught the earlier regressions): split the SAME way
   on BOTH sides. This file does the reference/topology-building half.
   build_trajectory_v2.py does the matching half - it applies the
   identical split to each raw frame before computing ordinals, so
   reference and frame residue counts stay identical throughout, and
   the AMBER<->CHARMM name aliases needed for the caps' own atoms
   (CH3/C/O/N vs CAY/CY/OY/NT/CAT) live there. Do not use either script
   without the other - a topology built here now has 2 more residues
   than v1, and only the paired build_trajectory_v2.py knows to expect
   that.
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


def split_terminal_caps(atoms):
    """
    Split CHARMM N-/C-terminal cap atoms fused onto the first/last residue
    of each chain into standalone ACE/NME residues (renamed to amber
    naming), matching amber99sb-ildn's own ACE/NME rtp entries. Must be
    paired with build_trajectory_v2.py, which applies the identical split
    to each raw frame - see this file's module docstring for why.
    """
    blocks, keys = _residue_blocks(atoms)
    first_idx_of_chain, last_idx_of_chain = {}, {}
    for idx, key in enumerate(keys):
        chain = key[0]
        if chain not in first_idx_of_chain:
            first_idx_of_chain[chain] = idx
        last_idx_of_chain[chain] = idx

    out = []
    n_ace, n_nme = 0, 0
    for idx, block in enumerate(blocks):
        chain = block[0]["chain"]
        nterm_atoms = [a for a in block if a["name"].strip().upper() in NTERM_CAP_ATOMS]
        cterm_atoms = [a for a in block if a["name"].strip().upper() in CTERM_CAP_ATOMS]

        if nterm_atoms and idx != first_idx_of_chain[chain]:
            raise ValueError(f"ACE-cap atoms found on a non-first residue of chain '{chain}' "
                              f"(resseq {block[0]['resseq']}) - refusing to split, investigate first.")
        if cterm_atoms and idx != last_idx_of_chain[chain]:
            raise ValueError(f"NME-cap atoms found on a non-last residue of chain '{chain}' "
                              f"(resseq {block[0]['resseq']}) - refusing to split, investigate first.")

        rest = [a for a in block if a not in nterm_atoms and a not in cterm_atoms]

        if nterm_atoms:
            try:
                ace_resseq = f"{int(block[0]['resseq']) - 1:>4d}"
            except ValueError:
                ace_resseq = block[0]["resseq"]
            for a in nterm_atoms:
                a["name"] = NTERM_CAP_ATOMS[a["name"].strip().upper()]
                a["resname"] = "ACE"
                a["resseq"] = ace_resseq
            out.extend(nterm_atoms)
            n_ace += 1

        out.extend(rest)

        if cterm_atoms:
            try:
                nme_resseq = f"{int(block[0]['resseq']) + 1:>4d}"
            except ValueError:
                nme_resseq = block[0]["resseq"]
            for a in cterm_atoms:
                a["name"] = CTERM_CAP_ATOMS[a["name"].strip().upper()]
                a["resname"] = "NME"
                a["resseq"] = nme_resseq
            out.extend(cterm_atoms)
            n_nme += 1

    return out, n_ace, n_nme


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

    atoms, n_ace, n_nme = split_terminal_caps(atoms)

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
        prev_chain = None
        for i, d in enumerate(atoms):
            cur_chain = d["chain"].strip() or "_"
            if prev_chain is not None and cur_chain != prev_chain:
                out.write("TER\n")
            serial_of[i] = serial
            out.write(fmt_atom(d, serial))
            serial += 1
            prev_chain = cur_chain
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
        print(f"  Disulfide bonds detected (SG-SG <= 2.2 A) and renamed CYS->CYX: {n_cyx} residue(s), "
              f"{len(conect_idx)} bond(s)")
        print(f"  Terminal caps split into standalone residues: {n_ace} ACE, {n_nme} NME "
              f"(must be matched by build_trajectory_v2.py, not v1's build_trajectory.py)")


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
