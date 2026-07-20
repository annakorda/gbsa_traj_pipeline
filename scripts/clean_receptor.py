#!/usr/bin/env python3
"""
Standalone receptor-only PDB cleaner for GBSA / GROMACS pdb2gmx input.
Takes a raw protein PDB (any source: crystal structure, MD frame, etc.)
and produces a GROMACS/amber-ff-compatible protein-only PDB.
"""
import re
import sys
import argparse

HIS_VARIANTS = {"HSD", "HSE", "HSP", "HID", "HIE", "HIP"}
DROP_WEIRD_ATOMS = {"CAY", "CY", "OY", "CAT", "NT"}

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
                d["resname"] = "HIS"
                rn = "HIS"

            if rn in drop_resnames:
                n_dropped_hetatm += 1
                dropped_resnames.add(rn)
                continue

            nm = d["name"].strip().upper()
            if nm in DROP_WEIRD_ATOMS:
                continue

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
    thr2 = 2.2 * 2.2
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            a, b = sg[i], sg[j]
            d2 = (a[3] - b[3]) ** 2 + (a[4] - b[4]) ** 2 + (a[5] - b[5]) ** 2
            if d2 <= thr2:
                conect_idx.add(tuple(sorted((a[0], b[0]))))

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
        print(f"  Disulfide bonds detected (SG-SG <= 2.2 A): {len(conect_idx)}")


def main():
    parser = argparse.ArgumentParser(description="Clean a receptor-only PDB for GROMACS/amber pdb2gmx.")
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
