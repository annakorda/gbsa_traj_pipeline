#!/usr/bin/env python3
"""
Build the receptor+ligand topology ONCE (no minimization, no box/solvent).
Uses unigbsa's own topology builder directly, matching the working
code path used internally by unigbsa-pipeline's mode=input / mode=em.

Input: a cleaned receptor-only PDB + a ligand mol2/mol/sdf file.
Output: complex.top, complex.pdb (reference structure, fixes the atom
order everything downstream must match), index.ndx.

If a "<ligand>.net_charge" sidecar exists (written by fix_ligand_naming.py),
forces that value through to acpype/antechamber instead of trusting
unigbsa's own net-charge auto-detection - see NOTES.md.
"""
import os
import shutil
import argparse

import unigbsa.simulation.topology as _unigbsa_topology
from unigbsa.simulation.topology import build_topol
from unigbsa.utils import generate_index_file


def main():
    parser = argparse.ArgumentParser(description="Build receptor+ligand topology once, no minimization.")
    parser.add_argument("-r", dest="receptor", required=True, help="Cleaned receptor-only PDB.")
    parser.add_argument("-l", dest="ligand", required=True, help="Ligand mol2/mol/sdf file.")
    parser.add_argument("-o", dest="outdir", default="topology", help="Output directory.")
    parser.add_argument("-pf", dest="proteinforce", default="amber99sb-ildn", help="Protein forcefield.")
    parser.add_argument("-lf", dest="ligandforce", default="gaff2", help="Ligand forcefield: gaff or gaff2.")
    parser.add_argument("-c", dest="charge_method", default="bcc", help="Ligand charge method (bcc or gas).")
    parser.add_argument("-nt", dest="threads", type=int, default=4, help="Threads for the ligand charge calc.")
    args = parser.parse_args()

    receptor = os.path.abspath(args.receptor)
    ligand = os.path.abspath(args.ligand)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    charge_file = ligand + ".net_charge"
    if os.path.exists(charge_file):
        with open(charge_file) as f:
            forced_charge = int(f.read().strip())
        print(f"Overriding unigbsa's net-charge auto-detection with {forced_charge:+d} from {charge_file}.")
        _unigbsa_topology.obtain_net_charge = lambda ligandfile, _c=forced_charge: _c

    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        indexfile = build_topol(
            receptor, ligand,
            outpdb="complex.pdb", outtop="complex.top",
            proteinforce=args.proteinforce, ligandforce=args.ligandforce,
            charge_method=args.charge_method, nt=args.threads,
        )
        if not os.path.exists(indexfile):
            indexfile = generate_index_file("complex.pdb")
        if os.path.abspath(indexfile) != os.path.join(outdir, "index.ndx"):
            shutil.copy2(indexfile, "index.ndx")
    finally:
        os.chdir(cwd)

    print("Built topology (this is the fixed atom order every frame must match):")
    print(f"  {outdir}/complex.top")
    print(f"  {outdir}/complex.pdb")
    print(f"  {outdir}/index.ndx")


if __name__ == "__main__":
    main()
