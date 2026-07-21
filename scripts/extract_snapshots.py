#!/usr/bin/env python3
"""
Create the pdb of only protein and ligand of a GPCRmd complex. Introduce
the dyn ID and the number of frames you want to extract to PDB. By the 
default the last 600 frameswill be extracted.

Output created in data/GPCRmd/dynID/XXX_fms
"""

import json
import argparse
import os
import subprocess
import tempfile

def get_traj_path(dynID):
    """Get the traj paths, psf, ligand name and sdf of a GPCRmd complex from the precomputed json file."""
    dynData={}
    jsonPath="/GPCRmd_files/Precomputed/compl_info_new.json" # read json with all receptor data including gpcrnum correspondence
    with open(jsonPath, "r") as f:
        jsonDic = json.load(f)
    for key, value in jsonDic.items():
        if value['dyn_id']==dynID:
            dynData['trajs'] = ['/GPCRmd_files/Dynamics/'+i for i in value['traj_fnames']]
            dynData['psf'] = '/GPCRmd_files' + value['topo_f']
            dynData['lig_sname'] = value['lig_sname']
            sdfs_all=value['dyn_comp']
            sdfs_lig=[i for i in sdfs_all if i['resname']==value['lig_sname']] 
            dynData['sdf']= '/GPCRmd_files' + sdfs_lig[0]['sdfile']
    print(f"""- Found {dynData['trajs']} \nand {dynData['psf']} for dynID {dynID}. \n- Ligand name: {dynData['lig_sname']}, SDF: {dynData['sdf']}""")
    return dynData


def run_vmd(traj_paths, psf_path, output_dir, ligand_resname=None, fms=600):
    """Open VMD, load PSF with each trajectory separately, and write PDBs
    for the last `fms` frames of each trajectory."""

    if not os.path.isfile(psf_path):
        raise FileNotFoundError(f"PSF file not found: {psf_path}")

    os.makedirs(output_dir, exist_ok=True)

    vmd_lines = []

    for traj_path in traj_paths:
        if not os.path.isfile(traj_path):
            raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

        tr = traj_path.split('/')[-1].split('.')[0]

        # Fresh molecule per trajectory so frame numbering isn't shared/cumulative
        vmd_lines.append(f'mol new "{psf_path}" type psf waitfor all')
        vmd_lines.append(f'mol addfile "{traj_path}" type xtc waitfor all')

        if ligand_resname:
            selection = f'(protein and chain P) or resname {ligand_resname}'
        else:
            print("Ligand resname not provided, extracting only protein and excluding water and ions.")
            selection = '(protein and chain P) or (not water and not ion and not protein and not lipid)'
        vmd_lines.append(f'set sel [atomselect top "{selection}"]')

        # Compute last `fms` frames of THIS trajectory
        vmd_lines.append('set nf [molinfo top get numframes]')
        vmd_lines.append(f'set start [expr {{$nf - {fms} < 0 ? 0 : $nf - {fms}}}]')
        vmd_lines.append('for {set i $start} {$i < $nf} {incr i} {')
        vmd_lines.append('    animate goto $i')
        vmd_lines.append(f'    set outfile [format "{output_dir}/{tr}_%d.pdb" [expr {{$i - $start + 1}}]]')
        vmd_lines.append('    $sel writepdb $outfile')
        vmd_lines.append('}')

        vmd_lines.append('mol delete top')

    vmd_lines.append('exit')

    script_text = "\n".join(vmd_lines)
    with tempfile.NamedTemporaryFile("w", suffix=".tcl", delete=False) as script_file:
        script_file.write(script_text)
        script_path = script_file.name

    try:
        subprocess.run(["vmd", "-dispdev", "text", "-e", script_path], check=True)
    finally:
        try:
            os.remove(script_path)
        except OSError:
            pass


def main():
    parser = argparse.ArgumentParser(description="Create protein-ligand PDB snapshots of a GPCRmd complex.")
    parser.add_argument("-dyn", dest="dynID", type=int, required=True, help="DynID of the GPCRmd complex.")
    parser.add_argument("-frames", dest="nframes", type=int, default=600, help="Number of frames to extract (default: last 600).")
    args = parser.parse_args()

    dynData= get_traj_path(args.dynID)
    output_dir = f"data/GPCRmd/{args.dynID}/{args.nframes}_fms"
    os.makedirs(output_dir, exist_ok=True)
    run_vmd(dynData['trajs'], dynData['psf'], output_dir, ligand_resname=dynData['lig_sname'])
    
    print(f"PDB files saved in {output_dir}")


if __name__ == "__main__":
    main()