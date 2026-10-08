#!/usr/bin/env python3
"""
Extract per-frame PDBs (protein+ligand) for one dyn_id from
config/miguel_b2ar_dyns_test.json - one shared structure.psf + one wrapped
trajectory per replicate. Mirrors extract_snapshots.py's output:
initial_frame.pdb + nvt_<N>_<frame>.pdb in the given output dir.
"""
import argparse
import json
import os
import subprocess
import tempfile


def get_dyn_data(dyn_id, input_json):
    with open(input_json) as f:
        data = json.load(f)
    for entry in data.values():
        if entry["dyn_id"] == dyn_id:
            return entry
    raise KeyError(f"dyn_id {dyn_id} not found in {input_json}")


def run_vmd(psf_path, traj_paths, output_dir, fms, ligand_resname):
    if not os.path.isfile(psf_path):
        raise FileNotFoundError(f"PSF file not found: {psf_path}")
    os.makedirs(output_dir, exist_ok=True)

    vmd_lines = []
    for idx, traj_path in enumerate(traj_paths):
        if not os.path.isfile(traj_path):
            raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

        label = f"nvt_{idx + 1}"
        traj_ext = os.path.splitext(traj_path)[1].lstrip(".")
        vmd_lines.append(f'mol new "{psf_path}" type psf waitfor all')
        vmd_lines.append(f'mol addfile "{traj_path}" type {traj_ext} waitfor all')
        vmd_lines.append(f'set sel [atomselect top "protein or resname {ligand_resname}"]')

        if idx == 0:
            vmd_lines.append('animate goto 0')
            vmd_lines.append(f'$sel writepdb "{output_dir}/initial_frame.pdb"')

        vmd_lines.append('set nf [molinfo top get numframes]')
        vmd_lines.append(f'set start [expr {{$nf - {fms} < 0 ? 0 : $nf - {fms}}}]')
        vmd_lines.append('for {set i $start} {$i < $nf} {incr i} {')
        vmd_lines.append('    animate goto $i')
        vmd_lines.append(f'    set outfile [format "{output_dir}/{label}_%d.pdb" [expr {{$i - $start + 1}}]]')
        vmd_lines.append('    $sel writepdb $outfile')
        vmd_lines.append('}')
        vmd_lines.append('mol delete top')
    vmd_lines.append('exit')

    script_text = "\n".join(vmd_lines)
    with tempfile.NamedTemporaryFile("w", suffix=".tcl", delete=False) as f:
        f.write(script_text)
        script_path = f.name
    try:
        subprocess.run(["vmd", "-dispdev", "text", "-e", script_path], check=True)
    finally:
        try:
            os.remove(script_path)
        except OSError:
            pass


def main():
    parser = argparse.ArgumentParser(description="Extract per-frame PDBs for a Miguel B2AR dyn (wrapped trajectories).")
    parser.add_argument("-dyn", dest="dyn_id", type=int, required=True)
    parser.add_argument("-input-json", dest="input_json", required=True)
    parser.add_argument("-frames", dest="nframes", type=int, default=None,
                         help="Frames per replicate (default: the JSON's own framenums).")
    parser.add_argument("-o", dest="output_dir", required=True)
    args = parser.parse_args()

    entry = get_dyn_data(args.dyn_id, args.input_json)
    fms = args.nframes if args.nframes is not None else max(entry["framenums"])
    run_vmd(entry["topo_f"], entry["traj_f"], args.output_dir, fms, entry["lig_sname"])
    print(f"PDB files saved in {args.output_dir}")


if __name__ == "__main__":
    main()
