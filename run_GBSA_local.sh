#!/bin/bash
# Non-SLURM smoke-test runner: the same pipeline as run_GBSA.sh (extract ->
# ligand prep -> clean -> topology -> match -> score), run directly and
# sequentially for all 3 replicates on any machine with /GPCRmd_files access -
# no sbatch, no array. Meant for quickly validating a dynID (few frames, fast)
# before committing to a full 600-frame sbatch production run.
#
# Usage: ./run_GBSA_local.sh [-dyn ID] [-frames N] [-nt THREADS]
# Defaults: -dyn 25 (matches run_GBSA.sh's current config), -frames 1, -nt 4
#
# Writes into the same data/GPCRmd/<dynID> and results/GPCRmd/<dynID>/fms_<N>_v2
# layout run_GBSA.sh uses, so ligand prep (mol2/net_charge/atom_map) done here
# is reused by a later full run of the same dynID, not redone.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_root}"

if [ -f "${repo_root}/local.env.sh" ]; then
  source "${repo_root}/local.env.sh"
fi

dynID=25
fmsIn=1
nt=4
while [ $# -gt 0 ]; do
  case "$1" in
    -dyn) dynID="$2"; shift 2 ;;
    -frames) fmsIn="$2"; shift 2 ;;
    -nt) nt="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

set +u
eval "$(conda shell.bash hook)"
gbsa_env_path="${gbsa_env_path:-/users/gpcr/mlopez/free_energies/gbsa_traj_pipeline/gbsa_env}"
conda activate "${gbsa_env_path}"
set -u

export OMPI_MCA_rmaps_base_oversubscribe=1

data_source="GPCRmd"
data_dyn_dir="data/${data_source}/${dynID}"
mkdir -p "${data_dyn_dir}"

# 0. Snapshots - only (re)extract if this exact frame count isn't already there.
if [ ! -f "${data_dyn_dir}/fms_${fmsIn}/initial_frame.pdb" ]; then
  echo "[$(date)] Extracting ${fmsIn} frame(s) for dyn ${dynID}"
  command -v module >/dev/null 2>&1 && module load VMD
  python3 scripts/extract_snapshots.py -dyn "${dynID}" -frames "${fmsIn}"
else
  echo "[$(date)] dyn ${dynID} already has fms_${fmsIn}, skipping extraction"
fi

dyn_config="${data_dyn_dir}/dyn_config.json"
lig_sdf=$(jq -r '.lig_sdf' "${dyn_config}")
lig_resname=$(jq -r '.lig_resname' "${dyn_config}")
mapfile -t trajIDs < <(jq -r '.trajIDs[]' "${dyn_config}")

frame_dir="${data_dyn_dir}/fms_${fmsIn}"
results_frame_dir="results/${data_source}/${dynID}/fms_${fmsIn}_v2"
mkdir -p "${results_frame_dir}"

# 1. Ligand mol2 - shared across trajIDs, built once.
ligand_mol2="${data_dyn_dir}/ligand_named.mol2"
if [ ! -f "${ligand_mol2}" ]; then
  echo "[$(date)] Building shared ligand mol2"
  python3 scripts/fix_ligand_naming.py \
    -i "${lig_sdf}" -o "${ligand_mol2}" -name "${lig_resname}" \
    -frame-pdb "${frame_dir}/initial_frame.pdb" -frame-resname "${lig_resname}"
fi

for trajID in "${trajIDs[@]}"; do
  echo "[$(date)] trajID=${trajID}"

  mapfile -t frame_files < <(
    for f in "${frame_dir}"/"${trajID}"_[0-9]*.pdb; do
      n=$(basename "$f" .pdb); n=${n##*_}
      echo "${n} ${f}"
    done | sort -n -k1,1 | cut -d' ' -f2-
  )
  frames=${#frame_files[@]}
  if [ "${frames}" -eq 0 ]; then
    echo "No frame PDB files found for ${trajID}" >&2
    exit 1
  fi

  # 2. Clean receptor.
  topo_dir="${data_dyn_dir}/topology_${trajID}_v2"
  mkdir -p "${topo_dir}"
  python3 scripts/clean_receptor_v2.py \
    -i "${frame_dir}/initial_frame.pdb" \
    -o "${topo_dir}/receptor_clean.pdb" \
    --drop-resname "${lig_resname}"

  # 3. Build topology.
  python3 scripts/build_topology_once.py \
    -r "${topo_dir}/receptor_clean.pdb" \
    -l "${ligand_mol2}" \
    -o "${topo_dir}" -nt "${nt}"

  # 4. Build the matched trajectory.
  gbsa_workdir="${results_frame_dir}/${trajID}_UNIGBSA"
  mkdir -p "${gbsa_workdir}"
  combined_pdb="${gbsa_workdir}/${trajID}_trajectory_${fmsIn}_frames.pdb"
  python3 scripts/build_trajectory_v2.py \
    -ref "${topo_dir}/complex.pdb" \
    -frames "${frame_files[@]}" \
    -ligand-resname "${lig_resname}" \
    -ligand-atom-map "${ligand_mol2}.atom_map" \
    -o "${combined_pdb}"

  # 5. Score.
  (
    cd "${gbsa_workdir}"
    unigbsa-traj \
      -i "${repo_root}/${topo_dir}/complex.pdb" \
      -p "${repo_root}/${topo_dir}/complex.top" \
      -ndx "${repo_root}/${topo_dir}/index.ndx" \
      -t "${repo_root}/${combined_pdb}" \
      -f "${repo_root}/config/mmpbsa_gb.in" \
      -nt "${nt}" \
      2>&1 | tee "${repo_root}/${gbsa_workdir}/${trajID}_gbsa_traj.log"
  )

  echo "[$(date)] trajID=${trajID} finished"
done

echo "[$(date)] All done. Check ${results_frame_dir}/*/Energy.csv"
