#!/bin/bash
#SBATCH --job-name=gbsa_miguel
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32000
#SBATCH --output=logs/slurm_miguel.%A_%a.out
#SBATCH --error=logs/slurm_miguel.%A_%a.err
#SBATCH --array=0-2

set -euo pipefail

repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${repo_root}"

if [ -f "${repo_root}/local.env.sh" ]; then
  source "${repo_root}/local.env.sh"
fi

set +u
module load VMD
module load Miniconda3/202601
eval "$(conda shell.bash hook)"
gbsa_env_path="${gbsa_env_path:-/users/gpcr/mlopez/free_energies/gbsa_traj_pipeline/gbsa_env}"
conda activate "${gbsa_env_path}"
set -u

export OMPI_MCA_rmaps_base_oversubscribe=1

# --- config -------------------------------------------------------------
dyn_json="config/miguel_b2ar_dyns_test.json"
dynID=7          # see dyn_json for the full list - 7 = cyanopindolol, inactive
fmsIn=500
drop_resname="CLA,SOD,POP,TIP"
# -------------------------------------------------------------------------

data_dyn_dir="data/miguel_b2ar/dyn${dynID}"
mkdir -p "${data_dyn_dir}"

# 0. Extract wrapped snapshots - shared across replicates.
snap_lockfile="${data_dyn_dir}/.snapshots.lock"
snap_donefile="${data_dyn_dir}/.snapshots.done"
frame_dir="${data_dyn_dir}/fms_${fmsIn}"

(
  flock -x 201
  if [ ! -f "${snap_donefile}" ]; then
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: extracting wrapped snapshots for dyn ${dynID}"
    python3 scripts/extract_miguel_snapshots_test.py \
      -dyn "${dynID}" -input-json "${dyn_json}" -frames "${fmsIn}" -o "${frame_dir}"
    touch "${snap_donefile}"
  else
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: snapshots already extracted, skipping"
  fi
) 201>"${snap_lockfile}"

lig_resname=$(jq -r --argjson id "${dynID}" '.[] | select(.dyn_id==$id) | .lig_sname' "${dyn_json}")
lig_sdf=$(jq -r --argjson id "${dynID}" '.[] | select(.dyn_id==$id) | .ligand_sdf' "${dyn_json}")

results_frame_dir="results/miguel_b2ar/dyn${dynID}/fms_${fmsIn}"
mkdir -p "${results_frame_dir}"

# 1. Ligand mol2 - shared across replicates.
ligand_mol2="${data_dyn_dir}/ligand_named.mol2"
lockfile="${data_dyn_dir}/.ligand_prep.lock"
donefile="${data_dyn_dir}/.ligand_prep.done"

(
  flock -x 200
  if [ ! -f "${donefile}" ]; then
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: building shared ligand mol2"
    python3 scripts/fix_ligand_naming.py \
      -i "${lig_sdf}" -o "${ligand_mol2}" -name "${lig_resname}" \
      -frame-pdb "${frame_dir}/initial_frame.pdb" -frame-resname "${lig_resname}"
    touch "${donefile}"
  else
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ligand mol2 already prepared, skipping"
  fi
) 200>"${lockfile}"

# --- per-replicate work (one SLURM array task per nvt_N) -----------------
replicate="nvt_$((SLURM_ARRAY_TASK_ID + 1))"
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: replicate=${replicate}"

mapfile -t frame_files < <(
  for f in "${frame_dir}/${replicate}"_[0-9]*.pdb; do
    n=$(basename "$f" .pdb); n=${n##*_}
    echo "${n} ${f}"
  done | sort -n -k1,1 | cut -d' ' -f2-
)
frames=${#frame_files[@]}
if [ "${frames}" -eq 0 ]; then
  echo "No frame PDB files found for ${replicate}" >&2
  exit 1
fi
echo "[$(date)] Found ${frames} frame PDBs for ${replicate}"

# 2. Clean receptor.
topo_dir="${data_dyn_dir}/topology_${replicate}"
mkdir -p "${topo_dir}"
python3 scripts/clean_receptor_v2.py \
  -i "${frame_dir}/initial_frame.pdb" \
  -o "${topo_dir}/receptor_clean.pdb" \
  --drop-resname "${lig_resname},${drop_resname}"

# 3. Build topology.
python3 scripts/build_topology_once.py \
  -r "${topo_dir}/receptor_clean.pdb" \
  -l "${ligand_mol2}" \
  -o "${topo_dir}" -nt "${SLURM_CPUS_PER_TASK}"

# 4. Build the matched trajectory.
gbsa_workdir="${results_frame_dir}/${replicate}_UNIGBSA"
mkdir -p "${gbsa_workdir}"
combined_pdb="${gbsa_workdir}/${replicate}_trajectory_${fmsIn}_frames.pdb"
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
    -nt "${SLURM_CPUS_PER_TASK}" \
    2>&1 | tee "${repo_root}/${gbsa_workdir}/${replicate}_gbsa_traj.log"
)

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} finished: replicate=${replicate}"
