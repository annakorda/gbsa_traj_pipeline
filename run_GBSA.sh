#!/bin/bash
#SBATCH --job-name=gbsa_6ZDR
#SBATCH --partition=normal
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --output=logs/slurm.%A_%a.out
#SBATCH --error=logs/slurm.%A_%a.err
#SBATCH --array=0-2

set -euo pipefail

cd /users/gpcr/mlopez/free_energies/gbsa_traj_pipeline
repo_root=/users/gpcr/mlopez/free_energies/gbsa_traj_pipeline

set +u
module load VMD
module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate ./gbsa_env
set -u

# --- config -------------------------------------------------------------
dynID=1033
fmsIn=600
data_dyn_dir="data/GPCRmd/${dynID}"
mkdir -p "${data_dyn_dir}"
# -------------------------------------------------------------------------

# 0. Create PDB snapshots of the three replicates of the GPCRmd trajectory
dyn_config="${data_dyn_dir}/dyn_config.json"
snap_lockfile="${data_dyn_dir}/.snapshots.lock"
snap_donefile="${data_dyn_dir}/.snapshots.done"

(
  flock -x 201
  if [ ! -f "${snap_donefile}" ]; then
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: extracting snapshots"
    python3 scripts/extract_snapshots.py -dyn "${dynID}" -frames "${fmsIn}"
    touch "${snap_donefile}"
  else
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: snapshots already extracted, skipping"
  fi
) 201>"${snap_lockfile}"

# Read config written by extract_snapshots.py
fms=$(jq -r '.fms' "${dyn_config}")
lig_sdf=$(jq -r '.lig_sdf' "${dyn_config}")
lig_resname=$(jq -r '.lig_resname' "${dyn_config}")
mapfile -t trajIDs < <(jq -r '.trajIDs[]' "${dyn_config}")

frame_dir="${data_dyn_dir}/fms_${fms}"
# _v2 suffix: keeps this run's topology/results separate from the original
# v1 run still sitting in fms_${fms}/ - nothing here overwrites that.
results_frame_dir="results/GPCRmd/${dynID}/fms_${fms}_v2"
mkdir -p "${results_frame_dir}"

# 1. Ligand: SDF -> mol2 with unique, indexed atom names. (shared, read-only after this)
ligand_mol2="${data_dyn_dir}/ligand_named.mol2"
lockfile="${data_dyn_dir}/.ligand_prep.lock"
donefile="${data_dyn_dir}/.ligand_prep.done"

(
  flock -x 200
  if [ ! -f "${donefile}" ]; then
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: building shared ligand mol2"
    python3 scripts/fix_ligand_naming.py \
      -i "${lig_sdf}" -o "${ligand_mol2}" -name "${lig_resname}"
    touch "${donefile}"
  else
    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ligand mol2 already prepared, skipping"
  fi
) 200>"${lockfile}"

# --- per-trajID work -----------------------------------------------------
trajID=${trajIDs[$SLURM_ARRAY_TASK_ID]}
echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: trajID=${trajID}"

# Check required snapshots (pdbs) match the introduced number of frames per traj id
# NOTE: a plain glob here sorts lexicographically (..._100.pdb before ..._2.pdb),
# scrambling frame order - confirmed by testing bash's actual glob expansion.
# MODEL N in the scored trajectory must be the N-th frame numerically, not the
# N-th glob match, so sort numerically by the trailing frame number instead.
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
echo "[$(date)] Found ${frames} frame PDBs for ${trajID}"

# 2. Clean receptor: each task uses the SAME shared initial_frame.pdb as input,
#    but writes its own output into its own per-trajID folder -> no write collision.
# _v2 topo_dir: separate from the original v1 topology_${trajID}/, which used
# clean_receptor.py/build_trajectory.py (v1) and is left untouched for comparison.
topo_dir="${data_dyn_dir}/topology_${trajID}_v2"
mkdir -p "${topo_dir}"

echo "2. python3 scripts/clean_receptor_v2.py \
  -i ${frame_dir}/initial_frame.pdb \
  -o ${topo_dir}/receptor_clean.pdb \
  --drop-resname ${lig_resname}"

python3 scripts/clean_receptor_v2.py \
  -i "${frame_dir}/initial_frame.pdb" \
  -o "${topo_dir}/receptor_clean.pdb" \
  --drop-resname "${lig_resname}"

# 3. Build topology for this trajID's own folder (duplicated, but identical
#    across all three since they're all built from the same initial_frame.pdb)
echo "3. python3 scripts/build_topology_once.py \
  -r ${topo_dir}/receptor_clean.pdb \
  -l ${ligand_mol2} \
  -o ${topo_dir} -nt ${SLURM_CPUS_PER_TASK}"

python3 scripts/build_topology_once.py \
  -r "${topo_dir}/receptor_clean.pdb" \
  -l "${ligand_mol2}" \
  -o "${topo_dir}" -nt "${SLURM_CPUS_PER_TASK}"

# 4. Build the multi-frame trajectory, matched to this trajID's own topology.
gbsa_workdir="${results_frame_dir}/${trajID}_UNIGBSA"
mkdir -p "${gbsa_workdir}"

combined_pdb="${gbsa_workdir}/${trajID}_trajectory_${fms}_frames.pdb"

echo "4. python3 scripts/build_trajectory_v2.py \
  -ref ${topo_dir}/complex.pdb \
  -frames ${frame_files[@]} \
  -ligand-resname ${lig_resname} \
  -o ${combined_pdb}"

python3 scripts/build_trajectory_v2.py \
  -ref "${topo_dir}/complex.pdb" \
  -frames "${frame_files[@]}" \
  -ligand-resname "${lig_resname}" \
  -o "${combined_pdb}"

# 5. Score all frames with GB.
echo "5. unigbsa-traj \
  -i ${topo_dir}/complex.pdb -p ${topo_dir}/complex.top -ndx ${topo_dir}/index.ndx \
  -t ${combined_pdb} -f config/mmpbsa_gb.in -nt ${SLURM_CPUS_PER_TASK}"

(
  cd "${gbsa_workdir}"
  /usr/bin/time -v unigbsa-traj \
    -i "${repo_root}/${topo_dir}/complex.pdb" \
    -p "${repo_root}/${topo_dir}/complex.top" \
    -ndx "${repo_root}/${topo_dir}/index.ndx" \
    -t "${repo_root}/${combined_pdb}" \
    -f "${repo_root}/config/mmpbsa_gb.in" \
    -nt "${SLURM_CPUS_PER_TASK}" \
    2>&1 | tee "${repo_root}/${gbsa_workdir}/${trajID}_gbsa_traj.log"
)

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} finished: trajID=${trajID}"