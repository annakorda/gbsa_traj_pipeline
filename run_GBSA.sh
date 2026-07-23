#!/bin/bash
#SBATCH --job-name=gbsa_6ZDR
#SBATCH --partition=long
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --output=results/GPCRmd/1033/fms_4/slurm.%A_%a.out
#SBATCH --error=results/GPCRmd/1033/fms_4/slurm.%A_%a.err
#SBATCH --array=0-2

set -euo pipefail

cd /users/gpcr/mlopez/free_energies/gbsa_traj_pipeline
repo_root=/users/gpcr/mlopez/free_energies/gbsa_traj_pipeline

set +u
module load Miniconda3/202601
eval "$(conda shell.bash hook)"
conda activate ./gbsa_env
set -u

# --- config -------------------------------------------------------------
dynID=1033
fms=4
trajIDs=(17356_trj_1033 17357_trj_1033 17358_trj_1033)
lig_sdf=data/13903_mol_5236.sdf
lig_resname=QGE
# -------------------------------------------------------------------------

data_dyn_dir="data/GPCRmd/${dynID}"            # dynID-level data dir: topology + ligand mol2 live here
frame_dir="${data_dyn_dir}/fms_${fms}"          # per-fms frame PDBs
results_frame_dir="results/GPCRmd/${dynID}/fms_${fms}"

# mkdir -p "${data_dyn_dir}"
mkdir -p "${results_frame_dir}"

ligand_mol2="${data_dyn_dir}/ligand_named.mol2"     # now dynID-level, not fms-level
lockfile="${data_dyn_dir}/.ligand_prep.lock"
donefile="${data_dyn_dir}/.ligand_prep.done"

# 1. Ligand: SDF -> mol2 with unique, indexed atom names.
#    Shared across all array tasks -> serialize with flock so only the
#    first task to grab the lock builds it; the rest wait, see donefile,
#    and skip straight through.
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

topo_dir="${data_dyn_dir}/topology_${trajID}"   # now nested under data/GPCRmd/${dynID}/
mkdir -p "${topo_dir}"

# Check required snapshots (pdbs) match the introduced number of frames per traj id
frame_files=(${frame_dir}/${trajID}_[0-9]*.pdb)
frames=${#frame_files[@]}
if [ "${frames}" -eq 0 ]; then
  echo "No frame PDB files found for ${trajID}" >&2
  exit 1
fi
echo "[$(date)] Found ${frames} frame PDBs for ${trajID}"

# 2. Receptor PDB: clean one frame, dropping the ligand atoms
echo "2. python3 scripts/clean_receptor.py \
  -i ${frame_dir}/${trajID}_1.pdb \
  -o ${frame_dir}/${trajID}_receptor_clean.pdb \
  --drop-resname ${lig_resname}"

python3 scripts/clean_receptor.py \
  -i "${frame_dir}/${trajID}_1.pdb" \
  -o "${frame_dir}/${trajID}_receptor_clean.pdb" \
  --drop-resname "${lig_resname}"

# 3. Build the gromacs topology for this trajID (pdb2gmx + acpype/antechamber AM1-BCC)
echo "3. python3 scripts/build_topology_once.py \
  -r ${frame_dir}/${trajID}_receptor_clean.pdb \
  -l ${ligand_mol2} \
  -o ${topo_dir} -nt ${SLURM_CPUS_PER_TASK}"

python3 scripts/build_topology_once.py \
  -r "${frame_dir}/${trajID}_receptor_clean.pdb" \
  -l "${ligand_mol2}" \
  -o "${topo_dir}" -nt "${SLURM_CPUS_PER_TASK}"

# 4. Build the multi-frame trajectory, matched to that topology's atom order.
#    Combined pdb now lives inside this trajID's UNIGBSA results folder.
gbsa_workdir="${results_frame_dir}/${trajID}_UNIGBSA"
mkdir -p "${gbsa_workdir}"

combined_pdb="${gbsa_workdir}/${trajID}_trajectory_${fms}_frames.pdb"

echo "4. python3 scripts/build_trajectory.py \
  -ref ${topo_dir}/complex.pdb \
  -frames ${frame_files[@]} \
  -ligand-resname ${lig_resname} \
  -o ${combined_pdb}"

python3 scripts/build_trajectory.py \
  -ref "${topo_dir}/complex.pdb" \
  -frames "${frame_files[@]}" \
  -ligand-resname "${lig_resname}" \
  -o "${combined_pdb}"

# 5. Score all frames with GB (igb=8, matching earlier production runs).
#    Run from inside gbsa_workdir so gmx_MMPBSA's scratch files (UNIGBSA/*)
#    and the log stay isolated per trajID, avoiding cross-task collisions.
echo "5. unigbsa-traj \
  -i ${topo_dir}/complex.pdb -p ${topo_dir}/complex.top -ndx ${topo_dir}/index.ndx \
  -t ${combined_pdb} -f config/mmpbsa_gb.in -nt ${SLURM_CPUS_PER_TASK}"

(
  cd "${gbsa_workdir}"
  unigbsa-traj \
    -i "${repo_root}/${topo_dir}/complex.pdb" \
    -p "${repo_root}/${topo_dir}/complex.top" \
    -ndx "${repo_root}/${topo_dir}/index.ndx" \
    -t "${repo_root}/${combined_pdb}" \
    -f "${repo_root}/config/mmpbsa_gb.in" \
    -nt "${SLURM_CPUS_PER_TASK}" \
    2>&1 | tee "${repo_root}/${gbsa_workdir}/${trajID}_gbsa_traj.log"
)

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} finished: trajID=${trajID}"