# gbsa_traj_pipeline

MM-GBSA scoring of MD snapshots against one built-once topology, no
minimization - single-point GB scoring of each frame exactly as given.

**Status: fixed and verified (v2).** The bugs that caused some frames to
clash after matching (frozen reference atoms overlapping real ones, plus
a frame-ordering bug and an MPI launch issue on the cluster) are fixed
and confirmed clean across a real 3×600-frame production run - see
`NOTES.md` for the full design reasoning and debugging history.

## Repo layout

```
gbsa_traj_pipeline/
├── data/                          <- not in git; extracted frames + ligand SDF live here
├── config/mmpbsa_gb.in            <- GB params (igb=8, intdiel=4.0, extdiel=78.5, saltcon=0.15)
├── env/gbsa_env.yml               <- environment spec
├── run_GBSA.sh                    <- cluster entry point (sbatch), runs the full pipeline below
├── scripts/
│   ├── extract_snapshots.py       <- VMD: GPCRmd dyn -> per-frame PDBs (once per dynID)
│   ├── fix_ligand_naming.py       <- SDF -> uniquely-named mol2 (once, for the ligand)
│   ├── clean_receptor_v2.py       <- strip + fix receptor for pdb2gmx (use this, not v1)
│   ├── build_topology_once.py     <- build complex.top/pdb/index.ndx one time
│   └── build_trajectory_v2.py     <- remap frames onto that fixed atom order (use this, not v1)
└── results/                       <- not in git; per-trajID Energy.csv/Dec.csv land here
```

`clean_receptor.py`/`build_trajectory.py` (no `_v2`) still exist for
reference but are the old, buggy versions - don't use them for a real
run. See `NOTES.md` if you're curious why they were replaced.

## Cluster usage (the normal way to run this)

```bash
git clone <this repo> ~/free_energies/gbsa_traj_pipeline
cd ~/free_energies/gbsa_traj_pipeline
mkdir -p logs
sbatch run_GBSA.sh
```

`run_GBSA.sh` already:
- activates Marta's existing built env (`gbsa_env`) directly - read-only
  use, doesn't interfere with her work, no need to build your own
- runs the full pipeline (extract -> clean -> topology -> match -> score)
  for all 3 GPCRmd replicates via `--array=0-2`
- requests 32 cores per replicate and includes the OpenMPI fix needed
  for `gmx_MMPBSA`'s internal `mpirun` to actually use them

**Only thing you need to change to run this as yourself:** the two
hardcoded path lines near the top of `run_GBSA.sh`:
```bash
cd /users/gpcr/<you>/free_energies/gbsa_traj_pipeline
repo_root=/users/gpcr/<you>/free_energies/gbsa_traj_pipeline
```
Everything else (Marta's env, the OpenMPI fix, core count) works as-is.

If you'd rather not reuse Marta's env (e.g. no `gpcr` group read access),
build your own instead and swap the `conda activate` line accordingly:
```bash
conda env create -f env/gbsa_env.yml -p ./gbsa_env
```

**Data isn't in git** (`data/`/`results/` are gitignored). `run_GBSA.sh`
extracts it automatically (step 0) the first time it runs for a given
`dynID` - nothing manual needed if you're just running `sbatch
run_GBSA.sh`. To pull frames for a **different** GPCRmd simulation
yourself (own use outside the sbatch flow), it's the same script
`run_GBSA.sh` calls, run directly:
```bash
module load VMD
python3 scripts/extract_snapshots.py -dyn <dynID> -frames 600
```
`<dynID>` is the GPCRmd dynamics ID for the simulation you want (`1033`
for the one this repo's fixes were verified against). This reads the
shared `/GPCRmd_files/Precomputed/compl_info_new.json` (not anyone's
personal directory) and writes `data/GPCRmd/<dynID>/fms_600/` plus
`dyn_config.json`, exactly like step 0 of the pipeline does - it's fully
self-contained, not dependent on anyone else's already-run copy.

## What "working" looks like

- `build_trajectory_v2.py`'s per-frame printout should show `0
  kept-from-reference` for both protein and ligand, every frame. Any
  nonzero count means an atom didn't match and is frozen at a stale
  position - don't trust the energies until that's 0 (see `NOTES.md`
  for how to debug it).
- The per-trajID log (`results/.../<trajID>_gbsa_traj.log`) should end
  with a clean results table and exit status 0, not a crash.
- `TOTAL` in `Energy.csv` should sit in a normal range (roughly -10 to
  -50 kcal/mol for a drug-like ligand), `complex`/`receptor`/`ligand`
  in the hundreds of thousands - not NaN, not overflowing, not billions.

## Adapting this to a different receptor/ligand

1. Whatever your ligand's residue name is in the raw frames, pass it as
   `-ligand-resname` - **not** as `-ref-ligand-resname`, which should
   stay `MOL` (`unigbsa` hardcodes the built ligand's resname to `MOL`
   regardless of what your mol2 called it).
2. `-nt` must be ≤ your frame count.
3. `fix_ligand_naming.py` assumes the SDF's atom order already matches
   the raw MD frame's atom order/naming - true for this ligand, not
   guaranteed for a new one. Check `build_trajectory_v2.py`'s ligand
   match count (should be 100%, not partial) before trusting a new
   ligand's results.

Read `NOTES.md` before running this on a receptor with a different
construct (e.g. different termini/caps) - the v2 fixes here are specific
to what this receptor actually has (ACE/NME caps, several histidine
tautomers, 4 disulfides); a different receptor could have different
gaps that need the same kind of direct verification, not just assumed
to be fine.
