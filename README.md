# gbsa_traj_pipeline

First test of `unigbsa-traj` (Uni-GBSA) on a handful of MD snapshots of one
receptor-ligand complex, scored as-is (no minimization, no box/solvent).

Build the ligand+receptor topology **once**, then only swap coordinates per
frame. Full reasoning is in `NOTES.md`.

## Repo layout

```
gbsa_traj_pipeline/
├── data/            <- put your inputs here (see below)
├── config/
│   └── mmpbsa_gb.in <- GB parameters matching the earlier production runs (igb=8)
├── env/
│   └── gbsa_env.yml <- conda/micromamba environment spec
├── scripts/
│   ├── clean_receptor.py       <- strip a receptor out of a complex frame
│   ├── build_topology_once.py  <- build complex.top/pdb/index.ndx one time
│   └── build_trajectory.py     <- remap N frames onto that fixed atom order
├── topology/        <- created by step 2
└── results/         <- created by step 4
```

## 0. Environment

If you're on this machine, the env already exists:

```bash
micromamba activate gbsa
```

On a fresh machine, rebuild it from the spec in this repo:

```bash
micromamba create -f env/gbsa_env.yml
micromamba activate gbsa
```

Sanity check the two tools this pipeline actually calls:

```bash
gmx --version
unigbsa-traj -v
```

## 1. Put your inputs in `data/`

Replace the four empty placeholder files with your real ones (same names):

- `data/frame1.pdb`, `data/frame2.pdb`, `data/frame3.pdb` — each one a full
  receptor+ligand complex snapshot (protein + ligand atoms only, no
  membrane/water/ions), same atoms in the same order across all three.
- `data/ligand.mol2` — the ligand with **correct bond orders** already
  assigned (from acpype/antechamber output, a template you fixed by hand,
  or a matching PubChem/ChEMBL structure - see the caveat in `NOTES.md`
  if you go that route).

You also need the ligand's residue name **as it appears in the frame PDBs**
(not necessarily the same name used inside the mol2) - substitute it for
`LIG` in step 4 below.

## 2. Clean the receptor (once, from frame1 only)

This just needs one frame - it fixes the topology, not the per-frame pose.

```bash
python3 scripts/clean_receptor.py \
  -i data/frame1.pdb \
  -o data/receptor_clean.pdb \
  --drop-resname LIG
```

Replace `LIG` with your ligand's actual residue name so it gets stripped out
of the receptor-only file. Check the printed summary - it should report 0
unexpected dropped chains/HETATMs.

## 3. Build the topology once

```bash
python3 scripts/build_topology_once.py \
  -r data/receptor_clean.pdb \
  -l data/ligand.mol2 \
  -o topology
```

This runs `pdb2gmx` on the receptor and acpype/antechamber (GAFF2 + AM1-BCC)
on the ligand, once. Produces `topology/complex.top`, `topology/complex.pdb`
(the reference structure - this fixes the atom order everything else must
match) and `topology/index.ndx`.

## 4. Build the 3-frame trajectory, matched to that topology's atom order

```bash
python3 scripts/build_trajectory.py \
  -ref topology/complex.pdb \
  -frames data/frame1.pdb data/frame2.pdb data/frame3.pdb \
  -ligand-resname LIG \
  -o results/trajectory_3frames.pdb
```

Again replace `LIG`. Check the printed match counts per frame - the "kept
from reference" count should be small and roughly constant across frames
(that's expected: it's the receptor hydrogens `pdb2gmx` rebuilt in step 3,
which never existed in the raw MD frames). If the ligand match count is 0,
your mol2 and frame PDBs don't share atom names - fix that before trusting
any result past this point.

## 5. Score all 3 frames with GB

```bash
unigbsa-traj \
  -i topology/complex.pdb \
  -p topology/complex.top \
  -ndx topology/index.ndx \
  -t results/trajectory_3frames.pdb \
  -f config/mmpbsa_gb.in \
  -nt 4 \
  2>&1 | tee results/gbsa_traj_test.log
```

`-f config/mmpbsa_gb.in` is what pins this to igb=8 (matching the earlier
production runs) - `unigbsa-traj`'s own built-in default is igb=2, which is
silently different, so always pass `-f` here rather than relying on `-m`.

## 6. Read the result

`results/gbsa_traj_test.log` prints one line per frame:

```
Frames    mode    delta_G(kcal/mole)
     1      GB              -34.2681
     2      GB              -31.9053
     3      GB              -35.7420
```

That's your per-frame MM-GBSA binding energy - no minimization involved,
straight single-point GB scoring of each snapshot as given.
