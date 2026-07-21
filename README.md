# gbsa_traj_pipeline

A verified, working MM-GBSA test: scores 3 real MD snapshots of one
receptor-ligand complex, no minimization - single-point GB scoring of
each frame exactly as given. Everything below uses the actual files
already in `data/` and produces the actual numbers we got running it.

Design reasoning and known limitations are in `NOTES.md` - read that
before running this on a different receptor/ligand.

## Repo layout

```
gbsa_traj_pipeline/
├── data/
│   ├── frame_100_GPCR_QGE.pdb    <- MD snapshot 1 (protein + ligand QGE only)
│   ├── frame_1000_GPCR_QGE.pdb   <- MD snapshot 2
│   ├── frame_2000_GPCR_QGE.pdb   <- MD snapshot 3
│   └── 13903_mol_5236.sdf        <- QGE ligand, correct bond orders (C18H17NO4S)
├── config/mmpbsa_gb.in           <- GB params matching earlier production runs (igb=8)
├── env/gbsa_env.yml              <- environment spec
├── scripts/
│   ├── fix_ligand_naming.py      <- SDF -> uniquely-named mol2 (once, for the ligand)
│   ├── clean_receptor.py         <- strip receptor out of a complex frame
│   ├── build_topology_once.py    <- build complex.top/pdb/index.ndx one time
│   └── build_trajectory.py       <- remap the 3 frames onto that fixed atom order
├── topology/                     <- created by step 3
└── results/                      <- created by step 4-5
```

## Install

Already exists on this machine:

```bash
micromamba activate gbsa
```

Fresh machine:

```bash
micromamba create -f env/gbsa_env.yml
micromamba activate gbsa
```

Sanity check:

```bash
gmx --version
unigbsa-traj -v
```

## Run the test

```bash
cd gbsa_traj_pipeline

# 1. Ligand: SDF -> mol2 with unique, indexed atom names (once)
python3 scripts/fix_ligand_naming.py \
  -i data/13903_mol_5236.sdf -o data/ligand_named.mol2 -name QGE

# 2. Receptor: clean one frame, dropping the QGE ligand atoms
python3 scripts/clean_receptor.py \
  -i data/frame_100_GPCR_QGE.pdb -o data/receptor_clean.pdb --drop-resname QGE

# 3. Build the topology once (pdb2gmx + acpype/antechamber AM1-BCC - takes a bit)
python3 scripts/build_topology_once.py \
  -r data/receptor_clean.pdb -l data/ligand_named.mol2 -o topology -nt 4

# 4. Build the 3-frame trajectory, matched to that topology's atom order
python3 scripts/build_trajectory.py \
  -ref topology/complex.pdb \
  -frames data/frame_100_GPCR_QGE.pdb data/frame_1000_GPCR_QGE.pdb data/frame_2000_GPCR_QGE.pdb \
  -ligand-resname QGE \
  -o results/trajectory_3frames.pdb

# 5. Score all 3 frames with GB (igb=8, matching earlier production runs)
unigbsa-traj \
  -i topology/complex.pdb -p topology/complex.top -ndx topology/index.ndx \
  -t results/trajectory_3frames.pdb -f config/mmpbsa_gb.in -nt 3 \
  2>&1 | tee results/gbsa_traj_test.log
```

(All commands above assume `micromamba activate gbsa` from the Install step is already active in your shell.)

Note on `-nt 3` in step 5: `gmx_MMPBSA` requires processes ≤ frames -
`-nt` here must not exceed how many frames you're scoring.

## What to check as you go (don't skip these)

- **Step 2** prints a summary - 0 unexpected dropped chains/HETATMs.
- **Step 4** prints per-frame match counts. On this data you should see
  `5076 matched, 10 kept-from-reference` for every frame. The 10 kept
  are known, checked-safe leftovers (N-terminal cap hydrogens, 5 distant
  histidines - see `NOTES.md`), not a new problem. If your numbers are
  very different from this, don't trust the energies - see `NOTES.md`
  for how to debug it (that's exactly how a real bug was caught here:
  a naive version of this matched 0 ligand atoms and separately let a
  hydrogen clash so hard it overflowed the energy calculation).

## Expected result

`results/gbsa_traj_test.log` should end with a table like:

```
   Frames mode       complex  ...      Gas  Solvation      TOTAL
0       1   gb  9.943378e+05  ... -36.3712   9.487106 -26.884094
1       2   gb  1.878720e+06  ... -37.9014   8.091573 -29.809827
2       3   gb  1.146861e+06  ... -37.3585   8.818791 -28.539709
```

`TOTAL` is the per-frame MM-GBSA binding energy in kcal/mol: **-26.88,
-29.81, -28.54** for frames 100/1000/2000. The large `complex` column
values are the receptor's own internal self-energy and cancel out in
the complex-receptor-ligand subtraction - they're not a sign of a
problem as long as `TOTAL` comes out in this normal range (roughly
-10 to -50 kcal/mol for a drug-like ligand). If instead you see
`TOTAL` missing, NaN, or `gmx_MMPBSA` crashing with "Some energy terms
are undefined," that's the clash/overflow failure mode - go back to
step 4's match counts.

## Adapting this to a different receptor/ligand/snapshot count

Swap the file names in steps 1-5 for your own. Two things from
`NOTES.md` you will hit immediately:

1. Whatever your ligand's residue name is in the raw frames, pass it as
   `-ligand-resname` in step 4 - **not** as `-ref-ligand-resname`, which
   should stay `MOL` (`unigbsa` hardcodes the built ligand's resname to
   `MOL` regardless of what your mol2 called it).
2. `-nt` in step 5 must be ≤ your frame count.

For more than a handful of frames, or a receptor with histidines/capped
termini near the binding site, read `NOTES.md`'s caveats section first -
this repo's matching logic has known gaps that happened not to matter
for this specific test.
