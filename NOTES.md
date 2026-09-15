# Design notes / gotchas found in the unigbsa source

## Why build the topology once and swap coordinates, instead of rebuilding per frame

A GROMACS `.top` only encodes chemistry (atom types, charges, bonds) - it
has no coordinates. Coordinates always live in a separate file. So "build
once" fixes the ligand's GAFF2 atom types + AM1-BCC charges one time (the
expensive, QM-based step), and every frame after that only needs a
coordinates file in the same atom order - no re-parametrization, no risk of
antechamber flipping a bond order between frames.

Checked directly in `unigbsa/simulation/topology.py::build_topol`: if the
ligand file ends in `.pdb`, it's routed through `build_protein()`
(`gmx pdb2gmx`, i.e. treated as a peptide) - **not** through acpype/GAFF.
A small-molecule ligand must be mol2/mol/sdf or it never gets GAFF atom
types/AM1-BCC charges at all.

## Getting a ligand mol2 with correct bond orders (no SDF from your own sim)

Three options, in order of reliability:

1. **A matching PubChem/ChEMBL SDF** (if the exact compound + protonation
   state is known) - most reliable for bond orders, since it's not
   geometry-guessed. The catch: its atom *names* almost certainly won't
   match the names used in your MD system's ligand - this is expected and
   now handled automatically, see "Ligand atom identity" below. 3D
   coordinates in the SDF don't matter either way - only the 2D bond graph
   and formal charges do, since coordinates get overwritten by each
   frame's real pose downstream.
2. **Maestro** (or any tool you already trust for bond assignment) - fine,
   same rule applies: only the bonds/formal-charge assignment need to be
   right, not the specific conformer geometry, since that also gets
   replaced by each frame's real coordinates.
3. **OpenBabel geometry perception** from one extracted reference-frame
   ligand PDB (`obabel ligand.pdb -O ligand.mol2`) - available in the `gbsa`
   env already, no RDKit needed. Works fine for typical drug-like molecules
   but sanity-check aromatic rings and net formal charge once before
   trusting it for every frame.

Whichever you use, do it **once** against one reference geometry - the
bond graph doesn't change frame to frame, only the 3D pose does.

## GB parameters: don't rely on unigbsa-traj's own defaults

`unigbsa-traj`'s built-in default (when no `-f` is given) is `igb=2`,
`extdiel=78.5` - confirmed from `unigbsa/data/mmpbsa.in`. The earlier
c0/c3/c5/c9/exp0 production runs used `igb=8, indi=4.0, exdi=80.0,
saltcon=0.15` (from `unigbsa_em.ini`). These are different GB models -
always pass `config/mmpbsa_gb.in` via `-f` to pin it to the same settings
used before, rather than trusting the `-m` flag's defaults.

Also found: `unigbsa/gbsa/parameters.py`'s template-writer only includes a
mode's section if the mode string matches **case-sensitively** against the
lowercase tags in the template file (`&gb`, `&pb`, ...), but the CLI's `-m`
flag defaults to `['GB']` (uppercase) - passing `-m GB` without `-f` would
silently produce an input file with no `&gb` section at all. Using `-f`
(as this repo does) bypasses that code path entirely, so it's not just
a style preference - it avoids a real bug in the templating logic.

## Matching assumption in build_trajectory_v2.py

Protein atoms are matched frame-to-frame by `(residue ORDINAL POSITION,
atom name)` - **not** by literal chain/resnum. Confirmed: `pdb2gmx` drops
the chain ID and renumbers residues (the ligand occupies residue slot 1
in the built topology, shifting every protein residue by +1), so literal
`(chain, resnum)` never lines up. Ordinal position (1st residue
encountered, 2nd, ...) survives that renumbering as long as residue
order/count is unchanged between the reference topology and each raw
frame. The ligand is matched by atom name alone (there's one copy).

A small, explicit alias table (`AMBER_TO_CHARMM_ALIASES` in
`build_trajectory_v2.py`) covers protein atom names that legitimately
differ between AMBER's naming (used by the built topology) and CHARMM's
naming (used by GPCRmd's raw frames) - backbone amide H, N-/C-terminal
H's, a couple of sidechain H's. This list was built by actually diagnosing
real clashes one at a time, not derived from a spec - it may not be
exhaustive for a receptor with residue types this one didn't stress.
**Always check the printed match counts before trusting a run** - see
"The frozen-atom clash bug" below for what a nonzero count means and how
to debug it.

## The frozen-atom clash bug, and why it mattered

Any reference atom `build_trajectory_v2.py` can't match in a given frame
keeps its topology-build-time coordinate, frozen, for every output frame.
This is silent - no error, no warning - it just quietly substitutes a
stale position. It only shows up as a clash in whichever frames happen to
bring a real, moving neighbor close to that fixed point, which is why the
symptom is "some frames spike, most don't," not a clean crash every time.

Three concrete sources of frozen atoms were found and fixed by directly
diagnosing this construct's actual chemistry (not by guessing):

1. **Histidine tautomer.** Collapsing `HSD`/`HSE`/`HSP` (CHARMM tautomer
   names) to generic `HIS` before `pdb2gmx` ever sees them is a trap:
   `pdb2gmx` is called with no `-his` flag
   (`unigbsa/simulation/topology.py::build_protein`), so it silently
   applies one uniform default tautomer to every histidine - wrong for
   whichever ones weren't actually that state. Fix (`clean_receptor_v2.py`):
   translate to the literal AMBER names (`HID`/`HIE`/`HIP`) instead of
   collapsing them - confirmed `amber99sb-ildn.ff/aminoacids.rtp` defines
   these directly, so `pdb2gmx` builds the correct sidechain per residue
   automatically, no flag or prompt needed.
2. **Disulfides must be communicated to pdb2gmx by residue naming, not
   CONECT.** Detecting SG-SG distances ≤2.2 Å and writing `CONECT`
   records doesn't work - confirmed via `pdb2gmx`'s own log output that
   **`pdb2gmx` ignores CONECT entirely** (`WARNING: all CONECT records
   are ignored`). Relying on `pdb2gmx`'s own distance-based auto-detect
   (via `specbond.dat`) instead triggers an interactive per-pair prompt
   that a scripted, piped run can't reliably answer for every pair. Fix:
   rename both residues of a detected disulfide pair `CYS`→`CYX` directly
   - confirmed `amber99sb-ildn.ff` defines `CYX`, so the bond forms
   deterministically with zero ambiguity, verified against `pdb2gmx`'s
   own geometry table (the renamed pairs sit at 0.2-1.0 Å in its internal
   SG-SG distance matrix; everything left as plain `CYS` sits several Å
   apart - the renaming picked exactly the same pairs `pdb2gmx`'s own
   geometry agrees are real bonds).
3. **Terminal caps must be translated, not dropped.** This construct's raw
   frames have a CHARMM N-terminal ACE patch fused onto the first residue
   (`CAY`/`HY1-3`/`CY`/`OY`) and a C-terminal NME patch fused onto the
   last (`NT`/`HNT`/`CAT`/`HT1-3`). Simply dropping these atoms leaves
   `pdb2gmx` to build ordinary charged termini instead of the true
   neutral capped ones.

   **This one took two failed attempts to fix correctly - worth
   understanding why, since it's a real trap for anyone touching this
   again.** Splitting the cap atoms into a standalone `ACE`/`NME` residue
   on the **topology-building side only** (`clean_receptor_v2.py`) builds
   correctly - confirmed via `pdb2gmx`'s log: `Identified residue ACE0 as
   a starting terminus` / `Identified residue NME320 as an ending
   terminus`, zero ambiguity. But it silently broke matching:
   - Splitting at the **front** (ACE) changes the reference's residue
     *count*, which shifts every downstream residue's ordinal index by +1
     relative to the raw frame (which still has ACE fused into residue 1,
     unsplit). Since most residues share common atom names (`N`, `CA`,
     `C`, `O`, `CB`...), this doesn't just fail to match - it **silently
     assigns the wrong residue's real coordinates** almost everywhere
     downstream (confirmed: reference ordinal 1 = MET's real atoms ended
     up matched against frame ordinal 1 = PRO's real atoms, since both
     have atoms named `N`/`CA`/etc). Undetectable by match-count alone,
     and far worse than the original frozen-atom bug.
   - Splitting at the **back** (NME) doesn't shift anything downstream
     (nothing follows it), but the new ordinal slot has nothing to match
     against in the raw frame at all - the *entire* NME residue (6 atoms:
     `N`,`H`,`CH3`,`HH31-33`) froze instead of just the 2 heavy atoms that
     froze when it was simply dropped. Splitting made it worse, not
     better.

   **The fix that actually works**: split identically on *both* sides -
   `clean_receptor_v2.py` splits the reference topology, and
   `build_trajectory_v2.py` applies the same split to every raw frame
   before computing ordinals (`split_terminal_caps`/
   `split_frame_terminal_caps` in each file respectively), so residue
   counts and ordinal alignment stay in sync throughout. Verified by
   directly re-checking ordinal 0/1 and the last two ordinals against
   both reference and frame content after the fix - not just trusting a
   lower kept-count. **The two scripts must be used as a pair** - matching
   `clean_receptor_v2.py`'s output against a mismatched frame-side split
   will silently misalign every residue.

   A real implementation bug was also caught this way, worth flagging as
   a general caution: an early version of the frame-side split computed
   the new residue number by reading `block[0]['resseq']` *inside* the
   atom-renaming loop - but `block[0]` was itself one of the atoms being
   renamed, so later iterations read their own already-mutated sibling's
   value and compounded it, silently splitting one residue into several
   spurious ones. Cache values you're about to mutate *before* the loop
   that mutates them, not inside it.

## Diagnosing a real clash: the actual playbook used

When a frame's `VDWAALS` overflows (prints as `*************` in sander's
mdout) or `TOTAL` looks wrong, this is how the real cause was tracked
down on this construct's production data - reusable for the next one:

1. **Get the real problem frame's exact identity.** A crash log's stated
   "frame N" refers to whatever position N was in the actual trajectory
   fed to `sander` - if that trajectory was built from a scrambled file
   list (see the `run_GBSA.sh` frame-ordering bug below), "frame N" does
   **not** mean the raw file literally named `_N.pdb`. Confirm by
   comparing a real, always-present atom's coordinates (e.g. a mid-chain
   residue's CA) between the raw file you assume it is and the actual
   `MODEL N` block of the trajectory that was really scored - if they
   don't match bit-for-bit, you have the wrong file.
2. **Extract that exact `MODEL N` block** from the real scored trajectory
   PDB and compute every pairwise interatomic distance between atoms in
   *different, non-adjacent* residues (skip same-residue and
   sequence-neighbor pairs - those are real bonds and will always be the
   closest, drowning out the actual clash). A non-bonded contact under
   ~1.2 Å between heavy atoms (or under ~1.0 Å involving an H) is a real
   clash - VDW repulsion scales as 1/r¹², so this is enough to explain
   an overflowed energy on its own.
3. **Confirm which side is frozen** by comparing that same atom's
   coordinates across two different `MODEL` blocks (e.g. `MODEL 1` vs the
   bad one) in the same trajectory - a frozen/unmatched reference atom
   will be bit-for-bit identical across frames; a real, correctly-matched
   atom will have moved.
4. Cross-check against `build_trajectory_v2.py`'s own printed
   `kept-from-reference` count for that frame - if it's the same known
   list of leftovers as every other frame (e.g. always exactly 5, always
   the same two termini), the clash is from a known, frozen atom class,
   not a new bug.

## run_GBSA.sh gotchas (cluster-specific, easy to silently reintroduce)

- **Frame ordering**: `frame_files=(${frame_dir}/${trajID}_[0-9]*.pdb)` -
  a plain bash glob - sorts **lexicographically**, not numerically.
  Confirmed directly (created 600 dummy files named `_1.pdb`..`_600.pdb`
  and ran the actual glob): the first 5 matches were `_100,_101,_102,
  _103,_104`, not `_1,_2,_3,_4,_5`. Since this array is passed straight to
  `-frames` and written out as `MODEL 1, 2, ... N` in that exact order,
  the scrambling means `MODEL N` in the scored trajectory is **not** the
  N-th real simulation frame - it's whatever frame happened to land at
  lexicographic position N. This silently corrupts frame *identity* for
  any per-frame analysis, independent of whether the atom-matching itself
  is correct. Fixed by sorting on the extracted numeric suffix instead of
  the raw glob order (see the `mapfile`/`sort -n` block in the script) -
  if anyone ever "simplifies" that back to a plain glob, this comes back.
- **OpenMPI "not enough slots"** at higher `-nt`: `gmx_MMPBSA`/`unigbsa`
  internally calls `mpirun -np <nt> --use-hwthread-cpus --allow-run-as-root`
  with no CLI way to add `--oversubscribe`. SLURM's `--cpus-per-task=N`
  correctly reserves N real cores for the job, but OpenMPI's own slot
  detection (when `mpirun` isn't launched via `srun`) can key off
  `SLURM_NTASKS` instead - which is `1` here, since the "3 replicates"
  come from `--array=0-2`, not from `--ntasks` (a single array task is
  one bash script running things as plain subprocesses, never spawning
  multiple SLURM tasks itself). OpenMPI sees "1 task" = "1 slot" and
  refuses to launch 32 processes, even though the 32 cores are genuinely
  allocated and idle. Fixed with `export OMPI_MCA_rmaps_base_oversubscribe=1`
  (the env-var equivalent of `--oversubscribe`) before the pipeline runs -
  this doesn't oversubscribe onto fewer real cores than requested, it
  just stops OpenMPI vetoing a launch based on a slot-count heuristic
  that was looking at the wrong SLURM variable for this task/cpus-per-task
  combination.
- **`--ntasks` vs `--array`**: don't set `--ntasks=3` thinking it
  represents "3 replicates" - that's what `--array=0-2` already does.
  `--ntasks` should be `1` (one process per array element); setting it
  higher multiplies the *per-array-task* resource request instead
  (`--ntasks=3 --cpus-per-task=32` asks for 96 cores per replicate, not
  32), which is why this surfaced as tasks failing to schedule /
  needlessly spanning multiple nodes once `--cpus-per-task` was raised
  from 1.
- **`repo_root`**: derived from `SLURM_SUBMIT_DIR` (falls back to `$(pwd)`
  outside SLURM), not hardcoded. `sbatch` sets `SLURM_SUBMIT_DIR` to
  wherever it was invoked from, so as long as you `cd` into your own
  clone before `sbatch run_GBSA.sh` (as the README says), this is always
  correct with nothing to edit. Previously this was a literal path,
  meaning every user had to hand-edit it back after every pull - easy to
  forget and silently run against someone else's checkout. The shared
  `gbsa_env` path is intentionally a separate variable
  (`gbsa_env_path`), not derived from `repo_root`, since the env lives in
  Marta's checkout regardless of whose `repo_root` this is. Per-user
  overrides (e.g. your own `gbsa_env_path` if you build your own env)
  live in `local.env.sh` (gitignored, see `local.env.sh.example`),
  sourced automatically if present.

## Multi-fragment receptors (a chain isn't always one continuous chain)

`clean_receptor_v2.py`/`build_trajectory_v2.py`'s cap-splitting originally
assumed exactly one ACE/NME pair per chain, at the chain's first/last
residue. Confirmed false for real GPCRmd dyns 15, 25 and 57: `clean_receptor_v2.py`
raised `NME-cap atoms found on a non-last residue`. Traced to the raw PDB's
own segid column (VMD writes it, though the PDB parser doesn't use it for
matching): chain `P` there is actually two separate PSF segments (`P0`,
`P1`) sharing one PDB chain letter, with a real residue-numbering gap
(resseq 250 -> 270) in between - almost certainly ICL3 replaced by a
crystallization fusion partner (T4-lysozyme/BRIL-style construct) that was
stripped out before the MD run, leaving two independently-capped
fragments. The other 4 dyns tested (16, 26, 55, 58) have only one segment.

Fix: fragment boundaries are now detected purely from where the cap atoms
themselves appear (an ACE-style atom is only valid immediately after a
fragment just closed by an NME-style atom, or at a chain's very first
residue) instead of a fixed first/last-of-chain assumption - handles any
number of fragments, fused anywhere in a chain, with no further changes
needed. `clean_receptor_v2.py` also now numbers synthesized ACE/NME
residues by reusing the adjacent real residue's resSeq with a fresh
insertion code (`250A`) instead of `resseq±1` - the old arithmetic would
silently collide if a future fusion break has no numbering gap - and
inserts an explicit `TER` at every fragment boundary. Validated directly
with real `gmx pdb2gmx` (not just the Python logic): all 7 dyns build with
exit 0 and zero interactive prompts; the two-fragment dyns correctly log
`Identified residue ACE33/NME250` and `ACE270/NME358` as two separate
starting/ending termini, producing two separate `.itp` files.

## Ligand net charge

`unigbsa`'s own net-charge auto-detection (`obtain_net_charge()`, which
round-trips the ligand file through OpenBabel) only recovers a nonzero
charge from an explicit SYBYL atom type (e.g. `N.4`) - confirmed it
silently returns 0 for a plain mol2 with generic element-symbol atom types,
which is exactly what `fix_ligand_naming.py` writes. A wrong net charge
doesn't fail loudly here - it hands `antechamber`/`sqm` an impossible odd
electron count for whatever the true charged/neutral split should be, and
`sqm` just crashes. Confirmed on real data: GPCRmd ligand Y00 (dyn 15/16)
has an SDF `M  CHG` record giving its amine nitrogen a +1 formal charge;
with that lost, `unigbsa` computed net charge 0, giving `sqm` 163
electrons (odd) for a claimed charge-0 molecule. All 7 ligands tested
across dyns 15/16/25/26/55/57/58 turned out to carry a +1 charge (all
protonated amines) - this isn't a rare edge case, it's most/all real
GPCRmd ligands.

Fix: `fix_ligand_naming.py` parses the SDF's own `M CHG` records directly
(authoritative, no guessing) and writes the summed net charge to a
`<ligand>.net_charge` sidecar. `build_topology_once.py` reads it and
monkeypatches `unigbsa.simulation.topology.obtain_net_charge` to return
that exact value, bypassing the lossy OpenBabel round-trip entirely.
Validated with real `antechamber`/`sqm`/`acpype`: all 7 dyns build their
full receptor+ligand topology with exit 0.

## Ligand atom identity (SDF naming vs the real MD frame's naming)

`fix_ligand_naming.py` used to number atoms sequentially by
element-appearance-order in the SDF (O1, C1, C2, ...) and trust that this
happened to match the raw MD frame's own atom naming - flagged as an
unverified assumption in this file's docstring from the start. Confirmed
false on every one of the 7 real GPCRmd ligands tested: the raw frame's
naming (CHARMM/PSF-style, e.g. hydrogens named by parent atom like
`H12X`/`H12Y`, or a PDB-CCD-style scheme like `CAI`/`HAIX`/`HAIY`/`HAIZ`)
has no relation to the SDF's naive ordering. `build_trajectory_v2.py`'s
by-name ligand matching silently paired up whichever atoms happened to
share a name and froze the rest at the reference build's pose - worse
than a crash, since it doesn't fail loudly.

Verified these are genuinely the same molecule, not a wrong-file problem,
before trusting a fix: element counts match exactly between the raw frame
and the SDF for both dyn 15 (Y00: C18 N1 O3 H24) and dyn 25 (5FW: C11 N1
O3 H18).

Fix, in `fix_ligand_naming.py`: try direct name matching against a raw
frame first (`-frame-pdb`); only if that isn't a full match, fall back to
real bond-graph matching - build the ligand's bond graph two ways (the
SDF's own authoritative bonds, and a geometry-perceived graph from the raw
frame's real 3D coordinates via RDKit's `rdDetermineBonds` - distance-based,
no OpenBabel) and graph-match them (element + connectivity only, ignoring
bond order/charge, which don't affect topological identity) via RDKit's
substructure matching. Raises rather than returning a partial match - a
wrong guess here is worse than a loud failure.

This mapping can't just be written into the output mol2's atom names,
though: confirmed `acpype`/`antechamber` ignores whatever name is written
in the input mol2 entirely and re-numbers every atom generically by
element + file order in its own output topology (checked directly - an
input mol2 written with real frame names and one written with plain
generic names produced byte-identical ligand atom ordering in the built
reference `complex.pdb`). So the mol2 body keeps the generic name scheme
acpype already reproduces reliably, and the real correspondence instead
goes into a `<ligand>.atom_map` sidecar (generic name -> real frame name).
`build_trajectory_v2.py -ligand-atom-map` translates through it before
the by-name lookup. All 7 dyns needed the graph-matching fallback (none
had directly-coinciding names) and all 7 resolved 100% of atoms - a sanity
check that held up independently: dyn 15's one nitrogen happens to be
named "N1" in both the SDF and the raw frame, and the graph match landed
it there on its own. Validated end-to-end: `build_trajectory_v2.py` now
reports `0 kept-from-reference` for both protein and ligand on all 7 dyns.
