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
   match the names used in your MD system's ligand, and `build_trajectory.py`
   matches ligand atoms by name. If you go this route, verify the match
   count printed by `build_trajectory.py` isn't 0/low - if names don't
   line up, you need a one-time atom-name remap (e.g. via a graph/substructure
   match) before this works, otherwise every frame silently reuses the
   PubChem structure's coordinates for the ligand instead of your actual
   simulated pose. 3D coordinates in the SDF don't matter either way -
   only the 2D bond graph and formal charges do, since coordinates get
   overwritten by each frame's real pose downstream.
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

## Matching assumption in build_trajectory.py

Protein atoms are matched frame-to-frame by `(chain, resnum, atom name)`;
the ligand is matched by atom name alone (there's one copy). This assumes
residue numbering is stable across your frames and against what
`pdb2gmx` produced in `topology/complex.pdb` (it generally preserves input
numbering, but always check the printed match counts before trusting a run
at scale - see step 4 in `README.md`).
