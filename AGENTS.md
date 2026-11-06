# AGENTS.md - RRG

RRG is the carbon-gas (CO, CO2, CH4, passive tracer `TR`) production-and-loss
gridded component of GEOS Chemistry. Module `RRG_GridCompMod`, file
`RRGMod.F90`.

**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `RRG` @ `786f7a1` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).

**Build status on this line:** not built, not ported. On `release/MAPL-v3`,
`@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)`
only, with the comment "Restore subdirs as components are ported to MAPL3.";
`RRG` appears only in the commented-out `alldirs` list, marked "not yet
ported to MAPL3 - restore when ported". The Chemistry top-level library
(`GEOS_ChemGridComp.F90`, which adds this component) is present but not
compiled on this line: the `esma_add_library` block that builds it is
commented out in `../CMakeLists.txt`. RRG does not depend on `Chem_Base` or
`Chem_Shared` (`CMakeLists.txt`). The code in this repository is an older
version of the MAPL2-line code (25 commits behind `main`), kept for a later
port; see "Differences from `main`" in `docs/ARCHITECTURE.md`.

## Placement
GEOSgcm fixture (`components.yaml` @ `f50622d`):
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@RRG`.
Parent: the Chemistry component (`../docs/ARCHITECTURE.md`, `../AGENTS.md`;
on this line the parent's documents are not written yet; the MAPL2-line
versions are on its `develop` branch), which adds RRG when `ENABLE_RRG:` is
true.

## Documents
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md): placement, file map and the
  nine component sections.
- [`docs/decisions.md`](docs/decisions.md): previous decisions.
- `DEVEL.org`: development notes (unverified).
- RRG has no `CHANGELOG.md`.

## Files
- `RRGMod.F90`: SetServices, Initialize, two Run phases, Finalize, instance
  and flux-table readers.
- `IMPORTS.h`, `EXPORTS.h`: fixed spec declarations (included into SetServices).
- `COchemMod.F90`, `CO2chemMod.F90`, `CH4chemMod.F90`: production/loss.
- `geos_SimplePhotolysisMod.F90`: photolysis rates.
- `surfaceMod.F90` (symlink to `posDef.surfaceMod.F90`; `std.` not built;
  `surfaceMod.F90X` is also present and not built),
  `integrationMod.F90`, `diagnosticsMod.F90`, `typesMod.F90`, `globalMod.F90`,
  `utilsMod.F90`.
- `RC/`: `RRG_GridComp.rc` (instances, flux tables, masks),
  `RRG_GridComp_ExtData.rc` and `.yaml`.

## Rules for agents
- Note bugs and oddities in `docs/ARCHITECTURE.md` section 9; do not fix them
  unless asked.
- Keep the documents current when code changes; cite files and routines,
  not line numbers; label inferences.
- The field set depends on `RC/RRG_GridComp.rc` at SetServices time; check it
  before reasoning about imports.

## Build and test
`CMakeLists.txt` defines library `RRG`, but on this line nothing adds this
directory to the build (see build status above). No stand-alone build or tests
here.
