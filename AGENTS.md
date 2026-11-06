# AGENTS.md - RRG

RRG is the carbon-gas (CO, CO2, CH4, passive tracer `TR`) production-and-loss
gridded component of GEOS Chemistry (MAPL2 line, branch `main`; verified
against `RRG` @ `a0779ba`). Module `RRG_GridCompMod`, file `RRGMod.F90`.

## Placement
GEOSgcm fixture (`components.yaml` @ `0f2e0fb`):
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@RRG`.
Parent: the Chemistry component (`../docs/ARCHITECTURE.md`, `../AGENTS.md`),
which adds RRG when `ENABLE_RRG:` is true.

## Documents
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md): placement, file map and the
  nine component sections.
- [`docs/decisions.md`](docs/decisions.md): previous decisions.
- Existing Sphinx sources (not edited by the above), index
  `docs/source/index.rst`: `docs/source/checkout.md` (checking out),
  `docs/source/build.md` (building), `docs/source/model_setup.md` (setting up
  a run).
- `DEVEL.org`: development notes (unverified).
- RRG has no `CHANGELOG.md`.

## Files
- `RRGMod.F90`: SetServices, Initialize, two Run phases, Finalize, instance
  and flux-table readers.
- `IMPORTS.h`, `EXPORTS.h`: fixed spec declarations (included into SetServices).
- `COchemMod.F90`, `CO2chemMod.F90`, `CH4chemMod.F90`: production/loss.
- `geos_SimplePhotolysisMod.F90`: photolysis rates.
- `surfaceMod.F90` (symlink to `posDef.surfaceMod.F90`; `std.` not built),
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
Built as library `RRG` by `CMakeLists.txt` inside the GEOSgcm build; no
stand-alone build or tests here.
