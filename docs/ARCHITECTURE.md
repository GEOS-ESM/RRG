# RRG architecture

**Audience:** coding agents and developers working on RRG.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `RRG` @ `786f7a1` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

**Build status on this line:** not built, not ported. On `release/MAPL-v3`,
`@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)`
only, with the comment "Restore subdirs as components are ported to MAPL3.";
`RRG` appears only in the commented-out `alldirs` list, marked "not yet
ported to MAPL3 - restore when ported". The Chemistry top-level library
(`GEOS_ChemGridComp.F90`, which adds this component) is present but not
compiled on this line: the `esma_add_library` block that builds it is
commented out in `../CMakeLists.txt`. RRG's own dependencies are MAPL, ESMF
and NetCDF Fortran (`CMakeLists.txt`); it does not use `Chem_Base` or
`Chem_Shared`. The code in this repository is an older version of the MAPL2-line code
(25 commits behind `main`; `MAPL_Add*Spec` specs, `MAPL_AddChild`), kept for a later port. Statements about the parent
describe `GEOS_ChemGridComp.F90` as it is on this line; the file is present
but not compiled.

**Differences from `main`.** `release/MAPL-v3` is 25 commits behind `main`
(0 ahead). This document describes the older code. Files that differ on
`main` (`git diff origin/release/MAPL-v3 origin/main`): `RRGMod.F90`
(`IMPORTS.h` and `EXPORTS.h` are included after `ProcessInstances`; oxidant
and aggregate pointers guarded; a `<species>.mw:` key), `IMPORTS.h` (oxidant
imports conditional), `EXPORTS.h` (units of some production/loss exports: `kg/kg/s` here; four are `kg/m3/s` on `main`), `RC/RRG_GridComp.rc`
(`TR.mw:`). `main` also has `.github/CODEOWNERS`, a docs workflow and Sphinx
sources under `docs/source/`; none is on this line. `surfaceMod.F90X` exists
on this line and not on `main`.

Paths below are relative to the RRG repository root. The parent is the
Chemistry gridded component; see
`@GEOSchem_GridComp/docs/ARCHITECTURE.md` (`../../docs/ARCHITECTURE.md` in the fixture)
(on this line the parent's documents are not written yet; the MAPL2-line
versions are on its `develop` branch) for how it adds and runs RRG, the
provider mechanism and the `ENABLE_*` flags. That document is not repeated
here.

## 1. Purpose

RRG is a small MAPL gridded component for the atmospheric part of the
carbon-gas cycle. For each species `C` it solves `dC/dt = P - L` with a
forward-Euler step (`RRGMod.F90`, module header; `integrationMod.F90`,
`integrate_forwardeuler`). Species are CO, CO2, CH4 and a passive tracer
`TR`. Each species is split into user-named *instances* (for example
emission-source tags) that are separate advected fields and are summed into
an aggregate (`RRGMod.F90`, `ProcessInstances`).

## 2. Place in the hierarchy

- Directory `@RRG` in the Chemistry repository checkout. Component name
  `RRG`, module `RRG_GridCompMod`, file `RRGMod.F90`, public entry
  `SetServices` (`RRGMod.F90`; the parent adds it in `../GEOS_ChemGridComp.F90`).
- Parent: `GEOSchem_GridComp`; it adds RRG with `MAPL_AddChild` only when
  `ENABLE_RRG:` is true (`../GEOS_ChemGridComp.F90`, SetServices).
- Children: none.
- Parent connections to RRG: `DELP`, `AIRDENS` from `CHEMENV`; `O3` from
  `PCHEM` when `PCHEM` is enabled (`../GEOS_ChemGridComp.F90`,
  SetServices, `MAPL_AddConnectivity` with `DST_ID = RRG`).
- Parent exports: `CO2` is re-exported from RRG when `CO2_PROVIDER: RRG`;
  names listed in `RATS_DIAGNOSTICS:` are exported with `CHILD_ID = RRG`
  (`../GEOS_ChemGridComp.F90`, SetServices).

## 3. Entry points

All set in `SetServices` with `MAPL_GridCompSetEntryPoint`, then
`MAPL_GenericSetServices` (`RRGMod.F90`, SetServices).

- `SetServices`: includes `IMPORTS.h` and `EXPORTS.h` first (on this line
  they come before the resource file is read; on `main` after
  `ProcessInstances`); opens `RRG_GridComp.rc` as a private `ESMF_Config`
  (`Initialize` loads the same file again, separately, for the tables); reads `strictMassBalance:`, `WellMixedSurfaceExchange:`,
  `UseResidual:`; calls `ProcessInstances` for CO, CO2, CH4, TR (adds
  internal specs from the instance lists); calls `RegisterFluxWithMAPL` for each species (adds one
  import per flux-table entry). Because the specs depend on the resource
  file, the field set is decided at SetServices time.
- `Initialize`: reads `RUN_DT:` and `GOCART_DT:`, calls
  `MAPL_GenericInitialize`, then `ReadFluxTable` and `ReadMasksTable` per
  species, sets `SetofHenryLawCts` to zeros on every instance field (so
  Moist does not scavenge them), and loads the photolysis tables from
  `photolysisFile:` (`readPhotTables`) (`RRGMod.F90`, Initialize).
- Run phase 1, `GridCompRun1`: aggregates instances, fills flux pointers
  (`fillFluxes`), applies surface fluxes (`surface_prodloss`), integrates
  (`integrate_forwardeuler`), optionally `util_accumulatenegatives`,
  aggregates again, and sets the emission exports `CO2_EM`, `CH4_EM`,
  `CO_EM`, `TR_EM`. On its first call it runs `ProcessExtdataMasks`
  (`RRGMod.F90`, GridCompRun1).
- Run phase 2, `GridCompRun2`: computes solar insolation
  (`MAPL_SunGetInsolation`), O2 and O3 columns and photolysis rates for CO2
  and CH4, then calls `CO_prodloss`, `CO2_prodloss`, `CH4_prodloss`
  (each when its species has instances), integrates, and sets diagnostics
  including the exports `CO2DRY`, `CH4DRY`, `CODRY` (`RRGMod.F90`, GridCompRun2).
- `Finalize`: deallocates instance and flux arrays, then
  `MAPL_GenericFinalize` (`RRGMod.F90`, Finalize).
- Which parent phase calls which RRG phase: the parent registers RRG with two
  phases (the parent's `SetServices`, `../GEOS_ChemGridComp.F90`); *Inferred:* phase 1 runs in
  the parent's Run1 and phase 2 in Run2, by registration order.

## 4. State

- Imports: `IMPORTS.h` holds the fixed specs: `T`, `Q`, `QCTOT`, `ZLE`,
  `DELP`, `PLE`, `AIRDENS`, `O3`, `AREA`, `ZPBL`, `PS`, `U10M`, `V10M`, and
  `RRG_OH`, `RRG_Cl`, `RRG_O1D` (unconditional on this line; on `main` only if
  CO or CH4 instances exist).
  On this line the aggregate `MAPL_GetPointer` calls in `GridCompRun1` and
  `GridCompRun2` are unconditional (on `main` they are guarded by
  `if (ispecies('X') .gt. 0)`).
  The `QTOT` import spec is commented out in `IMPORTS.h`, and its pointer line is commented in `GridCompRun2`.
  `CO_CH4` is added in `SetServices` when CO exists and CH4 does not.
  Per-flux imports come from `RegisterFluxWithMAPL`; `ProcessInstances`
  also adds an import when a `_masks::` table exists. `namerica_mask` is
  declared in `EXPORTS.h`.
- Sources of imports: `DELP`, `AIRDENS`, `O3` by connection from the
  parent (above); the oxidants `RRG_OH`, `RRG_Cl`, `RRG_O1D`, emissions
  and region masks by ExtData (`RC/RRG_GridComp_ExtData.yaml`;
  `RC/RRG_GridComp_ExtData.rc`). *Inferred:* the remaining meteorological
  imports are satisfied by the Chemistry parent or the AGCM export-import
  wiring, not by RRG files; not traced.
- Internal state: for each species, one aggregate field (`CO`, `CO2`,
  `CH4`, `TR`; `MAPL_RestartOptional`) and one field per instance named
  `<species>_<instance>` (units mol mol-1, friendly to `DYNAMICS:TURBULENCE:MOIST`
  unless the instance is passive, `add2export`). The per-instance restart
  attribute is commented out in `RegisterInstanceWithMAPL`
  (`RRGMod.F90`, ProcessInstances, RegisterInstanceWithMAPL).
- Exports (`EXPORTS.h`): `CO2DRY`, `CH4DRY`, `CODRY`, `CO2photj`,
  `CH4photj`, `O3col`, `O2col`, `CO2_ProdLoss`, `CO_Prod`, `CO_Loss`,
  `CH4_ProdLoss`, `CO2_EM`, `CH4_EM`, `TR_EM`, `CO_EM`, `CO2FDNL`,
  `namerica_mask`. The instance and aggregate fields are also exported.
- Module state: `RRGMod.F90` holds the instance arrays `CO`, `CO2`, `CH4`,
  `TR` as module variables of type `gas_instance`
  (`typesMod.F90`); `globalMod.F90` holds shared module-level data;
  `utilsMod.F90` keeps the registries. *Inferred:* one RRG instance per
  process is assumed because of these module variables.
- Restart: no restart-required instance fields; aggregates are
  restart-optional. Parent `AGCM.rc.tmpl` keeps RRG restart keys commented
  (*Unverified:* parent's restart configuration; not re-checked on this line).

## 5. Configuration that changes structure

From `RC/RRG_GridComp.rc` (read by `SetServices` and `Initialize`):

| Key | Effect |
|---|---|
| `CH4_instances:`, `CO_instances:`, `CO2_instances:` (and `TR_instances:`) | Instance names; define the internal fields. |
| `*_passive_instances:` | Instances that are not advected, mixed or convected. |
| `*_surface_flux_pairs::` | Instance to import-name table; options `d` (diurnal scaling), `p` (PBL mixing), a number (scale factor). |
| `*_masks::` | Instance region masks (integer or lat/lon box). |
| `WellMixedSurfaceExchange:`, `UseResidual:`, `strictMassBalance:` | Control switches in `cntrl`. |
| `photolysisFile:` | Photolysis table; read in `Initialize`. |

A `<species>.mw:` molecular-weight key (and `TR.mw:` in the rc file) exists on
`main` only; on this line `ProcessInstances` passes the built-in weights.

Parent keys: `ENABLE_RRG:`, `CO2_PROVIDER:`, `RATS_DIAGNOSTICS:`
(`../GEOS_ChemGridComp.F90`). The comments in `RC/RRG_GridComp.rc` call H2O, CH4 and CO2
feedbacks "not implemented yet".

## 6. Dependencies

- Libraries: MAPL, ESMF, NetCDF Fortran (`CMakeLists.txt`).
- Source modules (all in `CMakeLists.txt`, library name `RRG`):
  `RRGMod.F90` (the component); `COchemMod.F90`, `CO2chemMod.F90`,
  `CH4chemMod.F90` (`CO_prodloss`, `CO2_prodloss`, `CH4_prodloss`);
  `geos_SimplePhotolysisMod.F90` (photolysis rate tables and
  `CO2_photolysis_rate`, `CH4_photolysis_rate`); `surfaceMod.F90` (surface
  flux term, `surface_prodloss`); `integrationMod.F90`; `diagnosticsMod.F90`
  (`diag_prodloss`, `diag_sfcflux`); `typesMod.F90`; `globalMod.F90`;
  `utilsMod.F90` (`util_aggregate`, `util_addinstance`, ...).
- `surfaceMod.F90` is a symbolic link to `posDef.surfaceMod.F90`; the build
  compiles that name. `std.surfaceMod.F90` is not built; `surfaceMod.F90X`
  is also present and not built. The headers
  differ: the `posDef` version states a goal of no branching inside the
  nested loops and splits positive (source) and negative (loss) flux
  handling; the `std` version tests `inst%active` in the loop
  (`CMakeLists.txt`; `diff` of the three files).
- Installs (`RC/CMakeLists.txt`): `RRG_GridComp.rc`,
  `RRG_GridComp_ExtData.rc`, `RRG_GridComp_ExtData.yaml` to `etc`.
- Upstream: *Inferred:* photolysis code derives from the GOCART CH4
  component (comment in `GridCompRun2`).

## 7. Known traps

- The field set depends on `RRG_GridComp.rc` at SetServices time: adding an
  instance changes the internal state, the imports and (through ExtData) the
  required input files. An ExtData entry is needed for each flux import.
- Instance fields are friendly tracers; their transport is set by the
  `friendlyto` string at spec time, not by the Chemistry parent.
- `CO2_PROVIDER: RRG` exports `CO2` from RRG's aggregate; the aggregate
  exists only if `CO2_instances:` is non-empty (`ProcessInstances`).
- Surface-flux units and scaling: sign of flux matters (positive adds to
  `prod`, negative to `loss`); `DEVEL.org` says negative fluxes yield
  different outcomes (*Unverified:* development notes).
- Species-off cases: on this line the oxidant imports (`RRG_OH`, `RRG_Cl`,
  `RRG_O1D`) are declared and read unconditionally, also when no CO or CH4
  instance exists (on `main` they are conditional); CO without CH4 uses the
  `CO_CH4` import (`SetServices`).
- Module-level `CO/CO2/CH4/TR` pointers: two RRG instances in one executable
  would share them (*Inferred:*).

## 8. Decisions

See [`decisions.md`](decisions.md).

## 9. Observations (noted, not fixed) and open questions

1. The O2/O3 column code in `GridCompRun2` carries an author comment that
   the origin of the `0.50` and `1.e-4` factors is not understood
   (`RRGMod.F90`, GridCompRun2). The CO2 photolysis call passes a reversed
   level index `params%km-k+1`; a neutral question whether that is intended.
2. `GridCompRun2` fills `CO2DRY`, `CH4DRY`, `CODRY` (using `met%q` and
   `met%qctot`); `GridCompRun1` holds a commented-out earlier version of
   the same block that uses `met%qtot`. `QTOT` is not imported: its spec
   in `IMPORTS.h` and its pointer line in `GridCompRun2` are commented
   out (`RRGMod.F90`, GridCompRun1, GridCompRun2).
3. `RRG_GridComp.rc` `photolysisFile:` is an absolute path on a specific
   file system (`RC/RRG_GridComp.rc`); `RC/RRG_GridComp_ExtData.rc` also
   points at `deprecated` directories.
4. `RRGMod.F90` `Finalize` header says "simple ESMF wrapper" although it
   deallocates module data; `DEVEL.org` lists open tasks (mass balance,
   run with a species off, water-vapor feedback) (*Unverified:*).
5. `RC/RRG_GridComp_ExtData.rc` (classic) and `.yaml` are both installed;
   which one is read depends on the parent's ExtData configuration.
   Question for maintainers: is the `.rc` kept only for legacy runs?
6. `LONG_NAME` of the `CO_CH4` import is `source species`, unit `1`, while
   the CH4 aggregate is mol mol-1.
7. The parent's code (`../GEOS_ChemGridComp.F90`) agrees with the module
   name and the two-phase claim; no contradiction found.
