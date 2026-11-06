# RRG decisions

**Audience:** coding agents and developers working on RRG.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `RRG` @ `a0779ba` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).

### Species fields defined by the resource file
- Pattern: the internal state, imports and flux registrations are created in
  `SetServices` from `RRG_GridComp.rc` instance lists and flux tables
  (`RRGMod.F90`, SetServices, ProcessInstances, RegisterFluxWithMAPL).
- Likely reason: *Inferred:* add or remove emission-source tags without
  changing code (module comment: "users can add new species ... without
  having to add new code other than here").
- Question for maintainers: is the SetServices-time dependence on a private
  resource file intended to remain?

### Instances as separate advected fields plus an aggregate
- Pattern: each instance is its own friendly tracer; an aggregate field sums them
  (`RRGMod.F90`, RegisterInstanceWithMAPL, ProcessInstances; `utilsMod.F90`,
  `util_aggregate`).
- Likely reason: *Inferred:* track contribution by source/region.

### Two run phases
- Pattern: phase 1 applies surface fluxes and integrates; phase 2 computes
  photolysis, chemistry loss/production and integrates again
  (`RRGMod.F90`, GridCompRun1, GridCompRun2).
- Likely reason: *Inferred:* surface exchange and chemistry are placed at
  different points of the parent's sequence.

### Henry's law constants zero
- Pattern: all instance fields get `SetofHenryLawCts = 0` (`RRGMod.F90`, Initialize).
- Likely reason: the code comment says no scavenging should happen in Moist.

### Two surface-flux implementations
- Pattern: `surfaceMod.F90` is a symbolic link to `posDef.surfaceMod.F90`;
  `std.surfaceMod.F90` is kept but not built (`CMakeLists.txt`).
- Likely reason: *Inferred:* alternative scheme retained for comparison.
- Question for maintainers: is `std.surfaceMod.F90` still needed?
