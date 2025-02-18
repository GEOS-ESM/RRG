esma_set_this ()

set (srcs
    ${this}Mod.F90
    CH4chemMod.F90
    COchemMod.F90
    CO2chemMod.F90
    diagnosticsMod.F90
    geos_SimplePhotolysisMod.F90
    globalMod.F90
    integrationMod.F90
    surfaceMod.F90
    typesMod.F90
    utilsMod.F90
    )
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL ESMF::ESMF NetCDF::NetCDF_Fortran)

add_subdirectory (RC)
