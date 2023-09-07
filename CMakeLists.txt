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
  DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
