esma_set_this ()

set (srcs
    ${this}Mod.F90
    CH4chem_mod.F90
    COchem_mod.F90
    CO2chem_mod.F90
    diagnostics_mod.F90
    geos_SimplePhotolysisMod.F90
    global_mod.F90
    integration_mod.F90
    surface_mod.F90
    types_mod.F90
    utils_mod.F90
    )
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
