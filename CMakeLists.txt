esma_set_this ()

set (srcs
    ${this}Mod.F90
    types_mod.F90
    surface_mod.F90
    chemistry_mod.F90 
    global_mod.F90
    utils_mod.F90
    CO_mod.F90
    )
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
