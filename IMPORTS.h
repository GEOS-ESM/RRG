! 3-D
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'T',                                   &
       LONG_NAME  = 'air_temperature',                     &
       UNITS      = 'K',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'Q',                                   &
       LONG_NAME  = 'specific_humidity',                   &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)  

!!    call MAPL_AddImportSpec(GC,                            &
!!       SHORT_NAME = 'QTOT',                                &
!!       LONG_NAME  = 'mass_fraction_of_all_water',          &
!!       UNITS      = 'kg kg-1',                             &
!!       DIMS       = MAPL_DimsHorzVert,                     &
!!       VLOCATION  = MAPL_VLocationCenter,                  &
!!       RESTART    = MAPL_RestartSkip,     __RC__)  

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'QCTOT',                               &
       LONG_NAME  = 'mass_fraction_of_total_cloud_water',  &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter, __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'ZLE',                                 &
       LONG_NAME  = 'geopotential_height',                 &
       UNITS      = 'm',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'DELP',                                &
       LONG_NAME  = 'pressure_thickness',                  &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'PLE',                                 &
       LONG_NAME  = 'pressure_at_edges',                   &
       UNITS      = 'Pa',                                  &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationEdge,                    &
       RESTART    = MAPL_RestartSkip,     __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'AIRDENS',                             &
       LONG_NAME  = 'air_density',                         &
       UNITS      = 'kg m-3',                              &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'O3',                                  &
       LONG_NAME  = 'ozone',                               &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

    !  Add gas imports for CH4 & CO chemistry
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RRG_OH',                              &
       LONG_NAME  = 'hydroxyl',                            &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RRG_Cl',                              &
       LONG_NAME  = 'atomic Cl',                           &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'RRG_O1D',                             &
       LONG_NAME  = 'singlet O',                           &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

! 2-D
     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'AREA',                               &
        LONG_NAME  = '',                                   &
        UNITS      = 'm^2',                                &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'ZPBL',                               &
        LONG_NAME  = 'Planetary boundary layer height',    &
        UNITS      = 'm',                                  &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,    __RC__)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME = 'PS',                                 &
        LONG_NAME  = 'surface_pressure',                   &
        UNITS      = 'Pa',                                 &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)

    call MAPL_AddImportSpec(GC,                            &
        SHORT_NAME = 'U10M',                               &
        LONG_NAME  = '10-meter_eastward_wind',             &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)
    
    call MAPL_AddImportSpec(GC,                            &
        SHORT_NAME = 'V10M',                               &
        LONG_NAME  = '10-meter_northward_wind',            &
        UNITS      = 'm s-1',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
        RESTART    = MAPL_RestartSkip,   __RC__)
