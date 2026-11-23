    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'namerica_mask',                      &
        LONG_NAME  = '',            &
        UNITS      = '',                              &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2DRY',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4DRY',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CODRY',                              &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2photj',                           &
        LONG_NAME  = '',                                   &
        UNITS      = '',                                   &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4photj',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'O3col',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'O2col',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'mol mol-1',                          &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2_ProdLoss',                       &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO_Prod',                            &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO_Loss',                            &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4_ProdLoss',                       &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m3/s',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2_EM',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m3/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4_EM',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m3/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'TR_EM',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m3/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO_EM',                             &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2FDNL',                            &
        LONG_NAME  = '',                                   &
        UNITS      = '',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    do n = 1,NINSTANCES
       ! Loop over instances
       ! Get the species string
       name = trim(instances(n)%p%name)//'_COL'
       call MAPL_AddExportSpec(GC,                                                      &
            SHORT_NAME = trim(name),                                                    &
            LONG_NAME  = 'Column average mole fraction of '//trim(instances(n)%p%name), &
            UNITS      = 'mol mol-1 (dry air)',                                         &
            DIMS       = MAPL_DimsHorzOnly,                                             &
            VLOCATION  = MAPL_VLocationNone,                                            &
            __RC__)
    
       name = trim(instances(n)%p%name)//'_Prod'
       call MAPL_AddExportSpec(GC,                             &
            SHORT_NAME = trim(name),                           &
            LONG_NAME  = '',                                   &
            UNITS      = 'kg/m2/s',                            &
            DIMS       = MAPL_DimsHorzOnly,                    &
            VLOCATION  = MAPL_VLocationCenter,                 &
            __RC__)

       name = trim(instances(n)%p%name)//'_Loss'
       call MAPL_AddExportSpec(GC,                             &
            SHORT_NAME = trim(name),                           &
            LONG_NAME  = '',                                   &
            UNITS      = 'kg/m2/s',                            &
            DIMS       = MAPL_DimsHorzOnly,                    &
            VLOCATION  = MAPL_VLocationCenter,                 &
            __RC__)

    enddo

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CH4_COL',                                   &
        LONG_NAME  = 'Column average mole fraction of CH4_total', &
        UNITS      = 'mol mol-1 (dry air)',                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4_Prod',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4_Loss',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CO2_COL',                                   &
        LONG_NAME  = 'Column average mole fraction of CO2_total', &
        UNITS      = 'mol mol-1 (dry air)',                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2_Prod',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2_Loss',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CO_COL',                                   &
        LONG_NAME  = 'Column average mole fraction of CO_total', &
        UNITS      = 'mol mol-1 (dry air)',                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO_Prod',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO_Loss',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'TR_COL',                                   &
        LONG_NAME  = 'Column average mole fraction of TR_total', &
        UNITS      = 'mol mol-1 (dry air)',                       &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'TR_Prod',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'TR_Loss',                           &
        LONG_NAME  = '',                                   &
        UNITS      = 'kg/m2/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

