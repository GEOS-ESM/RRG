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
    UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzVert,                    &
        VLOCATION  = MAPL_VLocationCenter,                 &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CO2_EM',                             &
        LONG_NAME  = '',                                   &
    UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'CH4_EM',                             &
        LONG_NAME  = '',                                   &
    UNITS      = 'kg/kg/s',                            &
        DIMS       = MAPL_DimsHorzOnly,                    &
        VLOCATION  = MAPL_VLocationNone,                   &
    __RC__)

    call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME = 'TR_EM',                             &
        LONG_NAME  = '',                                   &
    UNITS      = 'kg/kg/s',                            &
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
