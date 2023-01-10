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
