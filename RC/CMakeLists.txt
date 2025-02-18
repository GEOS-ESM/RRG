set (resource_files
  RRG_GridComp.rc
  RRG_GridComp_ExtData.rc
  RRG_GridComp_ExtData.yaml
  )

install(
  FILES ${resource_files}
  DESTINATION etc
  )
