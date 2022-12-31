#include "MAPL_Generic.h"
!=============================================================================
!BOP

! !MODULE: GEOScarbon_GridCompMod

! !INTERFACE:
module GEOScarbon_GridCompMod
  
! !USES:
   use ESMF
   use MAPL
   use types_mod
   use global_mod
   use utils_mod
   use geos_simplephotolysismod
   
   implicit none
   private

   integer, parameter             :: DP=kind(1.0d0)
   integer                        :: status ! module-wide scope

   ! These are module-specific declarations that have a general interface
   ! with the gas module. Users can add new species to the system w/o
   ! having to add new code other than here in the main interface module.
   integer, save                  :: nCH4, nCO2, nCO ! number of instances
   type(gas_instance), pointer    :: CO2(:) => null()
   type(gas_instance), pointer    ::  CO(:) => null()
   type(gas_instance), pointer    :: CH4(:) => null()
   
! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

! !DESCRIPTION: This module implements the atmospheric component of the global
!               carbon cycle within the NASA GEOS modeling system using MAPL/ESMF
!
! NOTES:
! 1) TOTAL/AGGREGATE FIELDS: Totals of CO2, CO & CH4 are NOT independent instances.
!    They are used for diagnostics and coupling with other modules (e.g. grid comps)
!    But they are not directly modified by any process.   
!
! !REVISION HISTORY:
! 02Dec2022  M.S.Long   First pass

!EOP
!===========================================================================

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices ( GC, RC )

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!    DESCRIPTION: 

! !REVISION HISTORY:
! 02Dec2022  M.S.Long   First pass

!EOP
!============================================================================

!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (ESMF_Config)                   :: cfg
    type (ESMF_Config)                   :: universal_cfg

    character (len=255)                  :: name
    integer                              :: i,j,n
    logical                              :: present, found

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, config=universal_cfg, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

!==============================================================
!           S E T  T H E  I M P O R T  S T A T E

! 3-D
    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'T',                                   &
       LONG_NAME  = 'air_temperature',                     &
       UNITS      = 'K',                                   &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)

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

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'Q',                                   &
       LONG_NAME  = 'specific_humidity',                   &
       UNITS      = 'kg kg-1',                             &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       RESTART    = MAPL_RestartSkip,     __RC__)  

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
       SHORT_NAME = 'OH',                                  &
       LONG_NAME  = 'hydroxyl',                            &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'Cl',                                  &
       LONG_NAME  = 'atomic Cl',                           &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

    call MAPL_AddImportSpec(GC,                            &
       SHORT_NAME = 'O1D',                                 &
       LONG_NAME  = 'singlet O',                           &
       UNITS      = 'molecules cm-3',                      &
       DIMS       = MAPL_DimsHorzVert,                     &
       VLOCATION  = MAPL_VLocationCenter,                  &
       __RC__)

! 2-D
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

! ===============================================================
!      S E T  U P  T H E  E X P O R T  S T A T E
     ! Currently nothing
! ===============================================================

! ===============================================================
!  R E A D  C O N F I G  A N D  S E T U P  I N S T A N C E S
!
!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'GEOScarbon_GridComp.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'GEOScarbon_GridComp.rc does not exist!'
      VERIFY_(STATUS)
    end if

!   Read & set cntrl object
!   -----------------------
    call ESMF_ConfigFindLabel(cfg,label='strictMassBalance:',isPresent=present,rc=status)
    VERIFY_(STATUS) 
    if (present) then
       call ESMF_ConfigGetAttribute(cfg,cntrl%strictMassBalance,rc=status)
       VERIFY_(STATUS) 
    endif

!    Since this is all repeated for CO, CO2 & CH4, it could
!  be built into a routine or two. Depends on how clean we want
!  things

!  Carbon monoxide (CO)
!  CO: Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CO_instances:',rc=status)
   VERIFY_(STATUS)
   nCO = n

!  CO: define the total/aggregate field
   call MAPL_AddInternalSpec(gc,           &
        short_name ='CO_total',            &
        long_name  ='Aggregate CO field',  &
        units      ='kg kg-1',             &
        dims       =MAPL_DimsHorzVert,     &
        vlocation  =MAPL_VlocationCenter,  &
        restart    =MAPL_RestartOptional,  &
        add2export =.true., & !<-- is this what makes it available for HISTORY?
        __RC__)

!  CO: Get instances from RC file
!  ----------------------------------
   call ESMF_ConfigFindLabel(cfg,'CO_instances:',rc=status)
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
      name = TRIM(name)
      ! Register as tracer
      call RegisterInstanceWithMAPL( GC, 'CO', trim(name), rc=status )
      VERIFY_(STATUS)

      ! Add new active instance (assume active)
      call Util_AddInstance( CO, trim(name), 'CO', 28.0104, .true., status)
      VERIFY_(STATUS)
   end do

! If there are CO instances, then define an active residual by default
! ASSUMPTION: a species' residual is always the last instance 
   if (nCO .gt. 0) then
      call Util_AddInstance( CO, 'residual', 'CO', 28.0104, .true., status)
      VERIFY_(STATUS)
      call RegisterInstanceWithMAPL( GC, 'CO', 'residual', rc=status )
      VERIFY_(STATUS)
   endif

!  CO: Get passive instances and toggle them
   n = 0
   call ESMF_ConfigFindLabel(cfg,label='CO_passive_instances:',isPresent=present,rc=status)
   VERIFY_(STATUS)
   if (present) then
      n = ESMF_ConfigGetLen(cfg,label='CO_passive_instances:',rc=status)
      VERIFY_(STATUS)
      call ESMF_ConfigFindLabel(cfg,'CO_passive_instances:',rc=status)
      VERIFY_(STATUS)
      if (n .gt. 0) then
         do i = 1, n
            ! Get instance name
            call ESMF_ConfigGetAttribute(cfg,name,rc=status)
            VERIFY_(STATUS)

            found = .false.
            do j=1,size(CO)
            ! Search for name in list of instances
               if (index(CO(j)%name,trim(name)) .gt. 0) then
                  found = .true.
                  CO(j)%active = .false.
               endif
            enddo
            if (.not. found) then
               if (MAPL_am_I_root()) write(*,*) 'GEOScarbon::CO passive instance entry not in list'
               VERIFY_(-1)
            endif

         end do
      endif
   endif

!  Carbon dioxide (CO2)
!  CO2: define the total/aggregate field
   call MAPL_AddInternalSpec(gc,           &
        short_name ='CO2_total',           &
        long_name  ='Aggregate CO2 field', &
        units      ='kg kg-1',             &
        dims       =MAPL_DimsHorzVert,     &
       vlocation  =MAPL_VlocationCenter,  &
        restart    =MAPL_RestartOptional,  &
        add2export =.true., & !<-- is this what makes it available for HISTORY?
        __RC__)

!  CO2: Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CO2_instances:',rc=status)
   VERIFY_(STATUS)
   nCO2 = n

!  CO2: Get instances from RC file
!  ----------------------------------
   call ESMF_ConfigFindLabel(cfg,'CO2_instances:',rc=status)
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
      name = TRIM(name)
      ! Register as tracer
      call RegisterInstanceWithMAPL( GC, 'CO2', trim(name), rc=status )
      VERIFY_(STATUS)

      ! Add new active instance
      call Util_AddInstance( CO2, trim(name), 'CO2', 44.0, .true., status)
      VERIFY_(STATUS)
   end do
! If there are CO2 instances, then define an active residual by default   
! ASSUMPTION: a species' residual is always the last instance 
   if (nCO2 .gt. 0) then
      call Util_AddInstance( CO2, 'residual', 'CO2', 44.0, .true., status)
      VERIFY_(STATUS)
      call RegisterInstanceWithMAPL( GC, 'CO2', 'residual', rc=status )
      VERIFY_(STATUS)
   endif

!  CO2: Get passive instances and toggle them
   n = 0
   call ESMF_ConfigFindLabel(cfg,label='CO2_passive_instances:',isPresent=present,rc=status)
   VERIFY_(STATUS)
   if (present) then
      n = ESMF_ConfigGetLen(cfg,label='CO2_passive_instances:',rc=status)
      VERIFY_(STATUS)
      call ESMF_ConfigFindLabel(cfg,'CO2_passive_instances:',rc=status)
      VERIFY_(STATUS)
      if (n .gt. 0) then
         do i = 1, n
            ! Get instance name
            call ESMF_ConfigGetAttribute(cfg,name,rc=status)
            VERIFY_(STATUS)

            found = .false.
            do j=1,size(CO2)
            ! Search for name in list of instances
               if (index(CO2(j)%name,trim(name)) .gt. 0) then
                  found = .true.
                  CO2(j)%active = .false.
               endif
            enddo
            if (.not. found) then
               write(*,*) 'ERROR: GEOScarbon: CO2 passive instance entry not in list'
               VERIFY_(-1)
            endif

         end do
      endif
   endif

!  Methane (CH4)
!  CH4: define the total/aggregate field
   call MAPL_AddInternalSpec(gc,           &
        short_name ='CH4_total',           &
        long_name  ='Aggregate CH4 field', &
        units      ='kg kg-1',             &
        dims       =MAPL_DimsHorzVert,     &
        vlocation  =MAPL_VlocationCenter,  &
        restart    =MAPL_RestartOptional,  &
        add2export =.true., & !<-- is this what makes it available for HISTORY?
        __RC__)

!  CH4: Parse resource file
!  -------------------
   n = ESMF_ConfigGetLen(cfg,label='CH4_instances:',rc=status)
   VERIFY_(STATUS)
   nCH4 = n

!  CH4: Get instances from RC file
!  ----------------------------------
   call ESMF_ConfigFindLabel(cfg,'CH4_instances:',rc=status)
   VERIFY_(STATUS)
   do i = 1, n
      call ESMF_ConfigGetAttribute(cfg,name,rc=status)
      VERIFY_(STATUS)
      name = TRIM(name)
      ! Register as tracer
      call RegisterInstanceWithMAPL( GC, 'CH4', trim(name), rc=status )
      VERIFY_(STATUS)

      ! Add new active instance
      call Util_AddInstance( CH4, trim(name), 'CH4', 16.043, .true., status)
      VERIFY_(STATUS)
   end do
! If there are CH4 instances, then define an active residual by default   
! ASSUMPTION: a species' residual is always the last instance 
   if (nCH4 .gt. 0) then
      call Util_AddInstance( CH4, 'residual', 'CH4', 16.043, .true., status)
      VERIFY_(STATUS)
      call RegisterInstanceWithMAPL( GC, 'CH4', 'residual', rc=status )
      VERIFY_(STATUS)
   endif

!  CH4: Get passive instances and toggle them
   n = 0
   call ESMF_ConfigFindLabel(cfg,label='CH4_passive_instances:',isPresent=present,rc=status)
   VERIFY_(STATUS)
   if (present) then
      n = ESMF_ConfigGetLen(cfg,label='CH4_passive_instances:',rc=status)
      VERIFY_(STATUS)
      call ESMF_ConfigFindLabel(cfg,'CH4_passive_instances:',rc=status)
      VERIFY_(STATUS)
      if (n .gt. 0) then
         do i = 1, n
            ! Get instance name
            call ESMF_ConfigGetAttribute(cfg,name,rc=status)
            VERIFY_(STATUS)

            found = .false.
            do j=1,size(CH4)
            ! Search for name in list of instances
               if (index(CH4(j)%name,trim(name)) .gt. 0) then
                  found = .true.
                  CH4(j)%active = .false.
               endif
            enddo
            if (.not. found) then
               if (MAPL_am_I_root()) write(*,*) 'GEOScarbon::CH4 passive instance entry not in list'
               VERIFY_(-1)
            endif
         end do
      endif
   endif

! Fluxes
   call RegisterFluxWithMAPL( GC, cfg, 'CO' , RC )
   call RegisterFluxWithMAPL( GC, cfg, 'CO2', RC )
   call RegisterFluxWithMAPL( GC, cfg, 'CH4', RC )
! ===============================================================

!  Set entry points
!  ------------------------
   call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
   call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, GridCompRun, __RC__)
   
!   Store internal state in GC
!   --------------------------
!    call ESMF_UserCompSetInternalState ( GC, 'GEOScarbon_GridComp', wrap, STATUS )
!    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, IMPORT, EXPORT, CLOCK, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION:

! !REVISION HISTORY: 

!EOP
!============================================================================
!   !Locals 
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    character (len=255)                  :: name, photFile
    type (MAPL_MetaComp),      pointer   :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_Config)                   :: cfg, universal_cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type(ESMF_Field)                     :: field

!    real                                 :: CDT         ! chemistry timestep (secs)
!    real, pointer, dimension(:,:)        :: lats 
!    real, pointer, dimension(:,:)        :: lons
!    integer                              :: HDT         ! model     timestep (secs)
    integer                              :: i, n, dims(3)
    type (MAPL_VarSpec), pointer         :: INTERNALspec(:)  ! This is used to access GC information, e.g. field names, etc. (MSL)
    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=universal_cfg, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
!    call ESMF_UserCompGetInternalState(GC, 'GEOScarbon_GridComp', wrap, STATUS)
!    VERIFY_(STATUS)

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, localCellCountPerDim=dims, __RC__ )
    params%im = dims(1)
    params%jm = dims(2)
    params%km = dims(3)

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, params%HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, params%CDT, Label='GOCART_DT:', default=real(params%HDT), __RC__)

!   Set some quantities
!   ------------------- 
    params%AVO   = MAPL_AVOGAD
    params%AIRMW = MAPL_AIRMW

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'GEOScarbon_GridComp.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'GEOScarbon_GridComp.rc does not exist!'
      VERIFY_(STATUS)
    end if

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, &
                         LONS = params%LONS, &
                         LATS = params%LATS, &
                         INTERNALspec = INTERNALspec, &
                         __RC__ )

    if (nCO2 .gt. 0) call ReadFluxEntries( import, internal, cfg, 'CO2', __RC__ )
    if (nCO  .gt. 0) call ReadFluxEntries( import, internal, cfg, 'CO',  __RC__ )
    if (nCH4 .gt. 0) call ReadFluxEntries( import, internal, cfg, 'CH4', __RC__ )

!    if (MAPL_am_I_root()) call util_dumpinstances()

!  Set physical parameters
!  - Henry's Law constants. Only needed for GEOS
!    set to ZERO so that no scavenging happens in MOIST

    do i=1,nCO
       call ESMF_StateGet(internal, 'CO_'//trim(CO(i)%name), field, __RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', (/0,0,0,0/), __RC__)
    enddo
    do i=1,nCO2
       call ESMF_StateGet(internal, 'CO2_'//trim(CO2(i)%name), field, __RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', (/0,0,0,0/), __RC__)
    enddo
    do i=1,nCH4
       call ESMF_StateGet(internal, 'CH4_'//trim(CH4(i)%name), field, __RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', (/0,0,0,0/), __RC__)
    enddo

    call ESMF_ConfigFindLabel(cfg,label='photolysisFile:',rc=status)
    VERIFY_(STATUS) 
    call ESMF_ConfigGetAttribute(cfg,photfile,rc=status)
    VERIFY_(STATUS) 
    call readPhotTables( trim(photfile), params%km, status )
    VERIFY_(STATUS) 

!   Mask to prevent emissions from the Great Lakes and the Caspian Sea
!   ------------------------------------------------------------------
!    allocate(self%deep_lakes_mask(ubound(lons, 1),ubound(lons, 2)), __STAT__)
!    call deepLakesMask (lons, lats, real(MAPL_RADIANS_TO_DEGREES), self%deep_lakes_mask, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!============================================================================
!BOP
! !IROUTINE: GridCompRun

! !INTERFACE:
  subroutine GridCompRun (GC, import, export, clock, RC)
!   !USES:
    use surface_mod
    use global_mod
    use integration_mod

    ! Species chemistry modules
    use  COchem_mod
    use CO2chem_mod
    use CH4chem_mod

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: dC/dt = P - L
!               P & L are computed for all species and integrated simultaneously

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (MAPL_VarSpec), pointer      :: INTERNALspec(:)  ! This is used to access GC information, e.g. field names, etc. (MSL)

    integer i,j,k,n

    type(ESMF_Time)                   :: time
    integer                           :: iyr, imm, idd, ihr, imn, isc

    integer                           :: iCO2, iCH4, iCO
    real, pointer                     :: CO2_total(:,:,:), CH4_total(:,:,:), CO_total(:,:,:)
    logical :: first = .true.

    real, pointer, dimension(:,:,:) :: O3, OH, Cl, O1D 
    real(ESMF_KIND_R4), allocatable :: ZTH(:,:)
    real(ESMF_KIND_R4), allocatable :: SLR(:,:)
    type (MAPL_SunOrbit)            :: ORBIT
    type(ESMF_Alarm)                :: ALARM

    real                            :: r

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, INTERNALspec=INTERNALspec, __RC__)
    call MAPL_Get (mapl, ORBIT=ORBIT, RUNALARM=ALARM, __RC__)

!   Get current time
!   -----------------------------------
    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_PackTime(params%NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(params%NHMS,IHR,IMN,ISC)

! ===============================================================
!              G E T  D A T A  P O I N T E R S
!   Associate the met fields
!   -----------------------------------
    call MAPL_GetPointer(import,met%pblh, 'ZPBL', __RC__) ! pblh
    call MAPL_GetPointer(import,met%ps,   'PS',   __RC__) ! ps
    call MAPL_GetPointer(import,met%u10m, 'U10M', __RC__) ! u10
    call MAPL_GetPointer(import,met%v10m, 'V10M', __RC__) ! v10
    call MAPL_GetPointer(import,met%T,    'T',    __RC__) ! T
    call MAPL_GetPointer(import,met%zle,  'ZLE',  __RC__) ! zle
    call MAPL_GetPointer(import,met%ple,  'PLE',  __RC__) ! ple
    call MAPL_GetPointer(import,met%delp, 'DELP', __RC__) ! delp
    call MAPL_GetPointer(import,met%qctot,'QCTOT',__RC__) ! qtot
    call MAPL_GetPointer(import,met%rho,'AIRDENS',__RC__) ! rho
    CALL MAPL_GetPointer(import,     O3,    'O3', __RC__)
    CALL MAPL_GetPointer(import,     OH,    'OH', __RC__)
    CALL MAPL_GetPointer(import,     Cl,    'Cl', __RC__)
    CALL MAPL_GetPointer(import,    O1D,   'O1D', __RC__)

    allocate(  met%cosz(size(params%lats,1), size(params%lats,2)), __STAT__)
    allocate(  met%slr (size(params%lats,1), size(params%lats,2)), __STAT__)
    allocate(  met%o3col(params%im,params%jm,params%km), __STAT__)
    allocate(  met%photj, mold=met%o3col, __STAT__)

! Get the instance data pointers
    do i=1,NINSTANCES
       call MAPL_GetPointer(internal, instances(i)%p%data3d, trim(instances(i)%p%species)//'_'//trim(instances(i)%p%name), __RC__)
       allocate( instances(i)%p%prod, instances(i)%p%loss, mold=instances(i)%p%data3d, __STAT__ ) ! allocate the prod/loss arrays for each instance
       instances(i)%p%prod = 0.e0; instances(i)%p%loss = 0.e0
    enddo

!   Get pointers to the aggregates/totals
    call MAPL_GetPointer(internal, CO2_total, 'CO2_total',__RC__)
    aggregate(ispecies('CO2'))%q => CO2_total ! Aggregate is used under the hood

    call MAPL_GetPointer(internal,  CO_total,  'CO_total',__RC__) 
    aggregate(ispecies('CO'))%q  => CO_total  ! Aggregate is used under the hood

    call MAPL_GetPointer(internal, CH4_total, 'CH4_total',__RC__) 
    aggregate(ispecies('CH4'))%q => CH4_total ! Aggregate is used under the hood
! ===============================================================

! ===============================================================
!                S E T  U P  P H O T O L Y S I S
!  Update solar zenith angle
!  --------------------------
    call MAPL_SunGetInsolation(params%lons, params%lats, orbit, met%cosz, met%slr, clock=clock, __RC__)

!  Compute the O3 column
!  ---------------------
   r = 6.022e26*0.50/(28.97*9.8) ! r = Nsuba*0.50/(mwtAir*grav), copied from CFC_GridCompMod.F90
   met%O3col(:,:,1) = 1.1e15 + O3(:,:,1)*met%delp(:,:,1)*r
   DO k=2,params%km
      met%O3col(:,:,k) = met%O3col(:,:,k-1) + r * &
                        (O3(:,:,k-1) * met%delp(:,:,k-1) + O3(:,:,  k) * met%delp(:,:,  k))
   END DO
   
!  Compute the photolysis rate for CO2 + hv -> CO + O*
   met%photj = 0.e0
   do k=1,params%km
   do j=1,params%jm
   do i=1,params%im
      call jcalc4(i,j,params%km-k+1,met,met%photj(i,j,k))
   enddo
   enddo
   enddo

! ===============================================================

! ===============================================================
!                R U N  T H E  O P E R A T I O N S
!   Aggregate instances into the totals prior to operations
    call util_aggregate( RC )

!   Fill pointers for surface fluxes
    if (associated(sfc_flux)) call fillFluxes( import, sfc_flux, __RC__ )

!   Compute prod/loss & integrate
!   -- each species' chemistry
!   -- CURRENTLY: OH, O1D and Cl are in mcl/cm3
    call  CO_prodloss( CO, OH, O1D, Cl, CH4_total, CO2_total, RC )
    call CO2_prodloss(                                        RC ) ! Currently nothing in here. Just in case... 
    call CH4_prodloss( CH4, OH, O1D, Cl,                      RC )

!   -- surface fluxes for all instances
    call surface_prodloss( RC )

!   -- integration
    call integrate_forwardeuler( RC )

!   -- post processing
    if (cntrl%strictMassBalance) call util_accumulatenegatives( RC )

!   Aggregate instances into the totals after operations
    call util_aggregate( RC )
! ===============================================================

! ===============================================================
!      C O M P U T E  A N D  P A S S  D I A G N O S T I C S
! ===============================================================

! ===============================================================
!                            D O N E
!   Cleanup
    do i=1,NINSTANCES
       ! deallocate the prod/loss arrays for each instance
       nullify( instances(i)%p%prod, instances(i)%p%loss )
    enddo
    deallocate(met%cosz, met%slr, met%o3col, met%photj, __STAT__ )
    VERIFY_(status)

    RETURN_(ESMF_SUCCESS)

  end subroutine GridCompRun

!============================================================================

  subroutine fillFluxes( import, sfc_flux, RC )
    type (ESMF_State),              intent(inout) :: import     ! Import state
    type(surface_flux),    pointer, intent(inout) :: sfc_flux(:)
    integer,                        intent(out)   :: RC
    
    integer :: i

   __Iam__('GEOScarbon:fillFluxes')

    RC = 0

    ! This routine fills the flux fields from the ESMF import state
    ! -- loop over the registered fluxes
    do i=1,size(sfc_flux)
       ! Right now we are assuming all fluxes are 2-D
       CALL MAPL_GetPointer( import, sfc_flux(i)%flux, trim(sfc_flux(i)%shortname), notFoundOK=.FALSE., RC=RC )
       if (RC .eq. ESMF_RC_NOT_FOUND) then 
          write(*,*) 'Could not find flux, '//trim(sfc_flux(i)%shortname)//' in imports. Aborting'
          RETURN_(RC)
       endif
       if (RC .ne. ESMF_SUCCESS) then
          RETURN_(RC)
       endif
    enddo

    RETURN_(ESMF_SUCCESS)

  end subroutine fillFluxes

  subroutine RegisterFluxWithMAPL( GC, cfg, species, RC )

    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    type (ESMF_Config)                    :: cfg        ! Configuration
    character(*), intent(in)              :: species
    integer, optional                     :: RC
    
    character(len=255) :: string
    logical            :: tend

    __Iam__('RegisterFluxWithMAPL')

    ! Read the table in config
    call ESMF_ConfigFindLabel( cfg,trim(species)//'_surface_flux_pairs::',rc=RC )

    if (RC /= ESMF_SUCCESS) then ! Label not found
       if (MAPL_am_I_root()) write(*,*) '<<>> Fluxes not found for '//trim(species)//'. Exiting.'
       !_VERIFY(RC)
       RETURN_(ESMF_SUCCESS)
    endif

    ! Set up the read loop
    tend  = .false.
    do while (.not.tend)
       call ESMF_ConfigNextLine( cfg,tableEnd=tend,rc=RC )
       if (tend) cycle
       _VERIFY(RC)

       ! Get instance name
       call ESMF_ConfigGetAttribute ( cfg,value=string, default='',rc=RC) ! 1st field is the instance name
       if (string .ne. '') then ! If blank, cycle
          ! The the import pair
          call ESMF_ConfigGetAttribute ( cfg,value=string, default='',rc=RC) ! 2nd field is the flux import name
          if (string .ne. '') then
             ! Add to import state as a 2D flux
             call MAPL_AddImportSpec( GC,                            &
                  SHORT_NAME = trim(string),                         &
                  LONG_NAME  = trim(string),                         &
!                  UNITS      = '',                                   &
                  DIMS       = MAPL_DimsHorzOnly,                    &
                  VLOCATION  = MAPL_VLocationNone,                   &
                  RESTART    = MAPL_RestartSkip,   __RC__)
             if (MAPL_am_I_root()) write(*,*) 'GEOScarbon: added '//trim(string)//' to import state'
          else
             write(*,*) 'GEOScarbon:Unpaired flux in RC file for species ', trim(species) 
          endif
       else
          cycle
       endif
    end do
  end subroutine RegisterFluxWithMAPL

  subroutine RegisterInstanceWithMAPL( GC, species, name, DTM, RC )
    ! Named such because it is a MAPL-dependent routine
    ! written to keep a repetitive operation cleanly presentable.
    ! -- MSL
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    character(*)                          :: species, name
    logical, optional, intent(in)         :: DTM
    integer, optional                     :: RC

    character(len=32) :: friendlies

    __Iam__('RegisterInstanceWithMAPL')

    friendlies = 'DYNAMICS:TURBULENCE:MOIST'

    ! Toggle whether or not to advect/mix/convect
    if (present(DTM)) then
       if (.not. DTM) friendlies = ''
    endif

    call MAPL_AddInternalSpec(gc,&
         short_name =trim(species)//'_'//trim(name), &
         long_name  =trim(species)//' carbon field', &
         units      ='kg kg-1', &
         dims       =MAPL_DimsHorzVert, &
         vlocation  =MAPL_VlocationCenter, &
!         restart    =MAPL_RestartOptional, &
         friendlyto =trim(friendlies),     &
         add2export =.true., & !<-- is this what makes it available for HISTORY?
         __RC__)
    
    RETURN_(ESMF_SUCCESS)

  end subroutine RegisterInstanceWithMAPL

  subroutine ReadFluxEntries( import, internal, cfg, species, RC )
! Arguments
    integer,               optional             :: RC         ! Return code
    character(*),          intent(in)           :: species    ! Establishes the instance names
    type (ESMF_State),     intent(in)           :: import     ! Import state
    type (ESMF_State),     intent(in)           :: internal   ! Internal state
    type (ESMF_Config)                          :: cfg        ! Configuration

! Locals
    logical            :: tend, diurnal, pblmix
    integer            :: nelm, nlist, i
    character(len=255) :: string1, string2, string3, errmsg
    real               :: scalefactor

    __Iam__('ReadFluxEntries')

    if( MAPL_AM_I_ROOT() ) then
       print *, 'GEOScarbon_GridComp: Reading config file:'
       print *, '-----------------------------------------'
    endif

! Determine surface flux pairs
! ----------------------------
! This section links import fields representing surface fluxes to
! instance fields.
!
! the surface_flux object "sfc_flux" is global and associated through a
! use statement.
!

    ! Read the table in config
    call ESMF_ConfigFindLabel( cfg,trim(species)//'_surface_flux_pairs::',rc=RC )

    if (RC /= ESMF_SUCCESS) then ! Label not found
       if (MAPL_am_I_root()) write(*,*) '<<>> Flux pairing not found for '//trim(species)//'. Exiting.'
       !_VERIFY(RC)
       RETURN_(ESMF_SUCCESS)
    endif

    ! Set up the read loop
    tend  = .false.
    nlist = 0
    ! Make sure we can save the data already recorded in sfc_flux()
    if (associated(sfc_flux)) then
       nelm = size(sfc_flux)
    else
       nelm = 0
    endif

    do while (.not.tend)

!       scalefactor=1. ! Default to 1

       call ESMF_ConfigNextLine( cfg,tableEnd=tend,rc=RC )
       if (tend) cycle
       _VERIFY(RC)

       ! Get instance name
       call ESMF_ConfigGetAttribute ( cfg,value=string1, default='',rc=RC) ! 1st field is the instance name
       if (string1 /= '') then ! If not a blank line. Though, if a blank line, continue anyway.

          ! Find & verify instance
          RC = -1
          do i=1,size(instances)
             if (trim(string1) .eq. trim(instances(i)%p%name)) then
                 RC = 0
                exit
             endif
          enddo

          ! Error if not found
          if (RC /= ESMF_SUCCESS) then
             errmsg = 'GEOScarbon: Error registering flux pair. '//trim(string1)//' not an instance'
             _ASSERT(.false., errmsg)
          elseif (MAPL_am_I_root()) then
             write(*,*) 'GEOScarbon: Found instance '//trim(string1)
          endif
       else
          cycle
       endif
       ! The the import pair
       call ESMF_ConfigGetAttribute ( cfg,value=string2, default='',rc=RC) ! 1st field is the instance name
       if (string2 /= '') then ! If not blank...

          ! Find instance in internal state
          call ESMF_StateGet( import, itemSearch=trim(string2), RC=RC )

          ! Error if not found
          if (MAPL_am_I_root() .and. RC /= ESMF_SUCCESS) then
             write(*,*) 'GEOScarbon: Did not find '//trim(string2)//' in import state'
             VERIFY_(RC)
          elseif (MAPL_am_I_root()) then
             write(*,*) 'GEOScarbon: Found surface flux '//trim(string2)
          endif

       else ! string is empty, unpaired entry

          ! Exit with error
          _ASSERT(.FALSE.,'ERROR: GEOScarbon: Unpaired surface_flux entry in GEOScarbon_GridComp.rc')
          RETURN_(ESMF_SUCCESS)

       endif

       ! Diurnal or PBL mixing required?
       ! - there are three possible tags to read, though one or none may be used
       diurnal = .false.
       pblmix  = .false.
       call ESMF_ConfigGetAttribute ( cfg,value=string3, default='',rc=RC) ! 1st field is the instance name
       if (string3 /= '') then ! If not a blank line.

          string3 = ESMF_UtilStringUpperCase(string3, rc=status)

          if (trim(string3) .ne. 'P' .and. trim(string3) .ne. 'D' .and. .not. is_numeric(string3)) then
             errmsg = 'GEOScarbon::ConfigFile: invalid character in surface pairing table: '//trim(string3)
             _ASSERT(.false.,errmsg)
             RETURN_(ESMF_SUCCESS)
          endif
          
          if (trim(string3) .eq. 'P') pblmix =.true.
          if (trim(string3) .eq. 'D') diurnal=.true.
          if (is_numeric(string3)   ) read(string3,*) scalefactor

       endif
       call ESMF_ConfigGetAttribute ( cfg,value=string3, default='',rc=RC) ! 1st field is the instance name
       if (string3 /= '') then ! If not a blank line.

          string3 = ESMF_UtilStringUpperCase(string3, rc=status)

          if (trim(string3) .ne. 'P' .and. trim(string3) .ne. 'D' .and. .not. is_numeric(string3)) then
             errmsg = 'GEOScarbon::ConfigFile: invalid character in surface pairing table: '//trim(string3)
             _ASSERT(.false., errmsg)
             RETURN_(ESMF_SUCCESS)
          endif
          
          if (trim(string3) .eq. 'P') pblmix =.true.
          if (trim(string3) .eq. 'D') diurnal=.true.
          if (is_numeric(string3)   ) read(string3,*) scalefactor

       endif
       call ESMF_ConfigGetAttribute ( cfg,value=string3, default='',rc=RC) ! 1st field is the instance name
       if (string3 /= '') then ! If not a blank line.

          string3 = ESMF_UtilStringUpperCase(string3, rc=status)

          if (trim(string3) .ne. 'P' .and. trim(string3) .ne. 'D' .and. .not. is_numeric(string3)) then
             errmsg = 'GEOScarbon::ConfigFile: invalid character in surface pairing table: '//trim(string3)
             _ASSERT(.false., errmsg)
             RETURN_(ESMF_SUCCESS)
          endif
          
          if (trim(string3) .eq. 'P') pblmix =.true.
          if (trim(string3) .eq. 'D') diurnal=.true.
          if (is_numeric(string3)   ) read(string3,*) scalefactor

       endif

       ! If we got here, then we can proceed with registering the pair
       ! Increment list
       nlist = nlist + 1

       call util_addsurfaceflux(trim(string1), trim(string2), diurnal, pblmix, scalefactor, RC)
       VERIFY_(RC)
    enddo
    
    if (MAPL_am_I_root()) then
       write(*,*) '<<>> GEOScarbon::'//trim(species)//' paired surface fluxes <<>>', nlist
       do i=nelm+1,nlist
          if (sfc_flux(i)%diurnal .and. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%shortname)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with PBL mixing and diurnal scaling', sfc_flux(i)%scalefactor
          if (sfc_flux(i)%diurnal .and. .not. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%shortname)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with diurnal scaling', sfc_flux(i)%scalefactor
          if (.not. sfc_flux(i)%diurnal .and. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%shortname)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with PBL mixing', sfc_flux(i)%scalefactor
          if (.not. sfc_flux(i)%diurnal .and. .not. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%shortname)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with no scaling or mixing', sfc_flux(i)%scalefactor
       enddo
    endif

    if (nlist == 0) then
       RETURN_(ESMF_SUCCESS)
    end if

  end subroutine ReadFluxEntries
end module GEOScarbon_GridCompMod

