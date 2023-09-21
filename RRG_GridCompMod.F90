#include "MAPL_Generic.h"
!=============================================================================
!BOP

! !MODULE: RRG_GridCompMod

! !INTERFACE:
module RRG_GridCompMod
  
! !USES:
   use ESMF
   use MAPL
   use types_mod
   use global_mod
   use utils_mod
   use geos_simplephotolysismod
   use diagnostics

   implicit none
   private

   integer, parameter             :: DP=kind(1.0d0)
   integer                        :: status ! module-wide scope

   ! Users can add new species to the system w/o
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
!   For chem species C, the system simply solves dC/dt = P - L
!   P (prod) and L (loss) are computed and the system is integrated.
!   The only operator splitting is between surface fluxes and 3-D ops (chem). This
!   is done with two run methods. GridCompRun1 does surface fluxes, and GridCompRun2
!   does the chemistry.
!
! OPERATING ASSUMPTIONS:
! 1) Incoming surface fluxes are kg <species> m-2 s-1. e.g. kgCO2/m2/s
! 2) 

! NOTES:
! 1) TOTAL/AGGREGATE FIELDS: Totals of CO2, CO & CH4 are NOT independent instances.
!    They are used for diagnostics and coupling with other modules (e.g. grid comps)
!    They are not directly modified by any process.
!
! !REVISION HISTORY:
! 02Dec2022  M.S.Long   First pass
! 08Feb2023  M.S.Long   Preliminary tests against GOCART CH4, CO & CO2. Passed.
! 10Feb2023  M.S.Long   Convert operations from kg/kg to mol/mol. Species are 
!                       advected as mol/mol (<-assuming this is appropriate)
!
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

! using include is just so much cleaner
#include "IMPORTS.h" 

  call MAPL_AddImportSpec(GC,             &
       SHORT_NAME = 'CO_CH4',             &
       LONG_NAME  = 'source species',     &
       UNITS      = '1',                  &
       DIMS       = MAPL_DimsHorzVert,    &
       VLOCATION  = MAPL_VLocationCenter, &
       RESTART    = MAPL_RestartSkip,     &
       __RC__ )

! ===============================================================
!      S E T  U P  T H E  E X P O R T  S T A T E

! using include is just so much cleaner
#include "EXPORTS.h" 

! ===============================================================

! ===============================================================
!  R E A D  C O N F I G  A N D  S E T U P  I N S T A N C E S
!
!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'RRG_GridComp.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'RRG_GridComp.rc does not exist!'
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

!   Load instances    
!   --------------
    call ProcessInstances( GC, Cfg, CO,  'CO',  28.0104, nCO,  rc=status ); VERIFY_(STATUS)
    call ProcessInstances( GC, Cfg, CO2, 'CO2', 44.0098, nCO2, rc=status ); VERIFY_(STATUS)
    call ProcessInstances( GC, Cfg, CH4, 'CH4', 16.0422, nCH4, rc=status ); VERIFY_(STATUS)
    

!   Fluxes
!   ------
    call RegisterFluxWithMAPL( GC, cfg, 'CO' , status ); VERIFY_(status)
    call RegisterFluxWithMAPL( GC, cfg, 'CO2', status ); VERIFY_(status)
    call RegisterFluxWithMAPL( GC, cfg, 'CH4', status ); VERIFY_(status)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, GridCompRun1,       __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, GridCompRun2,       __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_FINALIZE, Finalize,      __RC__)
   
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

    integer                              :: i, n, dims(3)
    type (MAPL_VarSpec), pointer         :: INTERNALspec(:)  ! This is used to access GC information, e.g. field names, etc.
    
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
    params%AVO      = MAPL_AVOGAD*1e-3 ! Why is AVO in mcl/kmol?
    params%AIRMW    = MAPL_AIRMW
    grav            = MAPL_GRAV
    params%RadToDeg = 180./MAPL_PI

!   Load resource file
!   ------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'RRG_GridComp.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'RRG_GridComp.rc does not exist!'
      VERIFY_(STATUS)
    end if

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, &
                         LONS = params%LONS, &  ! Radians
                         LATS = params%LATS, &  ! Radians
                         INTERNALspec = INTERNALspec, &
                         __RC__ )

    if (nCO  .gt. 0) call ReadFluxTable( import, internal, cfg, 'CO',  __RC__ )
    if (nCO2 .gt. 0) call ReadFluxTable( import, internal, cfg, 'CO2', __RC__ )
    if (nCH4 .gt. 0) call ReadFluxTable( import, internal, cfg, 'CH4', __RC__ )

!   Read Masks
!   ----------
    if (nCO  .gt. 0) call ReadMasksTable( import, cfg, CO , __RC__ )
    if (nCO2 .gt. 0) call ReadMasksTable( import, cfg, CO2, __RC__ )
    if (nCH4 .gt. 0) call ReadMasksTable( import, cfg, CH4, __RC__ )

!DEBUG    if (MAPL_am_I_root()) call util_dumpinstances()
!DEBUG    if (MAPL_am_I_root()) call util_dumpfluxes()

!  Set physical parameters
!  - Henry's Law constants. Only needed for GEOS
!    set to ZERO so that no scavenging happens in MOIST

    do i=1,nCO
       call ESMF_StateGet(internal, trim(CO(i)%name), field, __RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', (/0,0,0,0/), __RC__)
    enddo
    do i=1,nCO2
       call ESMF_StateGet(internal, trim(CO2(i)%name), field, __RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', (/0,0,0,0/), __RC__)
    enddo
    do i=1,nCH4
       call ESMF_StateGet(internal, trim(CH4(i)%name), field, __RC__)
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
! !IROUTINE: GridCompRun1

! !INTERFACE:
  subroutine GridCompRun1 (GC, import, export, clock, RC)
!   !USES:
    use surface_mod
    use global_mod
    use integration_mod

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

    real, pointer                     :: CO2_total(:,:,:), CH4_total(:,:,:), CO_total(:,:,:)
    logical, save                     :: first = .true. ! I don't like using this but it has to happen. ExtData doesn't fill masks until run() and I don't want to repeat operations

    real, pointer                   :: ptr2d(:,:), ptr3d(:,:,:)
    type(ESMF_Alarm)                :: ALARM

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
    call MAPL_Get (mapl, RUNALARM=ALARM, __RC__)

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
    call MAPL_GetPointer(import,met%area, 'AREA', __RC__) ! Grid box area, m^2
    call MAPL_GetPointer(import,met%pblh, 'ZPBL', __RC__) ! pblh
    call MAPL_GetPointer(import,met%zle,  'ZLE',  __RC__) ! zle
    call MAPL_GetPointer(import,met%ple,  'PLE',  __RC__) ! ple
    call MAPL_GetPointer(import,met%delp, 'DELP', __RC__) ! delp
    call MAPL_GetPointer(import,met%qtot, 'QTOT', __RC__)

!XXXXXX! Get the instance data pointers
!XXXXXX    do i=1,NINSTANCES
!XXXXXX       call MAPL_GetPointer(internal, instances(i)%p%data3d, trim(instances(i)%p%name), __RC__)
!XXXXXX       allocate( instances(i)%p%prod, instances(i)%p%loss, mold=instances(i)%p%data3d, __STAT__ ) ! allocate the prod/loss arrays for each instance
!XXXXXX       instances(i)%p%prod = 0.e0; instances(i)%p%loss = 0.e0
!XXXXXX    enddo
!XXXXXX
!XXXXXX!   Get pointers to the aggregates/totals
!XXXXXX    call MAPL_GetPointer(internal, aggregate(ispecies('CO2'))%q, 'CO2', notFoundOK=.TRUE., __RC__)
!XXXXXX    if (associated(aggregate(ispecies('CO2'))%q)) CO2_total => aggregate(ispecies('CO2'))%q  ! Aggregate is used under the hood
!XXXXXX
!XXXXXX    call MAPL_GetPointer(internal, aggregate(ispecies('CO'))%q,  'CO' , notFoundOK=.TRUE., __RC__)
!XXXXXX    if (associated(aggregate(ispecies('CO'))%q))   CO_total => aggregate(ispecies('CO'))%q  ! Aggregate is used under the hood
!XXXXXX
!XXXXXX    call MAPL_GetPointer(internal, aggregate(ispecies('CH4'))%q, 'CH4', notFoundOK=.TRUE., __RC__)
!XXXXXX    if (associated(aggregate(ispecies('CH4'))%q)) CH4_total => aggregate(ispecies('CH4'))%q  ! Aggregate is used under the hood
!XXXXXX
!XXXXXX! ===============================================================
!XXXXXX
!XXXXXX! ===============================================================
!XXXXXX!                   P R O C E S S  M A S K S
!XXXXXX    if (first) then
!XXXXXX       call ProcessExtdataMasks( IMPORT, __RC__ ) 
!XXXXXX       first = .false.
!XXXXXX    endif
!XXXXXX
!XXXXXX! ===============================================================
!XXXXXX!                R U N  T H E  O P E R A T I O N S
!XXXXXX!   Aggregate instances into the totals prior to operations
!XXXXXX    call util_aggregate( RC )
!XXXXXX
!XXXXXX!   Fill pointers for surface fluxes
!XXXXXX    if (allocated(sfc_flux)) call fillFluxes( import, sfc_flux, __RC__ )
!XXXXXX
!XXXXXX!   -- surface fluxes for all instances
!XXXXXX    call surface_prodloss( RC )
!XXXXXX
!XXXXXX!   -- integration
!XXXXXX    call integrate_forwardeuler( RC )
!XXXXXX
!XXXXXX!   -- post processing
!XXXXXX    if (cntrl%strictMassBalance) call util_accumulatenegatives( RC )
!XXXXXX
!XXXXXX!   Aggregate instances into the totals after operations
!XXXXXX    call util_aggregate( RC )
!XXXXXX! ===============================================================
!XXXXXX
!XXXXXX! ===============================================================
!XXXXXX!      C O M P U T E  A N D  P A S S  D I A G N O S T I C S
!XXXXXX! ===============================================================
!XXXXXX    
!XXXXXX    ! Set emission diagnostic
!XXXXXX    call MAPL_GetPointer( export, Ptr2d, 'CO2_EM', __RC__ )
!XXXXXX    if (associated(Ptr2d)) then
!XXXXXX       call diag_sfcflux( 'CO2', Ptr2d, RC )
!XXXXXX       Ptr2d => null()
!XXXXXX    endif
!XXXXXX    ! Set emission diagnostic
!XXXXXX    call MAPL_GetPointer( export, Ptr2d, 'CH4_EM', __RC__ )
!XXXXXX    if (associated(Ptr2d)) then
!XXXXXX       call diag_sfcflux( 'CH4', Ptr2d, RC )
!XXXXXX       Ptr2d => null()
!XXXXXX    endif
!XXXXXX    ! Set emission diagnostic
!XXXXXX    call MAPL_GetPointer( export, Ptr2d, 'CO_EM', __RC__ )
!XXXXXX    if (associated(Ptr2d)) then
!XXXXXX       call diag_sfcflux( 'CO', Ptr2d, RC )
!XXXXXX       Ptr2d => null()
!XXXXXX    endif
!XXXXXX
!XXXXXX    !
!XXXXXX   !call MAPL_GetPointer( export, Ptr3D, 'CO2DRY', __RC__) 
!XXXXXX   !if (associated(Ptr3d)) then
!XXXXXX   !   Ptr3d = (CO2_total*MAPL_AIRMW/44.0098)/(1.e0 - met%qtot)
!XXXXXX   !   Ptr3d => null()
!XXXXXX   !endif
!XXXXXX   !
!XXXXXX   !call MAPL_GetPointer( export, Ptr3D, 'CH4DRY', __RC__) 
!XXXXXX   !if (associated(Ptr3d)) then
!XXXXXX   !   Ptr3d = (CH4_total*MAPL_AIRMW/16.0422)/(1.e0 - met%qtot)
!XXXXXX   !   Ptr3d => null()
!XXXXXX   !endif
!XXXXXX   !
!XXXXXX   !call MAPL_GetPointer( export, Ptr3D, 'CODRY', __RC__) 
!XXXXXX   !if (associated(Ptr3d)) then
!XXXXXX   !   Ptr3d = (CO_total*MAPL_AIRMW/28.0104)/(1.e0 - met%qtot)
!XXXXXX   !   Ptr3d => null()
!XXXXXX   ! endif
!XXXXXX
!XXXXXX! ===============================================================
!XXXXXX!                            D O N E
!XXXXXX!   Cleanup
!XXXXXX!   Lots of memory leak potential here. Need to be thorough
!XXXXXX    do i=1,NINSTANCES
!XXXXXX       ! deallocate the prod/loss arrays for each instance
!XXXXXX       deallocate( instances(i)%p%prod, instances(i)%p%loss, __STAT__ )
!XXXXXX       instances(i)%p%data3d => null()
!XXXXXX    enddo
!XXXXXX
!XXXXXX    CO2_total => null(); CO_total => null(); CH4_total => null()
!XXXXXX    do i=1,nspecies
!XXXXXX       if (associated(aggregate(i)%q)) aggregate(i)%q => null()
!XXXXXX    enddo
!XXXXXX
!XXXXXX    if (allocated(sfc_flux)) then
!XXXXXX       do i=1,size(sfc_flux)
!XXXXXX          if (associated(sfc_flux(i)%flux)) sfc_flux(i)%flux => null()
!XXXXXX       enddo
!XXXXXX    endif
    RETURN_(ESMF_SUCCESS)

  end subroutine GridCompRun1

!============================================================================
!BOP
! !IROUTINE: GridCompRun2

! !INTERFACE:
  subroutine GridCompRun2 (GC, import, export, clock, RC)
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
    logical, save                     :: first = .true. ! I don't like using this but it has to happen. ExtData doesn't fill masks until run() and I don't want to repeat operations

    real, pointer                   :: ptr2d(:,:), ptr3d(:,:,:), CO2ptr(:,:,:), CH4ptr(:,:,:), COptr(:,:,:)
    real, pointer, dimension(:,:,:) :: O3, OH, Cl, O1D 
    real(ESMF_KIND_R4), allocatable :: O3col(:,:,:), O2col(:,:,:), CO2photj(:,:,:), CH4photj(:,:,:)
    real(ESMF_KIND_R4), allocatable :: ZTH(:,:)
    real(ESMF_KIND_R4), allocatable :: SLR(:,:)
    type (MAPL_SunOrbit)            :: ORBIT
    type(ESMF_Alarm)                :: ALARM

    real                            :: r, m

   __Iam__('Run2')

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
    call MAPL_GetPointer(import,met%pblh,   'ZPBL', __RC__)
    call MAPL_GetPointer(import,met%T,      'T',    __RC__)
    call MAPL_GetPointer(import,met%zle,  'ZLE',  __RC__) ! zle
    call MAPL_GetPointer(import,met%ple,    'PLE',  __RC__)
    call MAPL_GetPointer(import,met%delp,   'DELP', __RC__)
    call MAPL_GetPointer(import,met%q,         'Q', __RC__)
!    call MAPL_GetPointer(import,met%qctot, 'QCTOT', __RC__)
    call MAPL_GetPointer(import,met%qtot,   'QTOT', __RC__)
    call MAPL_GetPointer(import,met%rho, 'AIRDENS', __RC__)
    CALL MAPL_GetPointer(import,     O3,      'O3', __RC__)
    CALL MAPL_GetPointer(import,     OH,  'RRG_OH', __RC__)
    CALL MAPL_GetPointer(import,     Cl,  'RRG_Cl', __RC__)
    CALL MAPL_GetPointer(import,    O1D, 'RRG_O1D', __RC__)

    allocate(  met%cosz(size(params%lats,1), size(params%lats,2)), __STAT__)
    allocate(  met%slr (size(params%lats,1), size(params%lats,2)), __STAT__)
    allocate(  O3col(params%im,params%jm,params%km), __STAT__)
    allocate(  O2col, CO2photj, CH4photj, mold=o3col, __STAT__)

    CH4photj = 0.

! Get the instance data pointers
    do i=1,NINSTANCES
       call MAPL_GetPointer(internal, instances(i)%p%data3d, trim(instances(i)%p%name), __RC__)
       allocate( instances(i)%p%prod, instances(i)%p%loss, mold=instances(i)%p%data3d, __STAT__ ) ! allocate the prod/loss arrays for each instance
       instances(i)%p%prod = 0.e0; instances(i)%p%loss = 0.e0
    enddo

!   Get pointers to the aggregates/totals
!   CO_total, CO2_total and CH4_total variables are made available for convenience.
    call MAPL_GetPointer(internal, aggregate(ispecies('CO2'))%q, 'CO2', notFoundOK=.TRUE.,__RC__)
    if (associated(aggregate(ispecies('CO2'))%q)) CO2_total => aggregate(ispecies('CO2'))%q  ! Aggregate is used under the hood

    call MAPL_GetPointer(internal, aggregate(ispecies('CO'))%q,  'CO' , notFoundOK=.TRUE.,__RC__)
    if (associated(aggregate(ispecies('CO'))%q))  CO_total  => aggregate(ispecies('CO'))%q  ! Aggregate is used under the hood

    call MAPL_GetPointer(internal, aggregate(ispecies('CH4'))%q, 'CH4', notFoundOK=.TRUE.,__RC__)
    if (associated(aggregate(ispecies('CH4'))%q)) CH4_total => aggregate(ispecies('CH4'))%q  ! Aggregate is used under the hood

! ===============================================================

! ===============================================================
!                S E T  U P  P H O T O L Y S I S
!  Update solar zenith angle
!  --------------------------
    call MAPL_SunGetInsolation(params%lons, params%lats, orbit, met%cosz, met%slr, clock=clock, __RC__)

!  Compute the O2 & O3 column
!  -- adapted from GOCART
!  --------------------------
    ! Below, I don't understand where '0.50' comes from (unless its for the average delp), and why a 1.e-4 factor is needed
    ! to make O2col work. This was transferred from CH4_GridCompMod.F90 in GOCART. -- MSL
    r = MAPL_AVOGAD*0.50/(MAPL_GRAV*MAPL_AIRMW) ! r = Nsuba*0.50/(mwtAir*grav), copied from CFC_GridCompMod.F90 6.022e26 = mcl/kmole
    m = MAPL_AIRMW/48.0e0  ! MW_air/MW_O3
    O3col(:,:,1) = 1.1e11       + m*O3(:,:,1)* met%delp(:,:,1)*r ! 1.1e11 = O3 above 80km (cm-2); O3 is mmr, we need vmr for this so MW conversion is done here
    O2col(:,:,1) = 7.072926E+19 +   0.20946  * met%delp(:,:,1)*r*1e-4 ! 7.072926E+19 = O2 above 80 km (cm-2); O2VMR = 0.20946
    DO k=2,params%km
       O3col(:,:,k) = O3col(:,:,k-1) + m*r*(O3(:,:,k-1)*met%delp(:,:,k-1) + O3(:,:,k)*met%delp(:,:,k))
       O2col(:,:,k) = O2col(:,:,k-1) + r*0.20946*(met%delp(:,:,k-1)+met%delp(:,:,k))*1e-4
    END DO
   
    O3col = 0.

!  Compute the photolysis rate for CO2 + hv -> CO + O*
    met%photj = 0.e0
    do k=1,params%km
    do j=1,params%jm
    do i=1,params%im
       call CO2_photolysis_rate(i,j,k, params%km-k+1, met, O3col(i,j,k),       CO2photj(i,j,k))
       call CH4_photolysis_rate(i,j,k,                met, O2col(i,j,k), 94.0, CH4photj(i,j,k))
    enddo
    enddo
    enddo

! ===============================================================

! ===============================================================
!                R U N  T H E  O P E R A T I O N S
!   Aggregate instances into the totals prior to operations
    call util_aggregate( RC )

!   Compute prod/loss & integrate
!   -- each species' chemistry
!   -- CURRENTLY: OH, O1D and Cl are in mcl/cm3

!   ==================================================
!   This section allows species operations
!   to use external sources for CO2 & CH4 in the 
!   case where a user disables a species (e.g. declares no instances of CH4)

!   At this point, local species (CO, CH4 & CO2) are mol/mol.
!   <<>> the commented code below is meant to enable CO to use 
!        RRG's CH4 in place of an offline field
!    if (nCH4 .gt. 0) then
!       CH4ptr => CH4_total
!       CH4ptr = CH4ptr !*params%AirMW/16.0422
!    else
       call MAPL_GetPointer(import, CH4ptr, 'CO_CH4',__RC__)
!    endif

!<<>>
!   Fill pointers for surface fluxes
    if (allocated(sfc_flux)) call fillFluxes( import, sfc_flux, __RC__ )

!   -- surface fluxes for all instances
    call surface_prodloss( RC )
!<<>>

    if (nCO2 .gt. 0) then
       CO2ptr => CO2_total
    else
       ! Need an external CO2 source. call MAPL_GetPointer(import, CO2ptr, 'CO_CH4',__RC__)
    endif

    if (nCO  .gt. 0) then
       COptr  => CO_total
    else
       ! Need an external CO source. call MAPL_GetPointer(import, COptr, 'CO_?',__RC__)
    endif
!   ==================================================

    call  CO_prodloss(  CO, OH, O1D, Cl, CO2photj, CH4photj, CH4ptr, CO2ptr,       RC )
    call CO2_prodloss( CO2, OH, COptr,                                             RC ) ! Currently nothing in here. Just in case... 
    call CH4_prodloss( CH4, OH, O1D, Cl, CH4photj,                                 RC )

    CH4ptr => null()
    CO2ptr => null()
    COptr  => null()

!   -- integration
    call integrate_forwardeuler( RC )

!   -- post processing
    if (cntrl%strictMassBalance) call util_accumulatenegatives( RC )

!   Aggregate instances into the totals after operations
    call util_aggregate( RC )
! ===============================================================

! ===============================================================
!      C O M P U T E  A N D  P A S S  D I A G N O S T I C S

!<<>>
    ! Set emission diagnostic
    call MAPL_GetPointer( export, Ptr2d, 'CO2_EM', __RC__ )
    if (associated(Ptr2d)) then
       call diag_sfcflux( 'CO2', Ptr2d, RC )
       Ptr2d => null()
    endif
!<<>>

    call MAPL_GetPointer( export, Ptr3d, 'CO2_ProdLoss', __RC__ )
    if (associated(Ptr3d)) then
       call diag_prodloss( CO2(:), Ptr3d, RC )
       Ptr3d = Ptr3d * met%delp/grav ! Convert to kg/m2/s
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'CO2DRY', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = (CO2_total                   )/(1.e0 - met%qtot)
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'CH4DRY', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = (CH4_total                   )/(1.e0 - met%qtot)
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'CODRY', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = (CO_total                   )/(1.e0 - met%qtot)
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'CO2photj', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = CO2photj
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'CH4photj', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = CH4photj
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'O3col', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = log10(O3col)
       Ptr3d => null()
    endif

    call MAPL_GetPointer( export, Ptr3D, 'O2col', __RC__) 
    if (associated(Ptr3d)) then
       Ptr3d = log10(O2col)
       Ptr3d => null()
    endif

! ===============================================================

! ===============================================================
!                            D O N E
!   Cleanup
!   Lots of memory leak potential here. Beware!
    do i=1,NINSTANCES
       ! deallocate the prod/loss arrays for each instance
       deallocate( instances(i)%p%prod, instances(i)%p%loss, __STAT__ )
       instances(i)%p%data3d => null()
    enddo
    CO2_total => null(); CO_total => null(); CH4_total => null()
    do i=1,nspecies
       aggregate(i)%q => null()
    enddo
    deallocate(met%cosz, met%slr, O3col, O2col, CO2photj, CH4photj, __STAT__ )

    RETURN_(ESMF_SUCCESS)

  end subroutine GridCompRun2

!============================================================================

   subroutine Finalize ( GC, IMPORT, EXPORT, clock, RC )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout) :: GC      ! Grid Component
   type(ESMF_State),    intent(inout) :: IMPORT      ! Import State
   type(ESMF_State),    intent(inout) :: EXPORT      ! Export State
   integer,             intent(out)   :: RC          ! Error return code:
                                                     !  0 - all is well
                                                     !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!  18Jan2023 M.Long - Adapted for GHG
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   integer :: i

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // '::' // 'Finalize_'

   ! Instances
   do i=1,size(CO)
      if (associated(CO(i)%data3d)) deallocate(CO(i)%data3d)
      if ( allocated(CO(i)%prod)  ) deallocate(CO(i)%prod)
      if ( allocated(CO(i)%loss)  ) deallocate(CO(i)%loss)
      if (associated(CO(i)%mask)  ) deallocate(CO(i)%mask)
   enddo
   if (associated(CO)) deallocate(CO)
   do i=1,size(CO2)
      if (associated(CO2(i)%data3d)) deallocate(CO2(i)%data3d)
      if ( allocated(CO2(i)%prod)  ) deallocate(CO2(i)%prod)
      if ( allocated(CO2(i)%loss)  ) deallocate(CO2(i)%loss)
      if (associated(CO2(i)%mask)  ) deallocate(CO2(i)%mask)
   enddo
   if (associated(CO2)) deallocate(CO2)
   do i=1,size(CH4)
      if (associated(CH4(i)%data3d)) deallocate(CH4(i)%data3d)
      if ( allocated(CH4(i)%prod)  ) deallocate(CH4(i)%prod)
      if ( allocated(CH4(i)%loss)  ) deallocate(CH4(i)%loss)
      if (associated(CH4(i)%mask)  ) deallocate(CH4(i)%mask)
   enddo
   if (associated(CH4)) deallocate(CH4)
   do i=1,size(instances)
      if (associated(instances(i)%p%data3d)) deallocate(instances(i)%p%data3d)
      if ( allocated(instances(i)%p%prod)  ) deallocate(instances(i)%p%prod)
      if ( allocated(instances(i)%p%loss)  ) deallocate(instances(i)%p%loss)
      if (associated(instances(i)%p%mask)  ) deallocate(instances(i)%p%mask)
   enddo
   if (allocated(instances)) deallocate(instances)

   ! Surface fluxes
   do i=1,size(sfc_flux)
      if (associated(sfc_flux(i)%flux)) deallocate(sfc_flux(i)%flux)
   enddo
   if (allocated(sfc_flux)) deallocate(sfc_flux)

   ! Aggregates
   do i=1,size(aggregate)
      if (associated(aggregate(i)%q)) deallocate(aggregate(i)%q)
   enddo
   if (allocated(aggregate)) deallocate(aggregate)

   ! Parameters
   if (associated(params%lats)) deallocate(params%lats)
   if (associated(params%lons)) deallocate(params%lons)
   if (allocated(species)) deallocate(species)

!  Finalize GEOS Generic
!  ---------------------
!ALT: do not deallocate "foreign objects"
   call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK, __RC__ )

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

!============================================================================
!      G E O S / M A P L / E S M F  S E R V I C E  R O U T I N E S
!============================================================================

  subroutine fillFluxes( import, sfc_flux, RC )
    type (ESMF_State),              intent(inout) :: import     ! Import state
    type(surface_flux),             intent(inout) :: sfc_flux(:)
    integer,                        intent(out)   :: RC
    
    integer :: i

   __Iam__('RRG:fillFluxes')

    RC = 0

    ! This routine fills the flux fields from the ESMF import state
    ! -- loop over the registered fluxes
    do i=1,size(sfc_flux)
       ! Right now we are assuming all fluxes are 2-D
       CALL MAPL_GetPointer( import, sfc_flux(i)%flux, trim(sfc_flux(i)%name), notFoundOK=.FALSE., RC=RC )
       if (RC .eq. ESMF_RC_NOT_FOUND) then 
          write(*,*) 'Could not find flux, '//trim(sfc_flux(i)%name)//' in imports. Aborting'
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
             if (MAPL_am_I_root()) write(*,*) 'RRG: added '//trim(string)//' to import state'
          else
             write(*,*) 'RRG:Unpaired flux in RC file for species ', trim(species) 
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
         units      ='mol mol-1', &
         dims       =MAPL_DimsHorzVert, &
         vlocation  =MAPL_VlocationCenter, &
!         restart    =MAPL_RestartOptional, &
         friendlyto =trim(friendlies),     &
         add2export =.true., & !<-- is this what makes it available for HISTORY?
         __RC__)
    
    RETURN_(ESMF_SUCCESS)

  end subroutine RegisterInstanceWithMAPL

  subroutine ReadFluxTable( import, internal, cfg, species, RC )
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

    __Iam__('Readfluxtable')

    if( MAPL_AM_I_ROOT() ) then
       print *, 'RRG_GridComp: Reading config file:'
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
    if (allocated(sfc_flux)) then
       nelm = size(sfc_flux)
    else
       nelm = 0
    endif

    do while (.not.tend)

       scalefactor=1. ! Default to 1

       call ESMF_ConfigNextLine( cfg,tableEnd=tend,rc=RC )
       if (tend) cycle
       _VERIFY(RC)

       ! Get instance name
       call ESMF_ConfigGetAttribute ( cfg,value=string1, default='',rc=RC) ! 1st field is the instance name
       if (string1 /= '') then ! If not a blank line. Though, if a blank line, continue anyway.

          ! Find & verify instance
          RC = -1
          do i=1,size(instances)
             if (trim(species)//'_'//trim(string1) .eq. trim(instances(i)%p%name)) then
                 RC = 0
                exit
             endif
          enddo

          ! Error if not found
          if (RC /= ESMF_SUCCESS) then
             errmsg = 'RRG: Error registering flux pair. '//trim(string1)//' not an instance of '//trim(species)
             _ASSERT(.false., errmsg)
          elseif (MAPL_am_I_root()) then
             write(*,*) 'RRG: Found instance '//trim(string1)//' for species '//trim(species)
          endif
          string1 = trim(species)//'_'//trim(string1)
       else
          cycle
       endif
       ! The the import pair
       call ESMF_ConfigGetAttribute ( cfg,value=string2, default='',rc=RC) ! 2nd field is the ExtData flux name
       if (string2 /= '') then ! If not blank...

          ! Find instance in internal state
          call ESMF_StateGet( import, itemSearch=trim(string2), RC=RC )

          ! Error if not found
          if (MAPL_am_I_root() .and. RC /= ESMF_SUCCESS) then
             write(*,*) 'RRG: Did not find '//trim(string2)//' in import state'
             VERIFY_(RC)
          elseif (MAPL_am_I_root()) then
             write(*,*) 'RRG: Found surface flux '//trim(string2)
          endif

       else ! string is empty, unpaired entry

          ! Exit with error
          _ASSERT(.FALSE.,'ERROR: RRG: Unpaired surface_flux entry in RRG_GridComp.rc')
          RETURN_(ESMF_SUCCESS)

       endif

       ! Diurnal or PBL mixing required?
       ! - there are three possible tags to read, though one or none may be used
       diurnal = .false.
       pblmix  = .false.
       call ESMF_ConfigGetAttribute ( cfg,value=string3, default='',rc=RC)
       if (string3 /= '') then ! If not a blank line.

          string3 = ESMF_UtilStringUpperCase(string3, rc=status)

          if (trim(string3) .ne. 'P' .and. trim(string3) .ne. 'D' .and. .not. is_numeric(string3)) then
             errmsg = 'RRG::ConfigFile: invalid character in surface pairing table: '//trim(string3)
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
             errmsg = 'RRG::ConfigFile: invalid character in surface pairing table: '//trim(string3)
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
             errmsg = 'RRG::ConfigFile: invalid character in surface pairing table: '//trim(string3)
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
       write(*,*) '<<>> RRG::'//trim(species)//' paired surface fluxes <<>>', nlist
       do i=nelm+1,nelm+nlist
          if (sfc_flux(i)%diurnal .and. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%name)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with PBL mixing and diurnal scaling', sfc_flux(i)%scalefactor
          if (sfc_flux(i)%diurnal .and. .not. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%name)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with diurnal scaling', sfc_flux(i)%scalefactor
          if (.not. sfc_flux(i)%diurnal .and. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%name)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with PBL mixing', sfc_flux(i)%scalefactor
          if (.not. sfc_flux(i)%diurnal .and. .not. sfc_flux(i)%pblmix) &
               write(*,'(a,f14.5)') '<<>> Flux '//trim(sfc_flux(i)%name)//' paired with '//trim(species)//'::'//trim(sfc_flux(i)%instance_pair)//' with no scaling or mixing', sfc_flux(i)%scalefactor
       enddo
    endif

    if (nlist == 0) then
       RETURN_(ESMF_SUCCESS)
    end if

  end subroutine ReadFluxTable

  subroutine ReadMasksTable( import, cfg, instance, RC )

    implicit none

    type (ESMF_State)                 :: import     ! Import state
    type(ESMF_Config)                 :: cfg     ! Configuration
    type(gas_instance), intent(inout) :: instance(:)
    integer, intent(out)              :: RC

    logical       :: tend, present, found
    character(32) :: name, string, lstr, species
    integer       :: i, imask, ierr, loc1, loc2, ii, jj, n, nterms
    real          :: lat1,lat2,lon1,lon2,ltmp
    real, pointer :: ptr2d(:,:)
    __Iam__('ReadMasksTable')

    RC    = 0 ! Assume success
    
    species = trim(instance(1)%species)

    ! Read the table in config
    call ESMF_ConfigFindLabel( cfg,trim(species)//'_masks::',isPresent=present,rc=status )
    VERIFY_(status)
    if (.not. present) return
    call ESMF_ConfigGetDim( cfg, lineCount=n, columnCount=nterms, rc=status) ! 'n' is dummy. lineCount isn't used
    VERIFY_(STATUS)
    call ESMF_ConfigFindLabel( cfg,trim(species)//'_masks::',isPresent=present,rc=status )
    VERIFY_(status)

    ! Set up the read loop
    tend  = .false.
    do while (.not.tend)
       ! Presets 
       imask = -999
       lon1  = -999.9
       lon2  = -999.9
       lat1  = -999.9
       lat2  = -999.9

       call ESMF_ConfigNextLine( cfg,tableEnd=tend,rc=status )
       VERIFY_(status)
 
       if (tend) cycle

       ! Get instance name
       ! 1st field is the instance name associated with the mask
!       nterms = ESMF_ConfigGetLen( cfg, rc=status ) ! TBD
       ! 1st entry is the associated instance
       call ESMF_ConfigGetAttribute ( cfg,value=name, default='',rc=status)
       if (name .ne. '') then ! If blank, cycle
          ! Loop over entries
          do n=2,nterms
             ! 1) Get term
             call ESMF_ConfigGetAttribute ( cfg,value=string,rc=status)
             if (string .eq. '') cycle

             ! 2) Find instance and allocate mask
             found = .false.
             do i=1,size(instance)
                if (trim(species)//'_'//trim(name) .eq. trim(instance(i)%name)) then
                   found = .true.
                   exit
                   ! i is now the instance index
                endif
             enddo
             if (.not. found) then
                write(*,*) 'ReadMasksTable: not able to find instance '//trim(species)//'_'//trim(name)//' associated with '//trim(instance(1)%species)
                RC = -1
                return
             endif
             instance(i)%hasmask = .true.
             if (.not. associated(instance(i)%mask)) then
                ! Only do this if not already associated. Otherwise, we can't
                ! superimpose additional masks. They'd get overwritten
                allocate(instance(i)%mask(params%im,params%jm),stat=status)
                VERIFY_(status)
                instance(i)%mask = 0
             endif

             ! 3) Is entry an integer reference or a lat lon box?
             if (index(string,',') .gt. 0) then
                ! 4) Process entry
                !    4b) If a lat/lon box
                call process_latlonbox( string, status )
             elseif (is_numeric(string)) then
                !    4a) If an integer reference. Test-for-integer happens in process_imask()
                call process_imask( instance(i), string, status )
             endif
             ! 5) Should be done. Dump info if needed
          enddo
       else
          cycle ! blank line. Continue
       endif
    end do

    return

  contains
    subroutine Process_LatLonBox( term, RC )

      implicit none
      
      character(*), intent(in)  :: term
      integer,      intent(out) :: RC

      logical :: xsubgrid = .false.
      logical :: ysubgrid = .false.

      real    :: dlat,dlon

      ! process lat1
      loc1 = index(term(:),',') ! Find the 1st delimiter
      if (loc1 .gt. 0) then
         read(term(1:loc1-1), *) lstr
         if (is_numeric(lstr)) then
            read(lstr,*) lat1
         else
            ! error
            write(*,*) 'ReadMasksTable: Error processing 1st latitude bound of mask '//trim(name)
            RC = -1
            return
         endif
      else
         ! No delimters? Error!
      endif

      ! process lat2
      loc2 = index(term(loc1+1:),',')+loc1 ! Find the 2nd delimiter
      if (loc2 .gt. 0) then
         read(term(loc1+1:loc2-1),*) lstr ! Read between the commas
         if (is_numeric(lstr)) then
            read(lstr,*) lat2
         else
            ! error
            write(*,*) 'ReadMasksTable: Error processing 2nd latitude bound of mask '//trim(name),trim(lstr)
            RC = -1
            return
         endif
         loc1 = loc2 ! stride forward
      else
         ! No delimters? Error!
      endif

      ! process lon1 & lon2 together
      loc2 = index(term(loc1+1:),',')+loc1 ! Find the 3rd delimiter
      if (loc2 .gt. 0) then
         read(term(loc1+1:loc2-1),*) lstr ! Read between the commas
         if (is_numeric(lstr)) then
            read(lstr,*) lon1
         else
            ! error
            write(*,*) 'ReadMasksTable: Error processing 1st longitude bound of mask '//trim(name),trim(lstr)
            RC = -1
            return
         endif
         read(term(loc2+1:),*) lstr ! Read last entry
         if (is_numeric(lstr)) then
            read(lstr,*) lon2
         else
            ! error
            write(*,*) 'ReadMasksTable: Error processing 2nd longitude bound of mask '//trim(name),trim(lstr)
            RC = -1
            return
         endif
         ! Fin
      else
         ! No delimters? Error!
         write(*,*) 'ReadMasksTable: Mask setting for '//trim(name)//' is neither an integer nor a proper lat/lon box'
         RC = -1
         return
      endif
      ! If we got here, then we have successfully processed a lat/lon box

      ! build a box from the lats/lons & superimpose on mask
      ! This is primitively done
      ! 1) convert box to radians
      lat1 = lat1 * MAPL_PI/180.0
      lat2 = lat2 * MAPL_PI/180.0
      if (lon1 .lt. 0.) then
         lon1 = (lon1 + 360.0) * MAPL_PI/180.0
      else
         lon1 = lon1 * MAPL_PI/180.0
      endif
      if (lon2 .lt. 0.) then
         lon2 = (lon2 + 360.0) * MAPL_PI/180.0
      else
         lon2 = lon2 * MAPL_PI/180.0
      endif
      ! 2) Populate mask
      write(*,'(a,8f12.2)') 'TESTING MASK: ', lat1,lat2,lon1,lon2, minval(params%lats), maxval(params%lats), minval(params%lons), maxval(params%lons)
      do jj=1,params%jm-1
         do ii=1,params%im-1
            if (lon2 .gt. lon1) then ! we don't straddle the prime meridian
               if ( params%lats(ii,jj) .ge. lat1 .and. &
                    params%lats(ii,jj+1) .le. lat2 .and. &
                    params%lons(ii,jj) .ge. lon1 .and. &
                    params%lons(ii+1,jj) .le. lon2 ) then
                  instance(i)%mask(ii,jj) = 1
                  write(*,'(a,4f12.2,a,i2,a,i2,a)') 'LAT/LON MASK: ', lat1,lat2,lon1,lon2,'(',ii,', ',jj,')'
               endif
            else ! We straddle the prime meridian
               if ( params%lats(ii,jj) .ge. lat1 .and. &
                    params%lats(ii,jj+1) .le. lat2 .and. &
                    params%lons(ii,jj) .ge. lon1 .and. &
                    params%lons(ii+1,jj)-2.*MAPL_PI .le. lon2 ) then
                  instance(i)%mask(ii,jj) = 1
               endif
            endif
         enddo
      enddo
      ii = params%im
      jj = params%jm
      dlat = params%lats(ii,jj)-params%lats(ii,jj-1)
      dlon = params%lons(ii,jj)-params%lons(ii,jj-1)
      if (lon2 .gt. lon1) then ! we don't straddle the prime meridian
         if ( params%lats(ii,jj) .ge. lat1 .and. &
              params%lats(ii,jj+1) .le. lat2 .and. &
              params%lons(ii,jj) .ge. lon1 .and. &
              params%lons(ii+1,jj) .le. lon2 ) then
            instance(i)%mask(ii,jj) = 1
            write(*,'(a,4f12.2,a,i2,a,i2,a)') 'LAT/LON MASK: ', lat1,lat2,lon1,lon2,'(',ii,', ',jj,')'
         endif
      else ! We straddle the prime meridian
         if ( params%lats(ii,jj) .ge. lat1 .and. &
              params%lats(ii,jj)+dlat .le. lat2 .and. &
              params%lons(ii,jj) .ge. lon1 .and. &
              (params%lons(ii,jj)-2.*MAPL_PI)+dlon .le. lon2 ) then
            instance(i)%mask(ii,jj) = 1
         endif
      endif

    end subroutine Process_LatLonBox

    subroutine Process_iMask( instance_, term, RC )

      implicit none

      type(gas_instance), intent(inout) :: instance_
      character(*),       intent(in)    :: term
      integer,            intent(out)   :: RC

      integer :: nmask

      nmask = findloc( instance_%imask, -999, DIM=1 ) ! Find the 1st blank entry

      ! Make sure we aren't over the number of masks
      if (nmask .eq. 0) then ! Looks like imask(:) is all full
         write(*,*) 'ReadMasksTable: Too many integer mask references. Increase MAXMASKS. Currently ',MAXMASKS
         RC = -1
         return
      endif

      ! If its a number, it better be an integer!
      read(term,'(i10)',iostat=ierr) instance_%imask(nmask) ! term should write into an integer w/o error, otherwise, not an integer
      if ( ierr .ne. 0 ) then
         ! not an integer
         write(*,'(a)') 'ReadMasksTable: Error reading mask named '//trim(instance_%name)//'. Expected an integer, got '//trim(term)//'. Exiting'
         RC = -1
         return
      endif

    end subroutine Process_iMask

  end subroutine ReadMasksTable

  subroutine ProcessExtdataMasks( import, RC ) 

    implicit none
    type (ESMF_State)           :: import     ! Import state
    integer, intent(out)        :: RC

    integer                     :: i,j,nst,n,nmasks
    real, pointer               :: Ptr2D(:,:)
    type(gas_instance), pointer :: instance
    character(32)               :: species

    RC = 0

    do nst=1,NINSTANCES
       instance => instances(nst)%p
       ! Loop over number of masks
       if (.not. instance%hasmask) cycle
       nmasks = findloc(instance%imask, -999, DIM=1)
       if (nmasks .eq. 0) nmasks = MAXMASKS ! This is because the imask vector is full (all .ne. -999)
       do n=1,nmasks
       if (instance%imask(n) .ne. -999) then
          species = instance%species
          call MAPL_GetPointer(import, Ptr2D, trim(species)//'_Mask', RC=RC)
          if (RC .ne. ESMF_SUCCESS) return
          where (Ptr2D .eq. instance%imask(n)) instance%mask = 1
          Ptr2D => null()
       endif
       enddo
       instance => null()
    enddo
    
  end subroutine ProcessExtdataMasks

  subroutine ProcessInstances( GC, cfg, GI, species, MW, nInst, RC )

    implicit none
    type (ESMF_GridComp),        intent(inout) :: GC      ! gridded component
    type (ESMF_Config),          intent(inout) :: cfg
    type(gas_instance), pointer, intent(inout) :: GI(:)
    character(*),                intent(in)    :: Species
    real,                        intent(in)    :: MW
    integer,                     intent(out)   :: nInst   ! Number of registered instances
    integer,                     intent(out)   :: RC      ! return code

    ! Local
    integer           :: i,j,n
    character(len=32) :: inst_name
    logical           :: isPresent, found
    logical           :: isActive = .true.

    __Iam__('ProcessInstances')

    RC = 0

    nInst = ESMF_ConfigGetLen(cfg,label=trim(species)//'_instances:',rc=status)
    VERIFY_(STATUS)

    if (nInst .eq. 0) return

    !  define the total/aggregate field
    call MAPL_AddInternalSpec(gc,                           &
         short_name =trim(species),                         &
         long_name  ='Aggregate '//trim(species)//' field', &
         units      ='mol mol-1',                           &
         dims       =MAPL_DimsHorzVert,                     &
         vlocation  =MAPL_VlocationCenter,                  &
         restart    =MAPL_RestartOptional,                  &
         add2export =.true.,                                & !<-- is this what makes it available for HISTORY?
         __RC__)

    !  Create a region mask import
    call MAPL_AddImportSpec(GC,                          &
         SHORT_NAME = trim(species)//'_Mask',      &
         LONG_NAME  = '',                                &
         UNITS      = '',                                &
         DIMS       = MAPL_DimsHorzOnly,                 &
         VLOCATION  = MAPL_VLocationNone,                &
         RESTART    = MAPL_RestartSkip,   __RC__)

    !  Get instances from RC file
    !  ----------------------------------
    call ESMF_ConfigFindLabel(cfg,trim(species)//'_instances:',rc=status)
    VERIFY_(STATUS)

    do i = 1, nInst
       call ESMF_ConfigGetAttribute(cfg,inst_name,rc=status)
       VERIFY_(STATUS)
       inst_name = TRIM(inst_name)
       ! Register as tracer
       call RegisterInstanceWithMAPL( GC, trim(species), trim(inst_name), rc=status )
       VERIFY_(STATUS)

       ! Add new active instance (assume active)
       call Util_AddInstance( GI, trim(inst_name), trim(species), MW, isActive, status)
       VERIFY_(STATUS)
    end do

    ! If there are instances, then define an active residual by default
    ! ASSUMPTION: a species' residual is always the last instance 
    call Util_AddInstance( GI, 'residual', trim(species), MW, isActive, status)
    VERIFY_(STATUS)
    nInst = nInst+1
    call RegisterInstanceWithMAPL( GC, trim(species), 'residual', rc=status )
    VERIFY_(STATUS)

    !  Get passive instances and toggle them
    n = 0
    call ESMF_ConfigFindLabel(cfg,label=trim(species)//'_passive_instances:',isPresent=ispresent,rc=status)
    VERIFY_(STATUS)
    if (ispresent) then
       n = ESMF_ConfigGetLen(cfg,label=trim(species)//'_passive_instances:',rc=status)
       VERIFY_(STATUS)
       call ESMF_ConfigFindLabel(cfg,trim(species)//'_passive_instances:',rc=status)
       VERIFY_(STATUS)
       if (n .gt. 0) then
          do i = 1, n
             ! Get instance name
             call ESMF_ConfigGetAttribute(cfg,inst_name,rc=status)
             VERIFY_(STATUS)

             found = .false.
             do j=1,size(GI)
                ! Search for name in list of active instances
                if (index(GI(j)%name,trim(inst_name)) .gt. 0) then
                   found = .true.
! ### What to do if found?
! ### outcome option. Die?
!                if (MAPL_am_I_root()) write(*,*) 'RRG: '//trim(species)//' passive instance entry already declared as an active instance'
!                VERIFY_(-1)
! ### outcome option. Toggle from active to passive
                   GI(j)%active = .false.
                endif
             enddo
             if (.not. found) then
! ### outcome option. Die if it's not already a declared instance.
!                if (MAPL_am_I_root()) write(*,*) 'RRG: '//trim(species)//' passive instance entry not in list of instances'
!                VERIFY_(-1)
! ### outcome option. Register as a new instance
                call Util_AddInstance( GI, trim(inst_name), trim(species), MW, .false., status)
                VERIFY_(STATUS)
                call RegisterInstanceWithMAPL( GC, trim(species), trim(inst_name), rc=status )
                VERIFY_(STATUS)
             endif

          end do
       endif
    endif

    return

  end subroutine ProcessInstances

end module RRG_GridCompMod
