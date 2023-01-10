!=============================================================================
!BOP

! !INTERFACE:
program simplebox

  ! !USES:
  use types_mod
  use global_mod
  use utils_mod
  use surface_mod
  use global_mod
  use integration_mod
  use M_memutils

  ! Species chemistry modules
  use  COchem_mod
  use CO2chem_mod
  use CH4chem_mod

  implicit none

  integer, parameter             :: DP=kind(1.0d0)
  integer                        :: status ! module-wide scope

  integer, save                  :: nCH4, nCO2, nCO ! number of instances
  type(gas_instance), pointer    :: CO2(:)
  type(gas_instance), pointer    ::  CO(:)
  type(gas_instance), pointer    :: CH4(:)

  integer                              :: i,j,k,n
  logical                              :: present, found
  character (len=255)                  :: COMP_NAME
  character (len=255)                  :: name, photFile
  integer                              :: dims(3)

  real, pointer                     :: CO2_total(:,:,:), CH4_total(:,:,:), CO_total(:,:,:)
  logical, save                     :: first = .true. ! I don't like using this but it has to happen. ExtData doesn't fill masks until run() and I don't want to repeat operations

  real, pointer                   :: ptr2d(:,:), ptr3d(:,:,:)

  real, pointer, dimension(:,:,:) :: testdata ! test data

  real, target, allocatable, dimension(:,:,:) :: CO2total, COtotal, CH4total

  integer                         :: timestep, nsteps

  ! Add new active instance (assume active)
  call Util_AddInstance( CO,  'test',     'CO', 28.0104, .true., status)
  call Util_AddInstance( CO,  'residual', 'CO', 28.0104, .true., status)

  !  Carbon dioxide (CO2)
  !  CO2: define the total/aggregate field
  call Util_AddInstance( CO2, 'test',     'CO2', 44.0,   .true., status)
  call Util_AddInstance( CO2, 'residual', 'CO2', 44.0,   .true., status)

  !  Methane (CH4)
  !  CH4: define the total/aggregate field
  call Util_AddInstance( CH4, 'test',     'CH4', 16.043, .true., status)
  call Util_AddInstance( CH4, 'residual', 'CH4', 16.043, .true., status)


  !****************************************************************************
  params%im = 1000
  params%jm = 1000
  params%km = 40

  params%CDT = 600
  params%HDT = 600

  params%AVO   = 6.022E23
  params%AIRMW = 28.97e0

  nsteps = 2

!  call util_dumpinstances()
  
  ! For testing, allocate the test data

  allocate(testdata(params%im,params%jm,params%km))
  allocate(CO2total,COtotal,CH4total, mold=testdata)

  testdata = 100.

  do timestep=1,nsteps
     call M_memused( memtotal, memused, percent_used, status)
     write(*,*) 'Mem total:', memtotal, ' Mem used: ', memused, ' % used: ', percent_used
  ! ===============================================================
  !              G E T  D A T A  P O I N T E R S
  !   Associate the met fields
  !   -----------------------------------
     do i=1,NINSTANCES
        instances(i)%p%data3d => testdata
        allocate( instances(i)%p%prod, instances(i)%p%loss, mold=instances(i)%p%data3d, stat=status ) ! allocate the prod/loss arrays for each instance
        instances(i)%p%prod = 0.e0; instances(i)%p%loss = 0.e0
     enddo

  !   Get pointers to the aggregates/totals
!  CO2_total => CO2total
!  aggregate(ispecies('CO2'))%q => CO2_total ! Aggregate is used under the hood

!  CO_total => COtotal
!  aggregate(ispecies('CO'))%q  => CO_total  ! Aggregate is used under the hood

!  CH4_total => CH4total
!  aggregate(ispecies('CH4'))%q => CH4_total ! Aggregate is used under the hood
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !                R U N  T H E  O P E R A T I O N S
!>>  !   Aggregate instances into the totals prior to operations
!>>  call util_aggregate( RC )
!>>
!>>  !   -- surface fluxes for all instances
!>>  call surface_prodloss( RC )
!>>
!>>  !   -- integration
!>>  call integrate_forwardeuler( RC )
!>>
!>>  !   -- post processing
!>>  if (cntrl%strictMassBalance) call util_accumulatenegatives( RC )
!>>
!>>  !   Aggregate instances into the totals after operations
!>>  call util_aggregate( RC )
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !      C O M P U T E  A N D  P A S S  D I A G N O S T I C S
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !                            D O N E
!>>  !   Cleanup
  do i=1,NINSTANCES
     ! deallocate the prod/loss arrays for each instance
     nullify( instances(i)%p%prod, instances(i)%p%loss )
  enddo
!>>
!>>  do i=1,NINSTANCES
!>>     call MAPL_GetPointer(internal, instances(i)%p%data3d, trim(instances(i)%p%species)//'_'//trim(instances(i)%p%name), __RC__)
!>>     allocate( instances(i)%p%prod, instances(i)%p%loss, mold=instances(i)%p%data3d, stat=status ) ! allocate the prod/loss arrays for each instance
!>>     instances(i)%p%prod = 0.e0; instances(i)%p%loss = 0.e0
!>>  enddo
!>>
!>>  !   Get pointers to the aggregates/totals
!>>  call MAPL_GetPointer(internal, CO2_total, 'CO2_total',__RC__)
!>>  aggregate(ispecies('CO2'))%q => CO2_total ! Aggregate is used under the hood
!>>
!>>  call MAPL_GetPointer(internal,  CO_total,  'CO_total',__RC__) 
!>>  aggregate(ispecies('CO'))%q  => CO_total  ! Aggregate is used under the hood
!>>
!>>  call MAPL_GetPointer(internal, CH4_total, 'CH4_total',__RC__) 
!>>  aggregate(ispecies('CH4'))%q => CH4_total ! Aggregate is used under the hood
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !                R U N  T H E  O P E R A T I O N S
!>>  !   Aggregate instances into the totals prior to operations
!>>  call util_aggregate( RC )
!>>
!>>  !   Compute prod/loss & integrate
!>>  !   -- each species' chemistry
!>>  !   -- CURRENTLY: OH, O1D and Cl are in mcl/cm3
!>>  call  CO_prodloss(  CO, OH, O1D, Cl, CO2photj, CH4photj, CH4_total, CO2_total, RC )
!>>  call CO2_prodloss(                                                             RC ) ! Currently nothing in here. Just in case... 
!>>  call CH4_prodloss( CH4, OH, O1D, Cl, CH4photj,                                 RC )
!>>
!>>  !   -- integration
!>>  call integrate_forwardeuler( RC )
!>>
!>>  !   -- post processing
!>>  if (cntrl%strictMassBalance) call util_accumulatenegatives( RC )
!>>
!>>  !   Aggregate instances into the totals after operations
!>>  call util_aggregate( RC )
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !      C O M P U T E  A N D  P A S S  D I A G N O S T I C S
!>>  ! ===============================================================
!>>
!>>  ! ===============================================================
!>>  !                            D O N E
!>>  !   Cleanup
!>>  do i=1,NINSTANCES
!>>     ! deallocate the prod/loss arrays for each instance
!>>     nullify( instances(i)%p%prod, instances(i)%p%loss )
!>>  enddo
!>>  deallocate(met%cosz, met%slr, O3col, O2col, CO2photj, CH4photj, stat=status )
  end do
  !============================================================================

contains

  !============================================================================

end program simplebox

