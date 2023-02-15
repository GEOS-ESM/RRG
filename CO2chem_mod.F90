module  CO2chem_mod

  implicit none
  
! !PUBLIC TYPES:
!
  PRIVATE
  PUBLIC  CO2_prodloss

CONTAINS

  subroutine CO2_prodloss( CO2inst, OH, CO, RC )

! !USES:
    use types_mod
    use global_mod
    use utils_mod
    
    implicit none

!   !INPUT PARAMETERS:
    type(gas_instance), pointer, intent(inout) :: CO2inst(:)

    ! Reactants
    real, pointer, intent(in)     ::  OH(:,:,:)
    real, pointer, intent(in)     ::  CO(:,:,:)

! !OUTPUT PARAMETERS:
    integer, intent(out) ::  rc
    INTEGER :: im, jm, km, ispc, nst
    INTEGER :: i, j, k, n

    !   Chem parameters
    real, allocatable  :: cvfac(:,:,:) ! Conversion factor
    real, allocatable  :: k_(:,:,:)    ! Rate constant
    real, pointer      :: prod(:,:,:), loss(:,:,:), CO2(:,:,:)

    integer :: STATUS

    !  Initialize local variables
    !  --------------------------
    rc = 0

    if (.not. associated(CO2inst) .or. size(CO2inst) .eq. 0) return ! Nothing to do

    im = params%im
    jm = params%jm
    km = params%km

    ispc = ispecies('CO2')

    !  Chemistry
    !  ASSUMPTION: all species units are input mol/mol
    !  --------------------------------------------------------

    allocate(k_(im,jm,km),    stat=RC)
    allocate(cvfac(im,jm,km), stat=RC)
    cvfac =  1e-3*met%rho*params%avo/params%airmw ! mol/mol <-> molec/cm3
!    cvfac =  1e-3*params%AVO*met%rho/28.0104e0 ! kg/kg <-> molec/cm3

    ! Process the prod terms into the residual
    nst = size(CO2inst)
    prod => CO2inst(nst)%prod(:,:,:)

    !  Production due to CO oxidation
    !  -------------------------------
    !            CO + OH -> CO2 + ... 
    do k=1,km
       k_(:,:,k) = 1.50E-13*(1.00+0.60E-05*(met%ple(:,:,k)+met%ple(:,:,k-1))*0.5e0) ! 2nd order (cm3/mcl/s): Where does this come from?
    enddo
    ! CO is in mol/mol.
!    prod = prod + k_*CO*OH*cvfac

    prod => null()

    !  Housekeeping
    !  ------------
    deallocate(cvfac, stat=RC)
    deallocate(k_, stat=RC)

    return
 end subroutine CO2_prodloss

end module CO2chem_mod
