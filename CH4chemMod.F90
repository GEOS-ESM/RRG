module  CH4chem_mod

  ! USES:

  implicit none

  private
  public  CH4_prodloss

  !
  ! !DESCRIPTION:
  !
  ! !REVISION HISTORY:
  !
  !
  !EOP
  !-------------------------------------------------------------------------

contains

  subroutine CH4_prodloss( CH4inst, OH, O1D, Cl, JV, RC )

    use global_mod
    use types_mod

    implicit none

    ! INPUT PARAMETERS:

    type(gas_instance), pointer, intent(inout)  :: CH4inst(:)

    ! Reactants
    real, pointer, intent(in)     ::  OH(:,:,:)
    real, pointer, intent(in)     :: O1D(:,:,:)
    real, pointer, intent(in)     ::  Cl(:,:,:)
    real,          intent(in)     ::  JV(:,:,:)

    ! OUTPUT PARAMETERS:

    integer, intent(out)                ::  rc

    ! DESCRIPTION: This routine implements the CH4 chemistry
    !
    ! !REVISION HISTORY:
    !
    !EOP
    !-------------------------------------------------------------------------

    real, pointer, dimension(:,:)   :: regionMask => null()

    INTEGER :: im, jm, km, nst, ios, idiag, iXj, ispc
    INTEGER :: i, j, k, n

    !   Chem parameters
    real, allocatable  :: cvfac(:,:,:) ! Conversion factor
    real, allocatable  :: k_(:,:,:)    ! Rate constant
    real, pointer      :: prod(:,:,:), loss(:,:,:), CH4(:,:,:)

    !  Initialize local variables
    !  --------------------------
    rc = 0

    if (.not. associated(CH4inst) .or. size(CH4inst) .eq. 0) return ! Nothing to do

    im = params%im
    jm = params%jm
    km = params%km

    !  Chemistry
    !  --------------------------------------------------------
    allocate(k_(im,jm,km), stat=RC)
    allocate(cvfac(im,jm,km), stat=RC)
    cvfac =  1e-3*met%rho*params%avo/params%airmw ! mol/mol <-> molec/cm3

    do nst = 1,size(CH4inst) ! cycle over instances
       loss  => CH4inst(nst)%loss(:,:,:) ! in mol/mol/s
       CH4   => CH4inst(nst)%data3d(:,:,:) ! CH4 pointer makes the code cleaner

       !  Loss rate [m^3 s^-1] for OH + CH4 => CO + products
       k_  = 2.45E-12*exp(-1775./met%T)
       loss = loss + k_*CH4*OH!*cvfac

       !  Loss rate [m^3 s^-1] for Cl + CH4 => CO + products
       k_  = 7.10E-12*exp(-1270./met%T)
       loss = loss + k_*CH4*Cl!*cvfac

       !  Loss rate [m^3 s^-1] for O(1D) + CH4 => CO + products
       k_ = 1.75E-10
       loss = loss + k_*CH4*O1D!*cvfac

       ! Loss rate from photolysis: CH4 + hv => ...
!       loss = loss + JV*CH4

       CH4  => null()
       loss => null()
    enddo

    !  Housekeeping
    !  ------------
    deallocate(cvfac, stat=RC)
    deallocate(k_, stat=RC)
    return

  end subroutine CH4_prodloss

end module CH4chem_mod
