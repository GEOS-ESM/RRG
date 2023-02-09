MODULE  COchem_mod

  ! !USES:

  IMPLICIT NONE

  ! !PUBLIC TYPES:
  !
  PRIVATE
  PUBLIC CO_prodloss

  !
  ! !DESCRIPTION:
  !
  !  This module implements the (pre-ESMF) CO Grid Component. 
  !
  ! !REVISION HISTORY:
  !
  !EOP
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
  !-------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:  CO_ops : ...
  !
  ! !INTERFACE:
  !

  SUBROUTINE CO_prodloss( COinst, OH, O1D, Cl, JV1, JV2, CH4, CO2, RC )

    !   !USES:
    use types_mod
    use global_mod
    use utils_mod

    IMPLICIT NONE

    !   !INPUT PARAMETERS:
    type(gas_instance), pointer, intent(inout) :: COinst(:)

    ! Reactants
    real, pointer, intent(in)     ::  OH(:,:,:)
    real, pointer, intent(in)     :: O1D(:,:,:)
    real, pointer, intent(in)     ::  Cl(:,:,:)
    real,          intent(in)     :: JV1(:,:,:)
    real,          intent(in)     :: JV2(:,:,:)
    real, pointer, intent(in)     :: CH4(:,:,:)
    real, pointer, intent(in)     :: CO2(:,:,:)

    !   !OUTPUT PARAMETERS:
    integer,          intent(out) :: RC

    ! !REVISION HISTORY:
    !
    !  Dec22, 2022 : M.S. Long - first crack. Adapted from GOCART's CO_GricCompMod.F90
    !                            CVS tag 'Tbw_Heracles-5_4_p3_SLES12' (BW = Brad Weir)
    !
    !EOP
    !-------------------------------------------------------------------------

    INTEGER :: im, jm, km, nst, ios, idiag, iXj, ispc
    INTEGER :: i, j, k, n

    !   Chem parameters
    real, allocatable  :: cvfac(:,:,:) ! Conversion factor from kg/kg -> mcl/cm3
    real, allocatable  :: k_(:,:,:)    ! Rate constant
    real, pointer      :: prod(:,:,:), loss(:,:,:), CO(:,:,:)

    integer :: STATUS

    !  Initialize local variables
    !  --------------------------
    rc = 0

    if (.not. associated(COinst) .or. size(COinst) .eq. 0) return ! Nothing to do

    im = params%im
    jm = params%jm
    km = params%km

    iXj = im * jm 

    ispc = ispecies('CO')

    !  Chemistry
    !  ASSUMPTION: all species units are input kg/kg
    !  CVFAC converts between units kg/kg <--> mcl/cm3
    !  --------------------------------------------------------

    allocate(k_(im,jm,km),    stat=RC)
    allocate(cvfac(im,jm,km), stat=RC)
    cvfac =  1e-3*met%rho*params%avo/params%airmw ! mol/mol <-> molec/cm3
!    cvfac =  1e-3*params%AVO*met%rho/28.0104e0 ! kg/kg <-> molec/cm3

    do nst = 1,size(COinst) ! cycle over instances
       loss  => COinst(nst)%loss(:,:,:) ! in kg/kg/s
       CO    => COinst(nst)%data3d(:,:,:) ! CO pointer makes the code cleaner
       
       !            CO + OH -> ...
       !  --------------
       ! the following uses the midpoint level pressure derived from PLE
       do k=1,km
          k_(:,:,k) = 1.50E-13!*(1.00+0.60E-05*(met%ple(:,:,k)+met%ple(:,:,k-1))*0.5e0) ! 2nd order (cm3/mcl/s): Where does this come from?
       enddo
       loss = loss + k_*CO*OH*cvfac
       
       CO   => null()
       loss => null()
    enddo

    ! Process the prod terms into the residual
    nst = size(COinst)
    prod => COinst(nst)%prod(:,:,:)

    !  Production due to CH4 oxidation
    !  -------------------------------
    !            CH4 + OH -> CO + ... 
    k_ = 2.45e-12*exp(-1775./met%t) ! 2nd order
    prod = prod + k_*CH4*OH*cvfac*28.0104/params%AirMW

    !            CH4 + Cl -> CO + ...
    k_ = 7.10e-12*exp(-1270./met%t) ! 2nd order
    prod = prod + k_*CH4*Cl*cvfac*28.0104/params%AirMW

    !            CH4 + O1D -> CO + ...
    k_ = 1.75e-10 ! 2nd order
    prod = prod + k_*CH4*O1D*cvfac*28.0104/params%AirMW

    !            CO2 + hv -> CO + O3P
    !            CH4 + hv -> 2H2O + CO + ... there is a bunch of branching in this. We're assuming 100% CO yield. Is this OK?
    !  ----------------------------------------------------------------------------
!    prod = prod + JV1*607.76522e-6*28.0104/44.0098!<<>>CO2    ! 1st order (1/s)
    prod = prod + JV1*CO2*28.0104/44.0098    ! 1st order (1/s)
!    prod = prod + JV2*CH4*28.0104/16.0422    ! 1st order (1/s)

    prod => null()

    ! Surface and imported fluxes are computed externally
    ! by routine, Surface_ProdLoss(), and are not species
    ! specific. MSL
    ! -------------------------------------------------
 
    !  Housekeeping
    !  ------------
    deallocate(cvfac, stat=RC)
    deallocate(k_, stat=RC)
    RETURN

!  CONTAINS
  END SUBROUTINE CO_prodloss

  !-------------------------------------------------------------------------

END MODULE COchem_mod
