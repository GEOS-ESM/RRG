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

  SUBROUTINE CO_prodloss( COinst, OH, O1D, Cl, CH4, CO2, RC )

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

    real, pointer, dimension(:,:)   :: regionMask => null()

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

    ! currently unused, since oxidants import as mcl/cm3 ... allocate(cvfac(im,jm,km),k_(im,jm,km), stat=RC)
    ! currently unused, since oxidants import as mcl/cm3 ... cvfac =  1e-3*params%AVO*met%rho/28.0104e0 ! kg/kg <-> molec/cm3
    allocate(k_(im,jm,km), stat=RC)
    do nst = 1,size(COinst) ! cycle over instances
       loss  => COinst(nst)%loss(:,:,:) ! in kg/kg/s
       CO    => COinst(nst)%data3d(:,:,:) ! CO pointer makes the code cleaner
       
       !            CO + OH -> ...
       !  --------------
       ! the following uses the midpoint level pressure derived from PLE
       do k=1,km
          k_(:,:,k) = 1.50E-13*(1.00+0.60E-05*(met%ple(:,:,k)+met%ple(:,:,k-1))*0.5e0) ! 2nd order (cm3/mcl/s): Where does this come from?
       enddo
       loss = loss + k_*CO*OH
       
       CO   => null()
       loss => null()
    enddo

    ! Process the prod terms into the residual
    nst = size(COinst)
    prod => COinst(nst)%prod(:,:,:)

    !  Production due to CH4 oxidation
    !  -------------------------------
    !            CH4 + OH -> CO + ... 
    k_ = 2.45E-12*1.00E-06*exp(-1775./met%t) ! 2nd order
    prod = prod + k_*CH4*OH 

    !            CH4 + Cl -> CO + ...
    k_ = 7.10E-12*1.00E-06*exp(-1270./met%t) ! 2nd order
    prod = prod + k_*CH4*Cl

    !            CH4 + O1D -> CO + ...
    k_ = 1.75e-10 ! 2nd orders
    prod = prod + k_*CH4*O1D

    !            CH4 + hv -> 2H2O + CO
    !            CO2 + hv -> CO + ???
    !     1st order (1/s); don't need CVFAC
    !  ----------------------------------------------------------------------------
    !prod = prod + met%photJ*CH4 <-- what to do about this? I don't think it's real <<>> MSL
    prod = prod + met%photJ*CO2

    prod => null()

    ! Surface and imported fluxes are computed externally
    ! by routine, Surface_ProdLoss(), and are not species
    ! specific. MSL
    ! -------------------------------------------------
    

    !<<>>!  Surface concentration [ppbv]
    !<<>>!  ----------------------------
    !<<>>   if (associated(CO_surface)) then
    !<<>>      CO_surface(1:im,1:jm) = w_c%qa(nbeg)%data3d(1:im,1:jm,km)*1.e9
    !<<>>   endif
    !<<>>
    !<<>>!  Column burden [kg m-2]
    !<<>>!  ----------------------
    !<<>>   if (associated(CO_column)) then
    !<<>>      CO_column(1:im,1:jm) = 0.
    !<<>>      do k = 1, km
    !<<>>         CO_column(1:im,1:jm) = CO_column(1:im,1:jm)                                &
    !<<>>                                + w_c%qa(nbeg)%data3d(1:im,1:jm,k)*mwtCO/MAPL_AIRMW * &
    !<<>>                                             w_c%delp(1:im,1:jm,k)/MAPL_GRAV
    !<<>>     enddo
    !<<>>   endif
    !<<>>
    !<<>>!  Dry-air mole fraction
    !<<>>!  ---------------------
    !<<>>   if (associated(CO_dry)) then
    !<<>>      CO_dry(1:im,1:jm,1:km) = w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) &
    !<<>>                                        / (1. - qtot(1:im,1:jm,1:km))
    !<<>>   endif
    !<<>>

    !  Housekeeping
    !  ------------
    ! currently unused, since oxidants import as mcl/cm3 ... deallocate(cvfac, stat=RC)
    deallocate(k_, stat=RC)
    RETURN

!  CONTAINS
  END SUBROUTINE CO_prodloss

  !-------------------------------------------------------------------------

END MODULE COchem_mod
