MODULE  CO_mod

  ! !USES:

  IMPLICIT NONE

  ! !PUBLIC TYPES:
  !
  PRIVATE

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

  SUBROUTINE CO_ops( params, met, OH, O1D, Cl, CH4, RC )

    !   !USES:
    use types_mod
    use global_mod, only : instances

    IMPLICIT NONE

    !   !INPUT PARAMETERS:
    type(parameters), intent(in)  :: params
    type(meteorology), intent(in) :: met

    real, pointer, intent(in)     ::  OH(:,:,:)
    real, pointer, intent(in)     :: O1D(:,:,:)
    real, pointer, intent(in)     ::  Cl(:,:,:)
    real, pointer, intent(in)     :: CH4(:,:,:)

    !   !OUTPUT PARAMETERS:
    integer,          intent(out) :: RC

    ! !REVISION HISTORY:
    !
    !  Dec22, 2022 : M.S. Long - first crack. Adapted from GOCART's CO_GricCompMod.F90
    !                            CVS tag 'Tbw_Heracles-5_4_p3_SLES12' (BW = Brad Weir)
    !
    !EOP
    !-------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: myname = 'CO_GridCompRun'
    CHARACTER(LEN=*), PARAMETER :: Iam = myname

    real, pointer, dimension(:,:)   :: regionMask => null()

    INTEGER :: im, jm, km, nst, ios, idiag, iXj
    INTEGER :: i, j, k, kReverse, n
    INTEGER :: nymd1, nhms1, ier(8)

    REAL    :: qmin, qmax

    !   Chem parameters
    integer, parameter :: nreacs = 4 ! number of reactions. Increment this if adding new rxns!
    real               :: rconst(nreacs) ! rate constant
    real, allocatable  :: pe(:,:,:), p(:,:,:), ndwet(:,:,:)
    real               :: cvfac ! Conversion factor from kg/kg -> mcl/cm3
    real               :: k_, tmp ! workspace
    real, pointer      :: prod, loss, c

    !  Photolysis (bweir: from StratChem, but aj is SINGLE)
    !  ----------
    REAL, ALLOCATABLE :: photJ(:,:,:), dCOPhot(:,:,:)
    REAL, ALLOCATABLE :: aj(:)
    REAL    :: szan

    real, pointer, dimension(:,:,:) :: ptr3d => null()
    real, pointer, dimension(:,:)   :: ptr2d => null()

    integer :: STATUS

    !  Initialize local variables
    !  --------------------------
    rc = 0
    im = params%im
    jm = params%jm
    km = params%km

    iXj = im * jm 

    !<<>>!   IF(gcCO%DBG) THEN
    !<<>>!      CALL pmaxmin('CO: eCO_bioburn', gcCO%eCO_bioburn, qmin, qmax, iXj,1, 1. )
    !<<>>!      CALL pmaxmin('CO: eCO_biofuel', gcCO%eCO_biofuel, qmin, qmax, iXj,1, 1. )
    !<<>>!      CALL pmaxmin('CO: eCO_fosfuel', gcCO%eCO_fosfuel, qmin, qmax, iXj,1, 1. )
    !<<>>!      CALL pmaxmin('CO: eCO_iso',     gcCO%eCO_iso,     qmin, qmax, iXj,1, 1. )
    !<<>>!      CALL pmaxmin('CO: eCO_mon',     gcCO%eCO_mon,     qmin, qmax, iXj,1, 1. )
    !<<>>!      CALL pmaxmin('CO: eCO_mtn',     gcCO%eCO_mtn,     qmin, qmax, iXj,1, 1. )
    !<<>>!   END IF

    !  Allocate temporary workspace
    !  ----------------------------
    !<<>>    if ( ios /= 0 ) then
    !<<>>       rc = 3
    !<<>>       return
    !<<>>    end if
    !  Wet-air number density
    !  ----------------------
    !<<>>    ndwet(1:im,1:jm,1:km) = rhowet(1:im,1:jm,1:km)*params%AVOGAD/params%AIRMW

    !  Chemistry
    !  We assume all units are in kg/kg
    !  --------------------------------------------------------

    !<<>>! Precompute the photolysis rate
    !<<>>!  photJ = 0.
    !<<>>!  call getJRates(status)
    !<<>>!  VERIFY_(status)

    ! Use do loop to minimize the need for 3d array allocation & deallocation
    do k=1,km
       do j=1,jm
          do i=1,im
             do nst = 1,size(instances) ! cycle over instances
                cvfac =  1e-3*params%AVO*met%rho(i,j,k)/28.0104e0 ! kg/kg -> molec/cm3
                prod  => instances(nst)%p%prod(i,j,k) !
                loss  => instances(nst)%p%loss(i,j,k) !
                c     => instances(nst)%p%data3d(i,j,k) ! c is used for cleanliness below

                ! Compute P-L
                ! ASSUMPTION: All species are in kg/kg.
                ! - CVFAC converts between units kg/kg <--> mcl/cm3

                !  Loss due to OH
                !  --------------
                ! the following uses the midpoint level pressure derived from PLE
                k_ = 1.50E-13*(1.00+0.60E-05*(met%ple(i,j,k)+met%ple(i,j,k-1))*0.5e0) ! 2nd order rate constant (cm3/mcl/s): Where does this come from?
                loss = loss + k_*c*OH(i,j,k)*cvfac ! rate

                !  Production due to CH4 oxidation
                !  -------------------------------
                !            CH4 + OH -> CO + ... 
                k_ = 2.45E-12*1.00E-06*exp(-1775./met%t(i,j,k))
                prod = prod + k_*CH4(i,j,k)*OH(i,j,k)*cvfac

                !            CH4 + Cl -> CO + ...
                k_ = 7.10E-12*1.00E-06*exp(-1270./met%t(i,j,k))
                prod = prod + k_*CH4(i,j,k)*Cl(i,j,k)*cvfac

                !            CH4 + O1D -> CO + ...
                k_ = 1.75e-10
                prod = prod + k_*CH4(i,j,k)*O1D(i,j,k)*cvfac

                !<<>>   if (associated(CO_prod)) then
                !<<>>      CO_prod(1:im,1:jm) = 0.
                !<<>>      do k = 1,km
                !<<>>         CO_prod(1:im,1:jm) = CO_prod(1:im,1:jm)                                    &
                !<<>>                              + (    rkch4_oh(1:im,1:jm,k)* gcCO%OHnd(1:im,1:jm,k)  &
                !<<>>                                  +  rkch4_cl(1:im,1:jm,k)* gcCO%Clnd(1:im,1:jm,k)  &
                !<<>>                                  + rkch4_o1d(1:im,1:jm,k)*gcCO%O1Dnd(1:im,1:jm,k)) &
                !<<>>                                * gcCO%CH4(1:im,1:jm,k)*mwtCO/MAPL_AIRMW              &
                !<<>>                                * w_c%delp(1:im,1:jm,k)/MAPL_GRAV
                !<<>>      enddo
                !<<>>   endif


                !  Calculate photolytic loss rates, J [s^-1] for
                !     CH4 + hv => 2H2O + CO
                !     CO2 + hv => CO + ???
                !  Notice that J and the losses are always computed. However, the setting 
                !  of the feedback switch(es) determines if the increments are actually applied
                !  ----------------------------------------------------------------------------

                !<<>>   ALLOCATE(photJ(1:im,1:jm,1:km), dCOPhot(1:im,1:jm,1:km), STAT=status)
                !<<>>   VERIFY_(STATUS)
                !<<>>   dCOPhot = 0.

                !  Change in CO number density [m^-3 s^-1] due to CH4 photolysis
                !<<>>!  -------------------------------------------------------------
                !<<>>!
                !<<>>!  dCOPhot = photJ*gcCO%CH4(1:im,1:jm,1:km)
                !<<>>
                !<<>>!  Change in CO number density [m^-3 s^-1] due to CO2 photolysis
                !<<>>!  -------------------------------------------------------------
                !<<>>   if (gcCO%numphoto > 0) then
                !<<>>      photJ = 0.
                !<<>>      allocate(aj(gcCO%numphoto), STAT=status)
                !<<>>      VERIFY_(STATUS)
                !<<>>
                !<<>>      do k = 1,km
                !<<>>         do j = j1,j2
                !<<>>            do i = i1,i2
                !<<>>               szan = 0.
                !<<>>               if (w_c%cosz(i,j) <= 1.) szan = acos(w_c%cosz(i,j))
                !<<>>!              bweir: Using 0 for O3 (FIXME)
                !<<>>               call jcalc4(km-k+1, szan, 0., p(i,j,k), t(i,j,k), aj, gcCO)
                !<<>>               photJ(i,j,k) = aj(12)
                !<<>>            enddo
                !<<>>         enddo
                !<<>>      enddo
                !<<>>   endif
                !<<>>
                !<<>>!  dCOPhot = photJ*gcCO%CO2(1:im,1:jm,1:km)
                !<<>>!  bweir: Using 400 ppm for CO2 (FIXME)
                !<<>>   dCOPhot = photJ*400.e-6
                !<<>>
                !<<>>!  Photolysis (bweir: why are we multiplying by ndwet?)
                !<<>>!  ----------
                !<<>>   if (associated(CO_phot)) THEN
                !<<>>      CO_phot(1:im,1:jm,1:km) = dCOPhot(1:im,1:jm,1:km)*ndwet(1:im,1:jm,1:km)
                !<<>>   endif
                !<<>>
                nullify(c,prod,loss)
             enddo
          enddo
       enddo
    enddo
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
    !<<>>!  CO Surface Emission Flux in kg m-2 s-1
    !<<>>!  --------------------------------------
    !<<>>   if (associated(CO_emis)) CO_emis(1:im,1:jm) = gcCO%COsfcFlux(1:im,1:jm)
    !<<>>
    !<<>>!   IF(gcCO%DBG) THEN
    !<<>>!     if (associated(CO_emis))    call pmaxmin('CO: emis',       CO_emis, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_loss))    call pmaxmin('CO: loss',       CO_loss, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_prod))    call pmaxmin('CO: prod',       CO_prod, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_phot))    call pmaxmin('CO: phot',       CO_phot, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_column))  call pmaxmin('CO: column',   CO_column, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_surface)) call pmaxmin('CO: surface', CO_surface, qmin, qmax, iXj,  1, 1. )
    !<<>>!     if (associated(CO_dry))     call pmaxmin('CO: dry',         CO_dry, qmin, qmax, iXj, km, 1. )
    !<<>>!   END IF
    !<<>>
    !<<>>!  Housekeeping
    !<<>>!  ------------
    DEALLOCATE(ndwet,p,pe,STAT=ier(1))
    DEALLOCATE(photJ, dCOPhot, aj, STAT=ier(1))

    RETURN

!  CONTAINS
  END SUBROUTINE CO_ops

  !-------------------------------------------------------------------------

END MODULE CO_mod
