MODULE  CO_mod

  ! !USES:

  !  bweir: for photolysis
  !<<>>   USE ESMF_CFIOFileMOD
  !<<>>   USE MAPL_CFIOMOD

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

  SUBROUTINE CO_ops( params, RC )

!   !USES:
    use types_mod

    IMPLICIT NONE

!   !INPUT PARAMETERS:
    type(parameters), intent(in)  :: params

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

    INTEGER :: im, jm, km, ios, idiag, iXj
    INTEGER :: i, j, k, kReverse, n
    INTEGER :: nymd1, nhms1, ier(8)

    REAL    :: qmin, qmax

!   Chem parameters
    integer, parameter :: nreacs = 4 ! number of reactions. Increment this if adding new rxns!
    real               :: rconst(nreacs) ! rate constant
    real, allocatable  :: pe(:,:,:), p(:,:,:), ndwet(:,:,:)
    real, allocatable  :: dCdt(:,:,:)
    real, allocatable  :: cvfac(:,:,:) ! Conversion factor from kg/kg -> mcl/cm3

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
    allocate(pe(im,jm,km+1), p(im,jm,km), ndwet(im,jm,km), stat=ios)
    if ( ios /= 0 ) then
       rc = 3
       return
    end if
!  Wet-air number density
!  ----------------------
!<<>>    ndwet(1:im,1:jm,1:km) = rhowet(1:im,1:jm,1:km)*params%AVOGAD/params%AIRMW

!  Chemistry
!  We assume all units are in kg/kg
!  --------------------------------------------------------

!  Loss due to OH
!  --------------
!<<>>   rkoh(1:im,1:jm,1:km) = 1.50E-13*1.00E-06*(1.00+0.60E-05*p(1:im,1:jm,1:km))

!  Production due to CH4 oxidation
!  -------------------------------
!<<>>   rkch4_oh(1:im,1:jm,1:km)  = 2.45E-12*1.00E-06*exp(-1775./t(1:im,1:jm,1:km))
!<<>>   rkch4_cl(1:im,1:jm,1:km)  = 7.10E-12*1.00E-06*exp(-1270./t(1:im,1:jm,1:km))
!<<>>   rkch4_o1d(1:im,1:jm,1:km) = 1.75E-10*1.00E-06
!<<>>
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
!<<>>
!<<>>!  Calculate photolytic loss rates, J [s^-1] for
!<<>>!     CH4 + hv => 2H2O + CO
!<<>>!     CO2 + hv => CO + ???
!<<>>!  Notice that J and the losses are always computed. However, the setting 
!<<>>!  of the feedback switch(es) determines if the increments are actually applied
!<<>>!  ----------------------------------------------------------------------------
!<<>>   ALLOCATE(photJ(1:im,1:jm,1:km), dCOPhot(1:im,1:jm,1:km), STAT=status)
!<<>>   VERIFY_(STATUS)
!<<>>   dCOPhot = 0.
!<<>>
!<<>>!  Change in CO number density [m^-3 s^-1] due to CH4 photolysis
!<<>>!  -------------------------------------------------------------
!<<>>!  photJ = 0.
!<<>>!  call getJRates(status)
!<<>>!  VERIFY_(status)
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
!<<>>!  Decrement the CO mole fraction due to oxidation 
!<<>>!  -----------------------------------------------
!<<>>   w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) = w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) &
!<<>>         - cdt * rkoh(1:im,1:jm,1:km)*gcCO%OHnd(1:im,1:jm,1:km)              &
!<<>>               * w_c%qa(nbeg)%data3d(1:im,1:jm,1:km)
!<<>>
!<<>>!  bweir: just checking
!<<>>   w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) = max(w_c%qa(nbeg)%data3d(1:im,1:jm,1:km), 0.)
!<<>>
!<<>>!  Compute and add surface emissions
!<<>>!  ---------------------------------
!<<>>   gcCO%COsfcFlux(1:im,1:jm) = 0.
!<<>>   call CO_Emission(rc)
!<<>>
!<<>>!  Increment the CO mole fraction due to photolysis
!<<>>!  ------------------------------------------------
!<<>>   w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) = w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) + cdt*dCOPhot(1:im,1:jm,1:km)
!<<>>
!<<>>!  Increment the CO mole fraction due to production
!<<>>!  ------------------------------------------------
!<<>>   w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) = w_c%qa(nbeg)%data3d(1:im,1:jm,1:km) &
!<<>>         + cdt * (    rkch4_oh(1:im,1:jm,1:km)* gcCO%OHnd(1:im,1:jm,1:km)    &
!<<>>                   +  rkch4_cl(1:im,1:jm,1:km)* gcCO%Clnd(1:im,1:jm,1:km)    &
!<<>>                   + rkch4_o1d(1:im,1:jm,1:km)*gcCO%O1Dnd(1:im,1:jm,1:km))   &
!<<>>               * gcCO%CH4(1:im,1:jm,1:km)
!<<>>
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

  CONTAINS

!-------------------------------------------------------------------------
!  Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
    SUBROUTINE getJRates(rc)

      IMPLICIT NONE

      INTEGER, INTENT(out) :: rc

      REAL, ALLOCATABLE :: o2Column(:,:,:)
      REAL, ALLOCATABLE :: SZARad(:,:)
      REAL, ALLOCATABLE :: SZADeg(:,:)
      REAL, ALLOCATABLE :: sinSZA(:,:)
      REAL, ALLOCATABLE :: zgrz(:,:)
      REAL, ALLOCATABLE :: sfaca(:,:)
      REAL, ALLOCATABLE :: arg(:,:)

      REAL, PARAMETER :: wavel = 1215.7
      REAL, PARAMETER :: O2xs  = 1.000E-20
      REAL, PARAMETER :: CH4xs = 2.000E-17
      REAL, PARAMETER :: sflux = 4.006E+11

  !  Constants for Chapman function at high solar zenith angle
  !  ---------------------------------------------------------
      REAL, PARAMETER :: hbar = 6.79
      REAL, PARAMETER :: zbar = 30.0
      REAL, PARAMETER :: r0   = 6.371E+03
      REAL, PARAMETER :: zp   = 60.0

      REAL, PARAMETER :: d1 = 1.060693
      REAL, PARAMETER :: d2 = 0.55643831
      REAL, PARAMETER :: d3 = 1.0619896
      REAL, PARAMETER :: d4 = 1.7245609
      REAL, PARAMETER :: d5 = 0.56498823
      REAL, PARAMETER :: d6 = 0.06651874

      REAL, PARAMETER :: O2ABV80KM = 7.072926E+19 ![cm^{-2}]
      REAL, PARAMETER :: O2VMR = 0.20946

  !  bweir: Hardcoding because I don't want to brick old RC files (FIXME)
      REAL, PARAMETER :: SZACUTOFF = 70.0

      REAL :: b, r, s
      REAL :: fO2VMR

      INTEGER :: status

      CHARACTER(LEN=255) :: Iam = "CO::getJRates"

      rc = 0
      b = SQRT(0.50*r0/hbar)

  !  O2 overhead number density profile [cm^{-2}]
  !  --------------------------------------------
      ALLOCATE(O2Column(1:im,1:jm,1:km), STAT=status)

      fO2VMR = O2VMR*5.00E-05
  !  O2Column(:,:,1) = O2ABV80KM+cellDepth(:,:,1)*ndwet(:,:,1)*fO2VMR
  !<<>>   O2Column(1:im,1:jm,1) = O2ABV80KM + w_c%delp(1:im,1:jm,1)*MAPL_AVOGAD/(MAPL_AIRMW*MAPL_GRAV)*fO2VMR

      DO k = 2,km
     !     O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
     !                                          cellDepth(:,:,  k)*ndwet(:,:,  k))*fO2VMR
     !<<>>      O2Column(1:im,1:jm,k) = O2Column(1:im,1:jm,k-1) &
     !<<>>                              + (w_c%delp(1:im,1:jm,k-1) + w_c%delp(1:im,1:jm,k)) &
     !<<>>                                *MAPL_AVOGAD/(MAPL_AIRMW*MAPL_GRAV)*fO2VMR
      END DO

  !   IF(gcCO%dbg) THEN
  !      CALL pmaxmin('CO: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
  !   END IF

  !  Grab some memory
  !  ----------------
  !<<>>   ALLOCATE(SZARad(1:im,1:jm), STAT=status)
  !<<>>   ALLOCATE(SZADeg(1:im,1:jm), STAT=status)
  !<<>>   ALLOCATE(sinSZA(1:im,1:jm), STAT=status)
  !<<>>
  !<<>>   WHERE(w_c%cosz(1:im,1:jm) > 1.00)
  !<<>>      SZARad(1:im,1:jm) = 0.00
  !<<>>   ELSEWHERE
  !<<>>      SZARad(1:im,1:jm) = ACOS(w_c%cosz(1:im,1:jm))
  !<<>>   ENDWHERE
  !<<>>   SZADeg(1:im,1:jm) = SZARad(1:im,1:jm)*radToDeg
  !<<>>   sinSZA(1:im,1:jm) = SIN(SZARad(1:im,1:jm))
  !<<>>
  !<<>>   ALLOCATE(zgrz(1:im,1:jm), STAT=status)
  !<<>>
  !<<>>   WHERE(SZADeg(1:im,1:jm) <= 90.00)
  !<<>>      zgrz(1:im,1:jm) = 1000.00
  !<<>>   ELSEWHERE
  !<<>>      zgrz(1:im,1:jm) = sinSZA(1:im,1:jm)*(zp+r0)-r0
  !<<>>   ENDWHERE
  !<<>>
  !<<>>!   IF(gcCO%dbg) THEN
  !<<>>!      CALL pmaxmin('CO: zgrz',     zgrz, qmin, qmax, iXj, 1,  1. )
  !<<>>!      CALL pmaxmin('CO: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
  !<<>>!   END IF
  !<<>>
  !<<>>   ALLOCATE(sfaca(1:im,1:jm), STAT=status)
  !<<>>   sfaca(1:im,1:jm) = 0.00
  !<<>>
  !<<>>!  Chapman function calculation from ACDB 2-D model
  !<<>>!  ------------------------------------------------
  !<<>>   DO j = j1,j2
  !<<>>      DO i = i1,i2
  !<<>>         Daytime: IF(SZADeg(i,j) < SZACUTOFF) THEN
  !<<>>            IF(SZADeg(i,j) < 70.00) THEN
  !<<>>               sfaca(i,j) = 1.00/w_c%cosz(i,j)
  !<<>>            ELSE IF(zgrz(i,j) > 0.00) THEN
  !<<>>               s = b*ABS(w_c%cosz(i,j))
  !<<>>               IF(s <= 8.00) THEN
  !<<>>                  s = (d1+d2*s)/(d3+d4*s+s**2)
  !<<>>               ELSE
  !<<>>                  s = d5/(d6+s)
  !<<>>               ENDIF
  !<<>>
  !<<>>               r = b*SQRT(MAPL_PI)
  !<<>>               sfaca(i,j) = r*s
  !<<>>
  !<<>>               IF(SZADeg(i,j) > 90.00) THEN
  !<<>>                  sfaca(i,j) = 2.00*r*EXP((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
  !<<>>               ENDIF
  !<<>>            ENDIF
  !<<>>         ENDIF Daytime
  !<<>>      ENDDO
  !<<>>   ENDDO
  !<<>>
  !<<>>   IF(gcCO%dbg) THEN
  !<<>>      CALL pmaxmin('CO: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
  !<<>>   END IF
  !<<>>
  !<<>>   ALLOCATE(arg(1:im,1:jm), STAT=status)
  !<<>>
  !<<>>!  At each layer, compute the rate constant, J [s^{-1}], if the sun is up
  !<<>>!  ----------------------------------------------------------------------
  !<<>>   DO k = 1,km
  !<<>>      WHERE(SZADeg(1:im,1:jm) < SZACUTOFF)
  !<<>>         arg(1:im,1:jm) = O2Column(1:im,1:jm,k)*O2xs*sfaca(1:im,1:jm)
  !<<>>         photJ(1:im,1:jm,k) = sflux*EXP(-arg(1:im,1:jm))*CH4xs
  !<<>>      ENDWHERE
  !<<>>   ENDDO
  !<<>>
  !<<>>   IF(gcCO%dbg) THEN
  !<<>>      CALL pmaxmin('CO: photJ', photJ, qmin, qmax, iXj, km,  1. )
  !<<>>   ENDIF
  !<<>>
  !<<>>   DEALLOCATE(SZARad, STAT=status)
  !<<>>   DEALLOCATE(SZADeg, STAT=status)
  !<<>>   DEALLOCATE(sinSZA, STAT=status)
  !<<>>   DEALLOCATE(zgrz, STAT=status)
  !<<>>   DEALLOCATE(sfaca, STAT=status)
  !<<>>   DEALLOCATE(arg, STAT=status)
  !<<>>   DEALLOCATE(O2Column, STAT=status)

      RETURN
    END SUBROUTINE getJRates

    SUBROUTINE interp_s()!k,sza,o3column,s,jo2,gcCO)
  ! ----------------------------------------------------------------------------
  ! NAME:
  !   interp_s
  !
  ! PURPOSE:
  !   Interpolate S values for each wavelength in table to specified O3
  !   column and zenith angle
  !
  ! INPUTS:
  !   k         Current layer number
  !   szaRad    Solar zenith angle [radians]
  !   o3column  Overhead o3 column value [cm^{-2}]
  !   gcCO      The GOCART::CO grid component, which contains
  !     sza_tab Solar zenith angle table
  !     o3_tab  Overhead O3 values table
  !     sdat    Radiative source function 
  !     o2jdat  Table of J(O2) values
  !
  ! OUTPUTS:
  !   s         S value for each wavelength at current k, interpolated to
  !               the given o3column and sza
  !   jo2       J(O2) values interpolated as above
  !
  ! 
  ! PROCEDURE:
  !   Bi-linear interpolation, for sza > 94 s=0, for O3 out of range use min/max
  !
  ! MODIFICATION HISTORY: 
  !   25 Aug 1993  Kawa
  !   10 Jul 1996  Kawa    For 28 levels and to handle J(O2) separately
  !   11 May 2012  Nielsen Accomodation for GEOS-5 FV cubed release
  !   30 Jan 2021  Weir    Copied from StratChem
  ! ----------------------------------------------------------------------------

      IMPLICIT NONE

  !<<>>   TYPE(CO_GridComp1), INTENT(IN) :: gcCO   ! Grid Component
  !<<>>
  !<<>>   INTEGER, INTENT(IN) :: k
  !<<>>   REAL, INTENT(IN) :: sza, o3column 
  !<<>>   REAL, INTENT(OUT) :: s(gcCO%nlam), jo2

      INTEGER :: ijj, ik, ikk, ikkm, il, is
      REAL :: omt, omu, t, u
      REAL, PARAMETER :: PI = 3.14159265

  !<<>>! For each input solar zenith angle, find the first element of gcCO%sza_tab that 
  !<<>>! is greater.  Use this element and previous one to determine the interpolated value.
  !<<>>! -----------------------------------------------------------------------------------
  !<<>>   DO is = 1,gcCO%nsza
  !<<>>      ijj = is 
  !<<>>      IF(gcCO%sza_tab(is) > sza) EXIT 
  !<<>>   ENDDO
  !<<>>      
  !<<>>! Zenith angle test       
  !<<>>! -----------------
  !<<>>   IF(sza > gcCO%sza_tab(gcCO%nsza)) THEN
  !<<>>!     Cell is dark, set s and jo2=0        
  !<<>>!     -----------------------------
  !<<>>      s(1:gcCO%nlam) = 0.
  !<<>>      jo2 = 0.
  !<<>>   ELSE  
  !<<>>!     Cell is illuminated     
  !<<>>!     -------------------
  !<<>>      t = (sza-gcCO%sza_tab(ijj-1))/(gcCO%sza_tab(ijj)-gcCO%sza_tab(ijj-1))
  !<<>>      omt = 1.-t
  !<<>>         
  !<<>>! For each overhead O3 column, find the first element in gcCO%o3_tab that is
  !<<>>! greater. Use this element and previous one to determine the interpolated value.
  !<<>>! -------------------------------------------------------------------------------
  !<<>>      DO is = 1,gcCO%numo3
  !<<>>         ikk = is 
  !<<>>         IF(gcCO%o3_tab(is,k) > o3column) EXIT
  !<<>>      ENDDO
  !<<>>
  !<<>>      ikkm = ikk-1 
  !<<>>      IF(ikk > 1 .AND. o3column <= gcCO%o3_tab(gcCO%numo3,k)) THEN
  !<<>>         u = (o3column-gcCO%o3_tab(ikkm,k))/(gcCO%o3_tab(ikk,k)-gcCO%o3_tab(ikkm,k))
  !<<>>         omu = 1.-u
  !<<>>
  !<<>>! Do bilinear interpolation for each wavelength.
  !<<>>! ----------------------------------------------
  !<<>>         DO il = 1,gcCO%nlam       
  !<<>>            s(il) = omt*omu*gcCO%sdat(ijj-1,ikkm,k,il)+t*omu*gcCO%sdat(ijj,ikkm,k,il)+ &
  !<<>>                    t*u*gcCO%sdat(ijj,ikk,k,il)+omt*u*gcCO%sdat(ijj-1,ikk,k,il)
  !<<>>         ENDDO
  !<<>>         jo2 = omt*omu*gcCO%o2jdat(ijj-1,ikkm,k)+t*omu*gcCO%o2jdat(ijj,ikkm,k)+ &
  !<<>>               t*u*gcCO%o2jdat(ijj,ikk,k)+omt*u*gcCO%o2jdat(ijj-1,ikk,k)
  !<<>>    
  !<<>>! Extrapolate ahead of table
  !<<>>! --------------------------
  !<<>>      ELSE IF (ikk == 1) THEN
  !<<>>         DO il = 1,gcCO%nlam
  !<<>>            s(il) = omt*gcCO%sdat(ijj-1,1,k,il)+t*gcCO%sdat(ijj,1,k,il)
  !<<>>         ENDDO
  !<<>>         jo2 = omt*gcCO%o2jdat(ijj-1,1,k)+t*gcCO%o2jdat(ijj,1,k)
  !<<>>
  !<<>>! Extrapolate beyond table
  !<<>>! ------------------------
  !<<>>      ELSE
  !<<>>         DO il = 1,gcCO%nlam
  !<<>>            s(il) = omt*gcCO%sdat(ijj-1,gcCO%numo3,k,il)+t*gcCO%sdat(ijj,gcCO%numo3,k,il)
  !<<>>         END DO 
  !<<>>         jo2 = omt*gcCO%o2jdat(ijj-1,gcCO%numo3,k)+t*gcCO%o2jdat(ijj,gcCO%numo3,k)
  !<<>>      ENDIF  
  !<<>>   ENDIF
  !<<>>      
  !<<>>   RETURN
    END SUBROUTINE interp_s

    SUBROUTINE jcalc4()!k,szan,o3column,press,kel,aj,gcCO)
  ! ---------------------------------------------------------------------------------
  ! NAME: jcalc4
  ! PURPOSE:
  !   Calculate photolysis rates
  ! INPUT:
  !   k         Current layer number
  !   levels    Number of layers
  !   szan      Solar zenith angle (radians)
  !   o3column  Overhead O3 values
  !   press     Mid-layer pressure (hPa)
  !   kel       Mid-layer temperature (K)
  ! OUTPUT:
  !   aj        Array of photolysis rates
  ! RESTRICTIONS:
  !   Currently set up for 23-J set (see var gcCO%nxdo)
  ! REQUIRED ROUTINES:
  !   interp_s
  ! MODIFICATION HISTORY: 
  !   26 Aug 1993 Kawa    Created
  !   23 Nov 1993 Kawa    Remade xtab to do multiplication by solar flux beforehand 
  !                        and removed inputs.
  !   25 Feb 1994         Add 3 additional Js, incl N2O
  !   18 Sep 1995         Add 2 additional Js, up to 22, and do CH2O special
  !   13 May 1996 Crum    Removed fossils, move toward Fortran 90
  !   10 Jul 1996         Modified to handle J(O2) separately and use 28 levels
  !    1 Apr 2009 Nielsen GEOS-5 form with standardized SC_GridComp interface.
  !    1 Jun 2009 Nielsen Updated to JPL 2006
  !   12 Dec 2010 Nielsen Updated to JPL 2010 following Luke Oman's testing.
  !   11 May 2012 Nielsen Accomodation for GEOS-5 FV cubed release
  !    3 Jun 2015 Liang   Updated to the new 50-slot table with addition of halons,
  !                       HCFCs, and 5 VSLSs
  !                       numphoto is now updated to 52
  !   30 Jan 2021 Weir    Copied from StratChem
  !
  ! WARNING: Photolysis reaction rate numbers 38-42 are calculated in MESO_PHOT.
  ! ---------------------------------------------------------------------------------
      IMPLICIT NONE
  !<<>>   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)
  !<<>>
  !<<>>   TYPE(CO_GridComp1), INTENT(INOUT) :: gcCO   ! Grid Component
  !<<>>
  !<<>>   INTEGER, INTENT(IN) :: k
  !<<>>   REAL, INTENT(IN) :: szan, o3column, press, kel
  !<<>>!  REAL(KIND=DBL), INTENT(OUT) :: aj(gcCO%numphoto)
  !<<>>!  bweir: demoted to single
  !<<>>   REAL, INTENT(OUT) :: aj(gcCO%numphoto)
  !<<>>
  !<<>>   INTEGER :: ilam,indt,ix
  !<<>>
  !<<>>   REAL :: alpha300, alphat, jo2, rjs(gcCO%nxdo), q1, q2, r1mq1
  !<<>>   REAL :: s(gcCO%nlam), sx(2,gcCO%nlam), tfac, wvl
  !<<>>
  !<<>>! Start with a clean slate
  !<<>>! ------------------------
  !<<>>   aj(1:gcCO%numphoto) = 0.
  !<<>>
  !<<>>! Interpolate radiative flux function values to model conditions
  !<<>>! --------------------------------------------------------------
  !<<>>   CALL interp_s(k,szan,o3column,s,jo2,gcCO)
  !<<>>   indt = kel-148.5
  !<<>>   indt = MAX(1,indt)
  !<<>>   indt = MIN(indt,200)
  !<<>>
  !<<>>! Preliminaries for CH2O quantum yield dependence on m, T, wavelength
  !<<>>! -------------------------------------------------------------------
  !<<>>   tfac = (kel-80.0)/80.0
  !<<>>
  !<<>>   DO ilam=1,gcCO%nlam
  !<<>>      ZeroS: IF(s(ilam) == 0.) THEN
  !<<>>         sx(1,ilam) = 0.00
  !<<>>         sx(2,ilam) = 0.00
  !<<>>      ELSE 
  !<<>>
  !<<>>         wvl = gcCO%rlam(ilam)*0.10
  !<<>>
  !<<>>         IF(wvl < 250.00) THEN
  !<<>>            q1 = 0.24
  !<<>>         ELSE IF(wvl >= 339.00) THEN
  !<<>>            q1 = 0.00
  !<<>>         ELSE
  !<<>>            q1 = gcCO%CH2O_aq(1) + gcCO%CH2O_aq(2)*wvl         + &
  !<<>>                                   gcCO%CH2O_aq(3)*wvl*wvl     + &
  !<<>>                                   gcCO%CH2O_aq(4)*wvl*wvl*wvl + &
  !<<>>                                   gcCO%CH2O_aq(5)*wvl*wvl*wvl*wvl
  !<<>>         ENDIF
  !<<>>
  !<<>>         r1mq1 = 1./(1.-q1)
  !<<>>
  !<<>>         IF(wvl < 330.00) THEN
  !<<>>            q2 = gcCO%xtab(ilam,22,indt)
  !<<>>         ELSE IF(wvl > 360.00) THEN
  !<<>>            q2 = 0.00
  !<<>>         ELSE
  !<<>>            alpha300 = 1.00E-03*(1./gcCO%xtab(ilam,22,1)-r1mq1)
  !<<>>            alphat = alpha300*(1.+0.05*(wvl-329.)*((300.-kel)/80.))
  !<<>>            q2 = 1.00/(r1mq1+alphat*press)
  !<<>>         ENDIF
  !<<>>
  !<<>>         IF(wvl .LT. 250.00) q2=0.5
  !<<>>
  !<<>>         sx(2,ilam) = s(ilam)*gcCO%xtab(ilam,21,indt)*q2
  !<<>>         sx(1,ilam) = s(ilam)*gcCO%xtab(ilam,21,indt)*q1
  !<<>>      ENDIF ZeroS
  !<<>>   ENDDO
  !<<>>
  !<<>>! J(BrONO2) through J(OCLO)
  !<<>>! -------------------------
  !<<>>   DO ix=1,14
  !<<>>      rjs(ix) = 0.
  !<<>>
  !<<>>      DO ilam=1,gcCO%nlam
  !<<>>         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
  !<<>>      ENDDO
  !<<>>   ENDDO
  !<<>>
  !<<>>! J(O2)
  !<<>>! -----
  !<<>>   rjs(15) = jo2
  !<<>>
  !<<>>! J(O3_O1D) through J(N2O)
  !<<>>! ------------------------
  !<<>>   DO ix=16,20
  !<<>>      rjs(ix) = 0.
  !<<>>
  !<<>>      DO ilam=1,gcCO%nlam
  !<<>>         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
  !<<>>      ENDDO
  !<<>>   ENDDO
  !<<>>
  !<<>>! J(CH2O)
  !<<>>! -------
  !<<>>   rjs(21) = 0.
  !<<>>   rjs(22) = 0.
  !<<>>   DO ilam=1,gcCO%nlam
  !<<>>      rjs(21) = rjs(21)+sx(1,ilam)
  !<<>>      rjs(22) = rjs(22)+sx(2,ilam)
  !<<>>   ENDDO
  !<<>>
  !<<>>! J(CO2 -> CO + O) through xH1211
  !<<>>! -------------------------------
  !<<>>   DO ix=23,gcCO%nxdo
  !<<>>      rjs(ix) = 0.
  !<<>>
  !<<>>      DO ilam=1,gcCO%nlam
  !<<>>         rjs(ix) = rjs(ix)+s(ilam)*gcCO%xtab(ilam,ix,indt)
  !<<>>      ENDDO
  !<<>>   ENDDO
  !<<>>               
  !<<>>! ---------------------------------------------------------------
  !<<>>! Order photolysis rates to match order in full chemistry model.  
  !<<>>! Sort rjs into CTM photolysis rate array, aj.  Order of rjs:
  !<<>>!
  !<<>>!  1-J(BrONO2)
  !<<>>!  2-J(BrO)
  !<<>>!  3-J(Cl2O2)
  !<<>>!  4-J(ClONO2)
  !<<>>!  5-J(H2O2)
  !<<>>!  6-J(HCl)
  !<<>>!  7-J(HNO3)
  !<<>>!  8-J(HO2NO2)
  !<<>>!  9-J(HOCl)
  !<<>>! 10-J(N2O5)
  !<<>>! 11-J(NO2)
  !<<>>! 12-J(NO3_NO)
  !<<>>! 13-J(NO3_NO2)
  !<<>>! 14-J(OClO)
  !<<>>! 15-J(O2)
  !<<>>! 16-J(O3_O1D)
  !<<>>! 17-J(O3_3P)
  !<<>>! 18-J(HOBr)
  !<<>>! 19-J(CH3OOH)
  !<<>>! 20-J(N2O)
  !<<>>! 21-J(CH2O_HCO)
  !<<>>! 22-J(CH2O_CO)
  !<<>>! 23-J(CO2 -> CO + O)
  !<<>>! 24-xCFC-11
  !<<>>! 25-xCFC-12
  !<<>>! 26-xCCl4
  !<<>>! 27-xCH3CCl3
  !<<>>! 28-xHCFC-22
  !<<>>! 29-xCFC-113
  !<<>>! 30-xCH3Cl
  !<<>>! 31-xCH3Br
  !<<>>! 32-xH1301
  !<<>>! 33-xH1211 
  !<<>>! 34-xH1202
  !<<>>! 35-xH2402
  !<<>>! 36-xCHBr3
  !<<>>! 37-xCH2Br2
  !<<>>! 38-xCH2ClBr
  !<<>>! 39-xCHClBr2
  !<<>>! 40-xCHCl2Br
  !<<>>! 41-xHCFC-141b
  !<<>>! 42-xHCFC-142b
  !<<>>! 43-xCFC-114 
  !<<>>! 44-xCFC-115
  !<<>>! 45-xOCS
  !<<>>! 46-
  !<<>>! 47-
  !<<>>! 48-
  !<<>>! 49-
  !<<>>! 50-
  !<<>>! ---------------------------------------------------------------
  !<<>>! Solar cycle goes here when ready  
  !<<>>!     aj( 1) = rjs(15)*gcCO%s_cycle(3,gcCO%iscyr)
  !<<>>! ----------------------------------------------------------------
  !<<>>   aj( 1) = rjs(15)
  !<<>>   aj( 2) = rjs(16)
  !<<>>   aj( 3) = rjs(17)
  !<<>>! H2O
  !<<>>! ---
  !<<>>   aj( 4) = 0.
  !<<>>   aj( 5) = rjs(13)
  !<<>>   aj( 6) = rjs(7)
  !<<>>   aj( 7) = rjs(11)
  !<<>>   aj( 8) = rjs(5)
  !<<>>   aj( 9) = rjs(10)
  !<<>>   aj(10) = rjs(21)
  !<<>>   aj(11) = rjs(22)
  !<<>>   aj(12) = rjs(23)
  !<<>>   aj(13) = rjs(19)
  !<<>>   aj(14) = rjs(20)
  !<<>>   aj(15) = rjs(4)
  !<<>>   aj(16) = 0.
  !<<>>   aj(17) = rjs(12)
  !<<>>   aj(18) = rjs(6)
  !<<>>   aj(19) = 0.
  !<<>>
  !<<>>! CH3Br(20) H1301(21) H12_24(22)
  !<<>>! ------------------------------
  !<<>>   aj(20) = rjs(31)
  !<<>>   aj(21) = rjs(32)
  !<<>>   aj(22) = rjs(33)
  !<<>>   aj(23) = rjs(9)
  !<<>>   aj(24) = rjs(8)
  !<<>>   aj(25) = rjs(18)
  !<<>>   aj(26) = 0.
  !<<>>   aj(27) = rjs(2)
  !<<>>   aj(28) = rjs(1)
  !<<>>
  !<<>>! F11(29) F12(30) CCl4(31) CHCCl3(32) HCFC(33) F113(34) CH3Cl(35)
  !<<>>! ---------------------------------------------------------------
  !<<>>   aj(29) = rjs(24)
  !<<>>   aj(30) = rjs(25)
  !<<>>   aj(31) = rjs(26)
  !<<>>   aj(32) = rjs(27)
  !<<>>   aj(33) = rjs(28)
  !<<>>   aj(34) = rjs(29)
  !<<>>   aj(35) = rjs(30)
  !<<>>   aj(36) = rjs(3)
  !<<>>   aj(37) = rjs(14)
  !<<>>
  !<<>>! ------------------------------------------
  !<<>>! WARNING: Photolysis reaction rate
  !<<>>! numbers 38-42 are calculated in MESO_PHOT.
  !<<>>! ------------------------------------------
  !<<>>! Add aj(43) which is J(Cl2O2) for partitioning but not Ox loss 
  !<<>>! which is aj(36). In lookup table J(Cl2O2) is J*qy where qy is 0.8 
  !<<>>! so multiply by 1.25 to equal J and used in part.F and partest.F
  !<<>>
  !<<>>   aj(43) = rjs(3)*1.25
  !<<>>
  !<<>>! QingLiang -- 06/03/2015
  !<<>>! CHBr3(44) CH2Br2(45) CH2BrCl(46) CHBrCl2(47) CHBr2Cl(48)
  !<<>>   aj(44) = rjs(36)
  !<<>>   aj(45) = rjs(37)
  !<<>>   aj(46) = rjs(38)
  !<<>>   aj(47) = rjs(39)
  !<<>>   aj(48) = rjs(40)
  !<<>>
  !<<>>! QingLiang -- 06/03/2015
  !<<>>! Add two new halons: H-1202 (49) H2402 (50) 
  !<<>>! and two new HCFCs: HCFC-141b (51) HCFC-142b (52) 
  !<<>>   aj(49) = rjs(34)
  !<<>>   aj(50) = rjs(35)
  !<<>>   aj(51) = rjs(41)
  !<<>>   aj(52) = rjs(42)
  !<<>>
  !<<>>! QingLiang -- 02/05/2016
  !<<>>! Add CFC-114 and CFC-115
  !<<>>! Add OCS for GOCART module
  !<<>>   aj(53) = rjs(43)
  !<<>>   aj(54) = rjs(44)
  !<<>>   aj(55) = rjs(45)
  !<<>>!  aj(53) = rjs(34)
  !<<>>!  aj(54) = rjs(34)
  !<<>>!  aj(55) = rjs(34)

      RETURN
    END SUBROUTINE jcalc4

    SUBROUTINE readPhotTables(fileName, rc)

      IMPLICIT NONE

  !  Read tables for photolysis in GOCART ... from a NetCDF file
  !
  !  Input parameters:
  !
      CHARACTER(LEN=*), INTENT(IN) :: fileName
  !
  !  Output parameters:
  !
      INTEGER, INTENT(OUT) :: rc
  !
  !  Restrictions:
  !  ASSERT that the number of pressure layers in the dataset equals km.
  !
  !  REVISION HISTORY:
  !  Nielsen     11 May 2012: First crack.
  !  Weir        29 Jan 2021: Pilferd from StratChem
  !-----------------------------------------------------------------------

      CHARACTER(LEN=255) :: Iam = "CO::readPhotTables"

  !  TYPE(ESMF_VM) :: vm

      INTEGER :: comm, info, unit, status
      INTEGER :: dimid, i, n

      INTEGER :: length

      INTEGER, PARAMETER :: nD = 7
      CHARACTER(LEN=255) :: dimName(nD)= (/"nsza  ", "numO3 ", "layers", &
           "nlam  ", "nts   ", "nxdo  ", "aqsize" /)

      INTEGER, PARAMETER :: nV = 7
      CHARACTER(LEN=255) :: varName(nV)= (/"sza    ", &
           "lambda ", "O3TAB  ",  "SDAT   ", &
           "O2JDAT ", "XTAB   ",  "CH2O_AQ" /)
      rc = 0

  ! Grab the virtual machine
  ! ------------------------
  !<<>>  CALL ESMF_VMGetCurrent(vm, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>#ifdef H5_HAVE_PARALLEL
  !<<>>
  !<<>>  CALL MPI_Info_create(info, status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>#ifdef NETCDF_NEED_NF_MPIIO
  !<<>>  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
  !<<>>#else
  !<<>>  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
  !<<>>#endif
  !<<>>
  !<<>>#else
  !<<>>
  !<<>>  IF(MAPL_AM_I_ROOT(vm)) THEN 
  !<<>>   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)
  !<<>>
  !<<>>#endif
  !<<>>
  !<<>>   IF(status /= NF_NOERR) THEN
  !<<>>    PRINT *,'Error opening file ',TRIM(fileName), status
  !<<>>    PRINT *, NF_STRERROR(status)
  !<<>>    VERIFY_(status)
  !<<>>   END IF
  !<<>>
  !<<>>   DO i = 1,nD
  !<<>>
  !<<>>    status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
  !<<>>    IF(status /= NF_NOERR) THEN
  !<<>>     PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
  !<<>>     PRINT *, NF_STRERROR(status)
  !<<>>     VERIFY_(status)
  !<<>>    END IF
  !<<>>
  !<<>>    status = NF_INQ_DIMLEN(unit, dimid, n)
  !<<>>    IF(status /= NF_NOERR) THEN
  !<<>>     PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
  !<<>>     PRINT *, NF_STRERROR(status)
  !<<>>    END IF
  !<<>>
  !<<>>    SELECT CASE (i)
  !<<>>     CASE (1)
  !<<>>      gcCO%nsza = n
  !<<>>     CASE (2)
  !<<>>      gcCO%numO3 = n
  !<<>>     CASE (3)
  !<<>>      ASSERT_(n == km)
  !<<>>     CASE (4)
  !<<>>      gcCO%nlam = n
  !<<>>     CASE (5)
  !<<>>      gcCO%nts = n
  !<<>>     CASE (6)
  !<<>>      gcCO%nxdo = n
  !<<>>     CASE (7)
  !<<>>      gcCO%aqsize = n
  !<<>>     CASE DEFAULT
  !<<>>    END SELECT
  !<<>>
  !<<>>   END DO
  !<<>>
  !<<>>#ifndef H5_HAVE_PARALLEL
  !<<>>
  !<<>>  END IF ! MAPL_AM_I_ROOT
  !<<>>
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%nsza, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%numO3, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%nlam, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%nts, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%nxdo, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%aqSize, 1, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>#endif
  !<<>>
  !<<>>  ALLOCATE(gcCO%sdat(gcCO%nsza,gcCO%numo3,km,gcCO%nlam), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%o2jdat(gcCO%nsza,gcCO%numo3,km), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%o3_tab(gcCO%numo3,km), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%xtab(gcCO%nlam,gcCO%nxdo,gcCO%nts), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%sza_tab(gcCO%nsza), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%CH2O_aq(gcCO%aqSize), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>  ALLOCATE(gcCO%rlam(gcCO%nlam), STAT=status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>#ifndef H5_HAVE_PARALLEL
  !<<>>
  !<<>>  IF(MAPL_AM_I_ROOT()) THEN
  !<<>>
  !<<>>#endif
  !<<>>
  !<<>>   DO i = 1,nV
  !<<>>
  !<<>>    status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
  !<<>>    IF(status /= NF_NOERR) THEN
  !<<>>     PRINT *,"Error getting varid for ", TRIM(varName(i)), status
  !<<>>     PRINT *, NF_STRERROR(status)
  !<<>>     VERIFY_(status)
  !<<>>    END IF
  !<<>>
  !<<>>    SELECT CASE (i)
  !<<>>     CASE (1)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%sza_tab)
  !<<>>     CASE (2)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%rlam)
  !<<>>     CASE (3)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%o3_tab)
  !<<>>     CASE (4)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%sdat)
  !<<>>     CASE (5)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%o2jdat)
  !<<>>     CASE (6)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%xtab)
  !<<>>     CASE (7)
  !<<>>      status = NF_GET_VAR_REAL(unit, n, gcCO%CH2O_aq)
  !<<>>     CASE DEFAULT
  !<<>>    END SELECT
  !<<>>
  !<<>>    IF(status /= NF_NOERR) THEN
  !<<>>     PRINT *,"Error getting values for ", TRIM(varName(i)), status
  !<<>>     PRINT *, NF_STRERROR(status)
  !<<>>     VERIFY_(status)
  !<<>>    END IF
  !<<>>
  !<<>>   END DO
  !<<>>
  !<<>>#ifdef H5_HAVE_PARALLEL
  !<<>>
  !<<>>   CALL MPI_Info_free(info, status)
  !<<>>   VERIFY_(status)
  !<<>>
  !<<>>#else
  !<<>>
  !<<>>  END IF ! MAPL_AM_I_ROOT
  !<<>>
  !<<>>  length = SIZE(gcCO%sza_tab)
  !<<>>  CALL MPI_Bcast(gcCO%sza_tab, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  length = SIZE(gcCO%rlam)
  !<<>>  CALL MPI_Bcast(gcCO%rlam, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  length = SIZE(gcCO%o3_tab)
  !<<>>  CALL MPI_Bcast(gcCO%o3_tab, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  length = SIZE(gcCO%sdat)
  !<<>>  CALL MPI_Bcast(gcCO%sdat, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  length = SIZE(gcCO%o2jdat)
  !<<>>  CALL MPI_Bcast(gcCO%o2jdat, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  length = SIZE(gcCO%xtab)
  !<<>>  CALL MPI_Bcast(gcCO%xtab, length, MPI_REAL, 0, comm, status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  CALL MAPL_CommsBcast(vm, gcCO%CH2O_aq, gcCO%aqsize, 0, RC=status)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>#endif
  !<<>>
  !<<>>  status = NF_CLOSE(unit)
  !<<>>  VERIFY_(status)
  !<<>>
  !<<>>  RETURN
    END SUBROUTINE readPhotTables

  END SUBROUTINE CO_ops

  !-------------------------------------------------------------------------

END MODULE CO_mod
