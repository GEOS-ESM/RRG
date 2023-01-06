#include "MAPL_Generic.h"
module geos_SimplePhotolysisMod
  use mapl
  use esmf
  use ESMF_CFIOFileMOD
  use MAPL_CFIOMOD

  use types_mod

  implicit none
  save

  INTEGER       :: numphoto
  INTEGER       :: nxdo
  INTEGER       :: nlam
  INTEGER       :: nsza
  INTEGER       :: numo3
  INTEGER       :: nts
  INTEGER       :: aqsize

  REAL, POINTER :: sdat(:,:,:,:)
  REAL, POINTER :: o2jdat(:,:,:)
  REAL, POINTER :: sza_tab(:)
  REAL, POINTER :: o3_tab(:,:)
  REAL, POINTER :: xtab(:,:,:)
  REAL, POINTER :: CH2O_aq(:)
  REAL, POINTER :: rlam(:)

contains

  subroutine CH4_photolysis_rate(i,j,k,met,O2col,SZAcutoff,JV)
!-------------------------------------------------------------------------
! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
! -- modified for GHG component
    implicit none

! Input variables
    integer,                    intent(in) :: i,j,k
    type(meteorology),          intent(in) :: met
    real,                       intent(in) :: O2col, SZAcutoff

! Output variables
!    integer, intent(out) :: rc
    real,    intent(out) :: JV

    real :: SZARad, SZADeg, sinSZA
    real :: zgrz, sfaca, arg

    real, parameter :: wavel = 1215.7
    real, parameter :: O2xs  = 1.000E-20
    real, parameter :: CH4xs = 2.000E-17
    real, parameter :: sflux = 4.006E+11
    
! Constants for Chapman function at high solar zenith angle
! ---------------------------------------------------------
    real, parameter :: hbar = 6.79
    real, parameter :: zbar = 30.0
    real, parameter :: r0   = 6.371E+03
    real, parameter :: zp   = 60.0
    
    real, parameter :: d1 = 1.060693
    real, parameter :: d2 = 0.55643831
    real, parameter :: d3 = 1.0619896
    real, parameter :: d4 = 1.7245609
    real, parameter :: d5 = 0.56498823
    real, parameter :: d6 = 0.06651874

    real :: b, f, r, s

    character(len=*), parameter :: Iam = "CH4::Calc_photolysis_rates"

    JV = 0.

    if (met%cosz(i,j) > 1.00) then
       SZARad = 0.00
    else
       SZARad = ACOS(met%cosz(i,j))
    endif
    SZADeg = SZARad*180./MAPL_PI

    if (SZADeg < szaCutoff) return ! no need to proceed, save the FLOPS

    ! We've determined the sun is (sufficiently) up, proceed

    b = sqrt(0.50*r0/hbar)
    sinSZA = sin(SZARad)
    if (SZADeg <= 90.00) then
       zgrz = 1000.00
    else
       zgrz = sinSZA*(zp+r0)-r0
    end if
    sfaca = 0.00

   ! Chapman function calculation from ACDB 2-D model
   ! ------------------------------------------------
    if (SZADeg < 70.00) then
       sfaca = 1.00/met%cosz(i,j)
    else if (zgrz > 0.00) then
       s = b*abs(met%cosz(i,j))
       if (s <= 8.00) then
          s = (d1+d2*s)/(d3+d4*s+s**2)
       else
          s = d5/(d6+s)
       end if
       r = b*sqrt(MAPL_PI)
       sfaca = r*s
       if (SZADeg > 90.00) then
          sfaca = 2.00*r*exp((r0+zbar)*(1.00-sinSZA)/hbar)-sfaca
       end if
    end if
       
   ! Compute the rate constant, J [s^{-1}]
   ! ----------------------------------------------------------------------

    arg = O2Col*O2xs*sfaca
    JV = sflux*exp(-arg)*CH4xs

 end subroutine CH4_photolysis_rate

  SUBROUTINE CO2_photolysis_rate(i,j,k,ki,met,o3col,aj)

    ! ---------------------------------------------------------------------------------
    ! NAME: was jcalc4
    ! PURPOSE:
    !   Calculate photolysis rates
    ! INPUT:
    !   k         Current layer number
    !   levels    Number of layers
    !   szan      Solar zenith angle (radians)
    !   o3column  Overhead O3 values
    !   press     Mid-layer pressure (hPa)
    ! OUTPUT:
    !   aj        Array of photolysis rates
    ! RESTRICTIONS:
    !   Currently set up for 23-J set (see var nxdo)
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
    implicit none
    integer,           intent(in)  :: i, j, k, ki
    type(meteorology), intent(in)  :: met
    real,              intent(in)  :: o3col
    real,              intent(out) :: aj

    ! Locals
    integer :: ilam,indt,ix
    real :: s(nlam)

    aj = 0.

    ! Interpolate radiative flux function values to model conditions
    ! --------------------------------------------------------------
    CALL interp_s(ki,acos(met%cosz(i,j)),o3col,s)
    indt = met%t(i,j,k)-148.5
    indt = MAX(1,indt)
    indt = MIN(indt,200)
    DO ilam=1,nlam
       aj = aj+s(ilam)*xtab(ilam,23,indt)
    ENDDO
    RETURN
  END SUBROUTINE CO2_photolysis_rate

  SUBROUTINE readPhotTables(fileName, km, rc)

    IMPLICIT NONE

    !  Read tables for photolysis in GOCART ... from a NetCDF file
    !
    !  Input parameters:
    !
    CHARACTER(LEN=*), INTENT(IN) :: fileName
    integer,          intent(in) :: km
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

    CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "CO::readPhotTables"

    TYPE(ESMF_VM) :: vm

    INTEGER :: comm, info, unit, status
    INTEGER :: dimid, i, n

    INTEGER :: length

    INTEGER, PARAMETER :: nD = 7
    CHARACTER(LEN=ESMF_MAXSTR) :: dimName(nD)= (/"nsza  ", "numO3 ", "layers", &
         "nlam  ", "nts   ", "nxdo  ", "aqsize" /)

    INTEGER, PARAMETER :: nV = 7
    CHARACTER(LEN=ESMF_MAXSTR) :: varName(nV)= (/"sza    ", &
         "lambda ", "O3TAB  ",  "SDAT   ", &
         "O2JDAT ", "XTAB   ",  "CH2O_AQ" /)
    rc = 0

    ! Grab the virtual machine
    ! ------------------------
    CALL ESMF_VMGetCurrent(vm, RC=status)
    VERIFY_(status)

    CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
    VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

    CALL MPI_Info_create(info, status)
    VERIFY_(status)
    CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
    VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
    status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
    status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

    IF(MAPL_AM_I_ROOT(vm)) THEN 
       status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

       IF(status /= NF_NOERR) THEN
          PRINT *,'Error opening file ',TRIM(fileName), status
          PRINT *, NF_STRERROR(status)
          VERIFY_(status)
       END IF

       DO i = 1,nD

          status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
          IF(status /= NF_NOERR) THEN
             PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
             PRINT *, NF_STRERROR(status)
             VERIFY_(status)
          END IF

          status = NF_INQ_DIMLEN(unit, dimid, n)
          IF(status /= NF_NOERR) THEN
             PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
             PRINT *, NF_STRERROR(status)
          END IF

          SELECT CASE (i)
          CASE (1)
             nsza = n
          CASE (2)
             numO3 = n
          CASE (3)
             ASSERT_(n == km)
          CASE (4)
             nlam = n
          CASE (5)
             nts = n
          CASE (6)
             nxdo = n
          CASE (7)
             aqsize = n
          CASE DEFAULT
          END SELECT

       END DO

#ifndef H5_HAVE_PARALLEL

    END IF ! MAPL_AM_I_ROOT

    CALL MAPL_CommsBcast(vm, nsza, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, numO3, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, nlam, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, nts, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, nxdo, 1, 0, RC=status)
    VERIFY_(status)
    CALL MAPL_CommsBcast(vm, aqSize, 1, 0, RC=status)
    VERIFY_(status)

#endif

    ALLOCATE(sdat(nsza,numo3,km,nlam), STAT=status)
    VERIFY_(status)
    ALLOCATE(o2jdat(nsza,numo3,km), STAT=status)
    VERIFY_(status)
    ALLOCATE(o3_tab(numo3,km), STAT=status)
    VERIFY_(status)
    ALLOCATE(xtab(nlam,nxdo,nts), STAT=status)
    VERIFY_(status)
    ALLOCATE(sza_tab(nsza), STAT=status)
    VERIFY_(status)
    ALLOCATE(CH2O_aq(aqSize), STAT=status)
    VERIFY_(status)
    ALLOCATE(rlam(nlam), STAT=status)
    VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

    IF(MAPL_AM_I_ROOT()) THEN

#endif

       DO i = 1,nV

          status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
          IF(status /= NF_NOERR) THEN
             PRINT *,"Error getting varid for ", TRIM(varName(i)), status
             PRINT *, NF_STRERROR(status)
             VERIFY_(status)
          END IF

          SELECT CASE (i)
          CASE (1)
             status = NF_GET_VAR_REAL(unit, n, sza_tab)
          CASE (2)
             status = NF_GET_VAR_REAL(unit, n, rlam)
          CASE (3)
             status = NF_GET_VAR_REAL(unit, n, o3_tab)
          CASE (4)
             status = NF_GET_VAR_REAL(unit, n, sdat)
          CASE (5)
             status = NF_GET_VAR_REAL(unit, n, o2jdat)
          CASE (6)
             status = NF_GET_VAR_REAL(unit, n, xtab)
          CASE (7)
             status = NF_GET_VAR_REAL(unit, n, CH2O_aq)
          CASE DEFAULT
          END SELECT

          IF(status /= NF_NOERR) THEN
             PRINT *,"Error getting values for ", TRIM(varName(i)), status
             PRINT *, NF_STRERROR(status)
             VERIFY_(status)
          END IF

       END DO

#ifdef H5_HAVE_PARALLEL

       CALL MPI_Info_free(info, status)
       VERIFY_(status)

#else

    END IF ! MAPL_AM_I_ROOT

    length = SIZE(sza_tab)
    CALL MPI_Bcast(sza_tab, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    length = SIZE(rlam)
    CALL MPI_Bcast(rlam, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    length = SIZE(o3_tab)
    CALL MPI_Bcast(o3_tab, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    length = SIZE(sdat)
    CALL MPI_Bcast(sdat, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    length = SIZE(o2jdat)
    CALL MPI_Bcast(o2jdat, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    length = SIZE(xtab)
    CALL MPI_Bcast(xtab, length, MPI_REAL, 0, comm, status)
    VERIFY_(status)

    CALL MAPL_CommsBcast(vm, CH2O_aq, aqsize, 0, RC=status)
    VERIFY_(status)

#endif

    status = NF_CLOSE(unit)
    VERIFY_(status)

    RETURN
  END SUBROUTINE readPhotTables
  
  !-------------------------------------------------------------------------
  SUBROUTINE interp_s(k,sza,o3column,s)
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
    !
    ! OUTPUTS:
    !   s         S value for each wavelength at current k, interpolated to
    !               the given o3column and sza
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

    INTEGER, INTENT(IN) :: k
    REAL, INTENT(IN) :: sza, o3column 
    REAL, INTENT(OUT) :: s(nlam)

    INTEGER :: ijj, ik, ikk, ikkm, il, is
    REAL :: omt, omu, t, u
    REAL, PARAMETER :: PI = 3.14159265

    ! For each input solar zenith angle, find the first element of sza_tab that 
    ! is greater.  Use this element and previous one to determine the interpolated value.
    ! -----------------------------------------------------------------------------------
    DO is = 1,nsza
       ijj = is 
       IF(sza_tab(is) > sza) EXIT 
    ENDDO

    ! Zenith angle test       
    ! -----------------
    IF(sza > sza_tab(nsza)) THEN
       !     Cell is dark, set s=0        
       !     -----------------------------
       s(1:nlam) = 0.
    ELSE  
       !     Cell is illuminated     
       !     -------------------
       t = (sza-sza_tab(ijj-1))/(sza_tab(ijj)-sza_tab(ijj-1))
       omt = 1.-t

       ! For each overhead O3 column, find the first element in o3_tab that is
       ! greater. Use this element and previous one to determine the interpolated value.
       ! -------------------------------------------------------------------------------
       DO is = 1,numo3
          ikk = is 
          IF(o3_tab(is,k) > o3column) EXIT
       ENDDO

       ikkm = ikk-1 
       IF(ikk > 1 .AND. o3column <= o3_tab(numo3,k)) THEN
          u = (o3column-o3_tab(ikkm,k))/(o3_tab(ikk,k)-o3_tab(ikkm,k))
          omu = 1.-u

          ! Do bilinear interpolation for each wavelength.
          ! ----------------------------------------------
          DO il = 1,nlam       
             s(il) = omt*omu*sdat(ijj-1,ikkm,k,il)+t*omu*sdat(ijj,ikkm,k,il)+ &
                  t*u*sdat(ijj,ikk,k,il)+omt*u*sdat(ijj-1,ikk,k,il)
          ENDDO
          ! Extrapolate ahead of table
          ! --------------------------
       ELSE IF (ikk == 1) THEN
          DO il = 1,nlam
             s(il) = omt*sdat(ijj-1,1,k,il)+t*sdat(ijj,1,k,il)
          ENDDO
          ! Extrapolate beyond table
          ! ------------------------
       ELSE
          DO il = 1,nlam
             s(il) = omt*sdat(ijj-1,numo3,k,il)+t*sdat(ijj,numo3,k,il)
          END DO
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE interp_s

end module geos_SimplePhotolysisMod
