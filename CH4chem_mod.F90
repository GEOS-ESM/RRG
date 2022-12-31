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

  subroutine Calc_photolysis_rates()!gcCH4, w_c, impChem, photJ, rc)
    !   !-------------------------------------------------------------------------
    !   ! Borrowed from meso_phot.F of StratChem, where number densities are cgs [cm^{-3}]
    !
    !   implicit none
    !
    !   ! Input variables
    !   type(CH4_GridComp1), intent(in)     :: gcCH4   ! Grid Component
    !   type(Chem_Bundle), intent(in)       :: w_c     ! Chemical tracer fields
    !   type(ESMF_State), intent(inout)     :: impChem ! Import State
    !
    !   ! Output variables
    !   integer, intent(out)                :: rc
    !   real, dimension(:,:,:), intent(out) :: photJ
    !
    !   real, allocatable :: o2Column(:,:,:), ndwet(:,:,:), cellDepth(:,:,:)
    !   real, allocatable :: SZARad(:,:), SZADeg(:,:), sinSZA(:,:)
    !   real, allocatable :: zgrz(:,:), sfaca(:,:), arg(:,:)
    !   real, pointer, dimension(:,:,:) ::  rhowet => null()
    !   real, pointer, dimension(:,:,:) ::  ple => null()
    !
    !   real, parameter :: wavel = 1215.7
    !   real, parameter :: O2xs  = 1.000E-20
    !   real, parameter :: CH4xs = 2.000E-17
    !   real, parameter :: sflux = 4.006E+11
    !
    !   ! Constants for Chapman function at high solar zenith angle
    !   ! ---------------------------------------------------------
    !   real, parameter :: hbar = 6.79
    !   real, parameter :: zbar = 30.0
    !   real, parameter :: r0   = 6.371E+03
    !   real, parameter :: zp   = 60.0
    !
    !   real, parameter :: d1 = 1.060693
    !   real, parameter :: d2 = 0.55643831
    !   real, parameter :: d3 = 1.0619896
    !   real, parameter :: d4 = 1.7245609
    !   real, parameter :: d5 = 0.56498823
    !   real, parameter :: d6 = 0.06651874
    !
    !   real, parameter :: O2Abv80km = 7.072926E+19 ![cm^{-2}]
    !   real, parameter :: O2VMR = 0.20946
    !
    !   real :: b, f, r, s
    !   integer :: status, i1, i2, j1, j2, km, i, j, k
    !
    !   character(len=*), parameter :: Iam = "CH4::Calc_photolysis_rates"
    !
    !   rc = 0
    !   photJ(:,:,:) = 0
    !
    !   i1 = w_c%grid%i1
    !   i2 = w_c%grid%i2
    !   j1 = w_c%grid%j1
    !   j2 = w_c%grid%j2
    !   km = w_c%grid%km
    !
    !   ! Get pointers
    !   call MAPL_GetPointer(impChem, rhowet, 'AIRDENS', __RC__)
    !   call MAPL_GetPointer(impChem, ple, 'PLE', __RC__)
    !
    !   b = sqrt(0.50*r0/hbar)
    !
    !   ! O2 overhead number density profile [cm^{-2}]
    !   ! --------------------------------------------
    !   allocate(O2Column(i1:i2,j1:j2,1:km), ndwet(i1:i2,j1:j2,km), cellDepth(i1:i2,j1:j2,km), STAT=status)
    !   VERIFY_(status)
    !
    !   !  Wet-air number density
    !   ndwet(i1:i2,j1:j2,1:km) = rhowet(i1:i2,j1:j2,1:km)*MAPL_AVOGAD/MAPL_AIRMW
    !   !  Cell depth
    !   do k=1,km
    !      cellDepth(i1:i2,j1:j2,k) = (ple(i1:i2,j1:j2,k)-ple(i1:i2,j1:j2,k-1)) / (rhowet(i1:i2,j1:j2,k)*MAPL_GRAV)
    !   end do
    !
    !   f = O2VMR*5.00E-05
    !   O2Column(:,:,1) = O2Abv80km+cellDepth(:,:,1)*ndwet(:,:,1)*f
    !
    !   do k = 2,km
    !      O2Column(:,:,k) = O2Column(:,:,k-1)+(cellDepth(:,:,k-1)*ndwet(:,:,k-1)+ &
    !                                     cellDepth(:,:,  k)*ndwet(:,:,  k))*f
    !   end do
    !
    !   !IF(gcCH4%DebugIsOn) THEN
    !      !CALL pmaxmin('CH4: O2Column', O2Column, qmin, qmax, iXj, km,  1. )
    !   !END IF
    !
    !   ! Grab some memory
    !   ! ----------------
    !   allocate(SZARad(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !   allocate(SZADeg(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !   allocate(sinSZA(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !
    !   where(w_c%cosz(i1:i2,j1:j2) > 1.00)
    !      SZARad(i1:i2,j1:j2) = 0.00
    !   elsewhere
    !      SZARad(i1:i2,j1:j2) = ACOS(w_c%cosz(i1:i2,j1:j2))
    !   endwhere
    !   SZADeg(i1:i2,j1:j2) = SZARad(i1:i2,j1:j2)*radToDeg
    !   sinSZA(i1:i2,j1:j2) = sin(SZARad(i1:i2,j1:j2))
    !
    !   allocate(zgrz(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !
    !   where(SZADeg(i1:i2,j1:j2) <= 90.00)
    !      zgrz(i1:i2,j1:j2) = 1000.00
    !   elsewhere
    !      zgrz(i1:i2,j1:j2) = sinSZA(i1:i2,j1:j2)*(zp+r0)-r0
    !   endwhere
    !
    !   !IF(gcCH4%DebugIsOn) THEN
    !      !CALL pmaxmin('CH4: zgrz', zgrz, qmin, qmax, iXj, 1,  1. )
    !      !CALL pmaxmin('CH4: cosz', w_c%cosz, qmin, qmax, iXj, 1,  1. )
    !   !END IF
    !
    !   allocate(sfaca(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !   sfaca(i1:i2,j1:j2) = 0.00
    !
    !   ! Chapman function calculation from ACDB 2-D model
    !   ! ------------------------------------------------
    !   do j = j1,j2
    !      do i = i1,i2
    !
    !         if (SZADeg(i,j) < gcCH4%szaCutoff) then ! Daytime
    !
    !            if (SZADeg(i,j) < 70.00) then
    !
    !               sfaca(i,j) = 1.00/w_c%cosz(i,j)
    !
    !            else if (zgrz(i,j) > 0.00) then
    !
    !               s = b*abs(w_c%cosz(i,j))
    !
    !               if (s <= 8.00) then
    !                  s = (d1+d2*s)/(d3+d4*s+s**2)
    !               else
    !                  s = d5/(d6+s)
    !               end if
    !
    !               r = b*sqrt(MAPL_PI)
    !               sfaca(i,j) = r*s
    !
    !               if (SZADeg(i,j) > 90.00) then
    !                  sfaca(i,j) = 2.00*r*exp((r0+zbar)*(1.00-sinSZA(i,j))/hbar)-sfaca(i,j)
    !               end if
    !
    !            end if
    !
    !         end if ! Daytime
    !
    !      end do
    !   end do
    !
    !   !IF(gcCH4%DebugIsOn) THEN
    !      !CALL pmaxmin('CH4: sfaca', sfaca, qmin, qmax, iXj, 1,  1. )
    !   !END IF
    !
    !   allocate(arg(i1:i2,j1:j2), STAT=status)
    !   VERIFY_(status)
    !
    !   ! At each layer, compute the rate constant, J [s^{-1}], if the sun is up
    !   ! ----------------------------------------------------------------------
    !   do k = 1,km
    !
    !      where(SZADeg(i1:i2,j1:j2) < gcCH4%szaCutoff)
    !         arg(i1:i2,j1:j2) = O2Column(i1:i2,j1:j2,k)*O2xs*sfaca(i1:i2,j1:j2)
    !         photJ(i1:i2,j1:j2,k) = sflux*exp(-arg(i1:i2,j1:j2))*CH4xs
    !      end where
    !
    !   end do
    !
    !   !IF(gcCH4%DebugIsOn) THEN
    !      !CALL pmaxmin('CH4: photJ', photJ, qmin, qmax, iXj, km,  1. )
    !   !END IF
    !
    !   deallocate(SZARad, STAT=status)
    !   VERIFY_(status)
    !   deallocate(SZADeg, STAT=status)
    !   VERIFY_(status)
    !   deallocate(sinSZA, STAT=status)
    !   VERIFY_(status)
    !   deallocate(zgrz, STAT=status)
    !   VERIFY_(status)
    !   deallocate(sfaca, STAT=status)
    !   VERIFY_(status)
    !   deallocate(arg, STAT=status)
    !   VERIFY_(status)
    !   deallocate(O2Column, ndwet, cellDepth, STAT=status)
    !   VERIFY_(status)
    !
  end subroutine Calc_photolysis_rates

  subroutine CH4_prodloss( CH4inst, OH, O1D, Cl, RC )

    use global_mod
    use types_mod

    implicit none

    ! INPUT PARAMETERS:

    type(gas_instance), pointer, intent(inout)  :: CH4inst(:)

    ! Reactants
    real, pointer, intent(in)     ::  OH(:,:,:)
    real, pointer, intent(in)     :: O1D(:,:,:)
    real, pointer, intent(in)     ::  Cl(:,:,:)

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
    real, allocatable  :: cvfac(:,:,:) ! Conversion factor from kg/kg -> mcl/cm3
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
    !  ASSUMPTION: all species units are input kg/kg
    !  --------------------------------------------------------
    allocate(k_(im,jm,km), stat=RC)
    do nst = 1,size(CH4inst) ! cycle over instances
       loss  => CH4inst(nst)%loss(:,:,:) ! in kg/kg/s
       CH4    => CH4inst(nst)%data3d(:,:,:) ! CH4 pointer makes the code cleaner
       !  Loss rate [m^3 s^-1] for OH + CH4 => CO + products
       k_  = 2.45E-12*exp(-1775./met%T)
       loss = loss + k_*CH4*OH

       !  Loss rate [m^3 s^-1] for Cl + CH4 => CO + products
       k_  = 7.10E-12*exp(-1270./met%T)
       loss = loss + k_*CH4*Cl

       !  Loss rate [m^3 s^-1] for O(1D) + CH4 => CO + products
       k_ = 1.75E-10
       loss = loss + k_*CH4*O1D

       CH4  => null()
       loss => null()
    enddo

    !   if (gcCH4%photolysis) then
    !  Calculate photolytic loss rates, J [s^-1], for CH4 + hv => 2H2O + CO
    !      allocate(photJ(i1:i2,j1:j2,1:km), dCH4Phot(i1:i2,j1:j2,1:km), STAT=status)
    !      call Calc_photolysis_rates(gcCH4, w_c, impChem, photJ, __RC__)
    !      dCH4Phot = -cdt * photJ(i1:i2,j1:j2,1:km) * w_c%qa(nbeg)%data3d(i1:i2,j1:j2,1:km) ! fwd Euler (1st order)
    !      deallocate(photJ, dCH4Phot, STAT=status)
    !   else
    !      dCH4Phot = 0.0
    !   end if

    !   if (gcCH4%CH4FeedBack .and. gcCH4%H2OFeedBack) then
    !      Q(i1:i2,j1:j2,1:km) = Q(i1:i2,j1:j2,1:km) + 2.00*cdt*dCH4Phot(i1:i2,j1:j2,1:km)*MAPL_H2OMW/MAPL_AIRMW
    !   end if

    !  Housekeeping
    !  ------------
    deallocate(k_, stat=RC)
    return

  end subroutine CH4_prodloss

end module CH4chem_mod
