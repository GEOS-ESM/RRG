MODULE Surface_Mod

! !USES:
   use ESMF
   use MAPL
   use types_mod

   implicit none
   public

   integer :: minpbl

contains

!-------------------------------------------------------------------------------
  SUBROUTINE surface_prodloss ( rc )

! !USES:

    use global_mod

    IMPLICIT NONE

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:
    integer, intent(out) :: rc ! Error return code

! !LOCAL VARIABLES
    ! Below, this is done for cleanliness

    integer :: i, j, k, n, nst, index, spc
    integer :: status
    real, allocatable :: fPBL(:,:,:), fDNL(:,:)
    real              :: fdC

    type(gas_instance), pointer :: inst

    if (.not. allocated(sfc_flux)) then ! no fluxes to compute.
       RC = 0
       return
    end if

!   1) IF NEEDED, calculate PBL height and diurnal factors
    if (any(sfc_flux(:)%pblmix)) then 
       allocate(fPBL(params%im,params%jm,params%km), stat=status)
       call surface_pblmix( fPBL, met, params, RC )
    endif

    if (any(sfc_flux(:)%diurnal)) then
       ! Ideally, this is only done once. BUT alas.
       ! Can't keep fDNL allocated locally like that (can we? Should we?)
       allocate(fDNL(params%im,params%jm), stat=status)
       call surface_DiurnalScaling ( fDNL, params )
    endif

!   2) Loop through and compute fluxes
!      -- this is written to keep the branching operations outside of nested do loops
!         and thus to a minimum
    do n=1,size(sfc_flux) ! should be > 0
       index = sfc_flux(n)%index
       inst => instances(index)%p ! For cleanliness & convenience

!      --------------------
       if (.not. inst%hasmask) then
!      -- N O  M A S K --
       if (sfc_flux(n)%pblmix) then
!      -- H A S  P B L M I X --
          ! ASSUMPTION: any flux that is mixed in the PBL will be positive 
          if (sfc_flux(n)%diurnal) then
             do k=minPBL,params%km
                inst%prod(:,:,k)   = inst%prod(:,:,k) + &
                                     sfc_flux(n)%flux(:,:) * fPBL(:,:,k) * fDNL(:,:)  * &
                                     sfc_flux(n)%scalefactor * grav / met%delp(:,:,k) * &
                                     params%airmw / inst%mw
             end do
          endif
          if (.not. sfc_flux(n)%diurnal) then
             do k=minPBL,params%km
                inst%prod(:,:,k)   = inst%prod(:,:,k) + &
                                     sfc_flux(n)%flux(:,:) * fPBL(:,:,k) * &
                                     sfc_flux(n)%scalefactor * grav / met%delp(:,:,k) * &
                                     params%airmw / inst%mw
             end do
          endif
       endif ! Non-PBL may be positive or negative (e.g. OCN, NEP)
       if (.not. sfc_flux(n)%pblmix) then
!      -- N O  P B L M I X --
          k = params%km ! at the surface
          if (sfc_flux(n)%diurnal) then
             do j=1,params%jm
             do i=1,params%im
                ! Source term?
                if (sfc_flux(n)%flux(i,j) .gt. 0.e0) then
                   inst%prod(i,j,k)   = inst%prod(i,j,k) +  &
                                                      sfc_flux(n)%flux(i,j) * fDNL(i,j) * &
                                                      sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                                      params%airmw / inst%mw
                   cycle
                endif
                ! Loss term?
                ! Don't branch. Just re-ask 'if'
                spc = inst%ispecies ! Species index
                if (sfc_flux(n)%flux(i,j) .lt. 0.e0 .and. aggregate(spc)%q(i,j,k).gt.0.e0) then ! Sink. Removes aggregate, not just one instance
                   if (cntrl%wellmixed_sfcexch) then
                      fdC = (sfc_flux(n)%flux(i,j) * fDNL(i,j) * sfc_flux(n)%scalefactor * grav / met%delp(i,j,k)) / aggregate(spc)%q(i,j,k) * params%airmw / inst%mw
                      ! Loop over all instances
                      if (inst%active) then
                         do nst=1,NINSTANCES
                            if (instances(nst)%p%ispecies .eq. spc .and. instances(nst)%p%active) &
                                 instances(nst)%p%loss(i,j,k) = instances(nst)%p%loss(i,j,k)-fdC*instances(nst)%p%data3d(i,j,k) ! Pay attention to the sign! Losses should still be stored as positive numbers!
                         enddo
                      endif
                      if (.not. inst%active) then ! not active, only remove the fraction of this instance (inst%q/aggregate%q)
                         inst%loss(i,j,k) = inst%loss(i,j,k)-fdC*inst%data3d(i,j,k)
                      endif
                   endif
                   if (.not. cntrl%wellmixed_sfcexch) then
                         inst%loss(i,j,k) = inst%loss(i,j,k) -  &
                                                      sfc_flux(n)%flux(i,j) * fDNL(i,j) * &
                                                      sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                                      params%airmw / inst%mw
                   endif
                endif
             enddo
             enddo
          endif
!         -- N O  P B L M I X -- N O  D I U R N A L --
          if (.not. sfc_flux(n)%diurnal) then
             do j=1,params%jm
             do i=1,params%im
                ! Source term?
                if (sfc_flux(n)%flux(i,j) .gt. 0.e0) then
                   inst%prod(i,j,k)   = inst%prod(i,j,k) +  &
                                        sfc_flux(n)%flux(i,j) * sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                        params%airmw / inst%mw
                   cycle
                endif
                ! Loss term?
                ! Don't branch. Just re-ask 'if'
                spc = inst%ispecies ! Species index
                if (sfc_flux(n)%flux(i,j) .lt. 0.e0 .and. aggregate(spc)%q(i,j,k).gt.0.e0) then ! Sink. Removes aggregate, not just one instance
                   if (cntrl%wellmixed_sfcexch) then
                      fdC = (sfc_flux(n)%flux(i,j) * sfc_flux(n)%scalefactor * grav / met%delp(i,j,k)) / aggregate(spc)%q(i,j,k) * params%airmw / inst%mw
                      if (inst%active) then
                         ! Active instance? Loop over all instances in the aggregate
                         do nst=1,NINSTANCES
                            if (instances(nst)%p%ispecies .eq. spc .and. instances(nst)%p%active) &
                                 instances(nst)%p%loss(i,j,k) = instances(nst)%p%loss(i,j,k)-fdC*instances(nst)%p%data3d(i,j,k) ! Pay attention to the sign! Losses should still be stored as positive numbers
                         enddo
                      endif ! not active, only remove the fraction of this instance (inst%q/aggregate%q)
                      if (.not. inst%active) then
                         inst%loss(i,j,k) = inst%loss(i,j,k)-fdC*inst%data3d(i,j,k)
                      endif
                   endif
                   if (.not. cntrl%wellmixed_sfcexch) then
                         inst%loss(i,j,k) = inst%loss(i,j,k) - &
                                        sfc_flux(n)%flux(i,j) * sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                        params%airmw / inst%mw
                   endif
                endif
             enddo
             enddo
          endif
       endif
       endif

       if (inst%hasmask) then
!      -- H A S  M A S K --
!  As of now, mask operations are not sparse. The full array is still processed
!  even if theres a multiply by zero. Implementing sparsity would, with only a few masks
!  promote load imbalance (if only slightly). When/if we find a lot of masks are being used
!  implementing a sparse version of this would not be hard and may be worth it.
       if (sfc_flux(n)%pblmix) then
          ! ASSUMPTION: any flux that is mixed in the PBL will be positive 
          if (sfc_flux(n)%diurnal) then
             do k=minPBL,params%km
                inst%prod(:,:,k)   = inst%prod(:,:,k) + inst%mask(:,:) * &
                                     sfc_flux(n)%flux(:,:) * fPBL(:,:,k) * fDNL(:,:) * &
                                     sfc_flux(n)%scalefactor * grav / met%delp(:,:,k)* &
                                     params%airmw / inst%mw
             end do
          endif
          if (.not. sfc_flux(n)%diurnal) then
             do k=minPBL,params%km
                inst%prod(:,:,k)   = inst%prod(:,:,k) + inst%mask(:,:) * &
                                     sfc_flux(n)%flux(:,:) * fPBL(:,:,k) * &
                                     sfc_flux(n)%scalefactor * grav / met%delp(:,:,k)* &
                                     params%airmw / inst%mw
             end do
          endif
       endif ! Non-PBL may be positive or negative (e.g. OCN, NEP)
       if (.not. sfc_flux(n)%pblmix) then
          k = params%km ! at the surface
          if (sfc_flux(n)%diurnal) then
             do j=1,params%jm
             do i=1,params%im
                ! Source term?
                if (sfc_flux(n)%flux(i,j) .gt. 0) then
                   inst%prod(i,j,k)   = inst%prod(i,j,k) + inst%mask(i,j) * &
                                                      sfc_flux(n)%flux(i,j) * fDNL(i,j) * &
                                                      sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                                      params%airmw / inst%mw
                   cycle
                endif
                ! Loss term?
                ! Don't branch. Just re-ask 'if'
                spc = inst%ispecies ! Species index
                if (sfc_flux(n)%flux(i,j) .lt. 0 .and. aggregate(spc)%q(i,j,k).gt.0.e0) then ! Sink. Removes aggregate, not just one instance
                   if (cntrl%wellmixed_sfcexch) then
                      fdC = inst%mask(i,j) * (sfc_flux(n)%flux(i,j) * fDNL(i,j) * sfc_flux(n)%scalefactor * grav / met%delp(i,j,k)) / aggregate(spc)%q(i,j,k) * params%airmw / inst%mw
                      if (inst%active) then
                         ! Active instance? Loop over all instances in the aggregate
                         do nst=1,NINSTANCES
                            if (instances(nst)%p%ispecies .eq. spc .and. instances(nst)%p%active) &
                                 instances(nst)%p%loss(i,j,k) = instances(nst)%p%loss(i,j,k)-fdC*instances(nst)%p%data3d(i,j,k) ! Pay attention to the sign! Losses should still be stored as positive numbers
                         enddo
                      endif ! not active, only remove the fraction of this instance (inst%q/aggregate%q)
                      if (.not. inst%active) then
                         inst%loss(i,j,k) = inst%loss(i,j,k)-fdC*inst%data3d(i,j,k)
                      endif
                   endif
                   if (.not. cntrl%wellmixed_sfcexch) then
                      inst%loss(i,j,k) = inst%loss(i,j,k) - inst%mask(i,j) * &
                                                  sfc_flux(n)%flux(i,j) * fDNL(i,j) * &
                                                  sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                                  params%airmw / inst%mw
                   endif
                endif
             enddo
             enddo
          endif
          if (.not. sfc_flux(n)%diurnal) then
             do j=1,params%jm
             do i=1,params%im
                ! Source term?
                if (sfc_flux(n)%flux(i,j) .gt. 0) then
                   inst%prod(i,j,k)   = inst%prod(i,j,k) + inst%mask(i,j) * &
                                                      sfc_flux(n)%flux(i,j) * &
                                                      sfc_flux(n)%scalefactor * &
                                                      grav / met%delp(i,j,k) * &
                                                      params%airmw / inst%mw
                   cycle
                endif
                ! Loss term?
                ! Don't branch. Just re-ask 'if'
                spc = inst%ispecies ! Species index
                if (sfc_flux(n)%flux(i,j) .lt. 0 .and. aggregate(spc)%q(i,j,k).gt.0.e0) then ! Sink. Removes aggregate, not just one instance
                   if (cntrl%wellmixed_sfcexch) then
                      fdC = inst%mask(i,j) * (sfc_flux(n)%flux(i,j) * sfc_flux(n)%scalefactor * grav /  met%delp(i,j,k)) / aggregate(spc)%q(i,j,k) * params%airmw / inst%mw
                      if (inst%active) then
                         ! Active instance? Loop over all instances in the aggregate
                         do nst=1,NINSTANCES
                            if (instances(nst)%p%ispecies .eq. spc .and. instances(nst)%p%active) &
                                 instances(nst)%p%loss(i,j,k) = instances(nst)%p%loss(i,j,k)-fdC*instances(nst)%p%data3d(i,j,k) ! Pay attention to the sign! Losses should still be stored as positive numbers
                         enddo
                      endif ! not active, only remove this instance
                      if (.not. inst%active) then
                         inst%loss(i,j,k) = inst%loss(i,j,k)-fdC*aggregate(spc)%q(i,j,k)!inst%data3d(i,j,k)
                      endif
                   endif
                   if (.not. cntrl%wellmixed_sfcexch) then
                      inst%loss(i,j,k) = inst%loss(i,j,k) - inst%mask(i,j) * &
                                                  sfc_flux(n)%flux(i,j) * &
                                                  sfc_flux(n)%scalefactor * grav / met%delp(i,j,k) * &
                                                  params%airmw / inst%mw
                   endif
                endif
             enddo
             enddo
          endif
       endif
!      ------------------       
    endif ! Mask/No-mask
    inst => null()
    end do

!   3) Cleanup
    if (allocated(fPBL)) deallocate(fPBL, stat=status)
    if (allocated(fDNL)) deallocate(fDNL, stat=status)

    rc = 0
    
    RETURN
  END SUBROUTINE surface_prodloss

!BOP
!
! !ROUTINE:  Chem_BiomassDiurnal - Applies diurnal cycle to biomass emissions.
!
! !INTERFACE:
  subroutine surface_DiurnalScaling ( Fout, params )

! !USES:

  IMPLICIT NONE

! !ARGUMENTS:

  real, intent(out)   :: Fout(:,:) ! Emissions scaling factor valid at NHMS
  type(parameters), pointer, intent(in) :: params

! !DESCRIPTION:
!
!      Applies diurnal cycle to biomass emissions.
!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  13nov2009  da Silva  First crack.
!  19Aug2020  E. Sherman - moved from Chem_UtilMod.F90 to process library
!
!EOP
!-------------------------------------------------------------------------
  
!      Hardwired diurnal cycle (multiplied by 100)
!      These numbers were derived from GOES-12
!      fire counts for 2003-2007.
!      -------------------------------------------
  integer, parameter :: N = 240
  real,    parameter :: DT = 86400. / N

!      Apply flat diurnal cycle for boreal forests as a
!      temporary solution to prevent very high aerosol
!      optical depth during the day
  real,    parameter :: Boreal(N) = 1.0
!      real,    parameter :: Boreal(N) = &
!      (/ 0.0277, 0.0292, 0.0306, 0.0318, 0.0327, 0.0335, &
!         0.0340, 0.0342, 0.0341, 0.0338, 0.0333, 0.0326, &
!         0.0316, 0.0305, 0.0292, 0.0278, 0.0263, 0.0248, &
!         0.0233, 0.0217, 0.0202, 0.0187, 0.0172, 0.0158, &
!         0.0145, 0.0133, 0.0121, 0.0110, 0.0100, 0.0091, &
!         0.0083, 0.0075, 0.0068, 0.0062, 0.0056, 0.0051, &
!         0.0046, 0.0042, 0.0038, 0.0035, 0.0032, 0.0030, &
!         0.0028, 0.0026, 0.0025, 0.0024, 0.0024, 0.0024, &
!         0.0024, 0.0026, 0.0027, 0.0030, 0.0033, 0.0036, &
!         0.0041, 0.0046, 0.0052, 0.0060, 0.0069, 0.0079, &
!         0.0090, 0.0104, 0.0119, 0.0137, 0.0157, 0.0180, &
!         0.0205, 0.0235, 0.0268, 0.0305, 0.0346, 0.0393, &
!         0.0444, 0.0502, 0.0565, 0.0634, 0.0711, 0.0794, &
!         0.0884, 0.0982, 0.1087, 0.1201, 0.1323, 0.1453, &
!         0.1593, 0.1742, 0.1900, 0.2069, 0.2249, 0.2439, &
!         0.2642, 0.2858, 0.3086, 0.3329, 0.3587, 0.3860, &
!         0.4149, 0.4455, 0.4776, 0.5115, 0.5470, 0.5840, &
!         0.6227, 0.6628, 0.7043, 0.7470, 0.7908, 0.8355, &
!         0.8810, 0.9271, 0.9735, 1.0200, 1.0665, 1.1126, &
!         1.1580, 1.2026, 1.2460, 1.2880, 1.3282, 1.3664, &
!         1.4023, 1.4356, 1.4660, 1.4933, 1.5174, 1.5379, &
!         1.5548, 1.5679, 1.5772, 1.5826, 1.5841, 1.5818, &
!         1.5758, 1.5661, 1.5529, 1.5365, 1.5169, 1.4944, &
!         1.4693, 1.4417, 1.4119, 1.3801, 1.3467, 1.3117, &
!         1.2755, 1.2383, 1.2003, 1.1616, 1.1225, 1.0832, &
!         1.0437, 1.0044, 0.9653, 0.9265, 0.8882, 0.8504, &
!         0.8134, 0.7771, 0.7416, 0.7070, 0.6734, 0.6407, &
!         0.6092, 0.5787, 0.5493, 0.5210, 0.4939, 0.4680, &
!         0.4433, 0.4197, 0.3974, 0.3763, 0.3565, 0.3380, &
!         0.3209, 0.3051, 0.2907, 0.2777, 0.2662, 0.2561, &
!         0.2476, 0.2407, 0.2352, 0.2313, 0.2289, 0.2279, &
!         0.2283, 0.2300, 0.2329, 0.2369, 0.2417, 0.2474, &
!         0.2536, 0.2602, 0.2670, 0.2738, 0.2805, 0.2869, &
!         0.2927, 0.2979, 0.3024, 0.3059, 0.3085, 0.3101, &
!         0.3107, 0.3102, 0.3087, 0.3061, 0.3026, 0.2983, &
!         0.2931, 0.2871, 0.2806, 0.2735, 0.2659, 0.2579, &
!         0.2497, 0.2412, 0.2326, 0.2240, 0.2153, 0.2066, &
!         0.1979, 0.1894, 0.1809, 0.1726, 0.1643, 0.1562, &
!         0.1482, 0.1404, 0.1326, 0.1250, 0.1175, 0.1101, &
!         0.1028, 0.0956, 0.0886, 0.0818, 0.0751, 0.0687 /)
  real,    parameter :: NonBoreal(N) = &
       (/ 0.0121, 0.0150, 0.0172, 0.0185, 0.0189, 0.0184, &
          0.0174, 0.0162, 0.0151, 0.0141, 0.0133, 0.0126, &
          0.0121, 0.0117, 0.0115, 0.0114, 0.0114, 0.0116, &
          0.0120, 0.0126, 0.0133, 0.0142, 0.0151, 0.0159, &
          0.0167, 0.0174, 0.0180, 0.0184, 0.0187, 0.0189, &
          0.0190, 0.0190, 0.0191, 0.0192, 0.0192, 0.0193, &
          0.0194, 0.0194, 0.0193, 0.0192, 0.0190, 0.0187, &
          0.0185, 0.0182, 0.0180, 0.0178, 0.0177, 0.0176, &
          0.0174, 0.0172, 0.0169, 0.0166, 0.0162, 0.0158, &
          0.0153, 0.0149, 0.0144, 0.0138, 0.0132, 0.0126, &
          0.0118, 0.0109, 0.0101, 0.0092, 0.0085, 0.0081, &
          0.0080, 0.0083, 0.0091, 0.0102, 0.0117, 0.0135, &
          0.0157, 0.0182, 0.0210, 0.0240, 0.0273, 0.0308, &
          0.0345, 0.0387, 0.0432, 0.0483, 0.0540, 0.0606, &
          0.0683, 0.0775, 0.0886, 0.1022, 0.1188, 0.1388, &
          0.1625, 0.1905, 0.2229, 0.2602, 0.3025, 0.3500, &
          0.4031, 0.4623, 0.5283, 0.6016, 0.6824, 0.7705, &
          0.8650, 0.9646, 1.0676, 1.1713, 1.2722, 1.3662, &
          1.4491, 1.5174, 1.5685, 1.6014, 1.6173, 1.6200, &
          1.6150, 1.6082, 1.6040, 1.6058, 1.6157, 1.6353, &
          1.6651, 1.7045, 1.7513, 1.8024, 1.8541, 1.9022, &
          1.9429, 1.9738, 1.9947, 2.0072, 2.0132, 2.0141, &
          2.0096, 1.9994, 1.9829, 1.9604, 1.9321, 1.8977, &
          1.8562, 1.8052, 1.7419, 1.6646, 1.5738, 1.4734, &
          1.3693, 1.2676, 1.1724, 1.0851, 1.0052, 0.9317, &
          0.8637, 0.8004, 0.7414, 0.6862, 0.6348, 0.5871, &
          0.5434, 0.5037, 0.4682, 0.4368, 0.4097, 0.3864, &
          0.3667, 0.3499, 0.3355, 0.3231, 0.3123, 0.3029, &
          0.2944, 0.2862, 0.2773, 0.2670, 0.2547, 0.2402, &
          0.2238, 0.2061, 0.1882, 0.1712, 0.1562, 0.1434, &
          0.1332, 0.1251, 0.1189, 0.1141, 0.1103, 0.1071, &
          0.1043, 0.1018, 0.0996, 0.0979, 0.0968, 0.0964, &
          0.0966, 0.0970, 0.0973, 0.0970, 0.0959, 0.0938, &
          0.0909, 0.0873, 0.0831, 0.0784, 0.0732, 0.0676, &
          0.0618, 0.0565, 0.0521, 0.0491, 0.0475, 0.0473, &
          0.0480, 0.0492, 0.0504, 0.0514, 0.0519, 0.0521, &
          0.0520, 0.0517, 0.0513, 0.0510, 0.0507, 0.0507, &
          0.0508, 0.0512, 0.0515, 0.0518, 0.0519, 0.0518, &
          0.0513, 0.0506, 0.0496, 0.0482, 0.0465, 0.0443, &
          0.0418, 0.0387, 0.0351, 0.0310, 0.0263, 0.0214 /)

!      Fixed normalization factors; a more accurate normalization would take
!      in consideration longitude and time step
!      ---------------------------------------------------------------------
  real*8, save :: fBoreal = -1., fNonBoreal = -1
  real,   save :: fDT=-1
  
  integer :: hh, mm, ss, ndt, i, j, k
  integer :: NN
  real :: secs, secs_local, aBoreal, aNonBoreal, alpha

  ! Below, this is done for cleanliness
  integer    :: im
  integer    :: jm
  integer    :: km
  real, pointer :: lons(:,:)
  real, pointer :: lats(:,:)
  integer    :: nhms
  real       :: cdt

  real       :: lat, lon ! local lat & lon in degrees
  
  im   = params%im
  jm   = params%jm
  km   = params%km
  cdt  = params%cdt
  nhms = params%nhms

!  lats => params%lats
!  lons => params%lons
!                              -----

!      Normalization factor depends on timestep
!      ----------------------------------------
!  if ( fDT /= cdt ) then
     fBoreal = 0.0
     fNonBoreal = 0.0
     NN = 0
     ndt = max(1,nint(cdt/DT))
     
     do k = 1, N, ndt
        NN = NN + 1
        fBoreal    = fBoreal    + Boreal(k)
        fNonBoreal = fNonBoreal + NonBoreal(k)
     end do
     
     fBoreal    = fBoreal / NN
     fnonBoreal = fnonBoreal / NN
!     fDT = cdt ! so it recalculates only if necessary
!  end if
  

!      Find number of secs since begining of the day (GMT)
!      ---------------------------------------------------
  hh = nhms/10000
  mm = (nhms - 10000*hh) / 100
  ss = nhms - 10000*hh - 100*mm
  secs = 3600.*hh + 60.*mm + ss
  
!      Apply factors depending on latitude
!      -----------------------------------
  do j = 1,jm
     do i = 1,im

        lat = params%lats(i,j)*params%radtodeg
        lon = params%lons(i,j)*params%radtodeg

!            Find corresponding index in hardwired diurnal cycle
!            240 = 24 * 60 * 60 secs / 360 deg
!            ---------------------------------------------------
        secs_local = secs + 240. * lon
        k = 1 + mod(nint(secs_local/DT),N)
        if ( k < 1 ) k = N + k

!            Apply diurnal cycle
!            -------------------
        aBoreal = Boreal(k) / fBoreal
        aNonBoreal = NonBoreal(k) / fNonBoreal
        
        if ( lat >= 50. ) then
           Fout(i,j) = aBoreal
        else if ( lat >= 30. ) then
           alpha = (lat - 30. ) / 20.
           Fout(i,j) = (1-alpha) * aNonBoreal + &
                alpha  * aBoreal
        else
           Fout(i,j) = aNonBoreal
        end if
     end do
  end do

end subroutine surface_DiurnalScaling
!==================================================================================

subroutine surface_pblmix( fPBL, met, params,  RC )
  ! Adapted from E. Nielsen's routine in GOCART

  implicit none

  ! Args
  type(meteorology), intent(in)  :: met
  type(parameters),  intent(in)  :: params
  integer,           intent(out) :: RC
  real,              intent(out) :: fPBL(:,:,:)

  ! Local
  real, allocatable              :: index(:), pblLayer(:,:)
  integer                        :: ios, i, j, k, kt

  ! Below, this is done for cleanliness
  integer    :: im
  integer    :: jm
  integer    :: km

  im = params%im
  jm = params%jm
  km = params%km
  
  rc    = 0

  ! Find the layer that contains the PBL.
  ! Layer thicknesses are ZLE(:,:,0:km).
  ! -------------------------------------
  ALLOCATE(index(0:km),STAT=ios)
  ALLOCATE(pblLayer(im,jm),STAT=ios)
  DO j=1,jm
     DO i=1,im
        index(0:km)=0
        WHERE(met%zle(i,j,0:km)-met%zle(i,j,km) > met%pblh(i,j)) index(0:km)=1
        pblLayer(i,j)=SUM(index)
     END DO
  END DO
  DEALLOCATE(index,STAT=ios)
  minPBL=MINVAL(pblLayer)

  ! Determine partitioning fraction based on layer thicknesses
  ! ----------------------------------------------------------
  fPBL(:,:,:)=0.00
  DO j=1,jm
     DO i=1,im
        kt=pblLayer(i,j)
        DO k=kt,km
           fPBL(i,j,k)=(met%zle(i,j,k-1)-met%zle(i,j,k))/(met%zle(i,j,kt-1)-met%zle(i,j,km))
        END DO
     END DO
  END DO

  ! Release memory
  ! --------------
  DEALLOCATE(pblLayer,STAT=ios)

end subroutine surface_pblmix

END MODULE Surface_Mod
