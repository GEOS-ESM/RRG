module diagnostics

  use global_mod
  use types_mod
  use utils_mod

  implicit none
  public

  contains

    subroutine diag_prodloss( GI, PL, RC )

      ! Receives a gas_instance object and sums up its P-L into
      ! a diagnostic pointer
      implicit none
      
      type(gas_instance), pointer, intent(in)    :: GI(:)     ! gas_instance object
      real,               pointer, intent(inout) :: PL(:,:,:) ! diagnostic pointer
      integer,                     intent(out)   :: RC        ! return code

      integer :: i

      RC = 0   ! Assume success

      PL = 0.e0

      do i=1,size(GI)
         PL = PL + (GI(i)%Prod - GI(i)%Loss)
      enddo

      return

    end subroutine diag_prodloss
  
    subroutine diag_sfcflux( SPC, FLX, RC )

      implicit none

      character(*),  intent(in)    :: spc
      real, pointer, intent(inout) :: flx(:,:)
      integer,       intent(out)   :: RC

      integer :: ispc, idx, inst, i
      logical :: is_actv, is_spc, has_mask ! Convenience vars

      RC = 0 ! Assume success

      ispc = ispecies(trim(SPC))

      FLX = 0.e0

      do i=1,size(sfc_flux)
         inst     = sfc_flux(i)%index                  ! Which associated instance?
         is_spc   = instances(inst)%p%ispecies == ispc ! Is this instance the right species?
         if (.not. is_spc) cycle                       !    if not cycle
         is_actv  = instances(inst)%p%active           ! Is it an active instance?
         has_mask = instances(inst)%p%hasMask

         ! Currently only summing over active instances.
         if (      has_mask .and. is_actv .and. is_spc) FLX = FLX + sfc_flux(i)%flux*instances(inst)%p%mask
         if (.not. has_mask .and. is_actv .and. is_spc) FLX = FLX + sfc_flux(i)%flux
         
      enddo

    end subroutine diag_sfcflux
  
end module diagnostics
