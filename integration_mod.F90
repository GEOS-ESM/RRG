module integration_mod

  implicit none
  public

  contains

    subroutine integrate_forwardeuler( params, RC )

      use global_mod, only : instances, NINSTANCES
      use types_mod

      implicit none
      
      type(parameters), intent(in)  :: params
      integer,          intent(out) :: RC

      integer :: i

      RC = 0

      do i=1,NINSTANCES
         instances(i)%p%data3d = instances(i)%p%data3d + &
                                 (instances(i)%p%prod - instances(i)%p%loss) * params%cdt
      enddo

      RETURN

    end subroutine integrate_forwardeuler

end module integration_mod
