module integration_mod

  use MAPL

  implicit none
  public

  contains

    subroutine integrate_forwardeuler( RC )

      use global_mod
      use types_mod

      implicit none
      
      integer,          intent(out) :: RC

      integer :: i

      RC = 0

      do i=1,NINSTANCES
         instances(i)%p%data3d = max(instances(i)%p%data3d + &
                                 (instances(i)%p%prod - instances(i)%p%loss) * params%cdt, 0.e0)
      enddo

      RETURN

    end subroutine integrate_forwardeuler

end module integration_mod
