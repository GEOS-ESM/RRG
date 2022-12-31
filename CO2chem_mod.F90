module  CO2chem_mod

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  CO2_prodloss

CONTAINS

   subroutine CO2_prodloss( RC )

! !USES:

   implicit none

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   integer, intent(out) ::  rc

!EOP
!-------------------------------------------------------------------------

!  --------------------------
   rc = 0

   return
 end subroutine CO2_prodloss

end module CO2chem_mod
