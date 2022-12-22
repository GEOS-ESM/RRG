module global_mod
  ! Global scope arrays, variables, and constants
  ! used throughout the GEOScarbon component but not specific
  ! to any coupling environment (ESMF/MAPL, etc.)
  !

  use types_mod
  
!  implicit none
  public 
  save

  type(aggregate),       pointer :: instances(:) => null()! A repository for all instances

  type(meteorology)              :: met
  type(surface_flux),    pointer :: sfc_flux(:)
  type(surface_flux),    pointer :: tmp_flux(:) ! Used for building the flux list
  type(parameters),      target  :: params

  integer                        :: total_instances = 0 ! over all gases. Initialized to zero

  real                           :: grav = 9.80665 ! This can be set at runtime if a different value is desired

  character(32), allocatable     :: species(:)
  integer                        :: nspecies = 0

end module global_mod
