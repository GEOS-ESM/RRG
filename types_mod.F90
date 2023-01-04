module types_mod
  implicit none
  public
  
  type gas_instance
     real, pointer, DIMENSION(:,:,:) :: data3d ! Abundance
     real, pointer, dimension(:,:,:) :: prod, loss
     real                            :: mw ! Molecular weight
     logical                         :: active = .true.  ! Assume active by default
     character(len=255)              :: name, species
     integer                         :: index
     integer                         :: ispecies
     logical                         :: hasMask = .false.
     integer                         :: imask = -999
     integer, pointer                :: mask(:,:) => null()
     ! Sparse mask indices. For future use
     integer, pointer                :: mask_inonzero(:) => null() ! We want this to be saved
     integer, pointer                :: mask_jnonzero(:) => null() ! We want this to be saved
  end type gas_instance

  type inst_
     type(gas_instance), pointer     :: p
  end type inst_

  type aggr
     real, pointer                   :: q(:,:,:)
  end type aggr

  type meteorology ! met
     ! Imports
     real, pointer, dimension(:,:)   :: pblh  
     real, pointer, dimension(:,:)   :: pco2   ! for OCN. currenly unused
     real, pointer, dimension(:,:)   :: ps    
     real, pointer, dimension(:,:)   :: sss    ! for OCN. currenly unused
     real, pointer, dimension(:,:)   :: sst    ! for OCN. currenly unused
     real, pointer, dimension(:,:)   :: u10m  
     real, pointer, dimension(:,:)   :: v10m  
     real, pointer, dimension(:,:,:) :: t     
     real, pointer, dimension(:,:,:) :: delp
     real, pointer, dimension(:,:,:) :: ple
     real, pointer, dimension(:,:,:) :: zle   
     real, pointer, dimension(:,:,:) :: rho    ! air density
     real, pointer, dimension(:,:,:) :: q     
     real, pointer, dimension(:,:,:) :: qctot 
     real, pointer, dimension(:,:,:) :: qtot   ! Calculated locally. Not currently an available import.

     ! Computed quantities
     real, allocatable               :: cosz(:,:), slr(:,:) ! insolation params
     real, allocatable               :: O3col(:,:,:)
     real, allocatable               :: photJ(:,:,:)
  end type meteorology

  type surface_flux
     real, pointer, DIMENSION(:,:)   :: flux
     logical                         :: diurnal = .false. ! Assume no by default
     logical                         :: pblmix  = .false. ! Assume no by default
     character(len=255)              :: shortname
     character(len=255)              :: instance_pair
     integer                         :: index ! which gas_instance index is this flux associated with?
     real                            :: scalefactor = 1. ! Default to 1.
  end type surface_flux

  type parameters
     ! local parameters used all over the place
     ! declaring the integers at pointers permits cleanliness throughout
     integer                        :: im,jm,km
     real                           :: CDT         ! chemistry timestep (secs)
     real, pointer, dimension(:,:)  :: lats 
     real, pointer, dimension(:,:)  :: lons
     integer                        :: HDT         ! model timestep (secs)
     integer                        :: NYMD, NHMS
!     integer                        :: minPBL ! Min k index for PBL
     real                           :: AVO
     real                           :: AIRMW
  end type parameters
  
  type toggles
     logical                        :: strictMassBalance = .false. ! default false
  end type toggles

end module types_mod
