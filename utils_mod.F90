module utils_mod

  use types_mod

  implicit none
  public

  contains

    subroutine util_accumulatenegatives( status )

      ! INCOMPLETE
      use global_mod, only: NINSTANCES, instances
      integer,      intent(out)      :: status

      integer i
      status = 0
      
      do i=1,NINSTANCES
         where(instances(i)%p%data3d .lt. 0.e0)
            ! Add the negative to the mass tracker
            ! TBD
            ! Set negative to zero
            instances(i)%p%data3d(:,:,:) = 0.d0
         end where
      end do

      return
    end subroutine util_accumulatenegatives

    subroutine util_addinstance( instance, name, gas, mw, active, status )

      use global_mod, only : NINSTANCES
      use global_mod, only : instances, aggregate
      use global_mod, only : species, nspecies

      ! This routine 
      ! 1) adds an instance to a blank or existing
      !    gas_instance object.
      ! 2) increments and manages the total instances object 'instances'
      ! 3) increments and manages the chemical species objects
      !    associated with instances
      !
      ! HISTORY:
      ! Dec. 9, 2022: M. Long - first pass. 

      ! arguments
      character(*), intent(in)       :: name, gas
      logical,      intent(in)       :: active
      real,         intent(in)       :: mw
      integer,      intent(out)      :: status
      type(gas_instance), pointer    :: instance(:)
      ! local
      integer                        :: i,n
      type(gas_instance), pointer    :: tmp(:)
      type(inst_), pointer           :: atmp(:)
      integer, allocatable           :: tind(:)
      character(32), allocatable     :: tspc(:)
      logical                        :: found

      n = 0

      ! store the current instance data locally
      if (associated(instance)) then
         n = size(instance)
         allocate(tmp(n), stat=status)
         tmp = instance
      endif

      ! Save the current index values
      if (n .gt. 0) then
         allocate(tind(n))
         tind = instance(:)%index
      endif

      ! free the instance object
      if (associated(instance)) deallocate(instance, stat=status)

      ! reallocate the instance and increment
      n = n+1
      NINSTANCES = NINSTANCES+1
      allocate(instance(n), stat=status)

      ! repopulate the incremented instance with the old data
      if (associated(tmp)) instance(1:n-1)    = tmp ! Pass the data back
      
      ! add the new instance info
      instance(n)%name    = trim(name)
      instance(n)%species = trim(gas)
      instance(n)%active  = active
      instance(n)%index   = NINSTANCES
      instance(n)%mw      = mw

      ! Restore any broken associations with the total instance object
      if (n .gt. 1) then
         do i=1,n-1
            instances(tind(i))%p => instance(i)
         enddo
      endif

      ! Set the species index (ispecies)
      if (nspecies .eq. 0) then ! no species yet
         nspecies = 1
         allocate(species(nspecies)) ! 1st species
         allocate(aggregate(nspecies))
         species(1) = gas
         instance(n)%ispecies = nspecies
      else ! 
         found = .false.
         do i=1,nspecies
            if (trim(gas) .eq. trim(species(i))) then
               found = .true.
               instance(n)%ispecies = i
               exit
            endif
         enddo
         ! if it wasn't found, increment & add a new species
         if (.not. found) then
            allocate(tspc(nspecies))
            tspc = species
            deallocate(species)
            nspecies = nspecies+1
            allocate(species(nspecies))
            species(1:nspecies-1) = tspc
            deallocate(tspc)
            species(nspecies) = gas
            instance(n)%ispecies = nspecies
            deallocate(aggregate)
            allocate(aggregate(nspecies))
         endif
      endif

      ! cleanup: eliminate the temporary instance
      deallocate(tmp, tind, stat=status)

      ! now increment the total instances object
      ! -- set up the temporary object to store data
      if (NINSTANCES .eq. 1) then
         allocate(instances(NINSTANCES), stat=status)
         instances(NINSTANCES)%p => instance(n)
      else   
         allocate(atmp(NINSTANCES), stat=status)
         do i=1,NINSTANCES-1
            atmp(i) = instances(i)
         enddo
         atmp(NINSTANCES)%p => instance(n)
         
         ! --increment the global object
         if (allocated(instances)) deallocate(instances, stat=status)
         allocate(instances(NINSTANCES), stat=status)
         
         ! Pass the saved data back
         do i=1,NINSTANCES
            instances(i) = atmp(i)
         enddo
         
         ! -- cleanup
         deallocate(atmp, stat=status)
      endif
      
      return

    end subroutine util_addinstance

    subroutine util_addsurfaceflux( pair, name, d, p, scalefactor, status)

      use global_mod
       
      character(*),                intent(in)    :: pair
      character(*),                intent(in)    :: name
      logical,                     intent(in)    :: d, p
      real,                        intent(in)    :: scalefactor
      integer,                     intent(out)   :: status

      integer :: i, n
      logical :: found

      type(surface_flux), pointer                :: tmp_flux(:) ! Used for building the flux list

      status = 0

      if (allocated(sfc_flux)) then
         n = size(sfc_flux)+1
      else
         n = 1
      endif

      if (n .gt. 1) then
         ! sfc_flux is already allocated size n-1
         allocate(tmp_flux(size(sfc_flux)), stat=status)
         if (status .ne. 0) return

         tmp_flux = sfc_flux ! data storage

         ! sfc_flux 'should' be allocated here
         if (allocated(sfc_flux)) deallocate(sfc_flux, stat=status)
         if (status .ne. 0) return
         allocate(sfc_flux(n), stat=status)
         if (status .ne. 0) return

         sfc_flux(1:n-1) = tmp_flux
         deallocate(tmp_flux, stat=status)
         if (status .ne. 0) return

         ! populate last entry
         sfc_flux(n)%instance_pair = pair
         sfc_flux(n)%shortname     = name
         sfc_flux(n)%diurnal       = d
         sfc_flux(n)%pblmix        = p
         sfc_flux(n)%scalefactor   = scalefactor
      else ! n=1
         ! neither flux types are allocated yet
         ! only need to add to sfc_flux()

         ! ... just in case
         if (allocated(sfc_flux)) deallocate(sfc_flux, stat=status)
         if (status .ne. 0) return
         allocate(sfc_flux(1), stat=status)
         if (status .ne. 0) return

         ! populate entry
         sfc_flux%instance_pair = pair
         sfc_flux%shortname     = name
         sfc_flux%diurnal       = d
         sfc_flux%pblmix        = p
         sfc_flux%scalefactor   = scalefactor
      endif

      ! Match surface flux's instance_pair to instance
      ! and set index
      found = .false.
      do i=1,size(instances)
         if (trim(sfc_flux(n)%instance_pair) .eq. trim(instances(i)%p%name)) then
            sfc_flux(n)%index = instances(i)%p%index ! Each instance can have more than one surface flux
            found = .true.
!DBG            write(*,*) '<<>> FOUND: ', trim(sfc_flux(n)%instance_pair),sfc_flux(n)%index
            exit
         endif
      enddo
      if (.not. found) then ! the instance isn't in the list!
         write(*,'(a)') 'GEOScarbon::surface flux: flux added without a registered instance,'
         write(*,'(a)') '                          Flux: '//trim(sfc_flux(n)%shortname)
         write(*,'(a)') '                          Instance: '//trim(sfc_flux(n)%instance_pair)
         write(*,'(a)') 'Potential causes: 1) fluxes were added before registering instances'
         write(*,'(a)') '                  2) the paired instance name is not a registered instance'
         status = -1
      endif

      return

    end subroutine util_addsurfaceflux

    subroutine util_aggregate( RC )

      use global_mod, only : instances, aggregate, NINSTANCES, nspecies

      integer, intent(out) :: RC
      
      integer :: i, n

      RC = 0

      do n=1,nspecies
         aggregate(n)%q = 0.e0
         do i=1,NINSTANCES
            if (instances(i)%p%ispecies .eq. n .and. instances(i)%p%active) aggregate(n)%q = aggregate(n)%q + instances(i)%p%data3d !
         enddo
      enddo

    end subroutine util_aggregate

    subroutine util_dumpinstances()
      use global_mod, only : NINSTANCES
      use global_mod, only : in => instances

      logical :: yesno
      integer :: i, n
      character(20) :: active

      ! 1) Check that the NINSTANCES equals the total number of 
      !    instances in object
      n = size(in)
      if (NINSTANCES .eq. n) then
         yesno = .true.
      else
         yesno = .false.
      endif

      ! First, dump info about the total instances
      !write(*,*) 'Does the index of total instances match the number of instances?'
      !if (yesno) write(*,*) 'It does! ', NINSTANCES
      !if (.not. yesno) write(*,*) 'It does not! Something is wrong: ', &
      !     NINSTANCES, ' .ne. ', n

      ! Dump instance info
      write(*,'(a)') ' name | gas | active? | index'
      do i=1,n
         if (in(i)%p%active) active = 'active'
         if (.not. in(i)%p%active) active = 'passive'
         write(*,'(a,i3)') '<<>> '//trim(in(i)%p%name)//' | '//trim(in(i)%p%species)//' | '//trim(active)//' | ', in(i)%p%index
      enddo

    end subroutine util_dumpinstances

    subroutine util_masksparsity( mask, inonzero, jnonzero, status )

      implicit none

      integer, intent(in)  :: mask(:,:)
      integer, pointer     :: inonzero(:)
      integer, pointer     :: jnonzero(:)
      integer, intent(out) :: status

      integer :: i,j,im,jm, nz, idx

      status = 0 ! Assume success

      ! Get dims
      im = size(mask,dim=1)
      jm = size(mask,dim=2)
      nz = count(mask.ne.0)

      if (associated(inonzero)) nullify(inonzero)
      if (associated(jnonzero)) nullify(jnonzero)

      allocate(inonzero(nz), jnonzero(nz), stat=status)
      if (status .ne. 0 ) return

      idx = 0
      do j=1,jm
         do i=1,im
            if (mask(i,j) .ne. 0 ) then
               idx = idx+1
               if (idx .gt. nz) then
                  write(*,*) 'util_masksparsity: improper mask indexing. idx > number of nonzero elements'
                  status = -1
                  return
               endif
               inonzero(idx) = i
               jnonzero(idx) = j
            endif
         enddo
      enddo
      
    end subroutine util_masksparsity

    function is_numeric(string)
      implicit none
      character(len=*), intent(in) :: string
      logical :: is_numeric
      real :: x
      integer :: e
      read(string,*,iostat=e) x
      is_numeric = e == 0
    end function is_numeric

    function ispecies(string)
      use global_mod, only : nspecies, species
      implicit none
      character(len=*), intent(in) :: string
      integer :: ispecies, i
      ispecies = -1 ! Default to 'not found'
      do i=1,nspecies
         if (trim(string) .eq. trim(species(i))) then
            ispecies = i
            return
         endif
      enddo
    end function ispecies

    function iinstance(string)
      use global_mod
      implicit none
      character(len=*), intent(in) :: string
      integer :: iinstance, i
      iinstance = -1 ! Default to 'not found'
      do i=1,NINSTANCES
         if (trim(string) .eq. trim(instances(i)%p%name)) then
            iinstance = i
            return
         endif
      enddo
    end function iinstance

end module utils_mod
