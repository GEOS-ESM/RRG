!#######################################################################
module M_memutils

  implicit none
  public

  real :: memused, memtotal, percent_used

  contains

  subroutine M_MemUsed ( memtotal, used, percent_used, RC )
     real, intent(out) :: memtotal, used, percent_used
     integer, optional, intent(OUT  ) :: RC

     character(len=32) :: meminfo   = '/proc/meminfo'
     character(len=32) :: string
     integer :: mem_unit
     real    :: multiplier, available

     character(len=64), parameter :: IAm="_MemUtils:_MemUsed"
     integer :: status

     available = -1
     memtotal = -1

     mem_unit = 999
     open(UNIT=mem_unit,FILE=meminfo,FORM='formatted',IOSTAT=STATUS)
     if (STATUS /= 0) then
        memtotal = 0.0
	used = 0.0
        percent_used = 0.0
        RC = 0
        return
     end if

     do
        read (mem_unit,'(a)', end=20) string
        if ( index ( string, 'MemTotal:' ) == 1 ) then  ! High Water Mark                                                                                                                                                                     
           read (string(10:LEN_trim(string)-2),*) memtotal
           multiplier = 1.0
           if (trim(string(LEN_trim(string)-1:)) == "kB" ) &
                multiplier = 1.0/1024. ! Convert from kB to MB                                                                                                                                                                                
           memtotal = memtotal * multiplier
        endif
        if ( index ( string, 'MemAvailable:' ) == 1 ) then  ! Resident Memory                                                                                                                                                                 
           multiplier = 1.0
           read (string(14:LEN_trim(string)-2),*) available
           if (trim(string(LEN_trim(string)-1:)) == "kB" ) &
                multiplier = 1.0/1024. ! Convert from kB to MB                                                                                                                                                                                
           available = available * multiplier
        endif
     enddo
20   close(mem_unit)

     if (memtotal >= 0 .and. available >= 0) then
        used = memtotal-available
        percent_used = 100.0*(used/memtotal)
     else
        ! fail, but don't crash                                                                                                                                                                                                               
        used = -1
        percent_used = -1
     end if

     RC = 0
     return
  end subroutine M_MemUsed

subroutine M_MemCommited ( memtotal, committed_as, percent_committed, RC )

real, intent(out) :: memtotal, committed_as, percent_committed
integer, optional, intent(OUT  ) :: RC

character(len=32) :: meminfo   = '/proc/meminfo'
character(len=32) :: string
integer :: mem_unit
real    :: multiplier

character(len=64), parameter :: IAm="MAPL_MemUtils:MAPL_MemCommited"
integer :: status


  multiplier = 1.0

  mem_unit = 999
  open(UNIT=mem_unit,FILE=meminfo,FORM='formatted',IOSTAT=STATUS)
  !_VERIFY(STATUS)                                                                                                                                                                                                                           

  if (STATUS /= 0) then
     memtotal = 0.0
     committed_as = 0.0
     percent_committed = 0.0
     RC = 0
     return
  end if

  do; read (mem_unit,'(a)', end=20) string
    if ( INDEX ( string, 'MemTotal:' ) == 1 ) then  ! High Water Mark                                                                                                                                                                         
      read (string(10:LEN_TRIM(string)-2),*) memtotal
      if (TRIM(string(LEN_TRIM(string)-1:)) == "kB" ) &
        multiplier = 1.0/1024. ! Convert from kB to MB                                                                                                                                                                                        
      memtotal = memtotal * multiplier
    endif
    if ( INDEX ( string, 'Committed_AS:' ) == 1 ) then  ! Resident Memory                                                                                                                                                                     
      read (string(14:LEN_TRIM(string)-2),*) committed_as
      if (TRIM(string(LEN_TRIM(string)-1:)) == "kB" ) &
        multiplier = 1.0/1024. ! Convert from kB to MB                                                                                                                                                                                        
      committed_as = committed_as * multiplier
    endif
  enddo
20 close(mem_unit)

   percent_committed = 100.0*(committed_as/memtotal)

 end subroutine M_MemCommited
end module M_memutils
