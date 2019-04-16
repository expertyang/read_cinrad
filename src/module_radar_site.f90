module radar_site
use string
implicit none

type t_radar_site
   character(len=20) :: name
   character(len=6 ) :: id, type
   real :: lat, lon, alt
end type

type(t_radar_site), dimension(:), allocatable :: radar_info
integer :: nsite

integer, parameter :: radar_unit = 901
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_radar_site(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: i

    open(radar_unit, file=filename, status="old")
    write(*,"(2A)") "Reading:",trim(filename)
    call count_line(radar_unit,nsite)
    nsite=nsite-1 
    if(allocated(radar_info))then
       deallocate(radar_info)
    endif
    allocate(radar_info(nsite))

    read(radar_unit,*)
    do i=1, nsite
        read(radar_unit,*) radar_info(i)%id, radar_info(i)%name,&
                           radar_info(i)%lat,radar_info(i)%lon, &
                           radar_info(i)%alt,radar_info(i)%type
    enddo 
    close(radar_unit)
    !write(*,"(A)") "Read Complete!"
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_radar_site(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: i

    open(radar_unit, file=filename, status="unknown")
    write(*,"(2A)") "Writing:",trim(filename)

    write(radar_unit,"(A2,3X,', ',A4,16X,2(', ',9X,A3),', ',5X,A3,', ',A4,2X)") "ID","Name","Lat","Lon","Alt","Type"
    do i=1, nsite
        write(radar_unit,"(A5,', ',A20,2(', ',F12.4),', ',F8.1,', ',A6)") &
                         radar_info(i)%id, radar_info(i)%name,&
                         radar_info(i)%lat,radar_info(i)%lon, &
                         radar_info(i)%alt,radar_info(i)%type
    enddo 
    close(radar_unit)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function search_radar_site(id)
    implicit none
    character(len=*), intent(in) :: id
    integer :: i

    search_radar_site=-1
    do i=1, nsite
       if(id==radar_info(i)%id)then
          search_radar_site=i
          exit
       endif
    enddo
    end function
end module
