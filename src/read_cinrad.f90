program read_cinrad
use radar_data
use radar_site
implicit none

integer, parameter :: iunit=21
type(t_radar_data) :: rdat

character(len=800) :: filename 
character(len=6  ) :: radar_id, radar_type
character(len=14 ) :: radar_time
integer            :: n, filesize, iargc

! read in radar_filename 
if(iargc()<1)then
   write(*,*) "Usage: read_cinrad radar_file"
   stop
endif
call getarg(1, filename)

! read in radar site information
call read_radar_site ("radar_info.txt")

call get_radar_filename_info(filename, radar_id, radar_time, radar_type)
!  write(*,*) trim(filename)//",", trim(radar_id)//",", trim(radar_time)//",", trim(radar_type)

! search radar_id in radar site list
n =search_radar_site(radar_id)
if(n<0)then
   write(*,*) "ID:",radar_id,"Not found!"
   stop
endif

rdat%latitude =radar_info(n)%lat
rdat%longitude=radar_info(n)%lon
rdat%altitude =radar_info(n)%alt

! read radar file 
filesize=0
inquire(file=filename,size=filesize)
write(*,"(A)") "=============================================================================="
write(*,"(A,I10,2F12.4,F8.1)") trim(filename), filesize, &
                               radar_info(n)%lat, radar_info(n)%lon, radar_info(n)%alt
if(radar_type/= radar_info(n)%type)then
   write(*,*) "Error radar_type in radar_info.txt file:", trim(radar_info(n)%type),"!!!"
endif
call read_radar_data(filename,radar_type,rdat)

end program
