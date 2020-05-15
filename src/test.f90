program read_cinrad
use radar_data, only: get_radar_filename_info, read_radar_data, &
                      write_radar_grads_station, t_radar_data, value_invalid, &
                      calculate_radar_loc
use radar_site, only: read_radar_site, search_radar_site, radar_info
implicit none

integer, parameter :: iunit=21
type(t_radar_data) :: rdat

character(len=800) :: filename 
character(len=6  ) :: radar_id, radar_type
character(len=14 ) :: radar_time
integer            :: i, j, k, n, filesize, iargc
real :: r, h

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
!  stop
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

call calculate_radar_loc(rdat)
 
   do k=1, rdat%ntilt
      if(rdat%ifref(k))then
         do j=1, rdat%nazim   (k)
            do i=1, rdat%nrgate(j,k)
               r=i*rdat%rgatesp(k)
               h=rdat%ralt(i,j,k)-rdat%altitude
               if(rdat%ref(i,j,k)<(20-0.004*h))then
               if (k==1.and.j==1) write(91,*) i, h, 20-0.004*h
               !if(rdat%ref(i,j,k)<(-80.+20.*log(r)/log(10.)).and.)then
                  rdat%ref(i,j,k)=-value_invalid
               endif

            enddo
         enddo
      endif
   enddo

call write_radar_grads_station("test."//trim(radar_id)//"."//trim(radar_time),rdat)


end program
