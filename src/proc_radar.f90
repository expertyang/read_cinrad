program proc_radar
use config, only: n_radar, radar_file, info_file, wrf_file, missing_r, &
                  read_namelist, if_use_wrf, if_radar_qc, &
                  if_output_qc, if_output_raw, if_output_obs, &
                  if_terminal_velocity, if_output_model, if_debug
use radar_data, only: t_radar_data, read_radar_data, get_radar_filename_info, &
                      get_dist_3d, write_radar_grads_station, value_invalid
use radar_site, only: radar_info, search_radar_site, &
                      read_radar_site, write_radar_site
use radar_qc, only: qc_radar, build_qc_radar, &
                    nzsnd, zsnd, usnd, vsnd, rfrsnd
use da_tools, only: da_get_radar_rv
use wrf_data, only: get_wrf_data
use radar_3dv, only: write_radar_3dv
implicit none

integer, parameter :: iunit=21
type(t_radar_data) :: rdat

character(len=800) :: filename 
character(len=6  ) :: radar_id, radar_type
character(len=14 ) :: radar_time
integer            :: nfile, i, j,  k, l, n, filesize, ierr, system, iargc

integer, dimension(3) :: vsize
real, dimension(:,:,:), allocatable :: tmp1, tmp2, bkgvel
real :: min_v, max_v
integer, dimension(201) :: dist
real, dimension(201) ::  pdist
real :: sumd, minr1,minr2, maxr1, maxr2

!type(date) :: radar_time

! read in namelist configuration
call read_namelist("namelist.input")

! read in china radar site information
call read_radar_site (info_file)
call write_radar_site("radar_info_in.txt")

if(if_use_wrf)then
   call get_wrf_data(wrf_file, ierr)
endif

open(401,file="radar_ref_dist.txt",status="unknown")
open(402,file="radar_vel_dist.txt",status="unknown")
open(403,file="radar_spw_dist.txt",status="unknown")
open(501,file="radar_vel.used.txt",status="unknown")

write(401,"(A,201(',',I20))")"val:",(i,i=-100,100)
write(402,"(A,201(',',I20))")"val:",(i,i=-100,100)
write(403,"(A,201(',',I20))")"val:",(i,i=-100,100)

! for each radar file
do l=1, n_radar
   filename=radar_file(l)
   call get_radar_filename_info(filename, radar_id, radar_time, radar_type)
   !  write(*,*) trim(filename)//",", trim(radar_id)//",", trim(radar_time)//",", trim(radar_type)

   ! search radar_id in radar site list
   n =search_radar_site(radar_id)
   if(n<0)then
      write(*,*) "ID:",radar_id,"Not found!"
      cycle
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

   call get_dist_3d(rdat%ref,dist)
   write(401,"(A,201(',',I20))") "raw:",dist
   call get_dist_3d(rdat%vel,dist)
   write(402,"(A,201(',',I20))") "raw:",dist
   call get_dist_3d(rdat%spw,dist)
   write(403,"(A,201(',',I20))") "raw:",dist

   rdat%radar_name=radar_info(n)%name
   rdat%radar_id  =radar_info(n)%id
   
   !call write_radar_data(trim(radar_id)//"."//trim(radar_time)//".dat",rdat)
   !!call write_radar_grads("radar",rdat)
   !call write_radar_csv          ("noqc."// trim(radar_id)//"."//trim(radar_time), rdat)
   if(if_output_raw) call write_radar_grads_station("raw."// trim(radar_id)//"."//trim(radar_time), rdat)

   write(*,*) allocated(tmp1),allocated(tmp2),allocated(bkgvel)
   if(allocated(tmp1))then
      deallocate(tmp1)
      deallocate(tmp2)
   endif
   if(allocated(bkgvel))then
      deallocate(bkgvel)
   endif
   vsize=ubound(rdat%vlat)
   write(*,*) "allocate bkgvel", vsize(1), vsize(2), vsize(3), size(bkgvel)
   allocate(bkgvel(vsize(1), vsize(2), vsize(3)))
   !bkgvel=0.
   ! dealiasing radical velocity
   if(if_use_wrf)then
      allocate(tmp1(vsize(1), vsize(2), vsize(3)))
      allocate(tmp2(vsize(1), vsize(2), vsize(3)))
      tmp1=value_invalid
      call da_get_radar_rv(vsize(1),        vsize(2),         vsize(3),        &
                           rdat%vlat,    rdat%vlon,     rdat%valt   , &
                           rdat%latitude,rdat%longitude,rdat%altitude,&
                           bkgvel)
     
      tmp2=rdat%vel
      !where(bkgvel==missing_r.or.bkgvel==value_ranfold.or.bkgvel==value_invalid) 
      !  bkgvel=rdat%vel
      !endwhere
      rdat%vel=bkgvel
      call get_dist_3d(bkgvel,dist)
      write(402,"(A,201(',',I20))") "model:",dist
    
      !do k=1, rdat%ntilt
      !   call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
      !   write(*,*) "bkgvel",k,min_v,max_v
      !enddo

      if(if_output_model) call write_radar_grads_station("model."// trim(radar_id)//"."//trim(radar_time), rdat)
      rdat%vel=tmp2
!     call unfold_radar(rdat, bkgvel)
!     if(if_output_qc) call write_radar_grads_station("unfold."// trim(radar_id)//"."//trim(radar_time), rdat)
   else
      bkgvel=rdat%vel(1:vsize(1),1:vsize(2),1:vsize(3))
   endif

!stop

   ! radar quality control
   if(if_radar_qc)then
      nzsnd=20
      if(allocated(zsnd))then
         deallocate(zsnd)
         deallocate(usnd)
         deallocate(vsnd)
         deallocate(rfrsnd)
      endif
      allocate(zsnd(nzsnd))
      allocate(usnd(nzsnd))
      allocate(vsnd(nzsnd))
      allocate(rfrsnd(nzsnd))
      zsnd=0
      usnd=0
      vsnd=0
      call qc_radar(rdat,bkgvel)
      !call write_radar_csv          ("radar."//trim(radar_id)//"."//trim(radar_time), rdat)
      call get_dist_3d(rdat%ref,dist)
      write(401,"(A,201(',',I20))") "qced:",dist
      call get_dist_3d(rdat%spw,dist)
      write(403,"(A,201(',',I20))") "qced:",dist
      call get_dist_3d(rdat%vel,dist)
      write(402,"(A,201(',',I20))") "qced:",dist

      sumd=sum(dist)
      pdist=dist/sumd
      minr1=-100
      do j=101,1,-1
         if(dist(j)<10.or.pdist(j)<3.E-4)then
            minr1=-100+j-0.5
             write(*,*) j, minr1, dist(j),pdist(j)
            exit
         endif
      enddo
      maxr1= 100
      do j=101,201
         if(dist(j)<10.or.pdist(j)<3.E-4)then
            maxr1=-100+j+0.5
            exit
         endif
      enddo
      !write(*,*) "minr,maxr",minr, maxr
      !where(rdat%vel<minr.or.rdat%vel>maxr)
      !   rdat%vel=missing_r
      !endwhere
      minr2=-100
      do j=2,101
         if(sum(pdist(1:j))>1.E-3)then
            minr2=-100+j-1.5
             write(*,*) j, minr2, dist(j),pdist(j)
            exit
         endif
      enddo
      maxr2= 100
      do j=200,101,-1
         if(sum(pdist(j:201))>1.E-3)then
            maxr2=-100+j+1.5
            exit
         endif
      enddo
      write(*,*) "minr1,minr2,maxr1,maxr2",minr1, minr2, maxr1, maxr2
      
      where(rdat%vel<max(minr1,minr2).or.rdat%vel>min(maxr1,maxr2))
         rdat%vel=missing_r
      endwhere
      if(if_output_qc) call write_radar_grads_station("qced."//trim(radar_id)//"."//trim(radar_time), rdat)
   else
      call build_qc_radar(rdat)
   endif

   ! final check
   where(rdat%ref<(20-0.004*(rdat%ralt-rdat%altitude)))
        rdat%ref=missing_r
   endwhere 
   where(rdat%spw<1.5.or.rdat%spw>8.)
        rdat%vel=missing_r
        rdat%spw=missing_r
   endwhere 
   where((rdat%ralt-rdat%altitude)>15000)
        rdat%ref=missing_r
   endwhere 
   where((rdat%valt-rdat%altitude)>15000)
        rdat%vel=missing_r
        rdat%spw=missing_r
   endwhere 
   call get_dist_3d(rdat%ref,dist)
   write(401,"(A,201(',',I20))") "limit:",dist
   call get_dist_3d(rdat%spw,dist)
   write(403,"(A,201(',',I20))") "limit:",dist
   call get_dist_3d(rdat%vel,dist)
   write(402,"(A,201(',',I20))") "limit:",dist
   call write_radar_grads_station("limit."//trim(radar_id)//"."//trim(radar_time), rdat)


   ! wrfda radar obs
   if(if_output_obs) call write_radar_3dv(trim(radar_id)//"."//trim(radar_time)//".3dv",rdat)

  ! if(if_use_wrf.and.if_radar_qc)then
  !   do k=1,rdat%ntilt
  !      if(rdat%ifvel(k))then
  !         do j=1,rdat%nazim(k)
  !            do i=1,rdat%nvgate(j,k)
  !               if(abs(rdat%vel(i,j,k))<200..and.abs(bkgvel(i,j,k))<200)then
  !                  write(501,"(A5,1X,A8,5F12.2,I8,2F12.2)") trim(radar_id),"radar rv", &
  !                        rdat%vlat(i,j,k),rdat%vlon(i,j,k),rdat%altitude,rdat%valt(i,j,k),rdat%vel(i,j,k),0,rdat%vel(i,j,k)-bkgvel(i,j,k),bkgvel(i,j,k) 
  !               endif
  !            enddo
  !         enddo
  !      endif
  !   enddo
  ! endif

   ! laps
   !call write_radar_laps(rdat)

   ! update radar site list
   if(rdat%ntilt>0)then
      radar_info(n)%lat=rdat%latitude 
      radar_info(n)%lon=rdat%longitude
      radar_info(n)%alt=rdat%altitude 
      radar_info(n)%type=radar_type
   endif
enddo
call write_radar_site("radar_info_out.txt")
close(401)
close(402)
close(403)
close(501)
end program
