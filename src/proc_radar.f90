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
   real, dimension(:,:,:), allocatable :: tmp1, tmpvel, bkgvel

   ! stats var
   real :: min_v, max_v
   integer, dimension(201) ::  dist
   real,    dimension(201) :: pdist
   real :: sumd, minr1,minr2, maxr1, maxr2


   ! read in namelist configuration
   call read_namelist("namelist.input")
   
   ! read in china radar site information
   call read_radar_site (info_file)
   
   if(if_use_wrf)then
      call get_wrf_data(wrf_file, ierr)
   endif
   
   open(401,file="radar_ref_dist.txt",status="unknown")
   open(402,file="radar_vel_dist.txt",status="unknown")
   open(403,file="radar_spw_dist.txt",status="unknown")
   
   write(401,"(A,201(',',I20))")"val:",(i,i=-100,100)
   write(402,"(A,201(',',I20))")"val:",(i,i=-100,100)
   write(403,"(A,201(',',I20))")"val:",(i,i=-100,100)
   
   ! for each radar file
   do l=1, n_radar
      filename=radar_file(l)
      !filename like : Z_RADR_I_Z9417_20180806180100_O_DOR_SA_CAP.bin 
      call get_radar_filename_info(filename, radar_id, radar_time, radar_type)
   
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
         write(*,*) "Different radar_type in filename:",trim(radar_type)," and in radar_info.txt file:", trim(radar_info(n)%type),"!!!"
      endif
      call read_radar_data(filename,radar_type,rdat)
   
      ! raw distribution
      call get_dist_3d(rdat%ref,dist)
      write(401,"(A,201(',',I20))") "raw:",dist
      call get_dist_3d(rdat%vel,dist)
      write(402,"(A,201(',',I20))") "raw:",dist
      call get_dist_3d(rdat%spw,dist)
      write(403,"(A,201(',',I20))") "raw:",dist
   
      rdat%radar_name=radar_info(n)%name
      rdat%radar_id  =radar_info(n)%id
      
      ! output raw data for plot
      if(if_output_raw)then
         call write_radar_grads_station("raw."// trim(radar_id)//"."//trim(radar_time), rdat)
      endif
   
      vsize=ubound(rdat%vlat)
      write(*,*) "allocate bkgvel", vsize(1), vsize(2), vsize(3), size(bkgvel)
      allocate(bkgvel(vsize(1), vsize(2), vsize(3)))
      
      ! dealiasing radical velocity
      if(if_use_wrf)then
         allocate(tmpvel(vsize(1), vsize(2), vsize(3)))
         call da_get_radar_rv(vsize(1),        vsize(2),         vsize(3),        &
                              rdat%vlat,    rdat%vlon,     rdat%valt   , &
                              rdat%latitude,rdat%longitude,rdat%altitude,&
                              bkgvel)
        
         tmpvel  = rdat%vel
         rdat%vel= bkgvel

         ! model vel distribution
         call get_dist_3d(bkgvel,dist)
         write(402,"(A,201(',',I20))") "model:",dist
       
         if(if_output_model)then
            call write_radar_grads_station("model."// trim(radar_id)//"."//trim(radar_time), rdat)
         endif
         rdat%vel=tmpvel

         deallocate(tmpvel)
      else
         bkgvel=rdat%vel(1:vsize(1),1:vsize(2),1:vsize(3))
      endif
   
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

         ! qc radar
         call qc_radar(rdat,bkgvel)

         !qced distribution
         call get_dist_3d(rdat%ref,dist)
         write(401,"(A,201(',',I20))") "qced:",dist
         call get_dist_3d(rdat%spw,dist)
         write(403,"(A,201(',',I20))") "qced:",dist
         call get_dist_3d(rdat%vel,dist)
         write(402,"(A,201(',',I20))") "qced:",dist
   
         ! distribution qc vel 
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
         if(if_output_qc)then
            call write_radar_grads_station("qced."//trim(radar_id)//"."//trim(radar_time), rdat)
         endif
      else
         ! for wrfda radar obs 
         call build_qc_radar(rdat)
      endif
   
      ! final gross check
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
      if(if_output_obs)then
          call write_radar_3dv(trim(radar_id)//"."//trim(radar_time)//".3dv",rdat)
      endif
   
      ! update radar site list
      if(rdat%ntilt>0)then
         radar_info(n)%lat =rdat%latitude 
         radar_info(n)%lon =rdat%longitude
         radar_info(n)%alt =rdat%altitude 
         radar_info(n)%type=radar_type
      endif
      deallocate(bkgvel)
   enddo
   call write_radar_site("radar_info_out.txt")
   close(401)
   close(402)
   close(403)
end program
