module radar_3dv
use libradar
implicit none

type t_point
    integer                            :: nlev
    real                               :: lat, lon, alt, u, v
    character(len=19)                  :: time
    real,    dimension(:), allocatable :: hgt, ref, vel, ref_err, vel_err, &
                                          ref_oma, ref_omb, vel_oma, vel_omb
    integer, dimension(:), allocatable :: ref_qc, vel_qc
end type

type t_radar
    character(len=12)                        :: id
    integer                                  :: nlev, npoint
    character(len=19)                        :: time
    real                                     :: lat, lon, alt
    type(t_point), dimension(:), allocatable :: point 
end type

!integer                                  :: nradar
!type(t_radar), dimension(:), allocatable :: radar

integer, parameter :: radar_unit=901

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allocate_point(point_dat)
   implicit none
   type(t_point), intent(inout) :: point_dat

   if(point_dat%nlev<1)return

   if(allocated (point_dat%hgt))then
      deallocate(point_dat%hgt    )
      deallocate(point_dat%ref    )
      deallocate(point_dat%ref_qc )
      deallocate(point_dat%ref_err)
      deallocate(point_dat%ref_omb)
      deallocate(point_dat%ref_oma)
      deallocate(point_dat%vel    )
      deallocate(point_dat%vel_qc )
      deallocate(point_dat%vel_err)
      deallocate(point_dat%vel_omb)
      deallocate(point_dat%vel_oma)
   endif
   allocate(point_dat%hgt    (point_dat%nlev))
   allocate(point_dat%ref    (point_dat%nlev))
   allocate(point_dat%ref_qc (point_dat%nlev))
   allocate(point_dat%ref_err(point_dat%nlev))
   allocate(point_dat%ref_omb(point_dat%nlev))
   allocate(point_dat%ref_oma(point_dat%nlev))
   allocate(point_dat%vel    (point_dat%nlev))
   allocate(point_dat%vel_qc (point_dat%nlev))
   allocate(point_dat%vel_err(point_dat%nlev))
   allocate(point_dat%vel_omb(point_dat%nlev))
   allocate(point_dat%vel_oma(point_dat%nlev))

   point_dat%ref_omb=0.
   point_dat%ref_oma=0.
   point_dat%vel_omb=0.
   point_dat%vel_oma=0.

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allocate_radar(radar_dat)
   implicit none
   type(t_radar), intent(inout) :: radar_dat

   if(radar_dat%npoint<1)return

   if(allocated (radar_dat%point))then
      deallocate(radar_dat%point)
   endif
   allocate(radar_dat%point(radar_dat%npoint))

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_radar_omb(filename, radar)
   implicit none
   character(len=*), intent(in) :: filename
   type(t_radar),    intent(out):: radar

   character(len=80  ) :: dat_fmt = "(2I8,A5,2F9.2,F17.7,2(2F17.7,I8,2F17.7))"
   integer             :: i, n, k, m
   character(len=10)   :: word

   open(radar_unit,file=filename,status="old")
   write(*,*) " Reading file:", trim(filename)
   read(radar_unit,*) word, radar%npoint
   write(*,*) "Point number in omb:", radar%npoint
   call allocate_radar(radar)

   do i=1, radar%npoint
      read(radar_unit,*) radar%point(i)% nlev
      call allocate_point(radar%point(i))
      do k=1, radar%point(i)% nlev
         read(radar_unit,fmt=dat_fmt) n, m, word, & 
                                      radar%point(i)% lat , &
                                      radar%point(i)% lon , &
                                      radar%point(i)% hgt    (k), &
                                      radar%point(i)% vel    (k), &
                                      radar%point(i)% vel_omb(k), &
                                      radar%point(i)% vel_qc (k), &
                                      radar%point(i)% vel_err(k), &
                                      radar%point(i)% vel_oma(k), &
                                      radar%point(i)% ref    (k), &
                                      radar%point(i)% ref_omb(k), &
                                      radar%point(i)% ref_qc (k), &
                                      radar%point(i)% ref_err(k), &
                                      radar%point(i)% ref_oma(k)
      enddo
   enddo
   close(radar_unit)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_netcdf(filename,radar) !,ifpassqc)
   use netcdf
   implicit none
   character(len=*),  intent(in)  :: filename
   type(t_radar),     intent(in)  :: radar
   !logical, optional, intent(out) :: ifpassqc

   real, parameter :: miss = -888888.0
   integer         :: i, n, k, m, ncid, ierr

   real,    allocatable, dimension(:)   :: lat, lon, hgt
   real,    allocatable, dimension(:,:) :: vel
   real,    allocatable, dimension(:,:) :: ref
   integer, allocatable, dimension(:)   :: lev

   integer :: nvar=5
   integer :: npoint, dim_npoint, dim_nvar, id_lat, id_lon, id_hgt, id_lev, & 
                      id_ref, id_ref_oma, id_ref_omb, id_ref_err, id_ref_qc, &
                      id_vel, id_vel_oma, id_vel_omb, id_vel_err, id_vel_qc
   logical            :: test_normal
   character(len=800) :: file3dv

   ierr=nf90_create(filename,nf90_clobber,ncid)
   IF(ierr==0)THEN
      WRITE(*,*) "Writing Netcdf File:",TRIM(filename)
   ELSE
      WRITE(*,*) "ERROR Open Netcdf File:",TRIM(filename)
      RETURN
   ENDIF

   npoint=0
   do n=1, radar%npoint
      npoint=npoint+radar%point(n)%nlev
   enddo

   if(allocated(lat))then
      deallocate(lat)
      deallocate(lon)
      deallocate(hgt)
      deallocate(vel)
      deallocate(ref)
      deallocate(lev)
   endif
   allocate(lat(npoint))
   allocate(lon(npoint))
   allocate(hgt(npoint))
   allocate(vel(nvar,npoint))
   allocate(ref(nvar,npoint))
   allocate(lev(npoint))

   m=0
   do n=1, radar%npoint
      do k=1,radar%point(n)%nlev
         m=m+1
         lat(m  )=radar%point(n)%lat
         lon(m  )=radar%point(n)%lon
         hgt(m  )=radar%point(n)%hgt    (k)
         vel(1,m)=radar%point(n)%vel    (k)
         vel(2,m)=radar%point(n)%vel_omb(k)
         vel(3,m)=radar%point(n)%vel_qc (k)
         vel(4,m)=radar%point(n)%vel_err(k)
         vel(5,m)=radar%point(n)%vel_oma(k)
         ref(1,m)=radar%point(n)%ref    (k)
         ref(2,m)=radar%point(n)%ref_omb(k)
         ref(3,m)=radar%point(n)%ref_qc (k)
         ref(4,m)=radar%point(n)%ref_err(k)
         ref(5,m)=radar%point(n)%ref_oma(k)
         lev(m  )=k
      enddo
   enddo
   
   ! define dimensions
   ierr = nf90_def_dim(ncid , 'npoint' , NF90_UNLIMITED, dim_npoint)
   ierr = nf90_def_dim(ncid , 'nvar'   , nvar, dim_nvar)
   !write(*,*) ierr, "nvar:",dim_nvar

  ! define Variables
   ierr = nf90_def_var(ncid, 'lat' , NF90_FLOAT , (/         dim_npoint/) , id_lat )
   ierr = nf90_def_var(ncid, 'lon' , NF90_FLOAT , (/         dim_npoint/) , id_lon )
   ierr = nf90_def_var(ncid, 'hgt' , NF90_FLOAT , (/         dim_npoint/) , id_hgt )
   ierr = nf90_def_var(ncid, 'lev' , NF90_INT   , (/         dim_npoint/) , id_lev )
   ierr = nf90_def_var(ncid, 'vel' , NF90_FLOAT , (/dim_nvar,dim_npoint/) , id_vel )
   !write(*,*) ierr, "vel:",id_vel, nf90_strerror(ierr)
   ierr = nf90_def_var(ncid, 'ref' , NF90_FLOAT , (/dim_nvar,dim_npoint/) , id_ref )
   !write(*,*) ierr, "ref:",id_ref, nf90_strerror(ierr)

   ! Variable Attributes 
   ierr = nf90_put_att(ncid, id_lat ,"_FillValue", miss )
   ierr = nf90_put_att(ncid, id_lon ,"_FillValue", miss )
   ierr = nf90_put_att(ncid, id_hgt ,"_FillValue", miss )
   ierr = nf90_put_att(ncid, id_vel ,"_FillValue", miss )
   ierr = nf90_put_att(ncid, id_ref ,"_FillValue", miss )
   
   ierr = nf90_put_att(ncid, id_lat ,"units", "degrees_north" )
   ierr = nf90_put_att(ncid, id_lon ,"units", "degrees_east" )
   ierr = nf90_put_att(ncid, id_hgt ,"units", "m" )
   ierr = nf90_put_att(ncid, id_vel ,"units", "m/s" )
   ierr = nf90_put_att(ncid, id_ref ,"units", "dB" )

   ierr = nf90_put_att(ncid, id_vel ,"comments", "obs,omb,qc,err,oma" )
   ierr = nf90_put_att(ncid, id_ref ,"comments", "obs,omb,qc,err,oma" )

   ierr = nf90_enddef(ncid)

   ! Put Variables Data
   ierr = nf90_put_var(ncid, id_lat , lat )
   ierr = nf90_put_var(ncid, id_lon , lon )
   ierr = nf90_put_var(ncid, id_hgt , hgt )
   ierr = nf90_put_var(ncid, id_lev , lev )
   ierr = nf90_put_var(ncid, id_vel , vel )
   ierr = nf90_put_var(ncid, id_ref , ref )

   ierr = nf90_close(ncid)
   !! test radar(n) normal distribution?
   !if(present(ifpassqc))then
   !   ifpassqc=test_normal(nvar,npoint,vel,miss)
   !endif
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_radar_3dv(filename, nradar, radar)
   implicit none
   character(len=*), intent(in)  :: filename
   integer,          intent(out) :: nradar
   type(t_radar), dimension(:), allocatable, intent(out) :: radar

   character(len=1024) :: line, word
   character(len=80  ) :: bhd_fmt, shd_fmt, dat_fmt
   integer             :: i, n, k
   real                :: dist, head

   open(radar_unit,file=filename,status="old")
   write(*,*) " Reading file:", trim(filename)
   read(radar_unit,*) line,word,nradar
   write(*,*) "Radar Number:", nradar
   if(allocated(radar))then
      deallocate(radar)
   endif
   allocate(radar(nradar))
   read(radar_unit,"(A)") line

   do n=1, nradar
      read(radar_unit,"(A)") bhd_fmt
      !write(*,*) "Big Header Format:", trim(bhd_fmt)
      read(radar_unit,fmt=bhd_fmt) word, radar(n)% id    , &
                                         radar(n)% lon   , &
                                         radar(n)% lat   , &
                                         radar(n)% alt   , &
                                         radar(n)% time  , &
                                         radar(n)% npoint, &
                                         radar(n)% nlev
      write(*,fmt=bhd_fmt) word, radar(n)% id    , &
                           radar(n)% lon   , &
                           radar(n)% lat   , &
                           radar(n)% alt   , &
                           radar(n)% time  , &
                           radar(n)% npoint, &
                           radar(n)% nlev
      call allocate_radar(radar(n))
      read(radar_unit,"(A)") shd_fmt
      !write(*,*) "Sub Header Format:", trim(shd_fmt)
      read(radar_unit,"(A)") dat_fmt
      !write(*,*) "Data Format:", trim(dat_fmt)
      do i=1, radar(n)%npoint
         read(radar_unit,fmt=shd_fmt) word, radar(n)%point(i)% time, &
                                            radar(n)%point(i)% lat , &
                                            radar(n)%point(i)% lon , &
                                            radar(n)%point(i)% alt , &
                                            radar(n)%point(i)% nlev
         !write(*,fmt=shd_fmt) word, radar(n)%point(i)% time, &
         !                     radar(n)%point(i)% lat , &
         !                     radar(n)%point(i)% lon , &
         !                     radar(n)%point(i)% alt , &
         !                     radar(n)%point(i)% nlev

         call disthead(radar(n)%point(i)%lat, radar(n)%point(i)%lon, radar(n)%lat, radar(n)%lon, head, dist)
         call ds2uv(head,1.,radar(n)%point(i)%u,radar(n)%point(i)%v)
         
         call allocate_point(radar(n)%point(i))
         do k=1, radar(n)%point(i)% nlev
            read(radar_unit,fmt=dat_fmt)  radar(n)%point(i)% hgt(k)    , &
                                          radar(n)%point(i)% vel(k)    , &
                                          radar(n)%point(i)% vel_qc(k) , &
                                          radar(n)%point(i)% vel_err(k), &
                                          radar(n)%point(i)% ref(k)    , &
                                          radar(n)%point(i)% ref_qc(k) , &
                                          radar(n)%point(i)% ref_err(k)
            !write(*,fmt=dat_fmt)  radar(n)%point(i)% hgt(k)    , &
            !                      radar(n)%point(i)% vel(k)    , &
            !                      radar(n)%point(i)% vel_qc(k) , &
            !                      radar(n)%point(i)% vel_err(k), &
            !                      radar(n)%point(i)% ref(k)    , &
            !                      radar(n)%point(i)% ref_qc(k) , &
            !                      radar(n)%point(i)% ref_err(k)
         enddo
      enddo
   enddo
   close(radar_unit)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine write_radar_3dv_dat(filename,radar)
!   implicit none
!   character(len=*), intent(in) :: filename
!   type(t_radar), intent(in) :: radar
!
!   character(len=1024)            :: line, word
!   character(len=80  ), parameter :: bhd_fmt="(A5,2X,A12,2(F8.3,2X),F8.1,2X,A19,2I6)", &
!                                     shd_fmt="(A12,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)", &
!                                     dat_fmt="( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )"
!   integer :: i, n, k
!
!   open(radar_unit,file=filename,status="unknown")
!   write(*,*) " Writing file:", trim(filename)
!   write(radar_unit,"(A,I4)") "Total Radar:",1
!   write(radar_unit,"(A)") "================================================================================"
!
!   write(radar_unit,"(A)") bhd_fmt
!   write(radar_unit,fmt=bhd_fmt) "RADAR", radar% id    , &
!                                          radar% lon   , &
!                                          radar% lat   , &
!                                          radar% alt   , &
!                                          radar% time  , &
!                                          radar% npoint, &
!                                          radar% nlev
!   write(radar_unit,"(A)") shd_fmt
!   write(radar_unit,"(A)") dat_fmt
!   do i=1, radar%npoint
!      write(radar_unit,fmt=shd_fmt) "FM-128      ", radar%point(i)% time, &
!                                                    radar%point(i)% lat , &
!                                                    radar%point(i)% lon , &
!                                                    radar%point(i)% alt , &
!                                                    radar%point(i)% nlev
!
!      do k=1, radar%point(i)% nlev
!
!         write(radar_unit,fmt=dat_fmt) radar%point(i)% hgt(k)    , &
!                                       radar%point(i)% vel(k)    , &
!                                       radar%point(i)% vel_qc(k) , &
!                                       radar%point(i)% vel_err(k), &
!                                       radar%point(i)% ref(k)    , &
!                                       radar%point(i)% ref_qc(k) , &
!                                       radar%point(i)% ref_err(k)
!      enddo
!   enddo
!   close(radar_unit)
!   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_3dv(filename, nradar, radar)
   implicit none
   character(len=*), intent(in) :: filename
   integer,          intent(in) :: nradar
   type(t_radar), dimension(nradar), intent(in) :: radar

   character(len=1024)            :: line, word
   character(len=80  ), parameter :: bhd_fmt="(A5,2X,A12,2(F8.3,2X),F8.1,2X,A19,2I6)", &
                                     shd_fmt="(A12,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)", &
                                     dat_fmt="( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )"
   integer :: i, n, k

   open(radar_unit,file=filename,status="unknown")
   write(*,*) " Writing file:", trim(filename)
   write(radar_unit,"(A,I4)") "Total Radar:",nradar
   write(radar_unit,"(A)") "================================================================================"

   do n=1, nradar
      write(radar_unit,"(A)") bhd_fmt
      write(radar_unit,fmt=bhd_fmt) "RADAR", radar(n)% id    , &
                                             radar(n)% lon   , &
                                             radar(n)% lat   , &
                                             radar(n)% alt   , &
                                             radar(n)% time  , &
                                             radar(n)% npoint, &
                                             radar(n)% nlev
      write(radar_unit,"(A)") shd_fmt
      write(radar_unit,"(A)") dat_fmt
      do i=1, radar(n)%npoint
         write(radar_unit,fmt=shd_fmt) "FM-128      ", radar(n)%point(i)% time, &
                                                       radar(n)%point(i)% lat , &
                                                       radar(n)%point(i)% lon , &
                                                       radar(n)%point(i)% alt , &
                                                       radar(n)%point(i)% nlev

         do k=1, radar(n)%point(i)% nlev

            write(radar_unit,fmt=dat_fmt) radar(n)%point(i)% hgt(k)    , &
                                          radar(n)%point(i)% vel(k)    , &
                                          radar(n)%point(i)% vel_qc(k) , &
                                          radar(n)%point(i)% vel_err(k), &
                                          radar(n)%point(i)% ref(k)    , &
                                          radar(n)%point(i)% ref_qc(k) , &
                                          radar(n)%point(i)% ref_err(k)
         enddo
      enddo
   enddo
   close(radar_unit)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine trim_radar(nradar,dkm,radar,radar_out)
   implicit none
   integer,                          intent(in)  :: nradar
   type(t_radar), dimension(nradar), intent(in)  :: radar
   type(t_radar), dimension(:), allocatable, intent(out) :: radar_out
   real,                             intent(in)  :: dkm

   integer :: i, j, k, n, nx, ny, np, i1, j1

   integer, dimension(:,:), allocatable :: idx
   real, dimension(:), allocatable :: x, y
   real :: rng, azm, minx, maxx, miny, maxy, dis1, dis2, x1, y1, dx

   write(*,*) "Triming Radar 3dv..."
   dx=dkm*1000.
   allocate(radar_out(nradar))
   do n=1, nradar
      if(allocated(x))then
         deallocate(x)
         deallocate(y)
      endif
      allocate(x(radar(n)%npoint))
      allocate(y(radar(n)%npoint))
      do i=1, radar(n)%npoint
         call disthead(radar(n)%lat,radar(n)%lon,radar(n)%point(i)%lat,radar(n)%point(i)%lon,azm,rng)    
         call rd2xy(rng,azm,x(i),y(i))
      enddo
      minx=floor  (minval(x)/dx)*dx
      maxx=ceiling(maxval(x)/dx)*dx
      miny=floor  (minval(y)/dx)*dx
      maxy=ceiling(maxval(y)/dx)*dx
      nx=(maxx-minx)/dx+1
      ny=(maxy-miny)/dx+1
      !write(*,*) "min,max:",minx,maxx,miny,maxy,nx,ny,minval(x),maxval(x),minval(y),maxval(y)
      if(allocated(idx))then
         deallocate(idx)
      endif 
      allocate(idx(nx,ny))
      idx=0
      do i=1, radar(n)%npoint
         i1=int((x(i)-minx)/dx+0.5)+1
         j1=int((y(i)-miny)/dx+0.5)+1
         if(i1<1.or.j1<1.or.i1>nx.or.j1>ny)then
            write(*,"(A,4F10.0,2I6,4F13.3)")"minmax:",minx,maxx,miny,maxy,nx,ny,minval(x),maxval(x),minval(y),maxval(y)
            write(*,"(A,I6,5F13.3,5I6)")"nearest",i, x(i), y(i), minx, miny, dx, i1, j1, nx, ny!, idx(i1,j1)
            cycle         
         endif
         if(idx(i1,j1)==0)then
            idx(i1,j1)=i
         else
            x1=x(idx(i1,j1))-(i1*dx+minx)
            y1=y(idx(i1,j1))-(j1*dx+miny)
            dis1=sqrt(x1*x1+y1*y1)
            x1=x(i)-(i1*dx+minx)
            y1=y(i)-(j1*dx+miny)
            dis2=sqrt(x1*x1+y1*y1)
!           write(*,*) "nearest",i, i1, j1, idx(i1,j1), dis1, dis2
            if(dis2<dis1) idx(i1,j1)=i
         endif
      enddo
      np=0
      do i=1, nx
         do j=1, ny
            if(idx(i,j)>0) np=np+1
         enddo
      enddo
      radar_out(n)%npoint=np
      radar_out(n)% id   =radar(n)% id    
      radar_out(n)% lon  =radar(n)% lon   
      radar_out(n)% lat  =radar(n)% lat   
      radar_out(n)% alt  =radar(n)% alt   
      radar_out(n)% time =radar(n)% time  
      radar_out(n)% nlev =radar(n)% nlev   

      call allocate_radar(radar_out(n))
      i1=0
      do i=1, nx
         do j=1, ny
            if(idx(i,j)>0)then
               i1=i1+1      
               radar_out(n)%point(i1)% time=radar(n)%point(idx(i,j))% time
               radar_out(n)%point(i1)% lat =radar(n)%point(idx(i,j))% lat 
               radar_out(n)%point(i1)% lon =radar(n)%point(idx(i,j))% lon 
               radar_out(n)%point(i1)% alt =radar(n)%point(idx(i,j))% alt 
               radar_out(n)%point(i1)% nlev=radar(n)%point(idx(i,j))% nlev

               call allocate_point(radar_out(n)%point(i1))
               radar_out(n)%point(i1)% hgt    =radar(n)%point(idx(i,j))% hgt
               radar_out(n)%point(i1)% vel    =radar(n)%point(idx(i,j))% vel
               radar_out(n)%point(i1)% vel_qc =radar(n)%point(idx(i,j))% vel_qc
               radar_out(n)%point(i1)% vel_err=radar(n)%point(idx(i,j))% vel_err
               radar_out(n)%point(i1)% ref    =radar(n)%point(idx(i,j))% ref
               radar_out(n)%point(i1)% ref_qc =radar(n)%point(idx(i,j))% ref_qc
               radar_out(n)%point(i1)% ref_err=radar(n)%point(idx(i,j))% ref_err

            endif 
         enddo
      enddo

     write(*,*) "Trim radar", trim(radar(n)%id),radar(n)%npoint,'->',radar_out(n)%npoint

   enddo
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_max(filename,nradar,radar)
   use libradar, only: quick_sort
   implicit none
   character(len=*),   intent(in) :: filename
   integer,          intent(in) :: nradar
   type(t_radar), dimension(nradar), intent(in) :: radar

   character(len=800) :: fileout, ctlfile, datfile, mapfile
   integer            :: i, j, k, n, ntotal
   real, parameter    :: miss=-888888.0

   INTEGER          :: NFLAG, NLEV
   CHARACTER(LEN=8) :: STID
   REAL             :: TIM, level

   real, dimension(:), allocatable ::  lat, lon, ref, vel, hgt, tmp
   integer, dimension(:), allocatable :: top 
   
   ntotal=0
   do n=1,nradar
      do i=1, radar(n)%npoint
         !ntotal=ntotal+radar(n)%point(i)% nlev 
         do j=1, radar(n)%point(i)% nlev
            write(901,*) n,i,j,radar(n)%point(i)%ref_qc(j),radar(n)%point(i)%ref_qc(j)>=0
            if(radar(n)%point(i)%ref_qc(j)>=0)then
               ntotal=ntotal+1 
            endif
         enddo 
      enddo 
   enddo
   allocate(lat(ntotal))
   allocate(lon(ntotal))
   allocate(ref(ntotal))
   allocate(vel(ntotal))
   allocate(hgt(ntotal))
   allocate(tmp(ntotal))
   allocate(top(ntotal))
 
   ref=miss
   vel=miss
   k=0
   do n=1,nradar
      do i=1, radar(n)%npoint
         do j=1, radar(n)%point(i)% nlev
            if(radar(n)%point(i)%ref_qc(j)>=0)then
               k=k+1
               lat(k)=radar(n)%point(i)%lat
               lon(k)=radar(n)%point(i)%lon
               ref(k)=radar(n)%point(i)%ref(j)
               vel(k)=radar(n)%point(i)%vel(j)
               hgt(k)=radar(n)%point(i)%hgt(j)
            endif
         enddo 
      enddo 
   enddo
   write(*,*) "in write_radar_max: sorting", ntotal
   tmp=ref
   do i=1, ntotal
      top(i)=i
   enddo
   call quick_sort(tmp,top,ntotal,1,ntotal)
   !!do i=1, ntotal
   !!   write(98 ,*) top(i), ref(top(i))
   !!enddo
   !call calc_top(-ref,ntotal,top,ntotal)
   !do i=1, ntotal
   !   write(99 ,*) top(i), ref(top(i))
   !enddo

   ctlfile=trim(filename)//".data.ctl"
   datfile=trim(filename)//".data.dat"
   mapfile=trim(filename)//".data.map"

   open(radar_unit,file=ctlfile,status="unknown")
   write(*,*) "write grads ctl file:", trim(ctlfile)
   write(radar_unit,"(A       )") "DSET ^"//trim(datfile)
   write(radar_unit,"(A       )") "DTYPE station"
   write(radar_unit,"(A       )") "OPTIONS sequential big_endian"
   write(radar_unit,"(A       )") "STNMAP "//trim(mapfile)
   write(radar_unit,"(A,F12.1 )") "UNDEF ", miss
   write(radar_unit,"(A       )") "TITLE  Station Data Sample"
   write(radar_unit,"(A       )") "TDEF   1 linear 00z01jan2000 1hr"
   write(radar_unit,"(A       )") "VARS 3"
   write(radar_unit,"(A,I4,A  )") "h    0  99  Height           (m)  "
   write(radar_unit,"(A,I4,A  )") "z    0  99  Reflectivity     (dB) "
   write(radar_unit,"(A,I4,A  )") "v    0  99  Radical Velocity (m/s)"
   write(radar_unit,"(A       )") "ENDVARS"
   close(radar_unit)

   open(radar_unit,file=datfile,form="unformatted",status="unknown")
   write(*,*) "write grads dat file:", trim(datfile)
   TIM = 0.0
   NFLAG = 1 ! If surface variables present.
   NLEV=1

   do k=1, ntotal
      write(STID,'(I8)') k
      write(radar_unit) STID,lat(top(k)),lon(top(k)),TIM,NLEV,NFLAG
      write(radar_unit) hgt(top(k)),ref(top(k)),vel(top(k))
   enddo
   NLEV = 0
   WRITE(radar_unit) STID,0.,0.,TIM,NLEV,NFLAG
   close(radar_unit)
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_point(filename,nradar,radar)
   implicit none
   character(len=*),   intent(in) :: filename
   integer,          intent(in) :: nradar
   type(t_radar), dimension(nradar), intent(in) :: radar

   character(len=800) :: fileout, ctlfile, datfile, mapfile
   integer            :: i, j, k, n
   real, parameter    :: miss=-888888.0

   INTEGER          :: NFLAG, NLEV
   CHARACTER(LEN=8) :: STID
   REAL             :: TIM, level

   real, dimension(:), allocatable ::  ref, vel, hgt
   
   NLEV=radar(1)%point(1)%nlev
   do n=1,nradar
      NLEV=max(nlev,maxval(radar(n)%point(:)%nlev))
   enddo
   allocate(ref(NLEV))
   allocate(vel(NLEV))
   allocate(hgt(NLEV))
 
   ctlfile=trim(filename)//".data.ctl"
   datfile=trim(filename)//".data.dat"
   mapfile=trim(filename)//".data.map"

   open(radar_unit,file=ctlfile,status="unknown")
   write(radar_unit,"(A       )") "DSET ^"//trim(datfile)
   write(radar_unit,"(A       )") "DTYPE station"
   write(radar_unit,"(A       )") "OPTIONS sequential big_endian"
   write(radar_unit,"(A       )") "STNMAP "//trim(mapfile)
   write(radar_unit,"(A,F12.1 )") "UNDEF ", miss
   write(radar_unit,"(A       )") "TITLE  Station Data Sample"
   write(radar_unit,"(A       )") "TDEF   1 linear 00z01jan2000 1hr"
   write(radar_unit,"(A       )") "VARS 5"
   write(radar_unit,"(A,I4,A  )") "h    ",NLEV,"  99  Height           (m)  "
   write(radar_unit,"(A,I4,A  )") "z    ",NLEV,"  99  Reflectivity     (dB) "
   write(radar_unit,"(A,I4,A  )") "rv   ",NLEV,"  99  Radical Velocity (m/s)"
   write(radar_unit,"(A,I4,A  )") "u    ",NLEV,"  99  U-Wind           (m/s)"
   write(radar_unit,"(A,I4,A  )") "v    ",NLEV,"  99  V-Wind           (m/s)"
   write(radar_unit,"(A       )") "ENDVARS"
   close(radar_unit)

   open(radar_unit,file=datfile,form="unformatted",status="unknown")
   TIM = 0.0
   NFLAG = 0 ! If surface variables present.

   do n=1, nradar
      do i=1, radar(n)%npoint
         ref=miss
         vel=miss
         hgt=miss
         write(STID,'(I3,I5)') n,i
         write(radar_unit) STID,radar(n)%point(i)%lat, radar(n)%point(i)%lon,TIM,NLEV,NFLAG
         do k=1, radar(n)%point(i)%nlev
            ref(k)=radar(n)%point(i)%ref(k)
            vel(k)=radar(n)%point(i)%vel(k)
            hgt(k)=radar(n)%point(i)%hgt(k)
         enddo
         do k=1, NLEV
            level=k
            !write(*,"(F6.0,5F12.1)") level,hgt(k),ref(k),vel(k),radar(n)%point(i)%u,radar(n)%point(i)%v

            write(radar_unit) level,hgt(k),ref(k),vel(k), &
                              vel(k)*radar(n)%point(i)%u, &
                              vel(k)*radar(n)%point(i)%v
         enddo
      enddo
   enddo
   NLEV = 0
   WRITE(radar_unit) STID,0.,0.,TIM,NLEV,NFLAG
   close(radar_unit)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine write_radar_point(filename)
!   implicit none
!   character(len=*),   intent(in) :: filename
!
!   character(len=800) :: fileout, ctlfile, datfile, mapfile
!   integer            :: i, j, k, n
!   real, parameter    :: miss=-888888.0
!
!   INTEGER          :: NFLAG, NLEV
!   CHARACTER(LEN=8) :: STID
!   REAL             :: TIM, level
!
!   real, dimension(:), allocatable ::  ref, vel, hgt
!   
!   NLEV=radar(1)%point(1)%nlev
!   do n=1,nradar
!      NLEV=max(nlev,maxval(radar(n)%point(:)%nlev))
!   enddo
!   allocate(ref(NLEV))
!   allocate(vel(NLEV))
!   allocate(hgt(NLEV))
! 
!   ctlfile=trim(filename)//".data.ctl"
!   datfile=trim(filename)//".data.dat"
!   mapfile=trim(filename)//".data.map"
!
!   open(radar_unit,file=ctlfile,status="unknown")
!   write(radar_unit,"(A       )") "DSET ^"//trim(datfile)
!   write(radar_unit,"(A       )") "DTYPE station"
!   write(radar_unit,"(A       )") "OPTIONS big_endian"
!   write(radar_unit,"(A       )") "STNMAP "//trim(mapfile)
!   write(radar_unit,"(A,F12.1 )") "UNDEF ", miss
!   write(radar_unit,"(A       )") "TITLE  Station Data Sample"
!   write(radar_unit,"(A       )") "TDEF   1 linear 00z01jan2000 1hr"
!   write(radar_unit,"(A       )") "VARS 5"
!   write(radar_unit,"(A,I4,A  )") "h    ",NLEV,"  99  Height           (m)  "
!   write(radar_unit,"(A,I4,A  )") "z    ",NLEV,"  99  Reflectivity     (dB) "
!   write(radar_unit,"(A,I4,A  )") "rv   ",NLEV,"  99  Radical Velocity (m/s)"
!   write(radar_unit,"(A,I4,A  )") "u    ",NLEV,"  99  U-Wind           (m/s)"
!   write(radar_unit,"(A,I4,A  )") "v    ",NLEV,"  99  V-Wind           (m/s)"
!   write(radar_unit,"(A       )") "ENDVARS"
!   close(radar_unit)
!
!   open(radar_unit,file=datfile,recordtype="stream",form="unformatted",status="unknown")
!   TIM = 0.0
!   NFLAG = 0 ! If surface variables present.
!
!   do n=1, nradar
!      do i=1, radar(n)%npoint
!         ref=miss
!         vel=miss
!         hgt=miss
!         write(STID,'(I3,I5)') n,i
!         write(radar_unit) STID,radar(n)%point(i)%lat, radar(n)%point(i)%lon,TIM,NLEV,NFLAG
!         do k=1, radar(n)%point(i)%nlev
!            ref(k)=radar(n)%point(i)%ref(k)
!            vel(k)=radar(n)%point(i)%vel(k)
!            hgt(k)=radar(n)%point(i)%hgt(k)
!         enddo
!         do k=1, NLEV
!            level=k
!            !write(*,"(F6.0,5F12.1)") level,hgt(k),ref(k),vel(k),radar(n)%point(i)%u,radar(n)%point(i)%v
!
!            write(radar_unit) level,hgt(k),ref(k),vel(k), &
!                              vel(k)*radar(n)%point(i)%u, &
!                              vel(k)*radar(n)%point(i)%v
!         enddo
!      enddo
!   enddo
!   NLEV = 0
!   WRITE(radar_unit) STID,0.,0.,TIM,NLEV,NFLAG
!   close(radar_unit)
!
!   end subroutine

end module
