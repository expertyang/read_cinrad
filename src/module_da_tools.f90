module da_tools
use config, only: missing_r
implicit none

! real :: missing_r = -888888.
 type t_loc
   integer :: i, j, k
   real    :: x, y, z, dx, dy, dxm ,dym, dz, dzm
 end type
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine da_get_radar_rv(ni, nj, nk, lat, lon, hgt, radar_lat, radar_lon, radar_alt, radar_rv)
   use wrf_data, only: wrf
   implicit none
   
   integer,                     intent(in)  :: ni, nj, nk
   real, dimension(ni, nj, nk), intent(in)  :: lat, lon, hgt
   real,                        intent(in)  :: radar_lat, radar_lon, radar_alt
   real, dimension(ni, nj, nk), intent(inout) :: radar_rv

   integer                          :: i, j, k, k1 
   real, dimension(ni, nj, nk)      :: model_u, model_v, model_w, model_qrn, model_p, model_ps

   type(t_loc), dimension(ni,nj,nk) :: loc
   type(t_loc)                      :: stn_loc
   real, dimension(:), allocatable  :: v_h
   real                             :: xr, yr, zr

   write(*,*) "wrf proj:",wrf%proj
   allocate(v_h(wrf%nz))
   call da_llxy(radar_lat, radar_lon, stn_loc)
   radar_rv=missing_r
   write(*,*) "size of loc wrf_psfc:",ubound(loc),ubound(wrf%psfc)
   do k=1, nk 
      do j=1, nj 
         do i=1, ni 
            if(lon(i,j,k)<-360.or.lon(i,j,k)>360) cycle
            call da_llxy(lat(i,j,k), lon(i,j,k), loc(i,j,k))
  !          write(603,*) "da_llxy:", i,j,k, lat(i,j,k), lon(i,j,k)
  !          write(*,*) "da_llxy:", i,j,k, lat(i,j,k), lon(i,j,k), loc(i,j,k)%i, loc(i,j,k)%j
            if((loc(i,j,k)%i==0).or.(loc(i,j,k)%j==0)) cycle
            ! write(*,*) loc(i,j,k)%x, loc(i,j,k)%i, loc(i,j,k)%dx, loc(i,j,k)%dxm
            ! write(*,*) loc(i,j,k)%y, loc(i,j,k)%j, loc(i,j,k)%dy, loc(i,j,k)%dym

            model_ps(i,j,k) = loc(i,j,k)%dym *( loc(i,j,k)%dxm * wrf%psfc(loc(i,j,k)%i,   loc(i,j,k)%j  )  &
                                               +loc(i,j,k)%dx  * wrf%psfc(loc(i,j,k)%i+1, loc(i,j,k)%j  )) &
                            + loc(i,j,k)%dy  *( loc(i,j,k)%dxm * wrf%psfc(loc(i,j,k)%i,   loc(i,j,k)%j+1)  &
                                               +loc(i,j,k)%dx  * wrf%psfc(loc(i,j,k)%i+1, loc(i,j,k)%j+1))

            do k1=1,wrf%nz
               v_h(k1) = loc(i,j,k)%dym *( loc(i,j,k)%dxm*wrf%h(loc(i,j,k)%i,  loc(i,j,k)%j,  k1)   &
                                          +loc(i,j,k)%dx *wrf%h(loc(i,j,k)%i+1,loc(i,j,k)%j,  k1))  &
                       + loc(i,j,k)%dy  *( loc(i,j,k)%dxm*wrf%h(loc(i,j,k)%i,  loc(i,j,k)%j+1,k1)   &
                                          +loc(i,j,k)%dx *wrf%h(loc(i,j,k)%i+1,loc(i,j,k)%j+1,k1))
            end do

            call da_to_zk(wrf%nz, hgt(i,j,k), v_h, loc(i,j,k))
            !write(901,*) "zk:", hgt(i,j,k), v_h, ";", loc(i,j,k)%z, loc(i,j,k)%k, loc(i,j,k)%dz, loc(i,j,k)%dzm

            call da_interp_lin_3d(wrf%nx,wrf%ny,wrf%nz,wrf%u  , loc(i,j,k), model_u  (i,j,k))
            call da_interp_lin_3d(wrf%nx,wrf%ny,wrf%nz,wrf%v  , loc(i,j,k), model_v  (i,j,k))
            call da_interp_lin_3d(wrf%nx,wrf%ny,wrf%nz,wrf%w  , loc(i,j,k), model_w  (i,j,k))
            call da_interp_lin_3d(wrf%nx,wrf%ny,wrf%nz,wrf%qrn, loc(i,j,k), model_qrn(i,j,k))
            call da_interp_lin_3d(wrf%nx,wrf%ny,wrf%nz,wrf%p  , loc(i,j,k), model_p  (i,j,k))


            xr = wrf%ds *(loc(i,j,k)%x - stn_loc%x)
            yr = wrf%ds *(loc(i,j,k)%y - stn_loc%y)
            zr = hgt(i,j,k) - radar_alt 

            if(model_u  (i,j,k)/=missing_r.and.&
               model_v  (i,j,k)/=missing_r.and.&
               model_w  (i,j,k)/=missing_r.and.&
               model_p  (i,j,k)/=missing_r.and.&
               model_ps (i,j,k)/=missing_r.and.&
               model_qrn(i,j,k)/=missing_r)then
               call da_radial_velocity(                &
                    model_p  (i,j,k), model_u (i,j,k), &
                    model_v  (i,j,k), model_w (i,j,k), &
                    model_qrn(i,j,k), model_ps(i,j,k), &
                    xr, yr, zr,       radar_rv(i,j,k))
            endif
            !write(901,*) "loc:",loc(i,j,k)%x,loc(i,j,k)%y,loc(i,j,k)%z,loc(i,j,k)%i,loc(i,j,k)%j,loc(i,j,k)%k,
            !write(901,*) loc(i,j,k)%dx,loc(i,j,k)%dy,loc(i,j,k)%dz,loc(i,j,k)%dxm,loc(i,j,k)%dym,loc(i,j,k)%dzm
            !write(901,*) i,j,k,model_u (i,j,k),model_v  (i,j,k),xr, yr, zr,       radar_rv(i,j,k)
         enddo
      enddo
   enddo

   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine da_radial_velocity(p,u,v,w,qrn,ps,x,y,z,rv)
   use config, only: if_terminal_velocity
   implicit none

   real, intent(in)  :: x, y, z
   real, intent(in)  :: p, u, v, w, qrn, ps
   real, intent(out) :: rv

   real :: r, alpha, vt
   real :: qrrc

   qrrc = 1.0e-3
   vt = 0.0  ! terminal velocity

   r=sqrt(x*x+y*y+z*z)
   alpha=(ps/p)**0.4

   if (if_terminal_velocity)then
      if (qrn <= qrrc)then
         vt=0.0
      else
         vt=5.4*alpha*qrn**0.125
      end if
   end if
   rv=u*x+v*y+(w-vt)*z
   rv=rv/r

   if(abs(rv)>100)then
      write(701,*) "da_radial_velocity",rv,r,x,y,z,u,v,w,vt
   endif
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine da_interp_lin_3d(nx,ny,nz,fm3d,loc,fo3d)
   implicit none

   integer,                   intent(in)    :: nx, ny, nz
   real, dimension(nx,ny,nz), intent(in)    :: fm3d     ! Input variable
   type(t_loc),               intent(in)    :: loc
   real,                      intent(inout) :: fo3d     ! Output variable


   integer             :: n, k
   real, dimension(nz) :: fmz


   fo3d=missing_r
   if (loc%i<=0.or.loc%j<=0.or.loc%k<=0) then
      return
   endif
   do k = 1, nz
      fmz(k) = loc%dym* ( loc%dxm*fm3d(loc%i,  loc%j,   k)  & 
                         +loc%dx *fm3d(loc%i+1,loc%j,   k)) &
             + loc%dy * ( loc%dxm*fm3d(loc%i,  loc%j+1, k)  &
                         +loc%dx *fm3d(loc%i+1,loc%j+1, k))
            
   end do
   fo3d = loc%dzm*fmz(loc%k) + loc%dz*fmz(loc%k+1)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! converts (lat, lon) into (x,y) coordinates
   subroutine da_llxy(lat, lon, loc)
   use map_utils, only: latlon_to_ij
   use wrf_data, only: wrf
   implicit none
   
   real,        intent(in)    :: lat, lon
   type(t_loc), intent(inout) :: loc

   call latlon_to_ij(wrf%proj, lat, lon, loc%x, loc%y)

   if(loc%x<1.or.loc%x>wrf%nx.or.loc%y<1.or.loc%y>wrf%ny)then ! outside domain
      loc%i  =0
      loc%dx =0
      loc%dxm=0
      loc%j  =0
      loc%dy =0
      loc%dym=0
      return
   endif
   
   loc%i  =int(loc%x)
   loc%dx =loc%x-loc%i
   loc%dxm=1.-loc%dx

   loc%j  =int(loc%y)
   loc%dy =loc%y-loc%j
   loc%dym=1.-loc%dy

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine da_to_zk(nz, obs_v, mdl_v, loc)
   implicit none

   integer,             intent(in)    :: nz
   real,                intent(in)    :: obs_v
   real, dimension(nz), intent(in)    :: mdl_v
   type(t_loc),         intent(inout) :: loc

   real                   :: zk
   integer                :: k, kts, kte
 
   zk = missing_r

   if (obs_v >= mdl_v(1) .and. obs_v <= mdl_v(nz)) then
      do k = 1, nz 
         if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
            zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
            exit
         end if
      end do
   end if

   if(zk==missing_r)then
      loc%z  =0
      loc%k  =0
      loc%dz =0
      loc%dzm=0
      return
   endif

   loc%z  = zk
   loc%k  = int(zk)
   loc%dz = zk-loc%k
   loc%dzm= 1.-loc%dz

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module
