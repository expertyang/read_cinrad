module radar_region
implicit none

private

real(kind=8), parameter :: PI         =3.141592653589793238
real(kind=8), parameter :: DEG2RAD    =(PI/180.)

! point flag: parameter
real, parameter :: vel0=0.5
real, parameter :: spw0=0.
real, parameter :: spw1=4.
real, parameter :: ref0=-83.
real, parameter :: rng0=2500.
real, parameter :: rng1=129000.
real :: v0, w0, w1, z0, r0, r1, at
! point flag:
integer, parameter :: & 
         flag_pt_unknown     =0  , & ! default
         flag_pt_used        =1  , & ! used point
         flag_pt_low_vel     =2  , & ! vel <=v0
         flag_pt_low_spw     =4  , & ! spw==w0
         flag_pt_high_spw    =8  , & ! spw >w1 
         flag_pt_low_snr     =16 , & ! low snr: ref< z0+20*log10(r), r unit is m
         flag_pt_near_radar  =32 , & ! range < r0: close to radar point 
         flag_pt_far_radar   =64     ! range > r1: far from radar point 

real,    parameter :: cv0=0.8, cv1=1., cv2=1.5, cv3=1., cv4=1.5, cv5=0
real               :: c0, c1, c2, c3, c4, c5, c6, a0=10, a1=100. !, s1=3.5, s2=0.4
integer            :: n0=2
integer            :: min_level, max_level

! region flag
integer, parameter :: flag_rg_used =1, flag_rg_unknown=0, flag_rg_delete=-1, flag_rg_zero=-2

! neighbor point flag
integer, parameter :: flag_cont     =  0, &
                      flag_neg_alias= -1, &
                      flag_pos_alias=  1, &
                      flag_mix      = -2, &
                      flag_noise    = -4, &
                      flag_zero     = -8, &
                      flag_unknown  =-16

! cluster flag
integer, parameter :: flag_cl_unknown=0, flag_cl_zero=-2

! output unit
integer, parameter :: cluster_unit=101, &
                       region_unit=102, &
                         sort_unit=103, &
                      nyquist_unit=104 , &
                    reference_unit=105, &
                   check_grid_unit=106, &
                    check_std_unit=107, &
                       border_unit=108
type t_grid
   integer :: num_range, num_azim, num_valid
   real    :: vmax, dr, da, area_total
   real,    dimension(:),   allocatable :: azim
   real,    dimension(:,:), allocatable :: vel, bkg, ref, spw, area, ref0
   integer, dimension(:,:), allocatable :: region, cluster, flag, level
   integer, dimension(:)  , allocatable :: start_i, end_i
end type

type t_point
   integer :: i, j
end type

type t_region
   integer :: num_point, num_neighbor, flag, level
   real    :: area_total, minv, maxv
   integer :: nyquist_number, base_border, base_region, cluster
   !logical :: if_removed
   integer,       dimension(:), allocatable :: neighbor_list ! num_neighbor
   integer,       dimension(:), allocatable :: neighbor_flag ! num_neighbor
   type(t_point), dimension(:), allocatable :: point_list ! num_point
end type

type t_border_list
   integer :: num_point
   integer, dimension(:), allocatable :: i1, j1 , i2, j2
end type

type t_cluster
   real    :: area_total, mean_omb, std_omb, minv, maxv
   integer :: num_region, num_point
   integer :: num_neighbor, nyquist_number
   integer, dimension(:),   allocatable :: region_list ! (num_region)
   integer, dimension(:,:), allocatable :: neighbor_list ! (num_neighbor,2)
   integer, dimension(:),   allocatable :: sorted_border_index !  (num_neighbor), sort num_border_effect big -> little
   integer, dimension(:),   allocatable :: num_border_point, num_border_effect, neighbor_flag !(num_neighbor
end type

integer :: num_cluster, num_region


real, private :: miss = -888888.
integer, parameter :: imiss=0

public :: dealias_region

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine init_grid(num_range, num_azim, azim, dr, vel, ref, spw, vmax, grid)
   use grd, only: write_grd, write_test_data 
   implicit none

   integer,                                   intent(in) :: num_range, num_azim
   real,                                      intent(in) :: vmax, dr
   real,    dimension(num_azim),              intent(in) :: azim
   real   , dimension(num_range,num_azim),    intent(in) :: vel, ref, spw
   type(t_grid  ),                            intent(out):: grid

   integer :: i, j, k, np, n, i1, j1
   real :: minv, maxv, r, c, z0r

   write(*,*) "v0,w0,w1,z0,r0,r1=",v0,w0,w1,z0,r0,r1
   if(allocated(grid%vel))then
      deallocate(grid%vel    )
      deallocate(grid%ref    )
      deallocate(grid%ref0   )
      deallocate(grid%spw    )
      deallocate(grid%area   )
      deallocate(grid%azim   )
      deallocate(grid%flag   )
      deallocate(grid%level  )
      deallocate(grid%region )
      deallocate(grid%start_i)
      deallocate(grid%  end_i)
      deallocate(grid%cluster)
   endif
   allocate(grid%start_i (num_azim))
   allocate(grid%  end_i (num_azim))
   allocate(grid%azim    (num_azim))
   allocate(grid%vel     (num_range,num_azim))
   allocate(grid%ref     (num_range,num_azim))
   allocate(grid%ref0    (num_range,num_azim))
   allocate(grid%spw     (num_range,num_azim))
   allocate(grid%area    (num_range,num_azim))
   allocate(grid%region  (num_range,num_azim))
   allocate(grid%cluster (num_range,num_azim))
   allocate(grid%flag    (num_range,num_azim))
   allocate(grid%level   (num_range,num_azim))

   grid%num_range=num_range
   grid%num_azim =num_azim
   grid%azim     =azim
   grid%dr       =dr
   grid%vel      =vel
   grid%ref      =ref
   grid%spw      =spw
   grid%vmax     =vmax

   !c0            =ceiling(cv0*grid%vmax)
   !c0            =10
   c0            =cv0*grid%vmax
   c1            =cv1*grid%vmax
   c2            =cv2*grid%vmax
   c3            =cv3*grid%vmax
   c4            =cv4*grid%vmax
   !c5            =ceiling(cv5*vmax)
   c5            =c1+cv5

   minv=minval(azim)
   maxv=maxval(azim)
   grid%da=(maxv-minv)/(num_azim-1)

   grid%flag   =flag_pt_unknown

   grid%start_i=0
   grid%  end_i=0
   
   grid%area_total=0.

   write(*,*) "init",grid%dr, grid%da
   ! flag and area
   c=grid%dr*grid%dr/1000./1000.*grid%da*DEG2RAD
   do j=1, num_azim
      do i=1, num_range
         if(is_valid_vel(vel(i,j)))then !/=miss)then
            if(grid%start_i(j)==0) grid%start_i(j)=i
            grid%end_i(j)=i
            grid%flag (i,j)=flag_pt_used
            grid%area (i,j)=i*c
            r=i*dr
            !1.low snr
            !if(grid%ref(i,j)/=miss.and.grid%ref(i,j)<(z0+(20*log10(r))+(r*at/1000.)))then
            !if(grid%ref(i,j)<grid%ref0(i,j))then
            !if(r<=r0)then
            !   z0r=z0
            !elseif(r<=r1)then
               z0r=(-100-z0)/(r1-r0)*(r-r0)+z0
            !else
            !   z0r=-100
            !endif
            !z0r=z0
            grid%ref0(i,j)=z0r+(20*log10(r))+(r*at/1000.)
            !if(grid%ref0(i,j)<10) grid%ref0(i,j)=10
            !grid%ref0(i,j)=z0-0.0008*r+20*log10(r)
            !if(r<r0)then 
            grid%ref0(i,j)=(0-z0)/(r0)*(r-0)+z0
            !else
            !grid%ref0(i,j)=z0
            !endif
            !if(grid%ref(i,j)<grid%ref0(i,j))then !.and.r<=r0
            !   grid%flag(i,j)=IOR (grid%flag(i,j),flag_pt_low_snr)
            !   grid%flag(i,j)=IEOR(grid%flag(i,j),flag_pt_used)
            !endif
            !2.near radar
            !if(r<=r0)then
            !   grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_near_radar)
            !endif
            !3.far radar
            !if(r>r1)then
            !   grid%flag(i,j)= IOR(grid%flag(i,j),flag_pt_far_radar)
            !   grid%flag(i,j)=IEOR(grid%flag(i,j),flag_pt_used)
            !endif
            !4.low vel
            !if(IAND(grid%flag(i,j),flag_pt_low_snr)>0.and.abs(grid%vel(i,j))<=v0.and.r<=r0)then
            !if(abs(grid%vel(i,j))<=v0)then
            !   grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
            !endif
            !if(grid%vel(i,j)<=v0.and.grid%ref(i,j)<(z0+(20*log10(r))+(r*at/1000.)))then
            !if(grid%vel(i,j)<=v0)then
            !   grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
            !endif
            !5.low spw
            !if(grid%spw(i,j)<=w0)then
            !   grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_spw)
            !endif
            !6.high spw
            !if(grid%spw(i,j)>w1)then
            !   grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_high_spw)
            !endif
         endif
      enddo
   enddo
   call write_test_data('flag.dat',num_range,num_azim,azim,dr,grid%flag,flag_pt_unknown)
   call write_grd      ('flag.grd',num_range,num_azim,grid%flag,flag_pt_unknown)
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   logical function is_valid_vel(vel)
   implicit none
   real, intent(in) :: vel
   if(abs(vel)<200)then
     is_valid_vel=.true.
   else
     is_valid_vel=.false.
   endif
   end function

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine create_region(grid, num_region, region)
   use grd, only: write_grd, write_test_data
   implicit none

   type(t_grid  ),                            intent(inout) :: grid
   integer,                                   intent(out):: num_region
   type(t_region), dimension(:), allocatable, intent(out):: region

   integer, dimension(:,:), allocatable :: l, r, u, d
   integer :: i, j, k, np, n, i1, j1, l0

   real    :: min_vel, v
   integer :: num_valid, nlevel
   real, dimension(:), allocatable :: v_min, v_max
   logical :: if_l0=.false., low_vel

   allocate(  l(grid%num_range, grid%num_azim))
   allocate(  r(grid%num_range, grid%num_azim))
   allocate(  u(grid%num_range, grid%num_azim))
   allocate(  d(grid%num_range, grid%num_azim))
   

   write(*,"(A,8F8.3)") "c0,c1,c2,c3,c4,c5=", c0,c1,c2,c3,c4,c5,grid%dr,grid%da
   min_vel=floor(-grid%vmax/c0)*c0
   nlevel =floor((grid%vmax-min_vel)/c0)+1
   if(if_l0)then
      nlevel=nlevel+1
      n=nlevel-1
   else
      n=nlevel
   endif

   allocate(v_min(nlevel))
   allocate(v_max(nlevel))

   
   do k=1, n
      v_min(k)=min_vel+(k-1)*c0
      v_max(k)=v_min(k)+c0
      if(v_min(k)==0)then
         l0=k
      endif
   enddo
   if(if_l0)then
      v_max(l0+1:nlevel)=v_max(l0:nlevel-1)
      v_min(l0+1:nlevel)=v_min(l0:nlevel-1)
      v_max(l0-1)=-v0
      v_min(l0  )=-v0
      v_max(l0  )= v0
      v_min(l0+1)= v0
   endif

   grid%region =imiss
   grid%cluster=imiss

   num_valid=0
   grid%level=miss
   !min_level=1
   !max_level=nlevel
   grid%area_total=0.
   do j=1, grid%num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         if(IAND(grid%flag (i,j),flag_pt_used)>0)then
            grid%level(i,j)=floor((grid%vel(i,j)-min_vel)/c0)+1
            if(if_l0)then
            !add a -v0<=v<=v0 level
               if(grid%vel(i,j)>=-v0.and.grid%vel(i,j)<0)then
                  grid%level(i,j)=grid%level(i,j)+1
               endif
               if(grid%vel(i,j)>v0)then
                  grid%level(i,j)=grid%level(i,j)+1
               endif
            endif
            
            grid%area_total=grid%area_total+grid%area  (i,j)
            num_valid=num_valid+1
         endif
      enddo
   enddo
   
   call write_test_data('level.dat',grid%num_range,grid%num_azim,grid%azim,grid%dr,grid%level,int(miss))
   call write_grd      ('level.grd',grid%num_range,grid%num_azim,grid%level,int(miss))

!   write(*,*) "in create region:", c1
   write(*,*) "num valid:" , num_valid
   write(*,*) "area total:", grid%area_total
   grid%num_valid=num_valid
   if (grid%num_valid==0) return

   n=0
   do k=1, nlevel
      write(*,"(2(A,I3),A,2F8.3,A,I8)") "level=", k, "/", nlevel, " var=", v_min(k), v_max(k), " n=",  n

      l=0
      r=0
      u=0
      d=0
      do j=1, grid%num_azim
         if(grid%start_i(j)<1) cycle
         do i=grid%start_i(j), grid%end_i(j)
            if(grid%level(i,j)==k)then
               i1=i+1
               j1=j+1
               if(j1>grid%num_azim) j1=1

               if(j1<=grid%num_azim)then
                  if(grid%level(i,j1)==k)then
                     r(i,j)=j1
                     !if(iand(grid%flag(i,j),flag_pt_low_vel)/=iand(grid%flag(i,j1),flag_pt_low_vel))then
                     !   r(i,j)=0
                     !endif
                  endif
               endif
               if(i1<=grid%num_range)then
                  if(grid%level(i1,j)==k)then
                     u(i,j)=i1
                     !if(iand(grid%flag(i,j),flag_pt_low_vel)/=iand(grid%flag(i1,j),flag_pt_low_vel))then
                     !   u(i,j)=0
                     !endif
                  endif
               endif
               if(r(i,j)>0) l(i,   r(i,j))=j
               if(u(i,j)>0) d(   u(i,j),j)=i
               
            endif
         enddo
      enddo

      do j=1, grid%num_azim
         if(grid%start_i(j)<1) cycle
         do i=grid%start_i(j), grid%end_i(j)
            if(grid%level(i,j)/=k) cycle
            if(grid%region(i,j)==imiss)then
               n=n+1
               call include_point(i,j,n,l,r,u,d,grid%region)
            endif
         enddo
      enddo
   enddo
   num_region=n

   write(*,*) "create ", num_region, "regions."
   call update_region(grid, num_region, region)

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   recursive subroutine include_point(i, j, reg_idx, l, r, u, d, region)
   implicit none
   integer,                 intent(in)    :: i, j, reg_idx
   integer, dimension(:,:), intent(in)    :: l, r, u, d
   integer, dimension(:,:), intent(inout) :: region

   if(region(i,j)/=imiss) return
   !write(121,*) i, j, reg_idx
   region(i,j)=reg_idx
   if(    r(i,j)>0) call include_point(  i,     r(i,  j), reg_idx, l, r, u, d, region)
   if(    u(i,j)>0) call include_point(u(i,   j),     j , reg_idx, l, r, u, d, region)
   if(    l(i,j)>0) call include_point(  i,     l(i,  j), reg_idx, l, r, u, d, region)
   if(    d(i,j)>0) call include_point(d(i,   j),     j , reg_idx, l, r, u, d, region)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine update_region(grid, num_region, region)
   use libradar, only: quick_sort
   use grd, only: write_grd, write_test_data, region_name
   implicit none
   type(t_grid)  ,                        intent(inout) :: grid
   integer,                               intent(inout) :: num_region
   type(t_region), dimension(:), allocatable, intent(inout) :: region

   integer :: i, j, k, m, ierr

   real,    dimension(:), allocatable :: a, w
   integer, dimension(:), allocatable :: d, np
   integer :: n, s, e
   integer, dimension(num_region) :: flag


   n=num_region

   allocate(a (n)) 
   allocate(d (n)) 
   allocate(np(n)) 
   allocate(w (n)) 
   write(*,*) "update region:", n

   ! count area
   m=0
   do k=1, n 
      d (k)=k
      np(k)=0
      w (k) =0
      do j=1,  grid%num_azim
         if(   grid%start_i(j)<1) cycle
         do i= grid%start_i(j), grid%end_i(j)
            if(grid%region(i,j)==k)then
               np(k)=np(k)+1
               w (k)=w (k) +grid%area(i,j)
            endif
         enddo
      enddo
      a(k)=w(k)
   enddo
   !write(193,*) num_region, a

   !sort
   !a=-a
   s=1
   e=n

   ! sort area
   ! big first, integer sort
   !call quick_sort(a,d,n,s,e)
   call quick_sort(a,d,n,s,e)
   !a=-a
   do k=1, n 
      if(a(k)==0)then
         num_region=k-1
         exit
      endif
   enddo
!stop

   call set_region(grid,imiss,0)
!   call write_grd_int('reg2.grd',grid%num_range,grid%num_azim,grid%region,int(miss))
   do k=1, num_region
      call set_region(grid,d(k),-k)
   enddo
!   call write_grd_int('reg3.grd',grid%num_range,grid%num_azim,grid%region,int(miss))

   grid%region=-grid%region
!   call write_grd_int('reg4.grd',grid%num_range,grid%num_azim,grid%region,int(miss))
   call set_region(grid,0,imiss)
!   call write_grd_int('reg5.grd',grid%num_range,grid%num_azim,grid%region,int(miss))

!   where(grid%region>30)
!        grid%region=miss
!   endwhere
!   num_region=30

   if(allocated(region))then
      deallocate(region)
   endif
   allocate(region(num_region))
   do n=1, num_region
      if(allocated (region(n)%point_list))then
         deallocate(region(n)%point_list)
      endif
      allocate(region(n)%point_list(np(d(n))))
   enddo
!
   np=0
    w=0
   flag=1
region(:)%minv=miss
region(:)%maxv=miss
   do j=1,  grid%num_azim
      if(   grid%start_i(j)<1) cycle
      do i= grid%start_i(j), grid%end_i(j)
         if(grid%region(i,j)<1) cycle
         k=grid%region(i,j)
         np(k)=np(k)+1
         w (k)=w (k)+grid%area(i,j)
         region(k)%point_list(np(k))%i=i
         region(k)%point_list(np(k))%j=j
         if(region(k)%maxv==miss)then
            region(k)%maxv=grid%vel(i,j)
         else
            region(k)%maxv=max(grid%vel(i,j),region(k)%maxv)
         endif
         if(region(k)%minv==miss)then
            region(k)%minv=grid%vel(i,j)
         else
            region(k)%minv=min(grid%vel(i,j),region(k)%minv)
         endif
         if(grid%flag(i,j)>0.and.abs(grid%vel(i,j))>v0)then
            flag(k)=0
         endif
      enddo
   enddo
!      write(203,*) n, np 
   write(sort_unit,"(A)") "==== Region Sort ====" 
   region(:)%flag          =flag_rg_unknown
   do k=1, num_region
      if(flag(k)==1)then
         region(k)%flag=flag_rg_zero
      endif
      region(k)%num_point =np(k)
      region(k)%area_total=w (k)
      write(sort_unit,"(A,2(I10,F12.3))") region_name(k), d(k), a(k), np(k), w(k)
   enddo

   deallocate(a)
   deallocate(d)
   deallocate(np)
   deallocate(w)
   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine create_cluster(grid, num_region, region, num_cluster, cluster, write_border_list)
!   implicit none
!   type(t_grid),                               intent(inout) :: grid
!   integer,                                    intent(in)    :: num_region
!   type(t_region), dimension(num_region),      intent(inout) :: region
!   integer,                                    intent(out)   :: num_cluster
!   type(t_cluster), dimension(:), allocatable, intent(inout) :: cluster
!   logical,                                    intent(in)    :: write_border_list
!
!   integer :: i, j, n
!   integer, dimension(num_region, num_region) :: num_point
!
!!   write(*,*) "create cluster... "
!   call get_neighbor_num_point(grid, num_region, num_point) 
!
!   ! count
!   do i=1, num_region
!      region(i)%num_neighbor  =0
!      do j=1, num_region
!         if(num_point(i,j)>=n0)then
!            region(i)%num_neighbor=region(i)%num_neighbor+1
!         endif
!      enddo
!   enddo
!
!   do i=1, num_region
!      if(allocated (region(i)%neighbor_list))then
!         deallocate(region(i)%neighbor_list)
!      endif
!      allocate(region(i)%neighbor_list(region(i)%num_neighbor))
!      n=0
!      do j=1, num_region
!         if(num_point(i,j)>=n0)then
!            n=n+1
!            region(i)%neighbor_list(n)=j
!         endif
!      enddo
!   enddo
!
!   ! set cluster index of region
!   region(:)%cluster =imiss
!
!   num_cluster=0
!   do i=1, num_region
!      if(region(i)%cluster==imiss)then
!         num_cluster=num_cluster+1
!         call include_region(num_region, region, i, num_cluster)
!      endif
!   enddo
!
!!   write(*,*) "num cluster before update", num_cluster
!   ! region cluster set complete set other structure in cluster, region and grid
!   call update_cluster(grid, num_region, region, num_cluster, cluster, write_border_list)
!!   write(*,*) "num cluster", num_cluster, sum(cluster(:)%num_region)
!   end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine include_region(num_region, region, i, set_idx)
   implicit none
   integer,                               intent(in)    :: i, set_idx, num_region
   type(t_region), dimension(num_region), intent(inout) :: region

   integer :: j, n 

   if(region(i)%cluster/=imiss) return
   !write(121,*) i, j, reg_idx
   region(i)%cluster=set_idx
   do j=1, region(i)%num_neighbor
      n=region(i)%neighbor_list(j)
      call include_region(num_region, region, n, set_idx)
   enddo
   end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine get_neighbor_num_point(grid, num_region, num_point)
!   implicit none
!   type(t_grid),                              intent(in)  :: grid
!   integer,                                   intent(in)  :: num_region
!   integer, dimension(num_region,num_region), intent(out) :: num_point
!
!   integer :: i, j, n, num_range, num_azim, n1, n2, i1, j1, i2, j2
!
!   !write(*,*) "in get_neighbor num point"
!   num_range=grid%num_range
!   num_azim =grid%num_azim
!
!   !write(*,*) "count border point..."
!   num_point=0 
!   do j=1, num_azim
!      if(grid%start_i(j)<1) cycle
!      do i=grid%start_i(j), grid%end_i(j)-1
!         i1=i+1
!         j1=j+1
!         if(j1>num_azim) j1=1
!         if(grid%region(i,j)/=imiss       .and.grid%region(i1,j)/=imiss       .and.&
!            grid%flag  (i,j)==flag_pt_used.and.grid%flag  (i1,j)/=flag_pt_used.and.&
!           (grid%region(i,j)/=grid%region(i1,j)))then
!            i2=i1
!            j2=j
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!!write(191,*) i1,j1,i2,j2,n1,n2
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!         endif
!         i1=i+1
!         j1=j+1
!         if(j1>num_azim) j1=1
!         if(grid%region(i,j)/=imiss       .and.grid%region(i,j1)/=imiss       .and.&
!            grid%flag  (i,j)==flag_pt_used.and.grid%flag  (i,j1)==flag_pt_used.and.&
!           (grid%region(i,j)/=grid%region(i,j1)))then
!            i2=i
!            j2=j1
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!!write(191,*) i1,j1,i2,j2,n1,n2
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!         endif
!      enddo
!   enddo
!
!   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_stdmax(grid, num_region, region, num_cluster, cluster)
   use grd, only: cluster_name, region_name, point_name
   implicit none
   type(t_grid),                            intent(inout) :: grid
   integer,                                 intent(in)    :: num_region, num_cluster
   type(t_region),  dimension(num_region),  intent(inout) :: region
   type(t_cluster), dimension(num_cluster), intent(inout) :: cluster

   real,    dimension(:,:), allocatable :: vel
   
   integer :: i, j, k, num_azim, num_range, nyq, n
   real    :: std, mean
   integer :: num 
   
   integer :: num_area
   type(t_grid) :: agrid
   type(t_region), dimension(:), allocatable :: area
   
   write(*,*) "Check maxstd", c5
   num_azim =grid%num_azim
   num_range=grid%num_range
   allocate(vel(num_range,num_azim))

   allocate(agrid%start_i (num_azim))
   allocate(agrid%  end_i (num_azim))
   allocate(agrid%azim    (num_azim))
   allocate(agrid%vel     (num_range,num_azim))
   allocate(agrid%ref     (num_range,num_azim))
   allocate(agrid%ref0    (num_range,num_azim))
   allocate(agrid%spw     (num_range,num_azim))
   allocate(agrid%area    (num_range,num_azim))
   allocate(agrid%region  (num_range,num_azim))
   allocate(agrid%cluster (num_range,num_azim))
   allocate(agrid%flag    (num_range,num_azim))
   allocate(agrid%level   (num_range,num_azim))

   agrid%num_range=grid%num_range
   agrid%num_azim =grid%num_azim
   agrid%azim     =grid%azim
   agrid%dr       =grid%dr
   agrid%da       =grid%da
   !agrid%vel      =grid%vel
   agrid%ref      =grid%ref
   agrid%spw      =grid%spw
   agrid%vmax     =grid%vmax
   agrid%start_i  =grid%start_i
   agrid%  end_i  =grid%  end_i
   agrid%area     =grid%area
   agrid%area_total=0

   !agrid%flag   =flag_pt_unknown
   
     vel=miss    
   ! mark low vel point
   !do j=1, num_azim
   !   if(grid%start_i(j)<1) cycle
   !   do i=grid%start_i(j), grid%end_i(j)
   !      n=grid%region(i,j)
   !      if(n==imiss) cycle
   !      nyq=region(n)%nyquist_number
   !      if(nyq/=miss)then
   !         if(IAND(grid%flag(i,j),flag_pt_low_snr)>0.and.abs(grid%vel(i,j))<=v0.and.(grid%dr*i)<=r0)then
   !            write(*,*) "low_snr and low_vel and near radar",i,j,grid%flag(i,j), grid%vel(i,j), r0
   !            cycle
   !         endif
   !         vel(i,j)=grid%vel(i,j)+2*nyq*grid%vmax
   !         agrid%flag(i,j)=flag_pt_used
   !      endif
   !   enddo
   !enddo
   
   std =0
   num =0
   mean=0
   ! mark low vel point
   do j=1, num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         n=grid%region(i,j)
         if(n==imiss) cycle
         nyq=region(n)%nyquist_number
         if(nyq/=miss)then
            if(IAND(grid%flag(i,j),flag_pt_low_snr)>0.and.abs(grid%vel(i,j))<=v0.and.(grid%dr*i)<=r0)then
      !        write(*,*) "low_snr and low_vel and near radar",i,j,grid%flag(i,j), grid%vel(i,j), r0
               cycle 
            endif
            vel(i,j)=grid%vel(i,j)+2*nyq*grid%vmax
            agrid%flag(i,j)=flag_pt_used
            num =num +1
            mean=mean+vel(i,j)
            std =std +vel(i,j)*vel(i,j)
         endif
      enddo
   enddo

   if(num>1)then
      mean=mean/num
      std=sqrt(std/num-mean*mean)

! only for v> vmax
      std=max(std,grid%vmax/3)
      write(check_std_unit,"(A,2F12.3)") "======== Check maxstd, mean,std=", mean, std
      write(check_std_unit,"(A)") "Cluster  Region      Point   unfold     vel     spw     ref flag      z0    area"
      do  j=1, num_azim
         if(grid%start_i(j)<1) cycle
         do i=grid%start_i(j), grid%end_i(j)
            if(vel(i,j)/=miss)then
               if(abs(vel(i,j)-mean)>3*std)then
                  n=grid%region(i,j)
                  if(n==imiss) cycle
                  if(region(n)%area_total<a0)then
                     if(region(n)%nyquist_number/=miss)then
                     write(nyquist_unit,*) "check std max remove region ", region_name(n)
                     write(check_std_unit,*) "check std max remove region ", region_name(n)
                     call remove_region(n, num_region, region, .true.)
                     endif
                  else
                  !else
                     write(check_std_unit,"(3(A,' '),4F8.2,I5,2F8.2)") cluster_name(grid%cluster(i,j)), region_name(n), point_name(i,j), &
                               vel(i,j), grid%vel(i,j), grid%spw(i,j), grid%ref(i,j), grid%flag(i,j), grid%ref0(i,j), region(n)%area_total
                     if(grid%flag(i,j)>flag_pt_used)then
                     call delete_point_region(grid,i,j,IEOR(grid%flag(i,j),flag_pt_used))
                     endif
                  endif
               endif
            endif
         enddo
      enddo

   !   vel=miss
   !   agrid%flag   =flag_pt_unknown
   !   do j=1, num_azim
   !      if(grid%start_i(j)<1) cycle
   !      do i=grid%start_i(j), grid%end_i(j)
   !         n=grid%region(i,j)
   !         if(n==imiss) cycle
   !         nyq=region(n)%nyquist_number
   !         if(nyq/=miss)then
   !            vel(i,j)=grid%vel(i,j)+2*nyq*grid%vmax
   !            agrid%flag(i,j)=flag_pt_used
   !         endif
   !      enddo
   !   enddo
   endif
   !agrid%vel      =vel

   !call create_region(agrid,num_area,area)

   !if(num_area<1) return
   !do n=1, num_area
   !   if(area(n)%area_total<a0)then
   !      do k=1,area(n)%num_point
   !        i=area(n)%point_list(k)%i
   !        j=area(n)%point_list(k)%j
   !        grid%region(i,j)=imiss
   !      enddo
   !   endif
   !enddo
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_grid(grid, num_region, region)
   implicit none
   type(t_grid),                          intent(inout) :: grid
   integer,                               intent(in)    :: num_region
   type(t_region), dimension(num_region), intent(inout) :: region

   real,    dimension(:,:), allocatable :: vel
   integer, dimension(:,:), allocatable :: flag

   integer, dimension(num_region, num_region) :: np
   integer, dimension(num_region) :: nyq_region

   integer :: i, j, num_azim, num_range, nyq, n, n1, n2
   integer :: i1, j1, i2, j2, num_in, num_out
   real    :: area_in, area_out
   integer :: kep_idx, rem_idx

   character(len=80) :: desc

   write(*,*) "Check grid: c5=", c5


   num_azim =grid%num_azim
   num_range=grid%num_range
   allocate(vel (num_range,num_azim))
   allocate(flag(num_range,num_azim))

     vel=miss 

    num_in=0
   area_in=0
   ! mark low vel point
   !flag=1 ! 0- normal region, 1-low vel region
   !do j=1, num_azim
   !   if(grid%start_i(j)<1) cycle
   !   do i=grid%start_i(j), grid%end_i(j)
   !      n=grid%region(i,j)
   !      if(n==imiss) cycle
   !      if(region(n)%nyquist_number/=miss)then
   !         flag(n)=0
   !      endif
   !      if(grid%flag(i,j)>0.and.abs(grid%vel(i,j))>v0)then
   !         flag(n)=0
   !      endif
   !   enddo
   !enddo

   !write(check_grid_unit,"(A,F8.3,A)") "======== check grid low vel point: v0=", v0, " ======="

   !nyq_region=miss
   !do j=1, num_azim
   !   !write(*,*) "check azim:",j
   !   if(grid%start_i(j)<1) cycle
   !   do i=grid%start_i(j), grid%end_i(j)
   !      n1=grid%region(i,j)
   !      if(n1==imiss) cycle
   !      if(flag(n1)/=1) cycle
   !      i1=i+1
   !      if(i1<=grid%end_i(j))then
   !         n2=grid%region(i1,j)
   !         if(n1/=n2.and.n2/=imiss)then
   !            if(region(n2)%nyquist_number/=miss.and.flag(n2)/=1.and.abs(grid%vel(i,j)-grid%vel(i1,j))<c6.and.grid%vel(i1,j)/=miss)then
   !               write(check_grid_unit,"(A,A,1X,A,F8.2,1X,I5,1X,A,1X,A,1X,A,F8.2,1X,I5,I5)") &
   !                    "low vel",       region_name(n1),point_name(i ,j),grid%vel(i ,j),flag(n1),&
   !                    cluster_name(n2),region_name(n2),point_name(i1,j),grid%vel(i1,j),flag(n2),region(n2)%nyquist_number
   !               if(nyq_region(n1)==miss.or.nyq_region(n1)==region(n2)%nyquist_number)then
   !                  nyq_region(n1)=region(n2)%nyquist_number
   !                  region(n1)%cluster=region(n2)%cluster
   !               else
   !                  nyq_region(n1)=miss-1
   !                  region(n1)%cluster=imiss
   !               endif
   !            endif
   !         endif
   !      endif
   !      j1=j+1
   !      if(j1>num_azim) j1=1
   !         n2=grid%region(i,j1)
   !         if(n1/=n2.and.n2/=imiss)then
   !            if(region(n2)%nyquist_number/=miss.and.flag(n2)/=1.and.abs(grid%vel(i,j)-grid%vel(i,j1))<c6.and.grid%vel(i,j1)/=miss)then
   !               write(check_grid_unit,"(A,A,1X,A,F8.2,1X,I5,1X,A,1X,A,1X,A,F8.2,1X,I5,I5)") &
   !                    "low vel",       region_name(n1),point_name(i ,j),grid%vel(i ,j),flag(n1),&
   !                    cluster_name(n2),region_name(n2),point_name(i,j1),grid%vel(i,j1),flag(n2),region(n2)%nyquist_number
   !               if(nyq_region(n1)==miss.or.nyq_region(n1)==region(n2)%nyquist_number)then
   !                  nyq_region(n1)=region(n2)%nyquist_number
   !                  region(n1)%cluster=region(n2)%cluster
   !               else
   !                  nyq_region(n1)=miss-1
   !                  region(n1)%cluster=imiss
   !               endif
   !            endif
   !         endif
   !      i2=i-1
   !      if(i2>=1)then
   !         n2=grid%region(i2,j)
   !         if(n1/=n2.and.n2/=imiss)then
   !            if(region(n2)%nyquist_number/=miss.and.flag(n2)/=1.and.abs(grid%vel(i,j)-grid%vel(i2,j))<c6.and.grid%vel(i2,j)/=miss)then
   !               write(check_grid_unit,"(A,A,1X,A,F8.2,1X,I5,1X,A,1X,A,1X,A,F8.2,1X,I5,I5)") &
   !                    "low vel",       region_name(n1),point_name(i ,j),grid%vel(i ,j),flag(n1),&
   !                    cluster_name(n2),region_name(n2),point_name(i2,j),grid%vel(i2,j),flag(n2),region(n2)%nyquist_number
   !               if(nyq_region(n1)==miss.or.nyq_region(n1)==region(n2)%nyquist_number)then
   !                  nyq_region(n1)=region(n2)%nyquist_number
   !                  region(n1)%cluster=region(n2)%cluster
   !               else
   !                  nyq_region(n1)=miss-1
   !                  region(n1)%cluster=imiss
   !               endif
   !            endif
   !         endif
   !      endif
   !      j2=j-1
   !      if(j2<1) j2=num_azim
   !         n2=grid%region(i,j2)
   !         if(n1/=n2.and.n2/=imiss)then
   !            if(region(n2)%nyquist_number/=miss.and.flag(n2)/=1.and.abs(grid%vel(i,j)-grid%vel(i,j2))<c6.and.grid%vel(i,j2)/=miss)then
   !               write(check_grid_unit,"(A,A,1X,A,F8.2,1X,I5,1X,A,1X,A,1X,A,F8.2,1X,I5,I5)") &
   !                    "low vel",       region_name(n1),point_name(i ,j),grid%vel(i ,j),flag(n1),&
   !                    cluster_name(n2),region_name(n2),point_name(i,j2),grid%vel(i,j2),flag(n2),region(n2)%nyquist_number
   !               if(nyq_region(n1)==miss.or.nyq_region(n1)==region(n2)%nyquist_number)then
   !                  nyq_region(n1)=region(n2)%nyquist_number
   !                  region(n1)%cluster=region(n2)%cluster
   !               else
   !                  nyq_region(n1)=miss-1
   !                  region(n1)%cluster=imiss
   !               endif
   !            endif
   !         endif
   !   enddo
   !enddo
   !do n=1, num_region
   !   if(flag(n)==1.and.nyq_region(n)>miss)then
   !      region(n)%nyquist_number=nyq_region(n)
   !      write(nyquist_unit,"(3A,I5,A)") "check grid set low vel region ",region_name(n)," nyq=",region(n)%nyquist_number, cluster_name(region(n)%cluster)
   !   endif
   !enddo

   do j=1, num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         n=grid%region(i,j)
         if(n==imiss) cycle
         nyq=region(n)%nyquist_number
         if(nyq/=miss)then
            vel(i,j)=grid%vel(i,j)+2*nyq*grid%vmax
            !if(abs(grid%vel(i,j))<=v0)then
            !   flag(i,j)=1
            !endif
             num_in= num_in+1
            area_in=area_in+grid%area(i,j)
         endif
      enddo
   enddo

   !write(check_grid_unit,"(A,F8.3,A)") "======== check grid low vel point: v0=", v0, " ======="

   !do j=1, num_azim
   !   !write(*,*) "check azim:",j
   !   if(grid%start_i(j)<1) cycle
   !   do i=grid%start_i(j), grid%end_i(j)
   !      i1=i+1
   !      if(i1<=grid%end_i(j))then
   !         if(flag(i,j)==1.and.vel(i1,j)/=miss.and.abs(vel(i,j)-vel(i1,j))>c6)then
   !            grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
   !            write(check_grid_unit,"(2(A,F8.2,1X),I5)") point_name(i,j),vel(i,j),point_name(i1,j),vel(i1,j),grid%flag(i,j)
   !         endif
   !      endif
   !      j1=j+1
   !      if(j1>num_azim) j1=1
   !         if(flag(i,j)==1.and.vel(i,j1)/=miss.and.abs(vel(i,j)-vel(i,j1))>c6)then
   !            grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
   !            write(check_grid_unit,"(2(A,F8.2,1X),I5)") point_name(i,j),vel(i,j),point_name(i,j1),vel(i,j1),grid%flag(i,j)
   !         endif
   !      i2=i-1
   !      if(i2>=1)then
   !         if(flag(i,j)==1.and.vel(i2,j)/=miss.and.abs(vel(i,j)-vel(i2,j))>c6)then
   !            grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
   !            write(check_grid_unit,"(2(A,F8.2,1X),I5)") point_name(i,j),vel(i,j),point_name(i2,j),vel(i2,j),grid%flag(i,j)
   !         endif
   !      endif
   !      j2=j-1
   !      if(j2<1) j2=num_azim
   !         if(flag(i,j)==1.and.vel(i,j2)/=miss.and.abs(vel(i,j)-vel(i,j2))>c6)then
   !            grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
   !            write(check_grid_unit,"(2(A,F8.2,1X),I5)") point_name(i,j),vel(i,j),point_name(i,j2),vel(i,j2),grid%flag(i,j)
   !         endif
   !   enddo
   !enddo
  ! ! delete low vel point
  ! do j=1, num_azim
  !    if(grid%start_i(j)<1) cycle
  !    do i=grid%start_i(j), grid%end_i(j)
  !       if(grid%region(i,j)<1)              cycle 
  !       if(.NOT.noise(i,j))                 cycle 
  !       if(flag(i,j)/=flag_pt_low_vel) cycle 
  !       call delete_point_region(grid,i,j,noise)
  !    enddo
  ! enddo

   !! mark noise point
   !do j=1, num_azim
   !   if(grid%start_i(j)<1) cycle
   !   do i=grid%start_i(j), grid%end_i(j)
   !      n=grid%region(i,j)
   !      if(n==imiss) cycle
   !      if(vel(i,j)/=miss)then
   !         if(grid%flag(i,j)>flag_pt_used)then
   !            noise(i,j)=.true.
   !         endif
   !      endif
   !   enddo
   !enddo

   ! check grid
   write(nyquist_unit,"(A,F8.3,A)") "======== check grid : c5=", c5, " ========"
   write(nyquist_unit,"(A)") " Region#1 level Region#2 level    Point#1    Point#2    v1    v2   f1   f2 Delete "
   desc=""
   flag=0
   do j=1, num_azim
      !write(*,*) "check azim:",j
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         j1=j+1
         if(j1>num_azim) j1=1
         call check_point(grid, num_region, region, vel, flag, i, j, i, j1)

         i1=i+1
         if(i1>grid%end_i(j)) cycle
         call check_point(grid, num_region, region, vel, flag, i, j, i1, j)
      enddo
   enddo

   do j=1, num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         if(flag(i,j)==1)then
            call delete_point_region(grid,i,j,IEOR(grid%flag(i,j),flag_pt_used))
         endif
      enddo
   enddo
   
    num_out=0
   area_out=0
   do j=1, grid%num_azim
      if(grid%start_i(j)<1)cycle
      do i=grid%start_i(j), grid%end_i(j)
         n=grid%region(i,j)
         if(n<1) cycle
         !if(region(n)%if_removed)cycle
         if(region(n)%nyquist_number==miss)cycle
          num_out= num_out+1
         area_out=area_out+grid%area(i,j)
      enddo
   enddo
   write(*,*) "After check grid remains ", num_out, "/", num_in , "points"
   write(*,*) "                         ", area_out, "/", area_in , "area"

   deallocate(vel)
   deallocate(flag)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_point(grid, num_region, region, vel, flag, i1, j1, i2, j2)
   use grd, only: region_name, point_name
   implicit none
   type(t_grid),                          intent(inout) :: grid
   integer     ,                          intent(in)    :: num_region 
   type(t_region), dimension(num_region), intent(inout) :: region
   real,           dimension(:,:),        intent(in)    :: vel
   integer,        dimension(:,:),        intent(inout) :: flag
   integer     ,                          intent(in)    :: i1, j1, i2, j2
  
   integer :: n1, n2, irem, jrem, ikep, jkep, nrem, nkep
   character(len=80) :: desc

   if(vel(i1,j1)==miss.or.vel(i2,j2)==miss) return
   if(IAND(grid%flag(i1,j1),flag_pt_used)==0.or.IAND(grid%flag(i2,j2),flag_pt_used)==0) return

   n1=grid%region(i1,j1)
   n2=grid%region(i2,j2)

   if(n1<1.or.n2<1) return

   desc=""
   if((abs(vel(i1,j1)-vel(i2,j2))>c5))then 
      if(grid%flag(i1,j1)==flag_pt_used.and.grid%flag(i2,j2)==flag_pt_used)then
                  
         if(region(n1)%level<region(n2)%level)then
            desc="region level"
            ikep=i1
            jkep=j1
            irem=i2
            jrem=j2
            nkep=n1
            nrem=n2
         elseif(region(n2)%level<region(n1)%level)then
            desc="region level"
            ikep=i2
            jkep=j2
            irem=i1
            jrem=j1
            nkep=n2
            nrem=n1
         else
            desc="region area"
            if(n1>n2)then
               ikep=i2
               jkep=j2
               irem=i1
               jrem=j1
               nkep=n2
               nrem=n1
            else
               irem=i2
               jrem=j2
               ikep=i1
               jkep=j1
               nkep=n1
               nrem=n2
            endif
         endif
         write(nyquist_unit,"(2(A9,I6),2(1X,A10),2F6.1,2I5,A9,1X,A)") &
               region_name(nrem),region(nrem)%level,&
               region_name(nkep),region(nkep)%level,&
               point_name (irem,jrem), point_name(ikep,jkep), &
               vel        (irem,jrem),        vel(ikep,jkep), &
               grid%flag  (irem,jrem),  grid%flag(ikep,jkep), region_name(nrem), trim(desc)
         if(region(nrem)%area_total>a1) write(region_unit,"(A)") "==== check grid remove region ==== ", region_name(nrem)
         call remove_region(nrem, num_region,region, .true.)
      elseif(grid%flag(i1,j1)==flag_pt_used.or.grid%flag(i2,j2)==flag_pt_used)then
         desc="noise point"
         if(grid%flag(i1,j1)==flag_pt_used)then
            irem=i2
            jrem=j2
            nrem=n2
            ikep=i1
            jkep=j1
            nkep=n1
         else
            irem=i1
            jrem=j1
            nrem=n1
            ikep=i2
            jkep=j2
            nkep=n2
         endif
         !call delete_point_region(grid,irem,jrem,IEOR(grid%flag(irem,jrem),flag_pt_used))
         flag(irem,jrem)=1
         write(nyquist_unit,"(2(A9,I6),2(1X,A10),2F6.1,2I5,1X,A10,1X,A)") &
               region_name(nrem),region(nrem)%level,&
               region_name(nkep),region(nkep)%level,&
               point_name (irem,jrem), point_name(ikep,jkep), &
               vel        (irem,jrem),        vel(ikep,jkep), &
               grid%flag  (irem,jrem),  grid%flag(ikep,jkep), point_name (irem,jrem), trim(desc)
      else
         desc="noise point"
         !call delete_point_region(grid,i1,j1,IEOR(grid%flag(i1,j1),flag_pt_used))
         !call delete_point_region(grid,i2,j2,IEOR(grid%flag(i2,j2),flag_pt_used))
         flag(i1,j1)=1
         flag(i2,j2)=1
         write(nyquist_unit,"(2(A9,I6),2(1X,A10),2F6.1,2I5,1X,2(A10,1X),A)") &
               region_name(n1),region(n1)%level,&
               region_name(n2),region(n2)%level,&
               point_name (i1,j1), point_name(i2,j2), &
               vel        (i1,j1),        vel(i2,j2), &
               grid%flag  (i1,j1),  grid%flag(i2,j2), &
               point_name (i1,j1), point_name(i2,j2), trim(desc)
      endif
   endif
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine delete_point_region(grid, i, j, flag)
   use grd, only: region_name, point_name
   implicit none
   type(t_grid),            intent(inout) :: grid
   integer     ,            intent(in)    :: i, j
   integer     , optional,  intent(in)    :: flag

   integer :: i1, j1, i2, j2, r1
   
   if(grid%region(i,j)==miss) return
   if(grid%flag  (i,j)<=1) return
   r1=grid%region(i,j)
   grid%region(i,j)=imiss
   i1=i-1
   j1=j-1
   i2=i+1
   j2=j+1
   if(j2>grid%num_azim) j2=1
   if(j1<1            ) j1=grid%num_azim
   if(present(flag))then
       write(check_grid_unit,"(3(A,1X),I5,F7.1,I5)") "remove point:", region_name(r1), point_name(i,j), grid%flag(i,j), grid%vel(i,j), flag
       !if(i1>=1             .and.IAND(grid%flag(i1,j),flag)>0.and.grid%region(i1,j)/=imiss) call delete_point_region(grid, i1, j,flag)
       !if(i2<=grid%num_range.and.IAND(grid%flag(i2,j),flag)>0.and.grid%region(i2,j)/=imiss) call delete_point_region(grid, i2, j,flag)
       !if(j1>=1             .and.IAND(grid%flag(i,j1),flag)>0.and.grid%region(i,j1)/=imiss) call delete_point_region(grid, i, j1,flag)
       !if(j2<=grid%num_azim .and.IAND(grid%flag(i,j2),flag)>0.and.grid%region(i,j2)/=imiss) call delete_point_region(grid, i, j2,flag)
   else
       write(check_grid_unit,"(3(A,1X),I5,F7.1)") "remove point:", region_name(r1), point_name(i,j), grid%flag(i,j), grid%vel(i,j)
       !if(i1>=1             .and.grid%flag(i1,j)>1.and.grid%region(i1,j)/=imiss) call delete_point_region(grid, i1, j)
       !if(i2<=grid%num_range.and.grid%flag(i2,j)>1.and.grid%region(i2,j)/=imiss) call delete_point_region(grid, i2, j)
       !if(j1>=1             .and.grid%flag(i,j1)>1.and.grid%region(i,j1)/=imiss) call delete_point_region(grid, i, j1)
       !if(j2<=grid%num_azim .and.grid%flag(i,j2)>1.and.grid%region(i,j2)/=imiss) call delete_point_region(grid, i, j2)
   endif
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   logical function check_region_continuity(grid, reg_idx)
!   implicit none
!   type(t_grid),                          intent(in)    :: grid
!   integer,                               intent(in)    :: reg_idx
!
!   integer :: i, j, num_azim, num_range, nyq, n, n1, n2
!   integer :: i1, j1, num_in, num_out
!
!   write(121,*) "Check region continuity:", reg_idx, c1
!   check_region_continuity=.true.
!   num_azim =grid%num_azim
!   num_range=grid%num_range
!   do j=1, num_azim
!      if(grid%start_i(j)<1) cycle
!      do i=grid%start_i(j), grid%end_i(j)
!         i1=i+1
!         if(i1>grid%end_i(j)) cycle
!         j1=j+1
!         if(j1>num_azim) j1=1
!         if(is_valid_vel(grid%vel(i,j)).and.is_valid_vel(grid%vel(i,j1)) .and. &
!            grid%region(i,j)==reg_idx.and.grid%region(i,j1)==reg_idx)then
!            if(abs(grid%vel(i,j)-grid%vel(i,j1))>c1)then
!               check_region_continuity=.false.
!               write(121,*) "check region continuinty failed,", reg_idx, i,j,i,j1,grid%vel(i,j),grid%vel(i,j1)
!               return
!            endif
!         endif 
!         if(is_valid_vel(grid%vel(i,j)).and.is_valid_vel(grid%vel(i1,j)) .and. &
!            grid%region(i,j)==reg_idx.and.grid%region(i1,j)==reg_idx)then
!            if(abs(grid%vel(i,j)-grid%vel(i1,j))>c1)then
!               check_region_continuity=.false.
!               write(121,*) "check region continuinty failed,", reg_idx, i,j,i1,j,grid%vel(i,j),grid%vel(i1,j)
!               return
!            endif
!         endif 
!      enddo
!   enddo
!
!   end function
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    subroutine update_cluster(grid, num_region, region, num_cluster, cluster, write_border_list)
!    implicit none
!    type(t_grid),                               intent(inout) :: grid
!    integer,                                    intent(in)    :: num_region
!    type(t_region), dimension(num_region),      intent(inout) :: region
!    integer,                                    intent(in)    :: num_cluster
!    type(t_cluster), dimension(:), allocatable, intent(out)   :: cluster
!    logical,                                    intent(in)    :: write_border_list
! 
!    integer :: i, j, k, n, np, s, e, w
!    real,    dimension(num_cluster) :: a
!    integer, dimension(num_cluster) :: d
! 
!    !write(*,*) "update cluster... "
! 
!    if(allocated(cluster))then
!       deallocate(cluster)
!    endif
!    allocate(cluster(num_cluster))
! 
!    ! region(:)%cluster already defined
!    ! sort cluster by area
!    a=0
!    do i=1, num_region
!       k=region(i)%cluster
!       a(k)=a(k)+region(i)%area_total
!       !a(k)=a(k)+region(i)%num_point
!    enddo
! 
!    do k=1, num_cluster
!       d(k)=k
!       cluster(k)%num_region =0
!       cluster(k)%num_point  =0
!       cluster(k)%area_total =0
!    enddo
! 
!    !a=-a
!    s=1
!    e=num_cluster
!    
!    call quick_sort_real(a,d,num_cluster,s,e)
!    !a=-a
!    
!    !change region(:)%cluster to sorted 
!    do i=1, num_region
!       if(region(i)%cluster==imiss) region(i)%cluster=0
!    enddo
!    do i=1, num_region
!       k=region(i)%cluster
!       do n=1, num_cluster
!          if(k==d(n))then
!             region(i)%cluster=-n
!             exit
!          endif
!       enddo
!    enddo
! 
!    ! set cluster
!    do i=1, num_region
!       if(region(i)%cluster/=0)then
!          region(i)%cluster=-region(i)%cluster
!          k=region(i)%cluster
!          cluster(k)%num_region=cluster(k)%num_region+1
!          cluster(k)%num_point =cluster(k)%num_point +region(i)%num_point
!          cluster(k)%area_total=cluster(k)%area_total+region(i)%area_total
!       endif
!    enddo
!    do i=1, num_region
!       if(region(i)%cluster==0) region(i)%cluster=imiss
!    enddo
! 
!    write(cluster_unit,"(A))") "================ Cluster List ================"
!    write(cluster_unit,"(A))") "cluster, num region,  num point,  area total"
!    do i=1, num_cluster
!       !write(*,*) "cluster update",  i
!       write(cluster_unit,"(A7,2(',',I11),',',F12.3)") &
!             cluster_name(i), cluster(i)%num_region,cluster(i)%num_point,cluster(i)%area_total
!    enddo
! 
!    ! set grid%cluster
!    do n=1, num_region 
!       do k=1, region(n)%num_point
!          i=region(n)%point_list(k)%i
!          j=region(n)%point_list(k)%j
!          grid%cluster(i,j)=region(n)%cluster
!       enddo
!    enddo
!    
!    write(sort_unit,"(A)") "================ Cluster Sort ================" 
!    do k=1, num_cluster
!       !write(*,*) "cluster sort",  k
!       write(sort_unit,"(A,2(I10,F12.3))") cluster_name(k), d(k), a(k), cluster(k)%num_point, cluster(k)%area_total
!    enddo
! 
!    write(cluster_unit,"(A)") "================ Cluster Region List ================"
!    write(cluster_unit,"(A)") "cluster, region#, region name, num point,  area total"
!    d=0
!    do k=1, num_cluster
!       if(allocated (cluster(k)%region_list))then
!          deallocate(cluster(k)%region_list)
!       endif
!       allocate( cluster(k)%region_list(cluster(k)%num_region))
!       write(cluster_unit,"(A)") "-----------------------------------------------------"
!       do n=1,num_region
!          if(region(n)%cluster==k)then
!             d(k)=d(k)+1
!             cluster(k)%region_list(d(k))=n
!          endif
!       enddo
!       !write(111,*) "cluster region list:", k, cluster(k)%num_region, d(k),ubound(cluster(k)%region_list), cluster(k)%num_point, cluster(k)%area_total
!       do i=1, cluster(k)%num_region
!          !write(111,*) i,cluster(k)%region_list(i)
!          !write(111,*) region(cluster(k)%region_list(i))%num_point
!       !write(*,*) "cluster update2",  k
!          write(cluster_unit,"(A7,',',I8,',',A12,',',I10,',',F12.3)") cluster_name(k), i, region_name(cluster(k)%region_list(i)), &
!                                   region(cluster(k)%region_list(i))%num_point, region(cluster(k)%region_list(i))%area_total
!       enddo
!       !write(111,*) (cluster(k)%region_list(i), region(cluster(k)%region_list(i))%num_point,',',i=1,cluster(k)%num_region)
!    enddo 
!    
!    if(write_border_list)then
!       call get_neighbor_flag(grid, num_region, region, num_cluster, cluster, .true.)
!    else
!       call get_neighbor_flag(grid, num_region, region, num_cluster, cluster, .false.)
!    endif
!    end subroutine
! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_region_miss(miss_in)
   implicit none
   real, intent(in) :: miss_in
  
   miss=miss_in
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_region_miss(miss_out)
   implicit none
   real, intent(out) :: miss_out
  
   miss_out=miss
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine dealias_region(num_range,num_azim,vel,bkg_in,ref,spw,vmax,azim,elv,dr,miss_var,nyq_num,flag,vel_out,&
                             ci0,ci1,ci2,ci3,ci4,ci5,vi0,wi0,wi1,zi0,ri0,ri1,ai0,ai1,ait,ni0,ival,rval)
   use grd, only: write_grd, write_grid_csv, get_filename, write_test_data, cluster_name, region_name, point_name
   implicit none

   integer,                                intent(in)  :: num_range, num_azim
   real,    dimension(num_range,num_azim), intent(in)  :: vel, bkg_in, ref, spw ! radical velocity, reference value, reflectivity, spectral width
   real,                                   intent(in)  :: vmax
   real,    dimension(num_azim),           intent(in)  :: azim
   real                        ,           intent(in)  :: dr, elv
   real,                                   intent(in)  :: miss_var
   real,                                   intent(in)  :: ci0,ci1,ci2,ci3,ci4,ci5,vi0,wi0,wi1,zi0,ri0,ri1,ai0,ai1,ait
   integer,                                intent(in)  :: ni0
   integer, dimension(num_range,num_azim), intent(out) :: nyq_num, flag
   real,    dimension(num_range,num_azim), intent(out) :: vel_out 
   integer, dimension(20),                 intent(out) :: ival
   ! 1: num point total valid
   ! 2: num region created
   ! 3: num region after merge
   ! 4: num cluster
   ! 5: num point in largest region
   ! 6: num point in largest cluster
   ! 7: num region in largest cluster
   ! 8: num point dealiased
   ! 9: num region dealiased
   !10: num cluster dealiased
   !11: num point aliased
   !12: num region aliased
   !13: minimum nyquist number
   !14: maximum nyquist number
   real,    dimension(40),                 intent(out) :: rval
   ! 1: minimum value
   ! 2: maximum value
   ! 3: Mean OMB in largest cluster
   ! 4: contour interval 
   ! 5: percent of point dealiased
   ! 6: percent of point aliased
   ! 7: c0
   ! 8: c1
   ! 9: c2
   ! 10: c3
   ! 11: c4
   ! 12: c5

   integer, dimension(num_range,num_azim) :: point_cluster_index
   real,    dimension(num_range,num_azim) :: bkg

   integer :: i, j, k, l, n, m, np, nsum, nyq, n1, n2, i1, j1, l1, npoint
   integer :: num_left, old_left, k1, k2

   type(t_region),  dimension(:), allocatable :: region
   type(t_cluster), dimension(:), allocatable :: cluster
   type(t_grid)                               :: grid

   logical :: new_cluster

   real   , dimension(:), allocatable :: area_dealiased_in_region
   integer, dimension(:), allocatable :: num_point_dealiased_in_cluster
   integer :: num_region_dealiased, num_point_dealiased, num_point_std_check
   integer :: num_cluster_dealiased, num_region_aliased, num_point_aliased
   integer :: minimum_nyquist_number, maximum_nyquist_number
   real    :: minimum_vel, maximum_vel, area_total_dealiased, area_total_aliased, area_total, area_zero
   real    :: max_vel, min_vel, area_no_nyq
   !real, dimension(:), allocatable :: mean_omb_cluster, std_omb_cluster

   logical :: first_in_cluster
   real(kind=8) :: mean_vel, std_vel, mean_o, mean_b, std_o, std_b, cor_ob

   character(len=800) :: filename

   write(*,*) "Dealiasing velocity (region based dealiasing method)...", elv

   call get_filename("cluster",filename,"txt")
   open(cluster_unit   ,file=filename,status="unknown")
   call get_filename("region" ,filename,"txt")
   open(region_unit    ,file=filename,status="unknown")
   call get_filename("sort"   ,filename,"txt")
   open(sort_unit      ,file=filename,status="unknown")
   call get_filename("nyquist",filename,"txt")
   open(nyquist_unit   ,file=filename,status="unknown")
   call get_filename("ref_chk",filename,"txt")
   open(reference_unit ,file=filename,status="unknown")
   call get_filename("grd_chk",filename,"txt")
   open(check_grid_unit,file=filename,status="unknown")
   call get_filename("border" ,filename,"txt")
   open(border_unit    ,file=filename,status="unknown")
   call get_filename("stdmax",filename,"txt")
   open(check_std_unit ,file=filename,status="unknown")

   call set_region_miss(miss_var)
   ival=miss
   rval=miss

   !min_range=min_height/sin(elv*DEG2RAD)
   !if(elv<6.5)then
   !  min_range=7500.
   !else
   !  min_range=750.
   !endif

   !azim=0
   call write_grd('obs.grd',num_range,num_azim,vel,miss)

   v0=vel0
   w0=spw0
   !w1=spw1
   w1=vmax/3
   z0=ref0
   !r0=rng0
   !r0=9916.3*(elv**(-0.259))
   !if(ri0>=0) r0=ri0*(elv**(-0.259))
   !r0=-7819*log(elv)+26699
   !if(ri0>=0) r0=-7819*log(elv)+ri0
   if(ri0>=0.and.ri1/=miss)then
      r0=ri0*(elv**(ri1))
   elseif(ri0>=0)then
      r0=ri0*(elv**(-0.24))
   elseif(ri1/=miss)then
      r0=12964*(elv**(ri1))
   else
      r0=12964*(elv**(-0.24))
   endif

   !if(elv<2) r0=1000
   r1=rng1

   if(vi0>=0) v0=vi0
   if(wi0>=0) w0=wi0
!  if(wi1>=0) w1=wi1
   if(zi0/=miss) z0=zi0
   !if(ri1>=0) r1=ri1
   if(ai0>=0) a0=ai0
   if(ai1>=0) a1=ai1
   if(ni0>0 ) n0=ni0
   if(wi1>-vmax) w1=w1+wi1
   at=ait

   call init_grid(num_range,num_azim,azim,dr,vel,ref,spw,vmax,grid)
   if(ci0/=miss) c0=ci0
   if(ci1/=miss) c1=c1+ci1
   if(ci2/=miss) c2=c2+ci2
   if(ci3/=miss) c3=c3+ci3
   if(ci4/=miss) c4=c4+ci4
   if(ci5/=miss) c5=c1+ci5
   c6=c5/5

!  vel_out=miss
!  where(grid%flag==flag_pt_used)
!     vel_out=vel
!  endwhere
   !where(vel>v0)
   !   vel_out=vel
   !endwhere
!  call calculate_torus_rv(num_range, num_azim, azim, elv, vmax, dr, miss, vel_out, bkg)
   !call calculate_torus_rv(num_range, num_azim, azim, elv, vmax, dr, miss, vel, bkg)
   bkg=bkg_in
   call write_test_data('bkg.dat'     ,num_range,num_azim,azim,dr,bkg,miss)

!stop

   write(*,*) "create region..."
   call create_region(grid,num_region,region)
   call write_grd      ('region.grd'     ,num_range,num_azim,grid%region,imiss)
   call write_test_data('region.dat'     ,num_range,num_azim,azim,dr,grid%region,imiss)
   write(*,*) "num region after create:", num_region 

   ival( 2)=num_region
   if (grid%num_valid==0) then
      ival( 1)=0
      ival( 8)=0
      ival(11)=0
      return
   endif

   region(:)%base_region   =imiss
   region(:)%base_border   =0
   region(:)%level         =miss
   region(:)%nyquist_number=miss
   

   !do n=1, num_region
   !   if(region(n)%area_total<a0)then
   !      i=region(n)%point_list(1)%i
   !      j=region(n)%point_list(1)%j
   !      if(abs(grid%vel(i,j))<=v0)then
   !         region(n)%flag=flag_rg_zero
   !         write(region_unit,"(2A,I4)") region_name(n), " set region flag to flag_rg_zero:", flag_rg_zero
   !         do l=1, region(n)%num_point
   !            i=region(n)%point_list(l)%i
   !            j=region(n)%point_list(l)%j
   !            grid%flag(i,j)=IOR(grid%flag(i,j),flag_pt_low_vel)
   !         enddo
   !      endif
   !   endif
   !enddo
   !call create_cluster(grid, num_region, region, num_cluster, cluster, .false.)
   call create_cluster(grid, num_region, region, num_cluster, cluster)
!stop
!   ! merge and create cluster
   !write(*,*) "merge region..."
   !call merge_region (grid, num_region, region, num_cluster, cluster)
   !region(:)%base_region   =imiss
   !region(:)%nyquist_number= miss
   !call write_test_data('merge_region.dat',num_range,num_azim,azim,dr,grid%region,imiss)
   !call write_grd('merge_region.grd',num_range,num_azim,grid%region,imiss)
   !call write_grid_csv('grid_merge.csv', num_range,num_azim, grid%vel, grid%region, grid%cluster,miss)
   !call write_region_csv('region_merge.csv', num_region, region(:)%cluster, region(:)%num_neighbor, region(:)%base_region, region(:)%nyquist_number, miss)
   !write(*,*) "num region after merge:", num_region 

   !region(:)%if_removed=.false.

   write(*,*) "Number of cluster:", num_cluster
   write(*,"(A,30I6)") "Cluster num region(top30):",  cluster(1:min(num_cluster,30))%num_region
   write(*,"(A,30I6)") "Cluster  of region(top30):",  region (1:min(num_region ,30))%cluster
   write(*,"(A,30I6)") "Flag     of region(top30):",  region (1:min(num_region ,30))%flag
   write(*,"(A,30F6.0)") "Area     of region(top30):",  region (1:min(num_region ,30))%area_total

   write(*,*) "calculate region nyquist number..."
   ! set nyquist folding number of region
   allocate(num_point_dealiased_in_cluster(num_cluster   ))
   allocate(area_dealiased_in_region(num_region))

   num_point_dealiased_in_cluster=0
   cluster(:)%nyquist_number=miss
   
   call write_grid_csv ('obs.csv',num_range,num_azim,grid%dr,grid%azim,grid%vel, grid%ref, grid%ref0, grid%spw, grid%flag, grid%region, grid%cluster, miss)

   area_no_nyq=0

   do n=1, num_cluster
!cluster(n)%area_total<ccmin*grid%area_total.or.
      if(cluster(n)%area_total>=a1) then
         first_in_cluster=.true.
         num_point_dealiased_in_cluster(n)=0
         write(nyquist_unit,"(A)") "---- set region nyquist number in "//cluster_name(n)//" ----"
         write(nyquist_unit,"(A)") "   set       base      base      set        border    border base      set        neighbor"
         write(nyquist_unit,"(A)") "   region#   region#   nyquist   nyquist    npoint    effect region_pt region_pt      flag"

         cluster(n)%nyquist_number=0
         do j=1, cluster(n)%num_region
            i=cluster(n)%region_list(j)
            if(region(i)%cluster==n)then
                 if(first_in_cluster.and.region(i)%flag/=flag_rg_zero.and. region(i)%area_total>=a0)then
                    region(i)%nyquist_number=0
                    write(nyquist_unit,"(2A10,7I10)") region_name(i),region_name(0),0,0,0,0,0,region(i)%num_point,0
                    region(i)%base_region=0
                    region(i)%level      =0
                    region(i)%base_border=region(i)%num_point
                    call calc_neighbor_nyquist_number(num_region, region, cluster(n), i)
                    first_in_cluster=.false.
                 endif
               !endif
            endif
         enddo
         !if(region(1)%nyquist_number==miss)then
         !write(102,*) "cluster", n, cluster(n)%num_region, ',',cluster(n)%region_list
         !write(102,*) "C1,region 1 nyq", n, j, i, region(1)%nyquist_number, region(1)%num_point,cluster(1)%num_point
         !endif
         do j=1, cluster(n)%num_region
            i=cluster(n)%region_list(j)
            if(region(i)%cluster==n)then
               if(region(i)%nyquist_number==miss.or.abs(region(i)%nyquist_number)>10)then
                  if(region(i)%nyquist_number/=miss)then
                     write(*,*) "set miss: Error in ", region_name(i), " nyq=",region(i)%nyquist_number
                     region(i)%nyquist_number=miss
                  endif
                  area_no_nyq=area_no_nyq+region(i)%area_total
                  write(nyquist_unit,*) "Cant calulate region ", region_name(i)," nyquist number", region(i)%area_total, " area_no_nyq=",area_no_nyq
               else
                  num_point_dealiased_in_cluster(n)=num_point_dealiased_in_cluster(n)+region(i)%num_point
               !   write(nyquist_unit,"(A,A,I10,1X,A,2(A,I10))") cluster_name(n)," dealiased ",num_point_dealiased_in_cluster(n),&
               !    region_name(i), " npoint=", region(i)%num_point, " nyq=", region(i)%nyquist_number
               endif
            endif
         enddo
         !write(102,*) "C2,region 1 nyq", n, j, i, region(1)%nyquist_number, region(1)%num_point,cluster(1)%num_point

         !if(cluster(n)%area_total<ccmin*grid%area_total.or.cluster(n)%num_point<nrmin)then
            if(n<=30) write(*,*) "number of point in cluster(dealiased):", n, cluster(n)%num_point, "(",num_point_dealiased_in_cluster(n),")"
      else
      !if(cluster(n)%area_total<a1)then
         cluster(n)%nyquist_number=miss
         area_no_nyq=area_no_nyq+cluster(n)%area_total
         write(nyquist_unit,"(2A,2I10,F8.3,2F12.3,A,F12.3)") 'set miss: Area too small in ', &
                                   cluster_name(n),cluster(n)%num_point, grid%num_valid, &
                                   a1, cluster(n)%area_total, grid%area_total, " area_no_nyq=",area_no_nyq
         !write(*,*) 'num in cluster <100', n
         !do j=1, cluster(n)%num_region
         !   i=cluster(n)%region_list(j)
         !   if(region(i)%cluster==n)then
         !      write(nyquist_unit,"(2A,I6,2A,F8.3,3I10,F8.3,3F12.3)") 'set miss: Area too small in ', cluster_name(n),j,"#", region_name(i), &
         !            ccmin, region(i)%num_point, cluster(n)%num_point, grid%num_valid, a1, region(i)%area_total, cluster(n)%area_total, grid%area_total
         !      region(i)%nyquist_number=miss
         !   endif
         !enddo
         !exit
      endif
      !if(sum(num_point_dealiased_in_cluster(1:n))>cratio*reg%num_valid)then
      !   exit
      !endif
      !write(102,*) "C3,region 1 nyq", n, j, i, region(1)%nyquist_number, region(1)%num_point,cluster(1)%num_point
   enddo

   write(nyquist_unit,"(A)") "================ Region Nyquist Number Calculated ================"
   do i=1, num_region
      if(region(i)%nyquist_number/=miss.and.cluster(region(i)%cluster)%nyquist_number/=miss)then
         write(nyquist_unit,"(A,1X,2A,I10,A,F12.3,A,I5)") cluster_name(region(i)%cluster),region_name(i),&
              " npoint=",region(i)%num_point," area=", region(i)%area_total, " nyq=", region(i)%nyquist_number
      endif
   enddo
   write(nyquist_unit,"(A,F12.3)") "No nyquist number region area=: ",area_no_nyq


   !nyq_num=miss
   !num_point_dealiased =0
   !num_region_dealiased=0 
   !vel_out=miss
   !do n=1, reg%num_region
   !   if(nyq_region(n)/=miss)then
   !   !write(*,*) "Region ", n, " nyquist number=", nyq_region(n), num_point_in_region(n) !, mean_obs(n), mean_bak(n), mean_bkg(n)
   !      num_region_dealiased=num_region_dealiased+1
   !      num_point_dealiased =num_point_dealiased +num_point_in_region(n)
   !      do m=1,  reg%point_list(n)%num_point
   !         i=reg%point_list(n)%range_index(m)
   !         j=reg%point_list(n)%azim_index(m)
   !         vel_out(i,j)=vel(i,j)+2*nyq_region(n)*vmax
   !         nyq_num(i,j)=nyq_region(n)
   !      enddo
   !   endif
   !enddo
   !write(*,*) "Before Reference check, Dealiased ", num_region_dealiased," region, ", num_point_dealiased, &
   !           " ponits out of ",reg%num_valid,"(",100.*num_point_dealiased/reg%num_valid,"%)."
   !call write_test_data('before_bkg.dat',num_range,num_azim,azim,dr,vel_out,miss)

   !reference check
   !allocate(mean_omb_cluster(num_cluster))
   !allocate( std_omb_cluster(num_cluster))
   cluster(:)%mean_omb=miss
   cluster(:)% std_omb=miss
   num_cluster_dealiased=0 
   vel_out=miss
   do l=1, num_cluster
      mean_vel=0
      std_vel=0
      area_total  =0
      if(cluster(l)%nyquist_number==miss) cycle
      write(reference_unit,"(3A,I10,A,F12.3,A)") "====== ", cluster_name(l), &
                                          " num_region=", cluster(l)%num_region, &
                                                " area=", cluster(l)%area_total, " ======"
      do k=1, cluster(l)%num_region
         n=cluster(l)%region_list(k)
         write(reference_unit,"(3A,I5,A,I7,A,F12.3,A)") "------ ", region_name(n), &
                                             " nyq=", region(n)%nyquist_number, &
                                          " npoint=", region(n)%num_point, &
                                          "   area=", region(n)%area_total , " ------"
         if(region(n)%nyquist_number==miss) cycle
         do m=1, region(n)%num_point
            i=region(n)%point_list(m)%i
            j=region(n)%point_list(m)%j
            if(.NOT.is_valid_vel(bkg(i,j))) cycle
            if(abs(grid%vel(i,j))<=v0) cycle
            if(grid%region(i,j)<1) cycle
            if(grid%flag(i,j)<flag_pt_used) cycle
            write(reference_unit,"(2A,1X,A,I6,3A,I5,2F8.2)") "reference check:", cluster_name(region(n)%cluster), region_name(n), m,"# ", point_name(i, j), &
                  " nyq=", region(n)%nyquist_number, vel(i,j), bkg(i,j)
            vel_out(i,j)=vel(i,j)+(cluster(l)%nyquist_number+region(n)%nyquist_number)*2*vmax
            mean_vel =mean_vel+grid%area(i,j)*(vel_out(i,j)-bkg(i,j))
             std_vel = std_vel+grid%area(i,j)*(vel_out(i,j)-bkg(i,j))*(vel_out(i,j)-bkg(i,j))
            area_total=area_total+grid%area(i,j)
         enddo
      enddo
      !write(*,*) "reference check npoint=", npoint, ",cluster=", l
      if(area_total>a0)then
         mean_vel=mean_vel/area_total
         !write(*,*) "calc std:",  std_vel, area_total, mean_vel, std_vel/area_total-mean_vel*mean_vel
         std_vel=(std_vel/area_total)-(mean_vel*mean_vel)
         if(std_vel>0)then
            std_vel=sqrt(std_vel)
         else
            std_vel=0 !abs(mean_vel)
         endif
         cluster(l)%mean_omb=mean_vel
         cluster(l)% std_omb= std_vel

         write(*,"(A,2F12.3,2F6.1,1X,A,2F8.2)") cluster_name(l), area_total, cluster(l)%area_total, cluster(l)%minv, cluster(l)%maxv, "vel-bkg, mean, std: ", mean_vel, std_vel

         k=0
         if(abs(mean_vel)>c4)then
            k=-floor(mean_vel/(2*vmax)+0.5)
         endif

         if(abs(mean_vel)<c3)then
            num_cluster_dealiased=num_cluster_dealiased+1 
            cycle
         endif
         mean_vel=mean_vel+k*2*vmax
         
         if(abs(mean_vel)<c3.and.abs(k)<10)then
            write(*,*) "reference check: change ", cluster_name(l), " nyq=", k
            num_cluster_dealiased=num_cluster_dealiased+1 
            cluster(l)%nyquist_number=k
            !do i=1, cluster(l)%num_region
            !   n=cluster(l)%region_list(i)
            !   if(region(n)%cluster/=l) cycle
            !   if(region(n)%nyquist_number==miss) cycle
            !   if(k/=0) write(99,"(3A,2(A,I10),I3,A)") cluster_name(l), ":", region_name(n), &
            !            " npoint", region(n)%num_point, " area_total", region(n)%area_total, &
            !            " nyquist", region(n)%nyquist_number, k," changed by reference check"
            !   region(n)%nyquist_number=region(n)%nyquist_number+k
            !enddo
         else
            cluster(l)%nyquist_number=miss
         endif
      else
         write(*,*) "vel-bkg, ",cluster_name(l)," valid area=",area_total,"< a0=",a0
         cluster(l)%nyquist_number=miss
      endif

      if(cluster(l)%nyquist_number==miss)then
         write(*,*) "reference check set nyq=miss for ", cluster_name(l)
         num_point_dealiased_in_cluster(l)=0
         cluster(l)%nyquist_number=miss
         do i=1, cluster(l)%num_region
            n=cluster(l)%region_list(i)
            if(region(n)%area_total>a1) write(region_unit,"(3A)") region_name(n), " removed by reference check ", cluster_name(l)
         !   if(region(n)%cluster/=l) cycle
         !   region(n)%nyquist_number=miss
         enddo
      endif
   enddo
   
   write(nyquist_unit,"(A)") "========== Cluster Nyquist Number Changed After Reference Check =========="
   do i=1, num_region
      if(region(i)%nyquist_number/=miss.and.cluster(region(i)%cluster)%nyquist_number/=miss)then
         region(i)%nyquist_number=region(i)%nyquist_number+cluster(region(i)%cluster)%nyquist_number
         if(cluster(region(i)%cluster)%nyquist_number/=0)then
         write(nyquist_unit,"(A,1X,2A,I10,A,F12.3,A,I5)") cluster_name(region(i)%cluster),region_name(i),&
              " npoint=",region(i)%num_point," area=", region(i)%area_total, " nyq=", region(i)%nyquist_number
         endif
      else
         region(i)%nyquist_number=miss
         if(region(i)%nyquist_number/=miss.and.cluster(region(i)%cluster)%nyquist_number==miss)then
         write(nyquist_unit,"(3A,I10,A,F12.3,A,I5)") "set miss: ",cluster_name(region(i)%cluster)," nyq=miss, ",&
                                                                  region_name (i), " nyq=miss"
         endif
      endif
   enddo

   call check_grid(grid, num_region, region)
!stop

   write(nyquist_unit,"(A)") "========== Check low vel Region =========="
   write(check_grid_unit,"(A)") "==== Check cluster point omb > 3*std_omb ===="
   write(check_grid_unit,"(A)") " Region Cluster      Point    vel    bak   mean    std    omb  3*std"
! check max and min vel, delete low vel region, check omb < 5*std_omb
   num_point_std_check=0
   do n=1, num_region
      if(region(n)%nyquist_number/=miss)then
         min_vel=miss
         max_vel=miss
         mean_vel=cluster(region(n)%cluster)%mean_omb+cluster(region(n)%cluster)%nyquist_number*2*vmax
          std_vel=cluster(region(n)%cluster)% std_omb
         l=0
         do m=1,  region(n)%num_point
            i=region(n)%point_list(m)%i
            j=region(n)%point_list(m)%j
            if(grid%region(i,j)<1) cycle
            vel_out(i,j)=vel(i,j)+2*region(n)%nyquist_number*vmax
            !if(bkg(i,j)/=miss.and.abs(vel_out(i,j)-bkg(i,j)-mean_vel)>3*std_vel)then
            !   write(check_grid_unit,"(2(A7,1X),A10,6F7.1)") &
            !         region_name(n), cluster_name(region(n)%cluster), point_name(i,j),  &
            !         vel_out(i,j), bkg(i,j), mean_vel, std_vel, vel_out(i,j)-bkg(i,j)-mean_vel, 3*std_vel
            !   grid%region(i,j)=imiss
            !   num_point_std_check=num_point_std_check+1
            !   write(check_grid_unit,"(2A)") &
            !         "set miss: omb >3*std:",region_name(n)//" "//cluster_name(region(n)%cluster)//" "//point_name(i,j)
            !         
            !   cycle
            !endif
            
            l=l+1
            if(min_vel==miss)then
               min_vel=vel_out(i,j)
            else
               min_vel=min(min_vel,vel_out(i,j))
            endif
            if(max_vel==miss)then
               max_vel=vel_out(i,j)
            else
               max_vel=max(max_vel,vel_out(i,j))
            endif
         enddo
         
         if(l==0)then !((min_vel>=-v0.and.max_vel<=v0).or..and.region(n)%area_total<a0)then ! low vel region or no point region
            region(n)%nyquist_number=miss 
            write(nyquist_unit,"(3A,I4,A,I8,A,F10.1,2F10.2)") "set miss:low vel ",region_name(n)//" "//cluster_name(region(n)%cluster), &
                                                          " nyq=", region(n)%nyquist_number, " np=",region(n)%num_point, &
                                                          " area=",region(n)%area_total, min_vel, max_vel
         endif
      endif
   enddo

!  call check_stdmax(grid, num_region, region, num_cluster, cluster)

   ! ! cluster vel
   ! vel_out=miss
   ! do n=1, num_region
   !    if(region(n)%nyquist_number/=miss.and.region(n)%cluster==1)then
   !       do m=1,  region(n)%num_point
   !          i=region(n)%point_list(m)%i
   !          j=region(n)%point_list(m)%j
   !          if(grid%region(i,j)<1) cycle
   !          vel_out(i,j)=vel(i,j)+2*region(n)%nyquist_number*vmax
   !       enddo
   !    endif
   ! enddo

   ! call extrapolate_rv(num_range,num_azim,azim,dr,vel_out,miss,bkg)
   vel_out=miss
!   vel_out=miss
!   do l=1, num_cluster
!      do k=1, cluster(l)%num_region
!         n=cluster(l)%region_list(k)
!         if(region(n)%nyquist_number==miss) cycle
!         mean_o=0
!          std_o=0
!         mean_b=0
!          std_b=0
!         cor_ob=0
!         area_zero  =0
!         do m=1, region(n)%num_point
!            i=region(n)%point_list(m)%i
!            j=region(n)%point_list(m)%j
!            vel_out(i,j)=vel(i,j)+region(n)%nyquist_number*2*vmax
!            if(abs(vel_out(i,j))<s1)then
!               mean_o =mean_o+grid%area(i,j)*(vel_out(i,j))
!                std_o = std_o+grid%area(i,j)*(vel_out(i,j))*(vel_out(i,j))
!               mean_b =mean_b+grid%area(i,j)*(bkg(i,j))
!                std_b = std_b+grid%area(i,j)*(bkg(i,j))*(bkg(i,j))
!               cor_ob =cor_ob+grid%area(i,j)*(vel_out(i,j))*(bkg(i,j))
!               area_zero=area_zero+grid%area(i,j)
!            endif
!         enddo
!         if(area_zero>0)then
!            mean_o=mean_o/area_zero
!            mean_b=mean_b/area_zero
!             std_o=sqrt(std_o/area_zero-mean_o*mean_o)
!             std_b=sqrt(std_b/area_zero-mean_b*mean_b)
!            cor_ob=(cor_ob   /area_zero-mean_o*mean_b)/std_o/std_b
!            write(*,"(2A,1X,A,4F8.2,F10.5)") "obs(mean,std),bak(mean,std),cor:", cluster_name(l), region_name(n), mean_o, std_o, mean_b, std_b, cor_ob
!            if(std_o==0.or.cor_ob<s2.or.abs(std_o-std_b)>s1)then ! delete zero region
!               write(*,*) "delete zero point in ", cluster_name(l), region_name(n),&
!                          " cor_ob<s2", cor_ob, cor_ob<s2, "abs(std_o-std_b)>s1",abs(std_o-std_b), abs(std_o-std_b)>s1
!               do m=1, region(n)%num_point
!                  i=region(n)%point_list(m)%i
!                  j=region(n)%point_list(m)%j
!                  vel_out(i,j)=vel(i,j)+region(n)%nyquist_number*2*vmax
!                  if(abs(vel_out(i,j))<s1)then
!                     grid%region(i,j)=imiss
!                  endif
!               enddo 
!            endif
!         endif
!      enddo
!   enddo

!   call write_grd('bf_chkgrd.grd',num_range,num_azim,vel_out,miss)
!   write(*,*) "call check grid, vmax:", reg%vmax
!   reg%vmax=vmax

   !point_cluster_index=miss
   !do j=1, num_azim
   !   if(reg%start_i(j)<1) cycle
   !   do i=reg%start_i(j),reg%end_i(j)
   !      n=reg%region(i,j)
   !      if(n<1) cycle
   !      if(region_cluster_index(n)==1)then
   !        point_cluster_index(i,j)=1
   !      else
   !        point_cluster_index(i,j)=0
   !      endif
   !   enddo
   !enddo
   call write_test_data('set_num.dat',num_range,num_azim,azim,dr,grid%cluster,imiss)
   call write_grd      ('set_num.grd',num_range,num_azim,grid%cluster,imiss)

   write(nyquist_unit,"(A)") "========== Final Region Nyquist Number =========="
   num_point_dealiased   =0
   area_total_dealiased  =0
   num_region_dealiased  =0 
   num_point_aliased     =0
   area_total_aliased    =0
   num_region_aliased    =0
   minimum_nyquist_number=0
   maximum_nyquist_number=0

   nyq_num    =miss
   minimum_vel=miss
   maximum_vel=miss
   vel_out    =miss

   l=0
   area_dealiased_in_region=0
   do n=1, num_region
      if(region(n)%nyquist_number/=miss)then
         min_vel=miss
         max_vel=miss
         l=l+1
         num_region_dealiased=num_region_dealiased+1
         if(region(n)%nyquist_number/=0)then
            num_region_aliased=num_region_aliased+1
            minimum_nyquist_number=min(minimum_nyquist_number,region(n)%nyquist_number)
            maximum_nyquist_number=max(maximum_nyquist_number,region(n)%nyquist_number)
         endif
         do m=1,  region(n)%num_point
            i=region(n)%point_list(m)%i
            j=region(n)%point_list(m)%j
            if(grid%region(i,j)<1) cycle
            !if(grid%ref(i,j)<10) cycle
            !if(IAND(grid%flag(i,j),flag_pt_low_snr)>0.and.abs(grid%vel(i,j))<=v0.and.(grid%dr*i)<=r0)then
            !   write(*,*) "low_snr and low_vel and near radar",i,j,grid%flag(i,j), grid%vel(i,j), r0
            !   cycle 
            !endif
            vel_out(i,j)=vel(i,j)+2*region(n)%nyquist_number*vmax
            num_point_dealiased=num_point_dealiased+1
            area_total_dealiased=area_total_dealiased+grid%area(i,j)
            area_dealiased_in_region(n)=area_dealiased_in_region(n)+grid%area(i,j)
            nyq_num(i,j)=region(n)%nyquist_number
            if(region(n)%nyquist_number/=0)then
               num_point_aliased=num_point_aliased+1
               area_total_aliased=area_total_aliased+grid%area(i,j)
            endif
            if(minimum_vel==miss)then
               minimum_vel=vel_out(i,j)
            else
               minimum_vel=min(minimum_vel,vel_out(i,j))
            endif
            if(maximum_vel==miss)then
               maximum_vel=vel_out(i,j)
            else
               maximum_vel=max(maximum_vel,vel_out(i,j))
            endif
            
            if(min_vel==miss)then
               min_vel=vel_out(i,j)
            else
               min_vel=min(min_vel,vel_out(i,j))
            endif
            if(max_vel==miss)then
               max_vel=vel_out(i,j)
            else
               max_vel=max(max_vel,vel_out(i,j))
            endif
         enddo
         if(l<=30.or.(region(n)%nyquist_number/=0.and.region(n)%area_total>a0))then
            write(*,"(A,A,I3,A,I8,A,2F10.1,2F10.2)") region_name(n)//" "//cluster_name(region(n)%cluster),&
                                                        " nyq=", region(n)%nyquist_number, " np=",region(n)%num_point,&
                                                             " area=",region(n)%area_total,area_dealiased_in_region(n),  min_vel, max_vel
         endif
         write(nyquist_unit,"(A,A,I3,A,I8,A,F10.1,2F10.2)") region_name(n)//" "//cluster_name(region(n)%cluster), &
                                                          " nyq=", region(n)%nyquist_number, " np=",region(n)%num_point, &
                                                          " area=",region(n)%area_total, min_vel, max_vel
      endif
   enddo
   write(*,*) "Processed ", num_region_dealiased," region, "

   write(*,"(8X,A12,2(A12,A9))") "total", "dealiased","(percent)","aliased","(percent)"
   write(*,"(A8,I12,2(I12,F9.3))") "point",grid%num_valid, &
                                    num_point_dealiased, 100.*num_point_dealiased/grid%num_valid, &
                                    num_point_aliased  , 100.*num_point_aliased  /grid%num_valid
   write(*,"(A8,F12.3,2(F12.3,F9.3))") "area",grid%area_total, &
                                        area_total_dealiased, 100.*area_total_dealiased/grid%area_total, &
                                        area_total_aliased  , 100.*area_total_aliased  /grid%area_total



   call write_grid_csv ('vel.csv',num_range,num_azim,grid%dr,grid%azim,vel_out, grid%ref, grid%ref0, grid%spw, grid%flag, grid%region, grid%cluster, miss)
   ival( 1)=grid%num_valid
   ival( 3)=num_region
   ival( 4)=num_cluster
   ival( 5)= region(1)%num_point
   ival( 6)=cluster(1)%num_point
   ival( 7)=cluster(1)%num_region
   ival( 8)=num_point_dealiased
   ival( 9)=num_region_dealiased
   ival(10)=num_cluster_dealiased
   ival(11)=num_point_aliased
   ival(12)=num_region_aliased
   ival(13)=minimum_nyquist_number
   ival(14)=maximum_nyquist_number
   ival(15)=num_point_std_check
   ival(16)=n0

   rval( 1)=minimum_vel
   rval( 2)=maximum_vel
   rval( 3)=cluster(1)%mean_omb
   rval( 4)=cluster(1)% std_omb
   rval( 5)=area_total_dealiased
   rval( 6)=area_total_aliased  
   rval( 7)=area_no_nyq
   rval( 8)=grid%area_total
   rval( 9)= region(1)%area_total
   rval(10)=cluster(1)%area_total
   rval(11)=100.*area_total_dealiased/grid%area_total
   rval(12)=100.*area_total_aliased  /grid%area_total
   rval(13)=c0
   rval(14)=c1
   rval(15)=c2
   rval(16)=c3
   rval(17)=c4
   rval(18)=c5
   rval(19)=v0
   rval(20)=w0
   rval(21)=w1
   rval(22)=z0
   rval(23)=r0
   rval(24)=r1
   rval(25)=a0
   rval(26)=a1
   rval(27)=at
   
   flag=grid%flag
   call write_test_data('nyq_num.dat',num_range,num_azim,azim,dr,nyq_num,int(miss))
   !write(121) vel
   call write_grd      ('vel_new.grd',num_range,num_azim,vel_out,miss)
   !write(122) nyq_num 
   call write_grd      ('nyq_num.grd',num_range,num_azim,nyq_num,int(miss))
   
   close(cluster_unit   )
   close(region_unit    )
   close(sort_unit      )
   close(nyquist_unit   )
   close(reference_unit )
   close(check_grid_unit)
   close(check_std_unit )
   close(border_unit    )
   !stop
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine remove_region(rem_idx, num_region, region, if_delete)
   use grd, only: region_name
   implicit none
   
   integer,                               intent(in)    :: rem_idx, num_region
   type(t_region), dimension(num_region), intent(inout) :: region
   logical, intent(in) :: if_delete

   integer :: i, j

   !if(region(rem_idx)%if_removed) return
   !region(rem_idx)%if_removed     =.true.
   if(region(rem_idx)%nyquist_number ==miss) return
   region(rem_idx)%nyquist_number =miss
   !if(region(rem_idx)%flag==flag_rg_zero)then
   !   if(if_delete) region(rem_idx)%flag=flag_rg_delete
   !else
      if(if_delete) region(rem_idx)%flag=flag_rg_delete


      if(region(rem_idx)%area_total>a1)then
         write(region_unit, "(A,A)") region_name(rem_idx), " removed"
      endif
      write(nyquist_unit,"(A,A10,A)") "remove region ", region_name(rem_idx), " and regions based on it"
      do j=1, region(rem_idx)%num_neighbor
         i=region(rem_idx)%neighbor_list(j)
         !if(region(i)%if_removed) cycle
         if(region(i)%nyquist_number ==miss) cycle
         if(region(i)%base_region==rem_idx)then ! region based on remove region
            call remove_region(i, num_region, region, .false.)
         endif
      enddo
   !endif
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nyq_region(reg_idx)=nyq_region(base_idx)+nyq_diff
   subroutine get_neighbor_nyquist(num_region, region, cluster, neighbor_idx, base_idx, reg_idx, nyq_diff)
   implicit none
   type(t_cluster),                       intent(in)  :: cluster
   integer,                               intent(in)  :: neighbor_idx, reg_idx, base_idx, num_region
   type(t_region), dimension(num_region), intent(in)  :: region
   integer,                               intent(out) :: nyq_diff

   nyq_diff=miss

   !if(region(base_idx)%if_removed.or.region(reg_idx)%if_removed) return
   !if(region( reg_idx)%nyquist_number==miss) return
   if(region(base_idx)%nyquist_number==miss) return
   !if(region(base_idx)%flag ==flag_rg_zero) return
   
   if(cluster%neighbor_list(neighbor_idx,1)/=base_idx.and.cluster%neighbor_list(neighbor_idx,2)/= reg_idx) return
   if(cluster%num_border_point (neighbor_idx)<1) return
   if(cluster%num_border_effect(neighbor_idx)<1) return

   
   if(cluster%neighbor_flag(neighbor_idx)==flag_cont)then
      nyq_diff=0
   else !if(region(reg_idx)%num_point>nrmin.and.cluster%num_border_point(neighbor_idx)>3)then
      if(cluster%neighbor_flag(neighbor_idx)==flag_pos_alias)then
         nyq_diff=1
      elseif(cluster%neighbor_flag(neighbor_idx)==flag_neg_alias)then
         nyq_diff=-1
      endif
   endif
   

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_nyquist_number(num_region, region, cluster, neighbor_idx, reg_idx, chk_idx)
   use grd, only: region_name
   implicit none
   type(t_cluster),                       intent(in)  :: cluster
   integer,                               intent(in)  :: neighbor_idx, reg_idx, chk_idx, num_region
   type(t_region), dimension(num_region), intent(inout)  :: region

   logical :: passed
   integer :: delete_idx, k, i

   if(region(reg_idx)%nyquist_number==miss) return
   if(region(chk_idx)%nyquist_number==miss) return
   !if(region(reg_idx)%if_removed)then
   !   region(reg_idx)%nyquist_number=miss
   !   return
   !endif
   !write(166,*) "in check_nyquist_number:", reg_idx,chk_idx,reg%vmax

   delete_idx=0
   passed=.true.

   call get_neighbor_nyquist(num_region, region, cluster, neighbor_idx, reg_idx, chk_idx, k)
   if(k==miss) return
   if(region(chk_idx)%nyquist_number/=region(reg_idx)%nyquist_number+k)then
      passed=.false. 
   endif
   if(.NOT.passed)then
      if(chk_idx<region(reg_idx)%base_region.and.region(reg_idx)%base_border<cluster%num_border_effect(neighbor_idx))then ! base region maybe wrong
         delete_idx=region(reg_idx)%base_region
         region(reg_idx)%nyquist_number=miss
         region(reg_idx)%base_region   =imiss
         region(reg_idx)%base_border   =0
         write(nyquist_unit,*) "set miss: Failed check", region_name(reg_idx), "based on ", region_name(neighbor_idx)
!         call remove_region(reg_idx, num_region, region)
      else
         !if(region(reg_idx)%base_border>region(chk_idx)%base_border)then ! set small region to miss
         !   delete_idx=chk_idx
         !elseif(region(reg_idx)%base_border<region(chk_idx)%base_border)then
         !   delete_idx=reg_idx
         !else
            if(reg_idx<chk_idx)then ! set small region to miss
               delete_idx=chk_idx
            else
               delete_idx=reg_idx
            endif
         !endif
      endif

      write(nyquist_unit,"(2(A,A10))") "check region failed ", region_name(reg_idx), " ",region_name(chk_idx)
      write(nyquist_unit,"(A)") " Nyquist#1 Nyquist#2  Npoint#1  Npoint#2 Flag  Delete"
      write(nyquist_unit,"(4I10,I5,1X,A7)") region(reg_idx)%nyquist_number, region(chk_idx)%nyquist_number, &
                  region(reg_idx)%num_point,      region(chk_idx)%num_point, &
                  cluster%neighbor_flag(neighbor_idx) , region_name(delete_idx)

      if(delete_idx>0)then
         if(region(delete_idx)%area_total>a1) write(region_unit,"(2A)") "==== nyquist check remove region ==== ", region_name(delete_idx)
         call remove_region(delete_idx, num_region, region, .true.)
      endif
   else
      write(nyquist_unit,"(2(A,A10,I5))") "check region passed ", region_name(reg_idx), region(reg_idx)%level," ",region_name(chk_idx), region(chk_idx)%level
      if(region(reg_idx)%base_region/=0)then
         if(region(chk_idx)%level<region(region(reg_idx)%base_region)%level)then
            region(reg_idx)%base_region=chk_idx
            region(reg_idx)%level=region(chk_idx)%level+1 
         endif
      endif
      
   endif

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculate reg_idx's neigh  bor's nyquist number while nyq_region(reg_idx) is known
   recursive subroutine calc_neighbor_nyquist_number(num_region, region, cluster, reg_idx)
   use grd, only: region_name
   implicit none
   ! region(reg_idx)'s nyq_region calculates based on region(base_idx)
   integer,                              intent(in)    :: num_region, reg_idx
   type(t_region), dimension(num_region),intent(inout) :: region
   type(t_cluster),                      intent(in)    :: cluster
   
!   real :: c1, c2
   integer :: i, j, k, base_idx
  
   if(region(reg_idx)%nyquist_number==miss) return 
   ! check first
   base_idx=region(reg_idx)%base_region
   if(base_idx>0)then
      do j=1, cluster%num_neighbor
         if(region(reg_idx)%nyquist_number==miss) return 
         if(cluster%neighbor_list(cluster%sorted_border_index(j),1)/=reg_idx) cycle
         i= cluster%neighbor_list(cluster%sorted_border_index(j),2)
         if(i==base_idx) cycle
         if(i== reg_idx) cycle
         !if(region(i)%if_removed) cycle
         
         if(cluster%num_border_point (cluster%sorted_border_index(j))<1) cycle
         if(cluster%num_border_effect(cluster%sorted_border_index(j))<1) cycle
         
         if(region(i)%nyquist_number/=miss)then ! check
            call check_nyquist_number(num_region, region, cluster, cluster%sorted_border_index(j), reg_idx, i)
         endif
      enddo
   endif

   if(region(reg_idx)%nyquist_number==miss) return
   !if(region(reg_idx)%if_removed) return
   ! calc neighbor nyquist number based on reg_idx 
   do j=1, cluster%num_neighbor
      !write(*,*)  "num neighbor:", j, cluster%num_neighbor, ubound(cluster%neighbor_list)
      if(cluster%neighbor_list(cluster%sorted_border_index(j),1)/=reg_idx) cycle
      i= cluster%neighbor_list(cluster%sorted_border_index(j),2)
      if(i==base_idx) cycle
      if(i== reg_idx) cycle
      if(region(i)%flag== flag_rg_delete) cycle
      if(i<reg_idx.and.cluster%num_border_effect(cluster%sorted_border_index(j))<1) cycle
      if(cluster%num_border_effect(cluster%sorted_border_index(j))<1) cycle
      !if(reg%point_list(i)%num_point > reg%point_list(reg_idx)%num_point) cycle ! only big region can be little region base
      !if(region(i)%if_removed) cycle
      if(cluster%num_border_point (cluster%sorted_border_index(j))<1) cycle
      if(cluster%num_border_effect(cluster%sorted_border_index(j))<1) cycle
      if(region(i)%level/=miss.and.region(reg_idx)%level>=region(i)%level) cycle

      if(region(i)%nyquist_number==miss)then ! calculate
         call get_neighbor_nyquist(num_region, region, cluster,cluster%sorted_border_index(j),reg_idx,i,k)
         if(k/=miss)then
            region(i)%nyquist_number=region(reg_idx)%nyquist_number+k
            write(nyquist_unit,"(2A10,7I10)")  region_name(i), region_name(reg_idx), region(reg_idx)%nyquist_number, &
                                 region(i)%nyquist_number, cluster%num_border_point(cluster%sorted_border_index(j)),&
                                 cluster%num_border_effect(cluster%sorted_border_index(j)), &
                                 region(reg_idx)%num_point, region(i)%num_point, &
                                 cluster%neighbor_flag(cluster%sorted_border_index(j))
            region(i)%base_region=reg_idx
            region(i)%base_border=cluster%num_border_effect(cluster%sorted_border_index(j))
            region(i)%level      =region(reg_idx)%level+1
            call calc_neighbor_nyquist_number(num_region, region, cluster, i)
         else
            write(nyquist_unit,"((2A10,I10,A10,5I10))") region_name(i), region_name(reg_idx), region(reg_idx)%nyquist_number,&
                            " not set", cluster%num_border_point(cluster%sorted_border_index(j)), &
                            cluster%num_border_effect(cluster%sorted_border_index(j)), &
                            region(reg_idx)%num_point,region(i)%num_point, &
                            cluster%neighbor_flag(cluster%sorted_border_index(j))
         endif
         
      endif
   enddo
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_neighbor_point_flag(grid,i1,j1,i2,j2,flag,ar1,ar2)
   implicit none
   type(t_grid), intent(inout) :: grid
   integer, intent(in) :: i1, j1, i2, j2
   integer, intent(out) :: flag
   real, intent(in) :: ar1, ar2

   integer :: uflag1, uflag2

   uflag1=ior(flag_pt_used,flag_pt_low_vel)
   uflag2=ior(ior(flag_pt_used,flag_pt_low_vel),flag_pt_low_snr)
!  write(999,*) "uflag", uflag 
   flag=flag_unknown
   !if((abs(grid%vel(i1,j1))<=v0).or.(abs(grid%vel(i2,j2))<=v0))then !.and.ar2<a0.and.ar1<a0 low vel, abs(v1-v2)<=(c1/3) 
   if((grid%flag(i1,j1)==uflag1.and.grid%flag(i2,j2)==flag_pt_used).or.(grid%flag(i1,j1)==flag_pt_used.and.grid%flag(i2,j2)==uflag1))then 
      if(abs(grid%vel(i1,j1)-grid%vel(i2,j2))>c6)then
         flag=flag_noise
      else
         flag=flag_cont
      endif
      return
   endif

   if(grid%flag(i1,j1)/=flag_pt_used.or.grid%flag(i2,j2)/=flag_pt_used)then ! dont use this point, effect -1
      flag=flag_noise
      return
   endif

   if(abs(grid%vel(i1,j1)-grid%vel(i2,j2))<=c1)then ! conti .or.(abs(grid%level(i1,j1)-grid%level(i2,j2))<=1)
      flag=flag_cont
   elseif(grid%vel(i2,j2)-grid%vel(i1,j1)<=-c2)then !grid%level(i1,j1)==max_level.and.grid%level(i2,j2)==min_level)then ! ! nyq(n2)=nyq(n1)+1
      flag=flag_pos_alias
   elseif(grid%vel(i2,j2)-grid%vel(i1,j1)>= c2)then !grid%level(i2,j2)==max_level.and.grid%level(i1,j1)==min_level)then ! ! nyq(n2)=nyq(n1)-1
      flag=flag_neg_alias
   elseif(abs(grid%vel(i1,j1)-grid%vel(i2,j2))>c1.and.abs(grid%vel(i1,j1)-grid%vel(i2,j2))<c2)then ! mix
      flag=flag_mix
   endif
   
   if(flag==flag_unknown)then
      if(grid%vel(i2,j2)-grid%vel(i1,j1)<=-c2)then
         flag=flag_pos_alias
      elseif(grid%vel(i2,j2)-grid%vel(i1,j1)>= c2)then
         flag=flag_neg_alias
      endif
   endif

   if(abs(grid%vel(i1,j1)-grid%vel(i2,j2))>c1.and.abs(grid%vel(i1,j1)-grid%vel(i2,j2))<c2)then
     if(abs(grid%vel(i1,j1))<=v0.or.abs(grid%vel(i2,j2))<=v0)then
        flag=flag_noise
        if(abs(grid%vel(i1,j1))<=v0)then
           grid%flag(i1,j1)=IOR(grid%flag(i1,j1),flag_pt_low_vel)
        endif
        if(abs(grid%vel(i2,j2))<=v0)then
           grid%flag(i2,j2)=IOR(grid%flag(i2,j2),flag_pt_low_vel)
        endif
     endif
   endif

   if(abs(grid%vel(i1,j1))<=v0.or.abs(grid%vel(i2,j2))<=v0.and.flag==flag_cont)then !< c1/2
      if(abs(grid%vel(i1,j1)-grid%vel(i2,j2))>c1/2)then
        flag=flag_noise 
        if(abs(grid%vel(i1,j1))<=v0)then
           grid%flag(i1,j1)=IOR(grid%flag(i1,j1),flag_pt_low_vel)
        endif
        if(abs(grid%vel(i2,j2))<=v0)then
           grid%flag(i2,j2)=IOR(grid%flag(i2,j2),flag_pt_low_vel)
        endif
      endif
   endif
   !if(flag==flag_cont)then
   !   if(abs(grid%vel(i1,j1))<=v0.or.abs(grid%vel(i2,j2))<=v0)then
   !      flag=flag_zero
   !   endif
   !endif
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine update_neighbor_point_flag(flag, n1, n2, num_region, neighbor_flag, num_effect)
   implicit none
   integer,                                   intent(in)    :: flag, n1, n2, num_region
   integer, dimension(num_region,num_region), intent(inout) :: neighbor_flag, num_effect

   if(flag==flag_cont)then ! conti
      if(neighbor_flag(n1,n2)==flag_unknown)then
         neighbor_flag(n1,n2)=flag_cont
      elseif(neighbor_flag(n1,n2)==flag_cont)then
         neighbor_flag(n1,n2)=flag_cont
      elseif(neighbor_flag(n1,n2)==flag_zero)then
         neighbor_flag(n1,n2)=flag_cont
      else
         neighbor_flag(n1,n2)=flag_mix
      endif

      if(neighbor_flag(n2,n1)==flag_unknown)then
         neighbor_flag(n2,n1)=flag_cont
      elseif(neighbor_flag(n2,n1)==flag_cont)then
         neighbor_flag(n2,n1)=flag_cont
      elseif(neighbor_flag(n2,n1)==flag_zero)then
         neighbor_flag(n2,n1)=flag_cont
      else
         neighbor_flag(n2,n1)=flag_mix
      endif
   elseif(flag==flag_zero)then 

      if(neighbor_flag(n1,n2)==flag_unknown)then
         neighbor_flag(n1,n2)=flag_zero
      elseif(neighbor_flag(n1,n2)==flag_cont)then
         neighbor_flag(n1,n2)=flag_cont
      elseif(neighbor_flag(n1,n2)==flag_zero)then
         neighbor_flag(n1,n2)=flag_zero
      else
         neighbor_flag(n1,n2)=flag_mix
      endif
      if(neighbor_flag(n2,n1)==flag_unknown)then
         neighbor_flag(n2,n1)=flag_zero
      elseif(neighbor_flag(n2,n1)==flag_cont)then
         neighbor_flag(n2,n1)=flag_cont
      elseif(neighbor_flag(n2,n1)==flag_zero)then
         neighbor_flag(n2,n1)=flag_zero
      else
         neighbor_flag(n2,n1)=flag_mix
      endif

   elseif(flag==flag_pos_alias)then ! nyq(n2)=nyq(n1)+1
   !if(abs(grid%vel(i2,j2))>c5.or.abs(grid%vel(i1,j1))>c2)then
   
      if(neighbor_flag(n1,n2)==flag_unknown)then
         neighbor_flag(n1,n2)=flag_pos_alias
      elseif(neighbor_flag(n1,n2)==flag_pos_alias)then
         neighbor_flag(n1,n2)=flag_pos_alias
      elseif(neighbor_flag(n1,n2)==flag_zero)then
         neighbor_flag(n1,n2)=flag_pos_alias
      else
         neighbor_flag(n1,n2)=flag_mix
      endif
      if(neighbor_flag(n2,n1)==flag_unknown)then
         neighbor_flag(n2,n1)=flag_neg_alias
      elseif(neighbor_flag(n2,n1)==flag_neg_alias)then
         neighbor_flag(n2,n1)=flag_neg_alias
      elseif(neighbor_flag(n2,n1)==flag_zero)then
         neighbor_flag(n2,n1)=flag_neg_alias
      else
         neighbor_flag(n2,n1)=flag_mix
      endif
   elseif(flag==flag_neg_alias)then ! nyq(n2)=nyq(n1)-1
   
      if(neighbor_flag(n1,n2)==flag_unknown)then
         neighbor_flag(n1,n2)=flag_neg_alias
      elseif(neighbor_flag(n1,n2)==flag_neg_alias)then
         neighbor_flag(n1,n2)=flag_neg_alias
      elseif(neighbor_flag(n1,n2)==flag_zero)then
         neighbor_flag(n1,n2)=flag_neg_alias
      else
         neighbor_flag(n1,n2)=flag_mix
      endif
      if(neighbor_flag(n2,n1)==flag_unknown)then
         neighbor_flag(n2,n1)=flag_pos_alias
      elseif(neighbor_flag(n2,n1)==flag_pos_alias)then
         neighbor_flag(n2,n1)=flag_pos_alias
      elseif(neighbor_flag(n2,n1)==flag_zero)then
         neighbor_flag(n2,n1)=flag_pos_alias
      else
         neighbor_flag(n2,n1)=flag_mix
      endif
   elseif(flag==flag_mix)then ! mix
   !endif
      neighbor_flag(n1,n2)=flag_mix
      neighbor_flag(n2,n1)=flag_mix
   elseif(flag==flag_noise.or.flag==flag_zero)then
      num_effect(n1,n2)=num_effect(n1,n2)-1
      num_effect(n2,n1)=num_effect(n2,n1)-1
   endif

   if(neighbor_flag(n1,n2)==flag_unknown)then
      if(flag/=flag_noise.and.flag/=flag_zero) neighbor_flag(n1,n2)=flag_mix
   endif
   if(neighbor_flag(n2,n1)==flag_unknown)then
      if(flag/=flag_noise.and.flag/=flag_zero) neighbor_flag(n2,n1)=flag_mix
   endif
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine get_neighbor_flag(grid, num_region, region, num_cluster, cluster, write_border_list)
!   implicit none
!   type(t_grid),                          intent(inout) :: grid
!   integer,                               intent(in)    :: num_region, num_cluster
!   type(t_region), dimension(num_region), intent(in)    :: region
!   type(t_cluster),dimension(num_cluster),intent(inout) :: cluster
!   logical,                               intent(in) :: write_border_list
!  
!   integer :: i, j, n, num_range, num_azim, n1, n2, i1, j1, i2, j2, flag
!   integer, dimension(:), allocatable :: a, d
!
!   integer,             dimension(:,:), allocatable :: num_point, num_effect
!   integer,             dimension(:,:), allocatable :: neighbor_flag
!
!   integer, parameter :: nb_num=300
!   type(t_border_list), dimension(nb_num,nb_num) :: border_list 
!
!   allocate(num_point    (num_region,num_region))
!   allocate(num_effect   (num_region,num_region))
!   allocate(neighbor_flag(num_region,num_region))
!   !write(*,*) "in get_neighbor_flag"
!   num_range=grid%num_range
!   num_azim =grid%num_azim
!
!   !write(*,*) "count border point..."
!   num_point =0
!   num_effect=0
!   !write(border_unit,*) "+==== Border List 1 ============"
!   neighbor_flag=flag_unknown
!!   write(204,*) "start_i, end_i",grid%start_i, grid%end_i
!   do j=1, num_azim
!      if(grid%start_i(j)<1) cycle
!      do i=grid%start_i(j), grid%end_i(j)-1
!!         if(i==57.and.j==269)then
!!            write(*,*) "P-0057_269",grid%region(i,j), grid%region(i+1,j), grid%region(i,j+1)
!!         endif
!         i1=i+1
!         j1=j+1
!         if(j1>num_azim) j1=1
!         !if(abs(grid%vel(i,j))<=2) cycle
!         !if(abs(grid%vel(i1,j))>2) then
!         if(grid%region(i,j)/=imiss.and.grid%region(i1,j)/=imiss.and.(grid%region(i,j)/=grid%region(i1,j)))then
!            i2=i1
!            j2=j
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!            !if(write_border_list.and.n1<=nb_num.and.n2<=nb_num)then
!            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
!            !endif
!
!            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag)
!            call update_neighbor_point_flag(flag, n1, n2, num_region, neighbor_flag, num_effect)
!
!         endif
!         !endif
!         i1=i+1
!         j1=j+1
!         if(j1>num_azim) j1=1
!         !if(abs(grid%vel(i,j1))>2) then
!         if(grid%region(i,j)/=imiss.and.grid%region(i,j1)/=imiss.and.(grid%region(i,j)/=grid%region(i,j1)))then
!            i2=i
!            j2=j1
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!            !if(write_border_list.and.n1<=nb_num.and.n2<=nb_num)then
!            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
!            !endif
!            
!            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag)
!            call update_neighbor_point_flag(flag, n1, n2, num_region, neighbor_flag, num_effect)
!
!         endif
!         !endif
!      enddo
!   enddo
!
!   num_effect=num_point+num_effect
!   
!! count num of neighbor
!!  write(*,*) ubound(cluster)
!   cluster(:)%num_neighbor=0
!   do i=1, num_region
!      do j=1, num_region
!         if(num_effect(i,j)<n0) cycle
!         n=region(i)%cluster
!!        write(*,*) "count neighbor", n,  cluster(n)%num_neighbor
!         cluster(n)%num_neighbor=cluster(n)%num_neighbor+1
!      enddo
!   enddo
!
!   ! allocate
!   do n=1, num_cluster
!      if(allocated( cluster(n)%neighbor_flag))then
!         deallocate(cluster(n)%neighbor_flag)
!         deallocate(cluster(n)%neighbor_list)
!         deallocate(cluster(n)%num_border_point)
!         deallocate(cluster(n)%num_border_effect)
!         deallocate(cluster(n)%sorted_border_index)
!      endif
!      if(cluster(n)%num_neighbor>0)then
!         allocate(cluster(n)%neighbor_flag      (cluster(n)%num_neighbor))
!         allocate(cluster(n)%num_border_point   (cluster(n)%num_neighbor))
!         allocate(cluster(n)%num_border_effect  (cluster(n)%num_neighbor))
!         allocate(cluster(n)%neighbor_list      (cluster(n)%num_neighbor,2))
!         allocate(cluster(n)%sorted_border_index(cluster(n)%num_neighbor))
!      endif
!   enddo
!
!!   write(98,*) "cluster num_neighbor 1", (cluster(n)%num_region,cluster(n)%num_neighbor,cluster(n)%num_point,',',n=1,num_cluster)
!
!   cluster(:)%num_neighbor=0
!   do i=1, num_region
!      do j=1, num_region
!         if(num_effect(i,j)<n0) cycle
!         n=region(i)%cluster
!         cluster(n)%num_neighbor    =cluster(n)%num_neighbor+1
!         cluster(n)%num_border_point (cluster(n)%num_neighbor  )=num_point (i,j)
!         cluster(n)%num_border_effect(cluster(n)%num_neighbor  )=num_effect(i,j)
!         cluster(n)%neighbor_list   (cluster(n)%num_neighbor,1)=i
!         cluster(n)%neighbor_list   (cluster(n)%num_neighbor,2)=j
!         cluster(n)%neighbor_flag   (cluster(n)%num_neighbor  )=neighbor_flag(i,j)
!      enddo
!   enddo
!   write(cluster_unit,"(A)") "================ Cluster Information ==============="
!   write(cluster_unit,"(A)") "   cluster    region  neighbor   cluster     cluster"
!   write(cluster_unit,"(A)") "      name    number    number num point  area total"
!   do n=1, num_cluster
!      !write(*,*) "cluster info",  n
!      write(cluster_unit,"(A10,3I10,F12.3)") cluster_name(n),cluster(n)%num_region,&
!            cluster(n)%num_neighbor,cluster(n)%num_point, cluster(n)%area_total
!   enddo
!   write(cluster_unit,"(A)") "================ Cluster Neighbor ================"
!   do n=1, num_cluster
!      !write(*,*) "cluster list",  n
!      if(cluster(n)%num_neighbor<1) cycle
!      write(cluster_unit,"(A)") "---- region neighbor list for cluster:"//cluster_name(n)//" ----"
!      write(cluster_unit,"(A)") "  neighbor  neighbor border_pt border_pt  neighbor"
!      write(cluster_unit,"(A)") "  region#1  region#2    number    effect      flag"
!      allocate(a(cluster(n)%num_neighbor))
!      allocate(d(cluster(n)%num_neighbor))
!      do i=1, cluster(n)%num_neighbor
!         write(cluster_unit,"(2A10,2I10,I10)") &
!               region_name(cluster(n)%neighbor_list(i,1)), &
!               region_name(cluster(n)%neighbor_list(i,2)), &
!               cluster(n)%num_border_point(i),             &
!               cluster(n)%num_border_effect(i),            &
!               cluster(n)%neighbor_flag(i)
!
!         a(i)=cluster(n)%num_border_effect(i)
!         d(i)=i
!      enddo
!      call quick_sort(a,d,cluster(n)%num_neighbor,1,cluster(n)%num_neighbor)
!      cluster(n)%sorted_border_index=d
!      deallocate(a)
!      deallocate(d)
!   enddo
!   !write(98,*) "cluster num_neighbor 2", cluster(:)%num_neighbor
!    
!   if(.NOT.write_border_list) return 
!   !allocate(border_list(num_region,num_region))
!   !write(*,*) "creating border list..."
!   do i=1, min(num_region,nb_num)
!      do j=1, min(num_region,nb_num)
!         border_list(i,j)%num_point=num_point(i,j)
!         if(allocated (border_list(i,j)%i1))then
!            deallocate(border_list(i,j)%i1)
!            deallocate(border_list(i,j)%j1)
!            deallocate(border_list(i,j)%i2)
!            deallocate(border_list(i,j)%j2)
!         endif
!         allocate(border_list(i,j)%i1(border_list(i,j)%num_point))
!         allocate(border_list(i,j)%j1(border_list(i,j)%num_point))
!         allocate(border_list(i,j)%i2(border_list(i,j)%num_point))
!         allocate(border_list(i,j)%j2(border_list(i,j)%num_point))
!      enddo
!   enddo
!   num_point=0
!   !write(border_unit,*) "+==================================="
!   !write(border_unit,*) "+==== Border List 2 ============"
!   do j=1, num_azim
!      if(grid%start_i(j)<1) cycle
!      do i=grid%start_i(j), grid%end_i(j)-1
!         i1=i+1
!         j1=j+1
!         if(j1>num_azim) j1=1
!         !if(abs(grid%vel(i,j))<=2) cycle
!         !if(abs(grid%vel(i1,j))>2) then
!         if(grid%region(i,j)/=imiss.and.grid%region(i1,j)/=imiss.and.(grid%region(i,j)/=grid%region(i1,j)))then
!            i2=i1
!            j2=j
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!            !if(n1<=nb_num.and.n2<=nb_num)then
!            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
!            !endif
!            if(grid%region(i1,j1)>=1.and.grid%region(i1,j1)<=nb_num.and. &
!               grid%region(i2,j2)>=1.and.grid%region(i2,j2)<=nb_num) then
!               border_list(n1,n2)%i1(num_point(n1,n2))=i1
!               border_list(n1,n2)%i2(num_point(n1,n2))=i2
!               border_list(n1,n2)%j1(num_point(n1,n2))=j1
!               border_list(n1,n2)%j2(num_point(n1,n2))=j2
!               border_list(n2,n1)%i1(num_point(n2,n1))=i2
!               border_list(n2,n1)%i2(num_point(n2,n1))=i1
!               border_list(n2,n1)%j1(num_point(n2,n1))=j2
!               border_list(n2,n1)%j2(num_point(n2,n1))=j1
!            endif
!         endif
!         i1=i+1
!         !endif
!         j1=j+1
!         if(j1>num_azim) j1=1
!         !if(abs(grid%vel(i,j1))>2) then
!         !if(grid%region(i,j1)>=1.and.grid%region(i,j1)<=30)then
!         if(grid%region(i,j)/=imiss.and.grid%region(i,j1)/=imiss.and.(grid%region(i,j)/=grid%region(i,j1)))then
!            i2=i
!            j2=j1
!            i1=i
!            j1=j
!            n1=grid%region(i1,j1)
!            n2=grid%region(i2,j2)
!            num_point(n1,n2)=num_point(n1,n2)+1
!            num_point(n2,n1)=num_point(n2,n1)+1
!            !if(n1<=nb_num.and.n2<=nb_num)then
!            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
!            !endif
!            if(grid%region(i1,j1)>=1.and.grid%region(i1,j1)<=nb_num .and.&
!               grid%region(i2,j2)>=1.and.grid%region(i2,j2)<=nb_num) then
!
!            if(num_point(n1,n2)>border_list(n1,n2)%num_point)stop
!!write(*,*) "border_list i1",n1,n2,num_point(n1,n2),ubound(border_list(n1,n2)%i1)
!               border_list(n1,n2)%i1(num_point(n1,n2))=i1
!               border_list(n1,n2)%i2(num_point(n1,n2))=i2
!               border_list(n1,n2)%j1(num_point(n1,n2))=j1
!               border_list(n1,n2)%j2(num_point(n1,n2))=j2
!               border_list(n2,n1)%i1(num_point(n2,n1))=i2
!               border_list(n2,n1)%i2(num_point(n2,n1))=i1
!               border_list(n2,n1)%j1(num_point(n2,n1))=j2
!               border_list(n2,n1)%j2(num_point(n2,n1))=j1
!            endif
!         endif
!         !endif
!      enddo
!   enddo
!
!   !write(*,*) "wrting border list..."
!   write(border_unit,*) "======== Region Border: c1=", c1, " c2=",c2, "========"
!   do i=1, min(num_region,nb_num)
!      do j=i+1, min(num_region,nb_num)
!         if(border_list(i,j)%num_point<1)cycle
!         write(border_unit,"(2(A,1X),3(A,I5))") region_name(i), region_name(j) , &
!                                                " Border_Pt",num_point    (i,j), &
!                                                " Effect_Pt",num_effect   (i,j), &
!                                                " flag"     ,neighbor_flag(i,j)
!
!         write(border_unit,"(2(A10,1X),6A7,2A4,2A3,A5,A7,3A5)")"Point#1","Point#2","v1","v2","z1","z2","w1","w2","f1","f2","l1","l2","flag","v1-v2","<c1",">c2","<-c2"
!         do n=1, border_list(i,j)%num_point
!            i1=border_list(i,j)%i1(n)
!            j1=border_list(i,j)%j1(n)
!            i2=border_list(i,j)%i2(n)
!            j2=border_list(i,j)%j2(n)
!            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag)
!            write(border_unit,"(2(A10,1X),6F7.1,2I4,2I3,I5,F7.1,3L5)") &
!                  point_name(border_list(i,j)%i1(n),border_list(i,j)%j1(n)), &
!                  point_name(border_list(i,j)%i2(n),border_list(i,j)%j2(n)), &
!                  grid%vel(i1,j1),grid%vel(i2,j2),&
!                  grid%ref(i1,j1),grid%ref(i2,j2),&
!                  grid%spw(i1,j1),grid%spw(i2,j2),&
!                  grid%flag(i1,j1),grid%flag(i2,j2),&
!                  grid%level(i1,j1),grid%level(i2,j2),flag,&
!                  grid%vel(i1,j1)-grid%vel(i2,j2),&
!                  abs(grid%vel(i1,j1)-grid%vel(i2,j2))<=c1,&
!                  (grid%vel(i1,j1)-grid%vel(i2,j2))>= c2, &
!                  (grid%vel(i1,j1)-grid%vel(i2,j2))<=-c2
!            
!         enddo
!      enddo
!   enddo
!
!   !deallocate(border_list)
!   deallocate(num_point    )
!   deallocate(num_effect   )
!   deallocate(neighbor_flag)
!   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_cluster(grid,old,new)
   implicit none
   type(t_grid),         intent(inout) :: grid
   integer,              intent(in)    :: old, new

   integer :: i, j

   do j=1, grid%num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)
         if(grid%cluster(i,j)==old)then
            grid%cluster(i,j)=new
         endif
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_region(grid,old,new)
   implicit none
   type(t_grid),         intent(inout) :: grid
   integer,              intent(in)    :: old, new

   integer :: i, j

   do j=1,  grid%num_azim
      do i=1, grid%num_range
      !if(   grid%start_i(j)<1) cycle
      !do i= grid%start_i(j), grid%end_i(j)
         if(grid%region(i,j)==old)then
            grid%region(i,j)=new
         endif
      enddo
   enddo
   end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine create_cluster(grid, num_region, region, num_cluster, cluster)
   use libradar, only: quick_sort
   use grd, only: cluster_name, region_name, point_name
   implicit none
   type(t_grid),                               intent(inout) :: grid
   integer,                                    intent(in)    :: num_region
   type(t_region), dimension(num_region),      intent(inout) :: region
   integer,                                    intent(out)   :: num_cluster
   type(t_cluster), dimension(:), allocatable, intent(inout) :: cluster

   integer :: i, j, k, n, num_range, num_azim, n1, n2, i1, j1, i2, j2, flag, s, e
   integer, dimension(:), allocatable :: a, d

   integer,             dimension(:,:), allocatable :: num_point, num_effect
   integer,             dimension(:,:), allocatable :: neighbor_flag

   integer, parameter :: nb_num=300
   type(t_border_list), dimension(nb_num,nb_num) :: border_list 

   ! count neighbor point and flag
   allocate(num_point    (num_region,num_region))
   allocate(num_effect   (num_region,num_region))
   allocate(neighbor_flag(num_region,num_region))

   num_range=grid%num_range
   num_azim =grid%num_azim

   !write(*,*) "count border point..."
   num_point =0
   num_effect=0
   
   neighbor_flag=flag_unknown

   do j=1, num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)-1
!         if(i==57.and.j==269)then
!            write(*,*) "P-0057_269",grid%region(i,j), grid%region(i+1,j), grid%region(i,j+1)
!         endif
         i1=i+1
         j1=j+1
         if(j1>num_azim) j1=1
         !if(abs(grid%vel(i,j))<=2) cycle
         !if(abs(grid%vel(i1,j))>2) then
         if(grid%region(i,j)/=imiss.and.grid%region(i1,j)/=imiss.and.(grid%region(i,j)/=grid%region(i1,j)))then
            i2=i1
            j2=j
            i1=i
            j1=j
            n1=grid%region(i1,j1)
            n2=grid%region(i2,j2)
            num_point(n1,n2)=num_point(n1,n2)+1
            num_point(n2,n1)=num_point(n2,n1)+1

            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag,region(n1)%area_total,region(n2)%area_total)
            call update_neighbor_point_flag(flag, n1, n2, num_region, neighbor_flag, num_effect)

         endif
         !endif
         i1=i+1
         j1=j+1
         if(j1>num_azim) j1=1
         !if(abs(grid%vel(i,j1))>2) then
         if(grid%region(i,j)/=imiss.and.grid%region(i,j1)/=imiss.and.(grid%region(i,j)/=grid%region(i,j1)))then
            i2=i
            j2=j1
            i1=i
            j1=j
            n1=grid%region(i1,j1)
            n2=grid%region(i2,j2)
            num_point(n1,n2)=num_point(n1,n2)+1
            num_point(n2,n1)=num_point(n2,n1)+1
            
            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag,region(n1)%area_total,region(n2)%area_total)
            call update_neighbor_point_flag(flag, n1, n2, num_region, neighbor_flag, num_effect)

         endif
         !endif
      enddo
   enddo

   num_effect=num_point+num_effect

   ! build region neighbor relation
   do i=1, num_region
      region(i)%num_neighbor  =0
      do j=1, num_region
         if(abs(neighbor_flag(i,j))<=1.and.num_effect(i,j)>=n0)then
            region(i)%num_neighbor=region(i)%num_neighbor+1
         endif
      enddo
   enddo

   do i=1, num_region
      if(allocated (region(i)%neighbor_list))then
         deallocate(region(i)%neighbor_list)
      endif
      allocate(region(i)%neighbor_list(region(i)%num_neighbor))
      n=0
      do j=1, num_region
         if(abs(neighbor_flag(i,j))<=1.and.num_effect(i,j)>=n0)then
            n=n+1
            region(i)%neighbor_list(n)=j
         endif
      enddo
   enddo

   ! set cluster index of region
   region(:)%cluster =imiss

   num_cluster=0
   do i=1, num_region
      if(region(i)%cluster==imiss)then
         num_cluster=num_cluster+1
         call include_region(num_region, region, i, num_cluster)
      endif
   enddo


   ! create cluster

   !write(*,*) "update cluster... "

   if(allocated(cluster))then
      deallocate(cluster)
   endif
   allocate(cluster(num_cluster))

   ! region(:)%cluster already defined
   ! sort cluster by area 
   if(allocated(a))then
      deallocate(a)
      deallocate(d)
   endif
   allocate(a(num_cluster))
   allocate(d(num_cluster))

   a=0
   do i=1, num_region
      k=region(i)%cluster
      a(k)=a(k)+region(i)%area_total
      !a(k)=a(k)+region(i)%num_point
   enddo

   do k=1, num_cluster
      d(k)=k
      cluster(k)%num_region =0
      cluster(k)%num_point  =0
      cluster(k)%area_total =0
   enddo

   !a=-a
   s=1
   e=num_cluster
   
   call quick_sort(a,d,num_cluster,s,e)
   !a=-a
   
   !change region(:)%cluster to sorted 
   do i=1, num_region
      if(region(i)%cluster==imiss) region(i)%cluster=0
   enddo

   do i=1, num_region
      k=region(i)%cluster
      do n=1, num_cluster
         if(k==d(n))then
            region(i)%cluster=-n
            exit
         endif
      enddo
   enddo

   ! set cluster
cluster(:)%minv=miss
cluster(:)%maxv=miss
   do i=1, num_region
      if(region(i)%cluster/=0)then
         region(i)%cluster=-region(i)%cluster
         k=region(i)%cluster
         if(cluster(k)%maxv==miss)then
            cluster(k)%maxv=region(i)%maxv
         else
            cluster(k)%maxv=max(region(i)%maxv,cluster(k)%maxv)
         endif
         if(cluster(k)%minv==miss)then
            cluster(k)%minv=region(i)%minv
         else
            cluster(k)%minv=min(region(i)%minv,cluster(k)%minv)
         endif
         cluster(k)%num_region=cluster(k)%num_region+1
         cluster(k)%num_point =cluster(k)%num_point +region(i)%num_point
         cluster(k)%area_total=cluster(k)%area_total+region(i)%area_total
      endif
   enddo
   do i=1, num_region
      if(region(i)%cluster==0) region(i)%cluster=imiss
   enddo

   write(cluster_unit,"(A))") "================ Cluster List ================"
   write(cluster_unit,"(A))") "cluster, num region,  num point,  area total"
   do i=1, num_cluster
      !write(*,*) "cluster update",  i
      write(cluster_unit,"(A7,2(',',I11),',',F12.3)") &
            cluster_name(i), cluster(i)%num_region,cluster(i)%num_point,cluster(i)%area_total
   enddo

   ! set grid%cluster
   do n=1, num_region 
      do k=1, region(n)%num_point
         i=region(n)%point_list(k)%i
         j=region(n)%point_list(k)%j
         grid%cluster(i,j)=region(n)%cluster
      enddo
   enddo
   
   write(sort_unit,"(A)") "================ Cluster Sort ================" 
   do k=1, num_cluster
      !write(*,*) "cluster sort",  k
      write(sort_unit,"(A,2(I10,F12.3))") cluster_name(k), d(k), a(k), cluster(k)%num_point, cluster(k)%area_total
   enddo

   write(cluster_unit,"(A)") "================ Cluster Region List ================"
   write(cluster_unit,"(A)") "cluster, region#, region name, num point,  area total"
   d=0
   do k=1, num_cluster
      if(allocated (cluster(k)%region_list))then
         deallocate(cluster(k)%region_list)
      endif
      allocate( cluster(k)%region_list(cluster(k)%num_region))
      write(cluster_unit,"(A)") "-----------------------------------------------------"
      do n=1,num_region
         if(region(n)%cluster==k)then
            d(k)=d(k)+1
            cluster(k)%region_list(d(k))=n
         endif
      enddo
      !write(111,*) "cluster region list:", k, cluster(k)%num_region, d(k),ubound(cluster(k)%region_list), cluster(k)%num_point, cluster(k)%area_total
      do i=1, cluster(k)%num_region
         !write(111,*) i,cluster(k)%region_list(i)
         !write(111,*) region(cluster(k)%region_list(i))%num_point
      !write(*,*) "cluster update2",  k
         write(cluster_unit,"(A7,',',I8,',',A12,',',I10,',',F12.3)") cluster_name(k), i, region_name(cluster(k)%region_list(i)), &
                                  region(cluster(k)%region_list(i))%num_point, region(cluster(k)%region_list(i))%area_total
      enddo
      !write(111,*) (cluster(k)%region_list(i), region(cluster(k)%region_list(i))%num_point,',',i=1,cluster(k)%num_region)
   enddo 

!  count num of neighbor
!  write(*,*) ubound(cluster)
   cluster(:)%num_neighbor=0
   do i=1, num_region
      do j=1, num_region
         if(abs(neighbor_flag(i,j))<=1.and.num_effect(i,j)>=n0)then
            n=region(i)%cluster
!           write(*,*) "count neighbor", n,  cluster(n)%num_neighbor
            cluster(n)%num_neighbor=cluster(n)%num_neighbor+1
         endif
      enddo
   enddo

   ! allocate
   do n=1, num_cluster
      if(allocated( cluster(n)%neighbor_flag))then
         deallocate(cluster(n)%neighbor_flag)
         deallocate(cluster(n)%neighbor_list)
         deallocate(cluster(n)%num_border_point)
         deallocate(cluster(n)%num_border_effect)
         deallocate(cluster(n)%sorted_border_index)
      endif
      if(cluster(n)%num_neighbor>0)then
         allocate(cluster(n)%neighbor_flag      (cluster(n)%num_neighbor))
         allocate(cluster(n)%num_border_point   (cluster(n)%num_neighbor))
         allocate(cluster(n)%num_border_effect  (cluster(n)%num_neighbor))
         allocate(cluster(n)%neighbor_list      (cluster(n)%num_neighbor,2))
         allocate(cluster(n)%sorted_border_index(cluster(n)%num_neighbor))
      endif
   enddo

!   write(98,*) "cluster num_neighbor 1", (cluster(n)%num_region,cluster(n)%num_neighbor,cluster(n)%num_point,',',n=1,num_cluster)

   cluster(:)%num_neighbor=0
   do i=1, num_region
      do j=1, num_region
         if(abs(neighbor_flag(i,j))<=1.and.num_effect(i,j)>=n0)then
         n=region(i)%cluster
         cluster(n)%num_neighbor    =cluster(n)%num_neighbor+1
         cluster(n)%num_border_point (cluster(n)%num_neighbor  )=num_point (i,j)
         cluster(n)%num_border_effect(cluster(n)%num_neighbor  )=num_effect(i,j)
         cluster(n)%neighbor_list   (cluster(n)%num_neighbor,1)=i
         cluster(n)%neighbor_list   (cluster(n)%num_neighbor,2)=j
         cluster(n)%neighbor_flag   (cluster(n)%num_neighbor  )=neighbor_flag(i,j)
         endif
      enddo
   enddo
   write(cluster_unit,"(A)") "================ Cluster Information ==============="
   write(cluster_unit,"(A)") "   cluster    region  neighbor   cluster     cluster"
   write(cluster_unit,"(A)") "      name    number    number num point  area total"
   do n=1, num_cluster
      !write(*,*) "cluster info",  n
      write(cluster_unit,"(A10,3I10,F12.3)") cluster_name(n),cluster(n)%num_region,&
            cluster(n)%num_neighbor,cluster(n)%num_point, cluster(n)%area_total
   enddo
   write(cluster_unit,"(A)") "================ Cluster Neighbor ================"
   if(allocated(a))then
      deallocate(a)
      deallocate(d)
   endif
   do n=1, num_cluster
      !write(*,*) "cluster list",  n
      if(cluster(n)%num_neighbor<1) cycle
      write(cluster_unit,"(A)") "---- region neighbor list for cluster:"//cluster_name(n)//" ----"
      write(cluster_unit,"(A)") "  neighbor  neighbor border_pt border_pt  neighbor"
      write(cluster_unit,"(A)") "  region#1  region#2    number    effect      flag"
      allocate(a(cluster(n)%num_neighbor))
      allocate(d(cluster(n)%num_neighbor))
      do i=1, cluster(n)%num_neighbor
         write(cluster_unit,"(2A10,2I10,I10)") &
               region_name(cluster(n)%neighbor_list(i,1)), &
               region_name(cluster(n)%neighbor_list(i,2)), &
               cluster(n)%num_border_point(i),             &
               cluster(n)%num_border_effect(i),            &
               cluster(n)%neighbor_flag(i)

         a(i)=cluster(n)%num_border_effect(i)
         d(i)=i
      enddo
      call quick_sort(a,d,cluster(n)%num_neighbor,1,cluster(n)%num_neighbor)
      cluster(n)%sorted_border_index=d
      deallocate(a)
      deallocate(d)
   enddo
   !write(98,*) "cluster num_neighbor 2", cluster(:)%num_neighbor
    
!  if(.NOT.write_border_list) return 
   !allocate(border_list(num_region,num_region))
   !write(*,*) "creating border list..."
   do i=1, min(num_region,nb_num)
      do j=1, min(num_region,nb_num)
         border_list(i,j)%num_point=num_point(i,j)
         if(allocated (border_list(i,j)%i1))then
            deallocate(border_list(i,j)%i1)
            deallocate(border_list(i,j)%j1)
            deallocate(border_list(i,j)%i2)
            deallocate(border_list(i,j)%j2)
         endif
         allocate(border_list(i,j)%i1(border_list(i,j)%num_point))
         allocate(border_list(i,j)%j1(border_list(i,j)%num_point))
         allocate(border_list(i,j)%i2(border_list(i,j)%num_point))
         allocate(border_list(i,j)%j2(border_list(i,j)%num_point))
      enddo
   enddo
   num_point=0
   !write(border_unit,*) "+==================================="
   !write(border_unit,*) "+==== Border List 2 ============"
   do j=1, num_azim
      if(grid%start_i(j)<1) cycle
      do i=grid%start_i(j), grid%end_i(j)-1
         i1=i+1
         j1=j+1
         if(j1>num_azim) j1=1
         !if(abs(grid%vel(i,j))<=2) cycle
         !if(abs(grid%vel(i1,j))>2) then
         if(grid%region(i,j)/=imiss.and.grid%region(i1,j)/=imiss.and.(grid%region(i,j)/=grid%region(i1,j)))then
            i2=i1
            j2=j
            i1=i
            j1=j
            n1=grid%region(i1,j1)
            n2=grid%region(i2,j2)
            num_point(n1,n2)=num_point(n1,n2)+1
            num_point(n2,n1)=num_point(n2,n1)+1
            !if(n1<=nb_num.and.n2<=nb_num)then
            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
            !endif
            if(grid%region(i1,j1)>=1.and.grid%region(i1,j1)<=nb_num.and. &
               grid%region(i2,j2)>=1.and.grid%region(i2,j2)<=nb_num) then
               border_list(n1,n2)%i1(num_point(n1,n2))=i1
               border_list(n1,n2)%i2(num_point(n1,n2))=i2
               border_list(n1,n2)%j1(num_point(n1,n2))=j1
               border_list(n1,n2)%j2(num_point(n1,n2))=j2
               border_list(n2,n1)%i1(num_point(n2,n1))=i2
               border_list(n2,n1)%i2(num_point(n2,n1))=i1
               border_list(n2,n1)%j1(num_point(n2,n1))=j2
               border_list(n2,n1)%j2(num_point(n2,n1))=j1
            endif
         endif
         i1=i+1
         !endif
         j1=j+1
         if(j1>num_azim) j1=1
         !if(abs(grid%vel(i,j1))>2) then
         !if(grid%region(i,j1)>=1.and.grid%region(i,j1)<=30)then
         if(grid%region(i,j)/=imiss.and.grid%region(i,j1)/=imiss.and.(grid%region(i,j)/=grid%region(i,j1)))then
            i2=i
            j2=j1
            i1=i
            j1=j
            n1=grid%region(i1,j1)
            n2=grid%region(i2,j2)
            num_point(n1,n2)=num_point(n1,n2)+1
            num_point(n2,n1)=num_point(n2,n1)+1
            !if(n1<=nb_num.and.n2<=nb_num)then
            !   write(border_unit,"(A,1X,A,I10,A,1X,A)") region_name(n1),region_name(n2),num_point(n1,n2),point_name(i1,j1),point_name(i2,j2)
            !endif
            if(grid%region(i1,j1)>=1.and.grid%region(i1,j1)<=nb_num .and.&
               grid%region(i2,j2)>=1.and.grid%region(i2,j2)<=nb_num) then

            if(num_point(n1,n2)>border_list(n1,n2)%num_point)stop
!write(*,*) "border_list i1",n1,n2,num_point(n1,n2),ubound(border_list(n1,n2)%i1)
               border_list(n1,n2)%i1(num_point(n1,n2))=i1
               border_list(n1,n2)%i2(num_point(n1,n2))=i2
               border_list(n1,n2)%j1(num_point(n1,n2))=j1
               border_list(n1,n2)%j2(num_point(n1,n2))=j2
               border_list(n2,n1)%i1(num_point(n2,n1))=i2
               border_list(n2,n1)%i2(num_point(n2,n1))=i1
               border_list(n2,n1)%j1(num_point(n2,n1))=j2
               border_list(n2,n1)%j2(num_point(n2,n1))=j1
            endif
         endif
         !endif
      enddo
   enddo

   !write(*,*) "wrting border list..."
   write(border_unit,*) "======== Region Border: c1=", c1, " c2=",c2, "========"
   do i=1, min(num_region,nb_num)
      do j=i+1, min(num_region,nb_num)
         if(border_list(i,j)%num_point<1)cycle
         write(border_unit,"(2(A,1X),3(A,I5))") region_name(i), region_name(j) , &
                                                " Border_Pt",num_point    (i,j), &
                                                " Effect_Pt",num_effect   (i,j), &
                                                " flag"     ,neighbor_flag(i,j)

         write(border_unit,"(2(A10,1X),6A7,2A4,2A3,A5,A7,3A5)")"Point#1","Point#2","v1","v2","z1","z2","w1","w2","f1","f2","l1","l2","flag","v1-v2","<c1",">c2","<-c2"
         do n=1, border_list(i,j)%num_point
            i1=border_list(i,j)%i1(n)
            j1=border_list(i,j)%j1(n)
            i2=border_list(i,j)%i2(n)
            j2=border_list(i,j)%j2(n)
            n1=grid%region(i1,j1)
            n2=grid%region(i2,j2)
            call get_neighbor_point_flag(grid,i1,j1,i2,j2,flag,region(n1)%area_total,region(n2)%area_total)
            write(border_unit,"(2(A10,1X),6F7.1,2I4,2I3,I5,F7.1,3L5)") &
                  point_name(border_list(i,j)%i1(n),border_list(i,j)%j1(n)), &
                  point_name(border_list(i,j)%i2(n),border_list(i,j)%j2(n)), &
                  grid%vel(i1,j1),grid%vel(i2,j2),&
                  grid%ref(i1,j1),grid%ref(i2,j2),&
                  grid%spw(i1,j1),grid%spw(i2,j2),&
                  grid%flag(i1,j1),grid%flag(i2,j2),&
                  grid%level(i1,j1),grid%level(i2,j2),flag,&
                  grid%vel(i1,j1)-grid%vel(i2,j2),&
                  abs(grid%vel(i1,j1)-grid%vel(i2,j2))<=c1,&
                  (grid%vel(i1,j1)-grid%vel(i2,j2))>= c2, &
                  (grid%vel(i1,j1)-grid%vel(i2,j2))<=-c2
            
         enddo
      enddo
   enddo

   !deallocate(border_list)
   deallocate(num_point    )
   deallocate(num_effect   )
   deallocate(neighbor_flag)

end subroutine
end module
