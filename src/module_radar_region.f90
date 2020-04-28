module radar_region
use grd
implicit none

private

type t_point_list
   integer :: num_point
   integer, dimension(:), allocatable :: range_index, azim_index !(num_point)
end type

type t_region
   integer :: num_range, num_azim, num_region, num_valid
   integer, dimension(:), allocatable :: start_range_index, end_range_index 
   integer, dimension(:,:), allocatable :: region, level
   !integer, dimension(:),   allocatable :: num_point
   real   , dimension(:,:), allocatable :: dat
   type(t_point_list), dimension(:), allocatable :: point_list !(num_region)
   real :: min_dat, cint
end type

real, private :: miss = -888888.


public :: remove_noise, dealias_region

contains

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
   subroutine remove_noise(num_range,num_azim,ref,miss_var,ref_out)
   implicit none

   integer,                             intent(in)  :: num_range, num_azim
   real, dimension(num_range,num_azim), intent(in)  :: ref
   real, dimension(num_range,num_azim), intent(out) :: ref_out 
   real,                                intent(in)  :: miss_var

   integer :: i, j, k, n, m, np, nsum, n1, n2

   type(t_region)                       :: reg

   integer, dimension(num_range,num_azim) :: cluster_index

   real,    dimension(:,:), allocatable :: reg_diff
   integer, dimension(:,:), allocatable :: num_diff
   integer, dimension(:)  , allocatable :: region_cluster_index, num_region_in_cluster, num_point_in_region, num_point_in_cluster
   integer :: num_cluster, num_region_dealiased, num_point_dealiased

   real    :: min_dat, max_dat
   logical :: first_in_cluster
   integer :: cint

   call set_region_miss(miss_var)
   !call write_grd_real('ref.grd',num_range,num_azim,ref,miss)
!   write(96) vel
!stop


   ref_out=ref
   where(ref_out<0)
      ref_out=miss
   endwhere
   cint=10
   call create_region(num_range,num_azim,ref_out,cint,reg)

   call write_grd_int('rlev.grd',num_range,num_azim,reg%level,int(miss))
!stop

!   call clean_no_neighbor_region (reg)
!   call write_grd_int('reg0.grd',num_range,num_azim,reg%region,int(miss))

   call merge_region (reg)

   !call clean_region (reg,100)
   !call clean_no_neighbor_region (reg)

   call write_grd_int('rreg.grd',num_range,num_azim,reg%region,int(miss))

   call get_neighbor_diff(reg,reg_diff,num_diff)
         
   allocate(region_cluster_index(reg%num_region))
   allocate(num_region_in_cluster(reg%num_region))
   region_cluster_index=0
   num_region_in_cluster=0
   num_cluster=0
   do while(.true.) !every region have cluster
      do i=1, reg%num_region
         if(region_cluster_index(i)==0) then
            num_cluster=num_cluster+1
            region_cluster_index(i)=num_cluster
            num_region_in_cluster(num_cluster)=num_region_in_cluster(num_cluster)+1
         endif
         do j=i+1, reg%num_region 
            if(num_diff(i,j)>0)then
               region_cluster_index(j)=region_cluster_index(i)
               num_region_in_cluster(num_cluster)=num_region_in_cluster(num_cluster)+1
            endif
         enddo
      enddo
      if(all(region_cluster_index/=0)) exit
   enddo
   write(*,*) "num cluster:", num_cluster, num_region_in_cluster(1:num_cluster)

   allocate(num_point_in_cluster(num_cluster))
   allocate(num_point_in_region(reg%num_region))
   write(*,*) ubound(num_point_in_cluster)
   cluster_index=miss
   ref_out      =miss
   num_point_in_cluster=0
   do m=1, num_cluster
      do n=1, reg%num_region
         if(region_cluster_index(n)==m)then
            write(99,*) m,n,reg%point_list(n)%num_point
            num_point_in_cluster(m)=num_point_in_cluster(m)+reg%point_list(n)%num_point
            do k=1,  reg%point_list(n)%num_point
               i=reg%point_list(n)%range_index(k)
               j=reg%point_list(n)%azim_index (k)
               cluster_index(i,j)=region_cluster_index(n)
               ref_out      (i,j)=ref(i,j)
            enddo
         endif
      enddo
      if(num_point_in_cluster(m)<100)then
         do n=1, reg%num_region
            do k=1,  reg%point_list(n)%num_point
               if(region_cluster_index(n)/=m) cycle
               i=reg%point_list(n)%range_index(k)
               j=reg%point_list(n)%azim_index (k)
               cluster_index(i,j)=miss
               ref_out      (i,j)=miss
            enddo
         enddo
         num_point_in_cluster(m)=0
         exit
      endif
      write(*,*) "num point in cluster", m,"=",num_point_in_cluster(m)
      if(sum(num_point_in_cluster(1:m))>0.8*reg%num_valid)then
         exit
      endif
   enddo
   call write_grd_int('cluster.grd',num_range,num_azim,cluster_index,int(miss))
      
   write(*,*) "Remove ", reg%num_valid-sum(num_point_in_cluster), &
              " ponits out of ",reg%num_valid,"(",100*(1.-(sum(num_point_in_cluster))/real(reg%num_valid)),"%)."

   !write(121) vel
   call write_grd_real('ref_new.grd',num_range,num_azim,ref_out,miss)
   !write(122) nyq_num 
!   call write_grd_int('n.grd',num_range,num_azim,nyq_num,int(miss))
   
   !stop
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine dealias_region(num_range,num_azim,vel,ref,vmax,miss_var,nyq_num,vel_out)
   implicit none

   integer, intent(in) :: num_range, num_azim
   real, dimension(num_range,num_azim),   intent(in)  :: vel, ref
   real,                     intent(in)  :: vmax
   real,    dimension(num_range,num_azim),intent(out) :: vel_out 
   integer, dimension(num_range,num_azim),intent(out) :: nyq_num 
   real,                     intent(in)  :: miss_var

   integer :: i, j, k, l, n, m, np, nsum, nyq, n1, n2, i1, j1, l1, npoint

   type(t_region)                    :: reg
   integer, dimension(:), allocatable :: nyg_region

   real,    dimension(:,:), allocatable :: reg_diff
   integer, dimension(:,:), allocatable :: num_diff
   integer, dimension(:)  , allocatable :: region_cluster_index, num_region_in_cluster, num_point_in_region, num_point_in_cluster, sum_nyq_in_cluster
   integer :: num_cluster, num_region_dealiased, num_point_dealiased

   logical :: first_in_cluster
   integer :: cint

!   call write_grd_real('vel.grd',num_range,num_azim,vel,value_invalid)
!   write(96) vel
!stop

   call set_region_miss(miss_var)

   vel_out=miss
   cint=max(min(4,int(vmax/4+0.5)),2)
   call create_region(num_range,num_azim,vel,cint,reg)

!   write(97) reg%level
!   write(*,*) ubound(reg%level)
!   call write_grd_int('lev.grd',num_range,num_azim,reg%level,int(miss))
!stop

!   write(98) reg%region
!   call clean_no_neighbor_region (reg)
!   call write_grd_int('reg0.grd',num_range,num_azim,reg%region,int(miss))

   call merge_region (reg)

   !call clean_region (reg,100)
   !call clean_no_neighbor_region (reg)
   allocate(nyg_region(reg%num_region))
   nyg_region=miss
!   write(99) reg%region
!   call write_grd_int('reg1.grd',num_range,num_azim,reg%region,int(miss))
!stop

   call get_neighbor_diff(reg,reg_diff,num_diff)
         
   allocate(region_cluster_index(reg%num_region))
   allocate(num_region_in_cluster(reg%num_region))
   region_cluster_index=0
   num_region_in_cluster=0
   num_cluster=0
   do while(.true.) !every region have cluster index
      do i=1, reg%num_region
         if(region_cluster_index(i)==0) then
            num_cluster=num_cluster+1
            region_cluster_index(i)=num_cluster
            num_region_in_cluster(num_cluster)=num_region_in_cluster(num_cluster)+1
         endif
         do j=i+1, reg%num_region 
            if(num_diff(i,j)>0)then
               region_cluster_index(j)=region_cluster_index(i)
               num_region_in_cluster(num_cluster)=num_region_in_cluster(num_cluster)+1
            endif
         enddo
      enddo
      if(all(region_cluster_index/=0)) exit
   enddo
   write(*,*) "region cluster:", num_cluster, num_region_in_cluster(1:num_cluster)

   allocate(num_point_in_region (reg%num_region))
   allocate(num_point_in_cluster(num_cluster   ))
   allocate(sum_nyq_in_cluster  (num_cluster   ))
   sum_nyq_in_cluster=0
   num_point_in_cluster=0
   num_point_in_region(1:reg%num_region)=reg%point_list(1:reg%num_region)%num_point
   nyg_region=miss
   do n=1, num_cluster
      first_in_cluster=.true.
      num_point_in_cluster(n)=0
      sum_nyq_in_cluster(n)=0
      do i=1, reg%num_region
         if(region_cluster_index(i)==n)then
            !if(num_region_in_cluster(n)==1)then
            !   !write(*,*) "miss nyquist",i, n, num_point_in_region(i)
            !   nyg_region(i)=miss
            !else
              if(first_in_cluster)then
                 nyg_region(i)=0
                 !write(*,*) "cluster init:", i, 0
                 first_in_cluster=.false.
              else
                 do j=1,i-1
                    if(nyg_region(j)/=miss.and.num_diff(i,j)>0.and.region_cluster_index(j)==n)then
                       if(nyg_region(i)==miss)then
                          if(abs(reg_diff(i,j)/vmax)<0.5.or.abs(reg_diff(i,j)/vmax)>1.5)then
                             nyg_region(i)=-floor(reg_diff(i,j)/vmax/2+0.5)+nyg_region(j)
                             !write(*,*) "nyq set:", i, j,  nyg_region(j), reg_diff(i,j), -floor(reg_diff(i,j)/vmax/2+0.5), &
                             !            nyg_region(i),num_point_in_region(i),num_point_in_region(j)
                          endif
                       else ! check
                          !if(nyg_region(i)/=-floor(reg_diff(i,j)/vmax/2+0.5)+nyg_region(j))then
                          if(abs(2*vmax*nyg_region(i)+reg_diff(i,j)-2*vmax*nyg_region(j))>0.5*vmax)then
                             write(*,*) "check region nyq:", i, j, reg_diff(i,j), nyg_region(i)*2*vmax, nyg_region(j)*2*vmax, &
                                        abs(2*vmax*nyg_region(i)+reg_diff(i,j)-2*vmax*nyg_region(j))
                             nyg_region(i)=miss
                             !stop
                          endif
                       endif
                    endif
                 enddo
              endif
              if(nyg_region(i)/=miss)then
                 num_point_in_cluster(n)=num_point_in_cluster(n)+num_point_in_region(i) 
                 sum_nyq_in_cluster(n)=sum_nyq_in_cluster(n)+nyg_region(i)*num_point_in_region(i)
              endif
            !endif
         endif
      enddo
      do i=-3,3
         if(abs(sum_nyq_in_cluster(n)+i*num_point_in_cluster(n))<abs(sum_nyq_in_cluster(n)))then
            sum_nyq_in_cluster(n)=sum_nyq_in_cluster(n)+i*num_point_in_cluster(n)
            do j=1, reg%num_region
               if(region_cluster_index(j)==n.and.nyg_region(j)/=miss)then
                  nyg_region(j)=nyg_region(j)+i
               endif
            enddo
         endif
      enddo
      if(sum(num_point_in_cluster(1:n))>0.8*reg%num_valid)then
         exit
      endif
      if(num_point_in_cluster(n)<100)then
         write(*,*) 'num in cluster <100', n
         do j=1, reg%num_region
            if(region_cluster_index(j)==n)then
               write(*,*) 'set nyq =miss in region', j
               nyg_region(j)=miss
            endif
         enddo
      endif
   enddo


   !reference check
   num_point_dealiased=0
   num_region_dealiased=0 
   do j=1, reg%num_azim
      do i=reg%start_range_index(j), reg%end_range_index(j)
         n=reg%region(i,j)
         if(vel_out(i,j)==miss.or.nyg_region(n)==miss) cycle
         vel_out(i,j)=vel(i,j)+nyg_region(n)*2*vmax
         if(abs(ref(i,j)-vel_out(i,j))>2.*vmax.and.abs(ref(i,j)-vel_out(i,j))<3.*vmax)then
            write(*,*) "reference check", i, j, n, vel_out(i,j), ref(i,j), vmax
            k=(ref(i,j)-vel_out(i,j))/2/vmax
            m=region_cluster_index(n)
            npoint=0
            do l=1, reg%num_region
               if(region_cluster_index(l)==m)then
                 do l1=1, reg%point_list(l)%num_point
                    i1=reg%point_list(l)%range_index(l1)
                    j1=reg%point_list(l)%azim_index (l1)
                    n1=reg%region(i1,j1)
                    if(vel_out(i1,j1)==miss.or.nyg_region(n1)==miss) cycle
                    vel_out(i1,j1)=vel(i1,j1)+(nyg_region(n1)+k)*2*vmax
                    if(abs(ref(i1,j1)-vel_out(i1,j1))>2.*vmax)then
                       npoint=npoint+1
                       exit
                    endif
                 enddo
               endif
               if(npoint>0) exit
            enddo
            if(npoint>0)then
               do l=1, reg%num_region
                  nyg_region(l)=miss
               enddo
            else
               do l=1, reg%num_region
                  nyg_region(l)=nyg_region(l)+k
               enddo
            endif
         endif
      enddo
   enddo

   nyq_num=miss
   num_point_dealiased=0
   num_region_dealiased=0 
   do n=1, reg%num_region
      if(nyg_region(n)/=miss)then
         num_region_dealiased=num_region_dealiased+1
         num_point_dealiased=num_point_dealiased+num_point_in_region(n)
         !write(*,*) "Region ", n, " nyquist number=", nyg_region(n), num_point_in_region(n) !, mean_obs(n), mean_bak(n), mean_ref(n)
         do m=1,  reg%point_list(n)%num_point
            i=reg%point_list(n)%range_index(m)
            j=reg%point_list(n)%azim_index(m)
            vel_out(i,j)=vel(i,j)+2*nyg_region(n)*vmax
            nyq_num(i,j)=nyg_region(n)
         enddo
      endif
   enddo
   write(*,*) "Dealiased ", num_region_dealiased," region, ", num_point_dealiased, &
              " ponits out of ",reg%num_valid,"(",100.*num_point_dealiased/reg%num_valid,"%)."

   !write(121) vel
!   call write_grd_real('vel_new.grd',num_range,num_azim,vel,value_invalid)
   !write(122) nyq_num 
!   call write_grd_int('n.grd',num_range,num_azim,nyq_num,int(miss))
   
   !stop
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_neighbor_diff(reg,reg_diff,num_diff)
   implicit none
   type(t_region), intent(in) :: reg
   integer, dimension(:,:), allocatable, intent(out) :: num_diff
   real   , dimension(:,:), allocatable, intent(out) :: reg_diff
  
   integer :: i, j, num_range, num_azim, n1, n2, i1, j1, i2, j2

!   write(*,*) "in get_neighbor"
   num_range=reg%num_range
   num_azim =reg%num_azim

   if(allocated(reg_diff))then
      deallocate(reg_diff)
      deallocate(num_diff)
   endif
!   write(*,*) "allocate get_neighbor"
   allocate(reg_diff(reg%num_region,reg%num_region))
   allocate(num_diff(reg%num_region,reg%num_region))

   num_diff=0
   reg_diff=0
   do j=1, num_azim
         if(reg%start_range_index(j)<1) cycle
      do i=reg%start_range_index(j), reg%end_range_index(j)-1
         i1=i+1
         j1=j+1
         if(j1>num_azim) j1=1
         if(reg%region(i,j)/=miss.and.reg%region(i1,j)/=miss.and.(reg%region(i,j)/=reg%region(i1,j)))then
            i2=i1
            j2=j
            i1=i
            j1=j
            n1=reg%region(i1,j1)
            n2=reg%region(i2,j2)
            num_diff(n1,n2)=num_diff(n1,n2)+1
            reg_diff(n1,n2)=reg_diff(n1,n2)+reg%dat(i1,j1)-reg%dat(i2,j2)
            num_diff(n2,n1)=num_diff(n2,n1)+1
            reg_diff(n2,n1)=reg_diff(n2,n1)+reg%dat(i2,j2)-reg%dat(i1,j1)
         endif
         i1=i+1
         j1=j+1
         if(j1>num_azim) j1=1
         if(reg%region(i,j)/=miss.and.reg%region(i,j1)/=miss.and.(reg%region(i,j)/=reg%region(i,j1)))then
            i2=i
            j2=j1
            i1=i
            j1=j
            n1=reg%region(i1,j1)
            n2=reg%region(i2,j2)
            num_diff(n1,n2)=num_diff(n1,n2)+1
            reg_diff(n1,n2)=reg_diff(n1,n2)+reg%dat(i1,j1)-reg%dat(i2,j2)
            num_diff(n2,n1)=num_diff(n2,n1)+1
            reg_diff(n2,n1)=reg_diff(n2,n1)+reg%dat(i2,j2)-reg%dat(i1,j1)
         endif
      enddo
   enddo
   where(num_diff/=0)
      reg_diff=reg_diff/num_diff
   endwhere
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine create_region(num_range, num_azim, dat, cint, reg)
   implicit none

   integer,                                intent(in) :: num_range, num_azim, cint
   real   , dimension(num_range,num_azim), intent(in) :: dat
   type(t_region),                         intent(out):: reg

   integer :: i, j, k, np, n
   real    :: min_dat, max_dat
   integer, dimension(8) :: l8, r8
   
   integer :: num_valid, num_no_region_point, nlevel, old_no_region_point
   logical :: new_region, change_region

   reg%num_range=num_range
   reg%num_azim=num_azim
   if(allocated (reg%level ))then
      deallocate(reg%level )
      deallocate(reg%region)
      deallocate(reg%dat)
      deallocate(reg%start_range_index)
      deallocate(reg%end_range_index)
   endif
   allocate(reg%level   (num_range,num_azim))
   allocate(reg%region  (num_range,num_azim))
   allocate(reg%dat     (num_range,num_azim))
   allocate(reg%start_range_index (num_azim))
   allocate(reg%end_range_index   (num_azim))

   num_valid=0
   reg%start_range_index=0
   reg%end_range_index  =0
   min_dat=miss
   max_dat=miss
   do j=1, num_azim
      do i=1, num_range
         if(abs(dat(i,j))<200.)then !/=miss)then
            if(reg%start_range_index(j)==0) reg%start_range_index(j)=i
            reg%end_range_index(j)=i
            num_valid=num_valid+1
            if(min_dat==miss)then
               min_dat=dat(i,j)
            else
               min_dat=min(min_dat,dat(i,j))
            endif
            if(max_dat==miss)then
               max_dat=dat(i,j)
            else
               max_dat=max(max_dat,dat(i,j))
            endif
         endif
      enddo
   enddo
   write(*,*) "num valid:", num_valid
   reg%num_valid=num_valid

   min_dat=floor(min_dat/cint)*cint
   reg%min_dat=min_dat
   reg%cint=cint
   write(*,*) "min_dat,max_dat:", min_dat, max_dat, cint

   reg%dat=dat
     
   reg%level=int((dat-min_dat)/cint)+1
   where(reg%level<=0)
      reg%level=miss
   endwhere
   write(*,*) "generate level"

!   write(102) reg%level

   nlevel=maxval(reg%level)
   reg%region=miss
   reg%num_region=0
   do n=1, nlevel
      old_no_region_point=0
      write(*,*) "loop level=",n, nlevel
      do while(.true.) ! all dat must have region
         num_no_region_point=0
         new_region=.true.
         do j=1, num_azim
            if(reg%start_range_index(j)<1) cycle
            do i=reg%start_range_index(j), reg%end_range_index(j)
               if(reg%level(i,j)/=n) cycle 
               change_region=.false.
               if(reg%level(i,j)/=miss.and.reg%region(i,j)==miss)then
                   
                  if(new_region)then
                     !write(111,*) "new init", i,j,r8
                     reg%num_region =reg%num_region+1
                     reg%region(i,j)=reg%num_region
                     new_region   =.false.
                     change_region=.true.
                  else
                     call get_n8_int(num_range,num_azim,reg%region,i,j,r8)
                     call get_n8_int(num_range,num_azim,reg%level ,i,j,l8)
                     do k=1, 8
                        if(l8(k)==reg%level(i,j).and.r8(k)/=miss)then
                           !write(111,*) "neighbor", i,j,r8(k),l8(k),reg%level(i,j)
                           reg%region(i,j)=r8(k)
                           change_region=.true.
                           exit
                        endif
                     enddo
                  endif
               endif

               ! expand region
               if(change_region.and.reg%region(i,j)/=miss)then
                  call expand_region(i,j)
               endif
               
               if(reg%level(i,j)/=miss.and.reg%region(i,j)==miss)then
                  num_no_region_point=num_no_region_point+1
               endif
            enddo
         enddo
         if(num_no_region_point==0) exit
         if(num_no_region_point==old_no_region_point)then
            new_region=.true.
         endif
         old_no_region_point=num_no_region_point
      enddo
   enddo

!   write(*,*) reg%num_region
!   write(103) reg%region
!stop
   call update_region(reg)
   write(*,*) "create region number:", reg%num_region
   contains
      subroutine expand_region(i,j)
      implicit none
         integer, intent(in) :: i,j 
         integer :: i1, j1, i_start, i_end, j_start, j_end
         integer :: i2, j2
         logical :: loop_end

         integer :: idiff=0, l1, l2, d
!         write(*,*) i,j
         if(i>1)then
            do i1=i-1,reg%start_range_index(j),-1
               j1=j

               i2=i1+1
               j2=j

               l1=reg%level(i1,j1)
               l2=reg%level(i2,j2)
               d=l1-l2

               loop_end=.false.  
               if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
                  !write(*,*) l1,l2,l1-l2 
                  reg%region(i1,j1)=reg%region(i2,j2)
               else
                  loop_end=.true.  
               endif
               if(loop_end) exit

               j_start=1
               j_end  =num_azim
               if(j>1)then 
                  do j1=j-1, 1,-1

                     i2=i1
                     j2=j1+1

                     l1=reg%level(i1,j1)
                     l2=reg%level(i2,j2)
                     d=l1-l2

                     loop_end=.false.  
                     if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
                        !write(*,*) l1,l2,l1-l2
                        reg%region(i1,j1)=reg%region(i2,j2)
                     else
                        loop_end=.true.  
                        j_start=j2
                     endif
                     if(loop_end) exit

                  enddo
               endif

               do j1=j+1, num_azim

                  i2=i1
                  j2=j1-1

                  l1=reg%level(i1,j1)
                  l2=reg%level(i2,j2)
                  d=l1-l2

                  loop_end=.false.  
                  if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
                     !write(*,*) l1,l2,l1-l2
                     reg%region(i1,j1)=reg%region(i2,j2)
                  else
                     loop_end=.true.  
                     j_end=j2
                  endif
                  if(loop_end) exit

               enddo
               !write(*,*) "expand:",i1,j_start, j_end
            enddo
         endif
         do i1=i+1,reg%end_range_index(j)
            j1=j

            i2=i1-1
            j2=j

            l1=reg%level(i1,j1)
            l2=reg%level(i2,j2)
            d=l1-l2

            loop_end=.false.  
            if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
               !write(*,*) l1,l2,l1-l2
               reg%region(i1,j)=reg%region(i2,j2)
            else
               loop_end=.true.  
               i_start=i2
            endif
            if(loop_end) exit

            j_start=1
            j_end  =num_azim
            if(j>1)then 
               do j1=j-1, 1,-1
                  i2=i1
                  j2=j1+1

                  l1=reg%level(i1,j1)
                  l2=reg%level(i2,j2)
                  d=l1-l2

                  loop_end=.false.  
                  if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
                     !write(*,*) l1,l2,l1-l2
                     reg%region(i1,j1)=reg%region(i2,j2)
                  else
                     loop_end=.true.  
                     j_start=j2
                  endif
                  if(loop_end) exit
               enddo
            endif

            do j1=j+1, num_azim
               i2=i1
               j2=j1-1

               l1=reg%level(i1,j1)
               l2=reg%level(i2,j2)
               d=l1-l2

               loop_end=.false.  
               if(abs(l1-l2)<=idiff.and.reg%region(i2,j2)/=miss)then
                  !write(*,*) l1,l2,l1-l2
                  reg%region(i1,j1)=reg%region(i2,j2)
               else
                  loop_end=.true.  
                  j_end=j2
               endif
               if(loop_end) exit
            enddo
            !write(*,*) "expand:",i1,j_start, j_end
         enddo
      end subroutine
   end subroutine
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine clean_no_neighbor_region(reg)
   implicit none
   type(t_region),         intent(inout) :: reg

   integer :: k 

   real, dimension(:,:), allocatable :: reg_diff
   integer, dimension(:,:), allocatable :: num_diff

   call get_neighbor_diff(reg,reg_diff,num_diff)
   do k=1, reg%num_region
      if(sum(num_diff(k,:))==0.and.reg%point_list(k)%num_point<100)then
!         write(*,*) "delete no neighbor region",k,reg%point_list(k)%num_point
         call set_region(reg,k,int(miss))
         reg%point_list(k)%num_point=0
      endif
   enddo
   call update_region(reg)
   deallocate(reg_diff)
   deallocate(num_diff)
   write(*,*) "region number after clean no neighbor region:", reg%num_region
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine clean_region(reg,min_np)
   implicit none
   type(t_region),         intent(inout) :: reg
   integer    ,             intent(in) :: min_np

   integer :: k 

   do k=1, reg%num_region
      if(reg%point_list(k)%num_point<min_np)then
         call set_region(reg,k,int(miss))
         reg%point_list(k)%num_point=0
      endif
   enddo
   call update_region(reg)
   write(*,*) "region number after clean:", reg%num_region
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copy_region(reg,reg_out)
   implicit none
   type(t_region),         intent(in) :: reg
   type(t_region),         intent(out) :: reg_out

   integer :: k, np, i, j

   if(allocated(reg_out%region))then
      deallocate(reg_out%region)
      deallocate(reg_out%level)
      do i=1, reg%num_region
         deallocate(reg_out%point_list(i)%range_index)
         deallocate(reg_out%point_list(i)%azim_index)
      enddo
      deallocate(reg_out%point_list       )
      deallocate(reg_out%dat              )
      deallocate(reg_out%start_range_index)
      deallocate(reg_out%end_range_index  )
   endif
   reg_out%num_range  =reg%num_range      
   reg_out%num_azim   =reg%num_azim      
   reg_out%num_region =reg%num_region 
   reg_out%min_dat    =reg%min_dat    
   reg_out%cint       =reg%cint    

   write(*,*) "copy region", reg_out%num_range, reg_out%num_azim, reg_out%num_region

   allocate(reg_out%region(reg_out%num_range,reg_out%num_azim))
   allocate(reg_out%level (reg_out%num_range,reg_out%num_azim))
   allocate(reg_out%dat   (reg_out%num_range,reg_out%num_azim))
   allocate(reg_out%point_list       (reg_out%num_region))
   allocate(reg_out%start_range_index(reg_out%num_azim))
   allocate(reg_out%end_range_index  (reg_out%num_azim))

   do i=1, reg%num_region
      allocate(reg_out%point_list(i)%range_index(reg%point_list(i)%num_point))
      allocate(reg_out%point_list(i)%azim_index (reg%point_list(i)%num_point))
   enddo
   
   reg_out%start_range_index=reg%start_range_index
   reg_out%end_range_index  =reg%end_range_index
   do i=1, reg%num_range
      do j=1, reg%num_azim
         reg_out%region(i,j)=reg%region(i,j)
         reg_out%level (i,j)=reg%level (i,j)
         reg_out%dat   (i,j)=reg%dat   (i,j)
      enddo
   enddo
   do i=1, reg_out%num_region
      reg_out%point_list(i)%num_point=reg%point_list(i)%num_point
      do k=1, reg_out%point_list(i)%num_point
         reg_out%point_list(i)%range_index(k)= reg%point_list(i)%range_index(k)
         reg_out%point_list(i)%azim_index(k) = reg%point_list(i)%azim_index(k)
      enddo
   enddo
         
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_region(reg,old,new)
   implicit none
   type(t_region),         intent(inout) :: reg
   integer,                intent(in)    :: old, new

   integer :: i, j

   do j=1, reg%num_azim
      if(reg%start_range_index(j)<1) cycle
      do i=reg%start_range_index(j), reg%end_range_index(j)
         if(reg%region(i,j)==old)then
            reg%region(i,j)=new
         endif
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine merge_region(reg)
   implicit none
   type(t_region),         intent(inout) :: reg

   integer :: i, j, k, n, m
   !integer, dimension(8) :: l8, r8, i8, j8
   
   real, dimension(:,:), allocatable :: reg_diff
   integer, dimension(:,:), allocatable :: num_diff

!   write(*,*) "in merge"
   call get_neighbor_diff(reg,reg_diff,num_diff)
!   write(*,*) "get diff"

   do i=1, reg%num_region
      if(reg%point_list(i)%num_point<1) cycle
      do j=i+1, reg%num_region 
         if(num_diff(i,j)>0.and.abs(reg_diff(i,j))<=abs(reg%cint))then
            call set_region(reg,j,i)
            reg%point_list(j)%num_point=0
         endif
      enddo
   enddo

!   write(*,*) "update_region"
   call update_region(reg)
   deallocate(reg_diff)
   deallocate(num_diff)
   write(*,*) "region number after merge:", reg%num_region
   end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine update_region(reg)
   use libradar, only: quick_sort
   implicit none
   type(t_region),         intent(inout):: reg

   integer :: i, j, k, np, m, ierr

   real   , dimension(:), allocatable :: a
   integer, dimension(:), allocatable :: d
   integer :: n, s, e

   if(allocated(reg%point_list))then
      do i=1, reg%num_region
         deallocate(reg%point_list(i)%range_index)
         deallocate(reg%point_list(i)%azim_index )
      enddo
      deallocate(reg%point_list)
   endif
!   write(*,*) "allocate point_list", reg%num_region
   allocate(reg%point_list(reg%num_region))
   reg%point_list(:)%num_point=0
   n=reg%num_region
!   write(*,*) "allocate a d"
   allocate(a(n)) 
   allocate(d(n)) 

   m=0
   do k=1, n 
      d(k)=k
      np=0
      do j=1, reg%num_azim
         if(reg%start_range_index(j)<1) cycle
         do i=reg%start_range_index(j), reg%end_range_index(j)
            if(reg%region(i,j)==k)then
               np=np+1
            endif
         enddo
      enddo
      a(k)=np
      reg%point_list(k)%num_point=np
   enddo

   a=-a
   s=1
   e=n

   call quick_sort(a,d,n,s,e)
   a=-a
   do k=1, n 
      if(a(k)==0)then
         reg%num_region=k-1
         exit
      endif
      write(201,*) k, a(k)
   enddo

   do k=1, reg%num_region
      call set_region(reg,d(k),-k)
   enddo

   where(reg%region/=miss)
      reg%region=-reg%region
   endwhere
   

   do i=1, reg%num_region
      if(allocated (reg%point_list(i)%range_index))then
         deallocate(reg%point_list(i)%range_index)
         deallocate(reg%point_list(i)%azim_index)
      endif
      np=abs(a(i))
      allocate(reg%point_list(i)%range_index(np))
      allocate(reg%point_list(i)%azim_index (np))
   enddo
!
   do k=1, reg%num_region
      np=0
      do j=1, reg%num_azim
         if(reg%start_range_index(j)<1) cycle
         do i=reg%start_range_index(j), reg%end_range_index(j)
            if(reg%region(i,j)==k)then
               np=np+1
               reg%point_list(k)%range_index(np)=i
               reg%point_list(k)%azim_index (np)=j
            endif
         enddo
      enddo
      reg%point_list(k)%num_point=np
   enddo

   deallocate(a)
   deallocate(d)
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_n8_int(num_range,num_azim,dat,i,j,n8,i8,j8)
   implicit none
   integer,                                intent(in) :: num_range, num_azim
   integer, dimension(num_range,num_azim), intent(in) :: dat
   integer,                                intent(in) :: i, j
   integer, dimension(8),                  intent(out):: n8
   integer, dimension(8), optional,        intent(out):: i8, j8

   integer :: i1,i2,j1,j2
   
   i1=i-1
   i2=i+1
   j1=j-1
   j2=j+1 
   if(i1<1 ) i1=1
   if(i2>num_range) i2=num_range
   if(j1<1 ) j1=num_azim
   if(j2>num_azim) j2=1
   n8(1)=dat(i1,j1)
   n8(2)=dat(i1,j )
   n8(3)=dat(i1,j2)
   n8(4)=dat(i ,j1)
   n8(5)=dat(i ,j2)
   n8(6)=dat(i2,j1)
   n8(7)=dat(i2,j )
   n8(8)=dat(i2,j2)
   if(present(i8))then
      i8(1)=i1
      i8(2)=i1
      i8(3)=i1
      i8(4)=i
      i8(5)=i
      i8(6)=i2
      i8(7)=i2
      i8(8)=i2
   endif
   if(present(j8))then
      j8(1)=j1
      j8(2)=j
      j8(3)=j2
      j8(4)=j1
      j8(5)=j2
      j8(6)=j1
      j8(7)=j
      j8(8)=j2
   endif
   if(i==i1) n8(1:3)=miss
   if(i==i2) n8(6:8)=miss
   end subroutine
end module
