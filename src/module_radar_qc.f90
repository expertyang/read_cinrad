module radar_qc
use radar_data
use libradar
use config
use radar_region
implicit none

   real    :: spwthrrat =0.67
   integer :: medfilt   =1
   real    :: winszrad=1000.,winszazim=10.,gcvrlim=1.6
   integer :: gclopt=1, bkgopt=0, shropt=1, rfropt=1,iordunf=2
   real    :: radarx=0., radary=0., radalt=0., dsort=0.5, velmedl=15., veldazl=30.
   real    :: refcheck=-200.,  velcheck=-200.

   real, parameter :: dspmiss = -991., anrflag=-800., spwflag=-931.

   integer                         :: nzsnd
   real, dimension(:), allocatable :: zsnd,usnd,vsnd,rfrsnd

   integer                             :: maxrgate,maxvgate,maxazim,maxelev
   real, dimension(:),     allocatable :: elvmnvol
   real, dimension(:,:),   allocatable :: rngrvol, azmvol, elvvol, rngvvol
   real, dimension(:,:,:), allocatable :: rxrvol, ryrvol, rzrvol, rxvvol, ryvvol, rzvvol

   integer, parameter        :: nsort=601
   integer, dimension(nsort) :: kntbin

!   real :: missing_r = -888888
   real :: rngmin=10000
   real, dimension(:),allocatable :: rngmaxr, rngmaxv
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine build_qc_radar(rdat,bkgvel)
   implicit none
   type(t_radar_data), intent(inout) :: rdat
   real, dimension(:,:,:), optional, intent(in) :: bkgvel

   integer :: i, j, k, istatus
   integer :: dims(3)

   real :: zagl, sfcr, cosaz, sinaz, sina
   
   character(len=14) :: strTime

   REAL :: deg2rad
   deg2rad = atan(1.)/45.


   !do k=1, rdat%ntilt
   !   sina=sin(deg2rad*rdat%rtilt(1,k))
   !   i=300./sina/rdat%rgatesp(k)
   !   rdat%ref(1:i,:,k)=value_invalid
   !   j=300./sina/rdat%vgatesp(k)
   !   write(*,*) "Set Missing:",k,rdat%rtilt(1,k), i, j
   !   rdat%vel(1:j,:,k)=value_invalid
   !   rdat%spw(1:j,:,k)=value_invalid
   !enddo

   dims     = ubound(rdat%ref)
   maxrgate = dims(1)
   maxazim  = dims(2)
   maxelev  = rdat%ntilt
   dims     = ubound(rdat%vel)
   maxvgate = dims(1)

   if(allocated(rngrvol))then
      deallocate(rngrvol )
      deallocate(rngvvol )
      deallocate(azmvol  )
      deallocate(elvvol  )
      deallocate(rxrvol  )
      deallocate(ryrvol  )
      deallocate(rzrvol  )
      deallocate(rxvvol  )
      deallocate(ryvvol  )
      deallocate(rzvvol  )
      deallocate(elvmnvol)
      deallocate(rngmaxr )
      deallocate(rngmaxv )
   endif
   allocate(rngrvol (maxrgate,          maxelev))
   allocate(rngvvol (maxvgate,          maxelev))
   allocate(azmvol  (          maxazim, maxelev))
   allocate(elvvol  (          maxazim, maxelev))
   allocate(rxrvol  (maxrgate, maxazim, maxelev))
   allocate(ryrvol  (maxrgate, maxazim, maxelev))
   allocate(rzrvol  (maxrgate, maxazim, maxelev))
   allocate(rxvvol  (maxvgate, maxazim, maxelev))
   allocate(ryvvol  (maxvgate, maxazim, maxelev))
   allocate(rzvvol  (maxvgate, maxazim, maxelev))
   allocate(elvmnvol(                   maxelev))
   allocate(rngmaxr (                   maxelev))
   allocate(rngmaxv (                   maxelev))

   elvvol=0.
   rngmaxr=rngmin
   rngmaxv=rngmin

   rxrvol=0
   ryrvol=0
   rzrvol=0
   rxvvol=0
   ryvvol=0
   rzvvol=0

   do k=1, maxelev
      do i=1, maxrgate
         rngrvol(i,k)=i*rdat%rgatesp(k)
      enddo
      do i=1, maxvgate
         rngvvol(i,k)=i*rdat%vgatesp(k)
      enddo
      do j=1, maxazim
         azmvol(j,k)=rdat%razim(j,k)
         elvvol(j,k)=rdat%rtilt(j,k)
      enddo
      elvmnvol(k)=minval(elvvol(:,k))
      rngmaxr (k)=max(rngmaxr(k),rdat%rgatesp(k)*rdat%nrgate(1,k))
      rngmaxv (k)=max(rngmaxv(k),rdat%vgatesp(k)*rdat%nvgate(1,k))

      do j=1,rdat%nazim(k)
         cosaz=cos(deg2rad*rdat%razim(j,k))
         sinaz=sin(deg2rad*rdat%razim(j,k))
         do i=1,rdat%nrgate(j,k)
            CALL beamhgt(rdat%rtilt(j,k),rngrvol(i,k),zagl,sfcr)
            rxrvol(i,j,k)=sinaz*sfcr
            ryrvol(i,j,k)=cosaz*sfcr
            rzrvol(i,j,k)=rdat%altitude+zagl
         enddo
         do i=1,rdat%nvgate(j,k)
            CALL beamhgt(rdat%rtilt(j,k),rngvvol(i,k),zagl,sfcr)
            rxvvol(i,j,k)=sinaz*sfcr
            ryvvol(i,j,k)=cosaz*sfcr
            rzvvol(i,j,k)=rdat%altitude+zagl
         enddo
      enddo
   enddo
   write(strTime,"(I4.4,5I2.2)") rdat%year,rdat%month ,rdat%day   , &
                                 rdat%hour,rdat%minute,rdat%second

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine qc_radar(rdat,bkgvel)
   implicit none
   type(t_radar_data), intent(inout) :: rdat
   real, dimension(:,:,:), optional, intent(in) :: bkgvel

   integer :: i, j, k, istatus
   integer :: dims(3)

   real :: zagl, sfcr, cosaz, sinaz, sina
   real :: min_z, max_z, min_v, max_v, min_w, max_w
   
   character(len=14) :: strTime
   real, dimension(:,:,:), allocatable :: vel_old
   integer, dimension(:,:,:), allocatable :: vel_nn
   integer, dimension(:), allocatable :: nn
   integer :: np, nyq, m, n
   integer, dimension(201) :: dist
   real,    dimension(201) :: p_dist

   integer, dimension(:,:), allocatable :: nyq_num

   REAL :: deg2rad
   deg2rad = atan(1.)/45.


   call calculate_radar_loc(rdat)
   call build_qc_radar(rdat,bkgvel)
   call anomalous_radial_radar(rdat)
   !call write_radar_csv          ("anr."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("anr."// trim(rdat%radar_id)//"."//trim(strTime), rdat)

   ! write(*,*) "after anomalous_radial"
   ! write(*,"(A5,6A10)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, maxelev
   !    call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
   !    call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
   !    call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
   !    write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
   call despeckle_radar(rdat)
   !call write_radar_csv          ("dsp."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("dsp."// trim(rdat%radar_id)//"."//trim(strTime), rdat)

   ! write(*,*) "after despeckle"
   ! write(*,"(A5,6A10)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, maxelev
   !    call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
   !    call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
   !    call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
   !    write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
   call median_filter_radar(rdat)
   !call write_radar_csv          ("mfl."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("mfl."// trim(rdat%radar_id)//"."//trim(strTime), rdat)

   ! write(*,*) "after median_filter"
   ! write(*,"(A5,6A10)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, maxelev
   !    call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
   !    call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
   !    call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
   !    write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
   if(spwthrrat > 0.) then
      CALL spw_threshold_radar(rdat)
    write(*,*) "before unfold"
    write(*,"(A5,6A15)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
    do k=1, maxelev
       call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
       call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
       call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
       write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
    enddo
   endif

   allocate(vel_old(maxvgate, maxazim, maxelev))
   allocate(vel_nn (maxvgate, maxazim, maxelev))
   write(*,*) "ubound1",ubound(vel_old)
   write(*,*) "ubound2",ubound(rdat%vel)
   vel_old(1:maxvgate,1:maxazim,1:maxelev)=rdat%vel(1:maxvgate,1:maxazim,1:maxelev)
   if(present(bkgvel))then
      where(abs(rdat%vel)>200)
         rdat%vel=missing_r
      endwhere
      call unfold_radar(rdat, bkgvel)
!      call unfold_4dd(rdat, missing_r, 1, istatus, soundVolume=bkgvel)
      !if(istatus/=1) stop
   else
      call unfold_radar(rdat)
   endif
   !if(present(bkgvel))then
   !do k=1, maxelev
   !   if(.NOT.rdat%ifvel(k)) cycle
   !   write(*,*) "=========== Level :", k, "====================="
   !   if(allocated(nyq_num))then
   !      deallocate(nyq_num)
   !   endif
   !   allocate(nyq_num(maxvgate,rdat%nazim(k)))
   !   call dealias_region(maxvgate,rdat%nazim(k),vel_old(1:maxvgate,1:rdat%nazim(k),k),&
   !                                           bkgvel,rdat%vmax(k),value_invalid,nyq_num,&
   !                                rdat%vel(1:maxvgate,1:rdat%nazim(k),k))
   !enddo 
   !endif
   !call write_radar_csv          ("unf."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("unf."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
    write(*,*) "after unfold"
    write(*,"(A5,6A15)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
    do k=1, maxelev
       call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
       call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
       call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
       write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
    enddo

   !call quad_refill_unfold_radar(rdat)
   !call write_radar_csv          ("qfu."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("qfu."// trim(rdat%radar_id)//"."//trim(strTime), rdat)

   !call copy_radar_data(rdat_in,rdat)
   write(*, *) 'Detect Anomalous Propagation for Ref...'
   call apdetect(maxrgate,maxvgate,maxazim,maxelev,               &
                 rdat%nrgate(1:maxazim,1:maxelev) ,               &
                 rdat%nazim (1:maxelev), maxelev  ,               &
                 rdat%nvgate(1:maxazim,1:maxelev) ,               &
                 rdat%nazim (1:maxelev), maxelev  ,               &
                 refcheck,velcheck,                               &
                 int(rdat%rgatesp(1)),int(rdat%vgatesp(maxelev)), &
                 winszrad,winszazim,rdat%vcp,gclopt,gcvrlim,      &
                 rngrvol,azmvol,elvvol,                           &
                 rngvvol,azmvol,elvvol,                           &
                 rdat%ref(1:maxrgate,1:maxazim,1:maxelev),        &
                 rdat%vel(1:maxvgate,1:maxazim,1:maxelev),istatus)
   ! write(*,*) "after apdetect"
   ! write(*,"(A5,6A10)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, maxelev
   !    call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
   !    call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
   !    call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
   !    write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
   do i=1,5
      call despeckle_radar(rdat)
   enddo
   ! write(*,*) "after despeckle 5"
   ! write(*,"(A5,6A10)") "k","Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, maxelev
   !    call get_maxmin_2d(rdat%ref(:,:,k),min_z, max_z)
   !    call get_maxmin_2d(rdat%vel(:,:,k),min_v, max_v)
   !    call get_maxmin_2d(rdat%spw(:,:,k),min_w, max_w)
   !    write(*,"(I5,6F15.2)") k, min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
                  

   !call write_radar_csv          ("apd."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   !call write_radar_grads_station("apd."// trim(rdat%radar_id)//"."//trim(strTime), rdat)
   
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gjelim(n,a,rhs,sol,work,work1d,eps,istatus)
   !
   !-----------------------------------------------------------------------
   !
   !  PURPOSE:
   !
   !    Solve an N x N system of equations.       [a][sol]=[rhs]
   !    Using Gauss-Jordan with array pivoting.
   !
   !-----------------------------------------------------------------------
   !
   !  INPUT :
   !
   !    n        Matrix order
   !    a        Array
   !    rhs      Right-hand side vector of matrix equation
   !    eps      Error checking threshold
   !
   !  OUTPUT:
   !    sol      Solution vector
   !    istatus  Status indicator
   !
   !  WORK SPACE:
   !    work     Work Array
   !    work1d   Work vector
   !
   !-----------------------------------------------------------------------
   !
   
   !
   !-----------------------------------------------------------------------
   !
   !  Variable Declarations.
   !
   !-----------------------------------------------------------------------
   !
     IMPLICIT NONE
   !
     integer, INTENT(IN)  :: n
     real,    INTENT(IN)  :: a(n,n)
     real,    INTENT(IN)  :: rhs(n)
     real,    INTENT(OUT) :: sol(n)
     real,    INTENT(OUT) :: work(n,n+1)
     real,    INTENT(OUT) :: work1d(n+1)
     real,    INTENT(IN)  :: eps
     integer, INTENT(OUT) :: istatus
   !
   !-----------------------------------------------------------------------
   !
   ! Misc. Local Variables
   !
   !-----------------------------------------------------------------------
   !
     real    :: pivot,const
     integer :: np1
     integer :: i,j,k,m
   !
   !-----------------------------------------------------------------------
   !
   ! Initialize the work array
   ! First set all elements to zero.
   ! Fill nxn with elements from input a
   ! Fill last column with RHS vector.
   !
   !-----------------------------------------------------------------------
   !
     np1=n+1
   
     do j=1, np1
       do i=1, n
         work(i,j)=0.0
       enddo
     enddo
   
     do j=1, n
       do i=1, n
         work(i,j)=a(i,j)
       enddo
     enddo
   
     do i=1,n
       work(i,np1)=rhs(i)
     enddo
   
     do j=1, n
   !
   !-----------------------------------------------------------------------
   !
   ! Find largest element in column j
   !
   !-----------------------------------------------------------------------
   !
       m=j
       pivot=ABS(work(m,j))
       do i=j+1,n
         if(ABS(work(i,j)) > pivot ) then
           m=i
           pivot=ABS(work(m,j))
         endif
       enddo
   !
   !-----------------------------------------------------------------------
   !
   ! Error trapping
   !
   !-----------------------------------------------------------------------
   !
       if( pivot < eps ) then
         do i=1, n
           sol(i)=0.
         enddo
         istatus=-1
         RETURN
       endif
   !
   !-----------------------------------------------------------------------
   !
   ! Swap rows
   !
   !-----------------------------------------------------------------------
   !
       if(m /= j) then
         do k=1,np1
           work1d(k)=work(j,k)
         enddo
         do k=1,np1
           work(j,k)=work(m,k)
           work(m,k)=work1d(k)
         enddo
       endif
   !
   !-----------------------------------------------------------------------
   !
   ! Normalize Row
   !
   !-----------------------------------------------------------------------
   !
       const=1./work(j,j)
       do k=1,np1
         work(j,k)=const*work(j,k)
       enddo
       work(j,j)=1.0
   !
   !-----------------------------------------------------------------------
   !
   ! Elimination
   !
   !-----------------------------------------------------------------------
   !
       do i=1,n
         if ( i /= j ) then
           const=work(i,j)
           do k=1,np1
             work(i,k)=work(i,k)-const*work(j,k)
           enddo
         endif
       enddo
     enddo
   !
   !-----------------------------------------------------------------------
   !
   ! Transfer last column to sol vector
   !
   !-----------------------------------------------------------------------
   !
     do i=1,n
       sol(i)=work(i,n+1)
     enddo
     istatus = 1
   
     RETURN
   end subroutine gjelim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine despeckle(dat,varcheck)
   implicit none
   real, dimension(:,:,:), intent(inout) :: dat
   real,                   intent(in)    :: varcheck

   integer                                :: ngate,nazim,nelev
   integer                                :: i,j,k
   integer                                :: kntgd,kntdsp
   real                                   :: sum,pctdsp
   integer, dimension(3)                  :: dims
   real   , dimension(:,:,:), allocatable :: tmp
  
   dims=ubound(dat)
   ngate=dims(1)
   nazim=dims(2)
   nelev=dims(3)

   allocate(tmp(ngate,nazim,nelev))
   tmp=0.
   do k=1,nelev
      do j=1,nazim
         do i=1,ngate
            if( dat(i,j,k) > varcheck ) tmp(i,j,k)=1.
         enddo
      enddo
   enddo

   kntgd=0
   kntdsp=0
   do k=1, nelev 
      do j=2,nazim-1
         do i=2,ngate-1
            if(tmp(i,j,k) > 0. ) then
               kntgd=kntgd+1
               sum=tmp(i-1,j+1,k)+tmp(i,j+1,k)+tmp(i+1,j+1,k) &
                  +tmp(i-1,j  ,k)+             tmp(i+1,j  ,k) &
                  +tmp(i-1,j-1,k)+tmp(i,j-1,k)+tmp(i+1,j-1,k)
               if( sum < 3. ) then
                  kntdsp=kntdsp+1
                  dat(i,j,k) = dspmiss
               endif
            endif 
         enddo 
      enddo 
   enddo 

   if(kntgd > 0 )then
     pctdsp=100.*(float(kntdsp)/float(kntgd))
   else
     pctdsp=0.
   endif 
   
   write(*,'(a,i8,a,i8,a,f6.1,a)') ' Despeckled ',kntdsp,' of ',kntgd,' data =',pctdsp,' percent.'
   deallocate(tmp)

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine despeckle_radar(rdat)
   implicit none
   type(t_radar_data),      intent(inout) :: rdat

   write(*,*) "Despeckle Ref..."
   call despeckle(rdat%ref,refcheck)

   write(*,*) "Despeckle Vel..."
   call despeckle(rdat%vel,velcheck)

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_radar_stats(rdat)
   implicit none
   type(t_radar_data), intent(in) :: rdat

   integer :: i, j, k
   real    :: rmin, rmax

   write(*,*) "Print Basic Stats Infomation..."
   do k=1, rdat%ntilt !Tilt loop
      write(*,"(A,I6)")  ' Processing base data... ', k
      if(rdat%ifref(k)) then
         ! Ref Stats
         write(*,"(A)") 'Ref    j    dtime     azim     elev   refmin   refmax' 
         do j = 1, rdat%nazim(k), 60
            rmin=999.
            rmax=-999.
            do i=1,rdat%nrgate(j,k)
               if(rdat%ref(i,j,k) /= value_invalid) then
                  rmin=min(rdat%ref(i,j,k),rmin)
                  rmax=max(rdat%ref(i,j,k),rmax)
               end if
            end do
            write(*,"(3X,I5,I9,4F9.1)") j,rdat%stime(j,k),rdat%razim(j,k),rdat%rtilt(j,k),rmin,rmax
         end do
      endif

      if(rdat%ifvel(k) ) then
         write(*,"(A,F10.2)") ' Nyquist velocity: ',rdat%vmax(k)
         ! Vel Stats
         write(*,"(A)") 'Vel     j     dtime    azim     elev     vr min   vr max' 
         do j = 1, rdat%nazim(k), 60
            rmin=999.
            rmax=-999.
            do i=1,rdat%nvgate(j,k)
               if(rdat%vel(i,j,k) /= value_invalid) then
                  rmin=min(rdat%vel(i,j,k),rmin)
                  rmax=max(rdat%vel(i,j,k),rmax)
               end if
            end do
            write(*,"(3X,I5,I9,4F9.1)") j,rdat%stime(j,k),rdat%razim(j,k),rdat%rtilt(j,k),rmin,rmax
         end do
      endif
   enddo

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine spw_threshold_radar(rdat)
   implicit none
   type(t_radar_data), intent(inout) :: rdat
   
   real    :: spwthresh
   integer :: i,j,k

   write(*,*) "Spw Threshold for Radial Vel..."
   do k=1, rdat%ntilt
      do j=1,rdat%nazim(k)
         if(rdat%vmax(k) > 0.) then
            spwthresh=spwthrrat*rdat%vmax(k)
            do i=1,rdat%nvgate(j,k)
               if(rdat%vel(i,j,k) > velcheck .AND. rdat%spw(i,j,k) > spwthresh ) rdat%vel(i,j,k) = spwflag
            enddo
         endif
      enddo
   enddo
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine median_filter(dat,varcheck)
   implicit none
   real, dimension(:,:,:), intent(inout) :: dat
   real,                   intent(in)    :: varcheck

   integer                                :: i, j, k, m, n, loc, is, ns, ne, ngate,nazim,nelev
   integer, dimension(3)                  :: dims
   real   , dimension(:,:,:), allocatable :: tmp

   real, dimension(9) :: sortdata
   real               :: swp
   integer            :: NN
   integer, parameter :: nmin = 3
   
   dims=ubound(dat)
   ngate=dims(1)
   nazim=dims(2)
   nelev=dims(3)

   allocate(tmp(ngate,nazim,nelev))
   
   do k=1, nelev
      do j=1,nazim
         do i=1,ngate
            tmp(i,j,k)=dat(i,j,k)
         enddo
      enddo
   enddo 
   
   do k=1, nelev
      do j=1,nazim
         do i=2,ngate-1
            NN=0
            do m=-1,1
               ns=-1
               if((j+ns  )<1    ) ns=0
               if((j+ns+2)>nazim) ns=-2
               ne=ns+2
               do n=ns,ne
                  swp=tmp(i+m,j+n,k)
                  if(swp > varcheck) then
                     NN=NN+1
                     if(NN == 1) then ! first
                        sortdata(NN)=swp
                     else
                        do loc=1,NN-1
                           if(swp < sortdata(loc)) then ! little -> big
                              do is=NN,loc+1,-1
                                 sortdata(is)=sortdata(is-1)
                              enddo
                              EXIT
                           endif
                        enddo
                        sortdata(loc)=swp
                     endif
                  endif
              enddo
           enddo
   
           if (NN > nmin) then
              dat(i,j,k)=sortdata((NN+1)/2)
           endif
         enddo   
      enddo 
   enddo
   deallocate(tmp)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine median_filter_radar(rdat)
   implicit none
   type(t_radar_data), intent(inout) :: rdat
   
   write(*,*) "Median Filter Ref..."
   call median_filter(rdat%ref,refcheck)

   write(*,*) "Median Filter Vel..."
   call median_filter(rdat%vel,velcheck)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine anomalous_radial_radar(rdat)
   implicit none
   type(t_radar_data), intent(inout) :: rdat

   real, parameter :: eps=1.0E-20
   real, parameter :: rminlsq=50.0E03 ! 50Km
   real            :: reflin(2,2),refrhs(2),sol(2)
   real            :: work(2,3),work1d(3)
   real            :: refcst,drefdr,dref,sum2,rms
   integer         :: nref,istatus
   
   real   , parameter   :: refnull = 0.0
   
   integer, allocatable :: kntref(:)
   real   , allocatable :: refmean(:)
   real   , allocatable :: refstdv(:)
   real                 :: refrng,sumref,sum2ref
   real                 :: sumrefglb,sum2refglb
   real                 :: refmeanglb,refstdvglb,fngate,fngatem1,vldratio,vldratglb
   real                 :: refvar,refnull2,fngatetot,fngatetotm1
   integer              :: i,j,jj,k,kk,kntrefglb,istat
   integer              :: kntazmflg,ngatetot

   integer :: maxgate
   integer :: ngate, nazim, nelev, kntgate
   real, dimension(:), allocatable :: rngs, refs
   character(len=200) :: strFormat

   maxgate=maxrgate

   allocate(kntref (maxazim))
   allocate(refmean(maxazim))
   allocate(refstdv(maxazim))

   allocate(rngs(maxgate))
   allocate(refs(maxgate))
   
   refnull2    = refnull*refnull

   nelev=rdat%ntilt
   do k=1, nelev 
      if(.NOT.rdat%ifref(k)) cycle
      nazim       = rdat%nazim(k)
      ngate       = rdat%nrgate(1,k)
      kntrefglb   = 0
      sumrefglb   = 0.
      sum2refglb  = 0.
      fngate      = float(ngate)
      fngatem1    = float(ngate-1)
      ngatetot    = ngate*nazim
      fngatetot   = float(ngatetot)
      fngatetotm1 = float(ngatetot-1)
      
      do j=1,nazim
      
         kntref(j) = 0
         sumref    = 0.
         sum2ref   = 0.
      
         do i=1,ngate
            if(rdat%ref(i,j,k) > refnull ) then
               kntref(j) = kntref(j)+1
               sumref    = sumref +rdat%ref(i,j,k)
               sum2ref   = sum2ref+(rdat%ref(i,j,k)*rdat%ref(i,j,k))
            else
               sumref  = sumref+refnull
               sum2ref = sum2ref+refnull2
            endif
         enddo
         kntrefglb  = kntrefglb+kntref(j)
         sumrefglb  = sumrefglb+sumref
         sum2refglb = sum2refglb+sum2ref
         refmean(j) = sumref/fngate
         refvar     = (sum2ref-(sumref*sumref/fngate))/fngatem1
         refstdv(j) = SQRT(refvar)
      enddo
     
      kntazmflg=0
      vldratglb=float(kntrefglb)/float(nazim*ngate)
      write(*,'(A,I4,A,F7.0,A)') 'Elevation',k,' ref>0 gates:',(vldratglb*100.),' percent.'
      if( kntrefglb > 100 ) then
         refmeanglb = sumrefglb/fngatetot
         refvar     = (sum2refglb-(sumrefglb*sumrefglb/fngatetot))/fngatetotm1
         refstdvglb = SQRT(refvar)
         write(*,'(2(A,F5.1))') 'Ref mean:',refmeanglb,' Std dev:',refstdvglb
      
         ! look for a mostly non-zero radial with positive linear trend
         do j=1,nazim
            vldratio=float(kntref(j))/fngate
            if(((refmeanglb > 0. .OR. vldratglb > 0.10) .and. vldratio > 0.5 ).or. & ! Precip situation mostly non-missing radial
               ( vldratio > 0.40 .AND. (refmean(j)-refmeanglb) > 5.0 )) then  ! clear situation, look for a high valid ratio
               !
               !  Least-squares linear trend
               !  refl = a + b*range
               !
               reflin=0.
               refrhs=0.
               kntgate=0
               do i=1,ngate
                  refrng=i*rdat%rgatesp(k)
                  if(refrng > rminlsq .AND. rdat%ref(i,j,k) > refcheck) then
                     kntgate=kntgate+1
                     rngs(kntgate)=refrng
                     refs(kntgate)=rdat%ref(i,j,k)
                     reflin(1,1) = reflin(1,1)+1.0
                     reflin(1,2) = reflin(1,2)+refrng
                     reflin(2,1) = reflin(2,1)+refrng
                     reflin(2,2) = reflin(2,2)+refrng*refrng
                     refrhs(1)   = refrhs(1)+rdat%ref(i,j,k)
                     refrhs(2)   = refrhs(2)+rdat%ref(i,j,k)*refrng
                  endif
               enddo
      
               CALL gjelim(2,reflin,refrhs,sol,work,work1d,eps,istatus)
      
               nref   = 0
               sum2   = 0.
               refcst = sol(1)
               drefdr = sol(2)
      
               do i=1,ngate
                  refrng=i*rdat%rgatesp(k)
                  if(refrng > rminlsq .AND. rdat%ref(i,j,k) > refcheck) then
                     nref = nref+1
                     dref = (refcst+drefdr*refrng) - rdat%ref(i,j,k)
                     sum2 = sum2+(dref*dref)
                  endif
               enddo
               rms=sqrt(sum2/float(nref))
               if( drefdr > 2.0E-05 .AND. drefdr < 9.0E-05 .AND. rms < 4.0) then
      
                  write(*,'(A,F7.1,A)') '*Anomalous radial detected (a).  Azim: ',rdat%razim(j,k),' degrees'
                  !
                  !   Flag reflectivity gates in this radial
                  !
                  write(strFormat,"(A,I4,A)") "(2(A,','),2(I4,','),F8.2,",kntgate,"(',',F10.1))"
                  write(int(ABS(anrflag)),FMT=strFormat) "rng",trim(rdat%radar_id),j,k,rdat%razim(j,k),(rngs(i),i=1,kntgate)
                  write(int(ABS(anrflag)),FMT=strFormat) "ref",trim(rdat%radar_id),j,k,rdat%razim(j,k),(refs(i),i=1,kntgate)

                  kntazmflg=kntazmflg+1
                  do i=1,ngate
                     if(rdat%ref(i,j,k) > refcheck) then
                        rdat%ref(i,j,k)=anrflag
                     endif
                  enddo
               endif
            endif
         enddo ! nazim
      endif 
      write(*,'(a,i4,a,i4,a)') ' Anomradial: ',kntazmflg, ' radials flagged of ',nazim,' radials this level'
   enddo
              
   deallocate(kntref)
   deallocate(refmean)
   deallocate(refstdv)
  
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine apdetect(nrefgates,nvelgates,maxazim,maxelev,                &
                       kntrgate,kntrazim,kntrelev,                         &
                       kntvgate,kntvazim,kntvelev,                         &
                       refcheck,velcheck,                                    &
                       irefgatsp,ivelgatsp,                                &
                       winszrad,winszazim,i_vcp,gcopt,gcvrlim,             &
                       rngrvol,azmrvol,elvrvol,                            &
                       rngvvol,azmvvol,elvvvol,                            &
                       refvol,velvol,istatus)
                       
   !
   !-----------------------------------------------------------------------
   !
   !  PURPOSE:
   !
   !  Anomalous propagation check and clutter check for radar data.
   !
   !-----------------------------------------------------------------------
   !
   !  INPUT :
   !
   !    maxgate   Maximum gates in a radial
   !    maxazim   Maximum radials per tilt
   !    maxelev   Maximum number of tilts
   !
   !    dazlim   Maximum value of azimuth difference (grid vs data) to accept
   !             Generally should be 30 degrees or less for velocity, 360 for refl
   !    rngmin   Minimum range (m) of data to use
   !            (10 000 m or more to eliminate near field ground targets).
   !    rngmax   Maximum range (m) of data to use
   !
   !    rngvvol  Range to gate in velocity 3-D volume
   !    azmvvol  Azimuth angle in velocity 3-D volume
   !    elvvvol  Elevation angle in velocity 3-D volume
   !    varvol   Radar data 3-D volume
   !
   !    xs       x coordinate of scalar grid points in physical/comp. space (m)
   !    ys       y coordinate of scalar grid points in physical/comp. space (m)
   !    zps      Vertical coordinate of scalar grid points in physical space(m)
   !
   !  OUTPUT:
   !    varvol   Radar data 3-D volume
   !    istatus  Status indicator
   !
   !-----------------------------------------------------------------------
   !
   
   !
   !-----------------------------------------------------------------------
   !
   !  Variable Declarations.
   !
   !-----------------------------------------------------------------------
   !
     IMPLICIT NONE
   !
     integer, INTENT(IN) :: nrefgates
     integer, INTENT(IN) :: nvelgates
     integer, INTENT(IN) :: maxazim
     integer, INTENT(IN) :: maxelev
   
     integer, INTENT(IN) :: kntrgate(maxazim,maxelev)
     integer, INTENT(IN) :: kntrazim(maxelev)
     integer, INTENT(IN) :: kntrelev
   
     integer, INTENT(IN) :: kntvgate(maxazim,maxelev)
     integer, INTENT(IN) :: kntvazim(maxelev)
     integer, INTENT(IN) :: kntvelev
   
     integer, INTENT(IN) :: i_vcp
     integer, INTENT(IN) :: gcopt
     real, INTENT(IN)    :: gcvrlim
   
     real, INTENT(IN)    :: winszrad
     real, INTENT(IN)    :: winszazim
     real, INTENT(IN)    :: refcheck
     real, INTENT(IN)    :: velcheck
     integer, INTENT(IN)    :: irefgatsp
     integer, INTENT(IN)    :: ivelgatsp
     real, INTENT(IN)    :: rngrvol(nrefgates,maxelev)
     real, INTENT(IN)    :: azmrvol(maxazim,maxelev)
     real, INTENT(IN)    :: elvrvol(maxazim,maxelev)
     real, INTENT(IN)    :: rngvvol(nvelgates,maxelev)
     real, INTENT(IN)    :: azmvvol(maxazim,maxelev)
     real, INTENT(IN)    :: elvvvol(maxazim,maxelev)
     real, INTENT(INOUT) :: refvol(nrefgates,maxazim,maxelev)
     real, INTENT(INOUT) :: velvol(nvelgates,maxazim,maxelev)
   
     integer, INTENT(OUT) :: istatus
   !
   !-----------------------------------------------------------------------
   !
   ! Misc. Local Variables
   !
   !-----------------------------------------------------------------------
   !
   
     real                :: tmp(MAX(nrefgates,nvelgates),maxazim)
     integer :: ii,jj,kk,kkv,kv,kr2,i,j,k,knt
     integer :: mvelok,spin,flags
     integer :: igate,ivgate,jazim,kelev,jazmin,jazmax
     integer :: iigate,jjazim,jazmref,jazmvel
     integer :: irngmin,irngmax
     integer :: igatebgn,igateend,jazmbgn,jazmend,indgate,indazm
     integer :: igspratio,iwinszrad,iwinszazm
     integer :: kntcheck,kntapref,kntapvel
   
     real :: azmdiff,maxdbz
     real :: summdbz,sumvel,sumvel2,sumtdbz,sumtvel
     real :: dbzdiff,sign,refprev,prevvel,veldiff
     real :: all_counts,spinchange_counts,signcnt
     real :: delev,avgelv,avgelvv,avgelvr,avgelvr2,appct
     real :: mdbz,tdbz,deltdbz,spinchange,stdvel,mvel,tvel
     real :: elvmax
     LOGICAL :: found
   
     real, parameter :: elvmincmp = 1.1
   ! real, parameter :: elvmincmp = 7.1
     real, parameter :: dbzthr = 10.0
     real, parameter :: apflag = -888.0
     real, parameter :: gcflag = -890.0
     real, parameter :: spinthr = 15.0
     real, parameter :: ddbzthr = -12.0
     real, parameter :: ddbzthr2 = -20.0
     real, parameter :: mvelthr = 2.3
     real, parameter :: apvelthr = 5.0
   !
   ! Parameters applied for gcopt=2
   !
     real, parameter :: dbzclutterL = -5
     real, parameter :: dbzclutterH = 200
     real, parameter :: velclutter = 1.0
     real :: drnglim
     integer ::  iii,jjj,kkk,ku
   !
   !-----------------------------------------------------------------------
   !
   !  Misc initializations
   !
   !-----------------------------------------------------------------------
   !
     istatus=0
   
     if( kntrelev > 1 ) then
   
       igspratio=max(int(float(irefgatsp)/float(2*ivelgatsp)),1)
       iwinszrad=max(int(winszrad/float(2*irefgatsp)),1)
       iwinszazm=max(int(0.5*winszazim),1)
   !
   !-----------------------------------------------------------------------
   !
   !  Establish the index of the reflectivity tilt above elvmincmp degrees.
   !  For NEXRAD the first two unique tilt angles should be 0.5
   !  and 1.5 degrees.
   !
   !-----------------------------------------------------------------------
   !
       found=.false.
       avgelvr2=0.
       do kk=2,kntrelev
         avgelvr2=0.
         do jazim=1,kntrazim(kk)
           avgelvr2=avgelvr2+elvrvol(jazim,kk)
         enddo
         avgelvr2=avgelvr2/float(kntrazim(kk))
         print *, ' kk= ',kk,'  avgelvr2=',avgelvr2
         if(avgelvr2 > elvmincmp) then
           found=.true.
           EXIT
         endif
       enddo
       if(found) then
       kr2=kk
       write(*,'(1x,a,f9.2,a,f9.2,a,i6)') 'Elev 1st ref tilt above ',     &
             elvmincmp,' deg: ',avgelvr2,' kr2=',kr2
   
   
   !
   !-----------------------------------------------------------------------
   !
   !  For now, AP detection is only done on any tilts below elvmincmp.
   !
   !-----------------------------------------------------------------------
   !
       do kk = 1, (kr2-1)
   
         tmp=0.
         kntcheck=0
         kntapref=0
         kntapvel=0
   
   !-----------------------------------------------------------------------
   !
   !  Establish the velocity tilt that most closely matches the
   !  this reflectivity elevation angle.
   !
   !-----------------------------------------------------------------------
   !
         avgelvr=0.
         do jazim=1,kntrazim(kk)
           avgelvr=avgelvr+elvrvol(jazim,kk)
         enddo
         avgelvr=avgelvr/float(kntrazim(kk))
   
         kv=1
         delev=999.
         avgelvv=-99.
         do kkv=1,kntvelev
           avgelv=0.
           do jazim=1,kntvazim(kkv)
             avgelv=avgelv+elvvvol(jazim,kkv)
           enddo
           avgelv=avgelv/float(kntvazim(kkv))
           if(abs(avgelv-avgelvr) < delev ) then
             delev=abs(avgelv-avgelvr)
             avgelvv=avgelv
             kv=kkv
           endif
         enddo
         write(*,'(a,f9.2,a,i6)') ' Elev of nearest velocity tilt: ',      &
           avgelvv,' kv=',kv
   
         do jjazim = 1,kntrazim(kk)
   !
   !  Find nearest azimuth in reflectivity data above elvmincmp
   !
           jazmref=1
           azmdiff = 999.0
           do jazim = 1,kntrazim(kr2)
             if(abs(azmrvol(jazim,kr2)-azmrvol(jjazim,kk))<azmdiff)then
               azmdiff = abs(azmrvol(jazim,kr2)-azmrvol(jjazim,kk))
               jazmref = jazim
             endif
           enddo
   !
   !-----------------------------------------------------------------------
   !
   !  Find nearest azimuth in velocity data
   !
   !-----------------------------------------------------------------------
   !
           jazmvel=1
           azmdiff = 999.0
           do jazim = 1,kntvazim(kv)
             if(abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk)) < azmdiff        &
                .and. kntvgate(jazim,kv)>0 ) then
               azmdiff = abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))
               jazmvel = jazim
             endif
           enddo
   !
   !-----------------------------------------------------------------------
   !
   !  Calculate mean ref, texture of ref, spin
   !
   !-----------------------------------------------------------------------
   !
           do iigate= 1,kntrgate(jjazim,kk)
             if(refvol(iigate,jjazim,kk) > refcheck )then
               kntcheck=kntcheck+1
               summdbz=0.
               sumtdbz=0.
               all_counts=0.
               spinchange_counts=0.
               spin = 1
   ! iwinszrad:  radial  window
   ! iwinszazm:  azimuth window
               irngmin = max((iigate - iwinszrad),1)
               irngmax = min((iigate + iwinszrad),kntrgate(jjazim,kk))
               jazmin = jjazim - iwinszazm
               jazmax = jjazim + iwinszazm
               if(jazmin <= 0) jazmin = jazmin + kntrazim(kk)-2
               if(jazmax > kntrazim(kk) )jazmax = jazmax-kntrazim(kk)+2
   !
   !-----------------------------------------------------------------------
   !
   !  Loop forward from jazmin to jazmax
   !
   !-----------------------------------------------------------------------
   !
               dbzdiff=0.
               refprev=0.
               signcnt=0.
               if(jazmax > jazmin)then
                 do jazim = jazmin,jazmax
                   do igate = irngmin,irngmax
                     if(refvol(igate,jazim,kk)>refcheck)then
                       summdbz=summdbz+refvol(igate,jazim,kk)
                       sumtdbz=sumtdbz+dbzdiff*dbzdiff
                       all_counts = all_counts + 1
                       dbzdiff = refvol(igate,jazim,kk)-refprev
                       if(dbzdiff > dbzthr)then
                         if(spin<0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = 1
                         endif
                       endif
                       if(dbzdiff < -dbzthr)then
                         if(spin>0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = -1
                         endif
                       endif
                       if(dbzdiff > 0.)then
                         signcnt = signcnt + 1.
                       else
                         signcnt = signcnt - 1.
                       endif
                       refprev = refvol(igate,jazim,kk)
                     endif
                   enddo ! igate
                 enddo ! jazim
   !
   !-----------------------------------------------------------------------
   !
   ! if jazmax<kazmin
   !
   !-----------------------------------------------------------------------
   !
               else
                 do jazim = jazmin,kntrazim(kk)
                   do igate = irngmin,irngmax
                     if(refvol(igate,jazim,kk)>refcheck)then
                       summdbz=summdbz+refvol(igate,jazim,kk)
                       sumtdbz=sumtdbz+dbzdiff*dbzdiff
                       all_counts = all_counts + 1
                       dbzdiff = refvol(igate,jazim,kk)-refprev
                       if(dbzdiff > dbzthr)then
                         if(spin<0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = 1
                         endif
                       endif
                       if(dbzdiff < -dbzthr)then
                         if(spin>0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = -1
                         endif
                       endif
                       if(dbzdiff > 0.)then
                         signcnt = signcnt + 1.
                       else
                         signcnt = signcnt - 1.
                       endif
                       refprev = refvol(igate,jazim,kk)
                     endif
                   enddo ! igate
                 enddo ! jazim
                 do jazim = 1,jazmax
                   do igate = irngmin,irngmax
                     if(refvol(igate,jazim,kk)>refcheck)then
                       summdbz=summdbz+refvol(igate,jazim,kk)
                       sumtdbz=sumtdbz+dbzdiff*dbzdiff
                       all_counts = all_counts + 1
                       dbzdiff = refvol(igate,jazim,kk)-refprev
                       if(dbzdiff > dbzthr)then
                         if(spin<0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = 1
                         endif
                       endif
                       if(dbzdiff < -dbzthr)then
                         if(spin>0)then
                           spinchange_counts = SPINchange_counts + 1
                           spin = -1
                         endif
                       endif
                       if(dbzdiff > 0.)then
                         signcnt = signcnt + 1.
                       else
                         signcnt = signcnt - 1.
                       endif
                       refprev = refvol(igate,jazim,kk)
                     endif
                   enddo ! igate
                 enddo ! jazim
               endif
   
               if(all_counts > 0.)then
                 mdbz= summdbz/all_counts
                 tdbz= sumtdbz/all_counts
                 spinchange= (spinchange_counts/all_counts)*100.
               endif
   !
   !  Calculate difference in reflectivity between this gate and the
   !  first level above elvmincmp.
   !
               if(iigate <= kntrgate(jazmref,kr2) .and.   &
                  refvol(iigate,jazmref,kr2) > refcheck ) then
                  deltdbz=                            &
                    refvol(iigate,jazmref,kr2)-refvol(iigate,jjazim,kk)
               else
                 if(i_vcp==31 .or. i_vcp==32) then
                   deltdbz= -32.0-refvol(iigate,jjazim,kk)
                 else
                   deltdbz= -5.0-refvol(iigate,jjazim,kk)
                 endif
               endif
   !
   !  If the delta-dBZ check fails when comparing to the nearest gate,
   !  account for tilt of echoes by finding max of gates in the neighborhood.
   !
               if (deltdbz < ddbzthr) then
                 jazmbgn = max((jazmref-1),1)
                 jazmend = min((jazmref+1),kntrazim(kr2))
                 maxdbz = refvol(iigate,jazmref,kr2)
                 do indazm = jazmbgn,jazmend
                   igatebgn = max((iigate-1),1)
                   igateend = min((iigate+1),kntrgate(indazm,kr2))
                   do indgate = igatebgn,igateend
                     maxdbz = max(maxdbz,refvol(indgate,indazm,kr2))
                   enddo
                 enddo
                 if( maxdbz > refcheck )                                    &
                   deltdbz = maxdbz - refvol(iigate,jjazim,kk)
                 endif
   !
   !-----------------------------------------------------------------------
   !
   !  Find std of vel, texture of vel, mean vel
   !
   !  First, find nearest velocity gate
   !
   !  Nearest velocity tilt (kv) and velocity azimuth (jazmvel) were
   !  found earlier
   !
   !-----------------------------------------------------------------------
   !
                 do igate=2,kntvgate(jazmvel,kv)-1
                   if(rngvvol(igate,kv) >= rngrvol(iigate,kk)) EXIT
                 enddo
                 if(abs(rngvvol(igate,kv)-rngrvol(iigate,kk)) <            &
                    abs(rngvvol(igate-1,kv)-rngrvol(iigate,kk)) ) then
                   ivgate=igate
                 else
                   ivgate=igate-1
                 endif
   
                 if(abs(rngvvol(ivgate,kv)-rngrvol(iigate,kk))<            &
                    4*max(irefgatsp,ivelgatsp) .and.                       &
                    abs(azmvvol(jazmvel,kv)-azmrvol(jjazim,kk))<           &
                    float(iwinszazm))then
                   sumvel2=0.
                   sumvel =0.
                   sumtvel=0.
                   mvelok=0
                   irngmin = max((ivgate-igspratio),1)
                   irngmax = min((ivgate+igspratio),kntvgate(jazmvel,kv))
                   jazmin = jazmvel - iwinszazm
                   jazmax = jazmvel + iwinszazm
                   if( jazmin < 1 ) jazmin = jazmin+kntvazim(kv)
                   if( jazmax > kntvazim(kv) ) jazmax = jazmax-kntvazim(kv)
   !-----------------------------------------------------------------------
   !
   !  Loop forward from jazmin to jazmax
   !
   !-----------------------------------------------------------------------
   !
                   if(jazmax > jazmin)then
                     do jazim = jazmin,jazmax
                       prevvel = velvol(irngmin,jazim,kv)
                       do igate = irngmin+1,irngmax
                         if(velvol(igate,jazim,kv)>velcheck)then
                           sumvel=sumvel+velvol(igate,jazim,kv)
                           sumvel2=sumvel2+                                &
                             velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                           mvelok=mvelok+1
                           veldiff=velvol(igate,jazim,kv) - prevvel
                           sumtvel=sumtvel+veldiff*veldiff
                           prevvel=velvol(igate,jazim,kv)
                         endif
                     enddo ! igate
                   enddo ! jazim
   !
   !-----------------------------------------------------------------------
   !
   ! if jazmax<kazmin
   !
   !-----------------------------------------------------------------------
   !
                 else
                   do jazim = jazmin,kntvazim(kv)
                     prevvel = velvol(irngmin,jazim,kv)
                     do igate = irngmin+1,irngmax
                       if(velvol(igate,jazim,kv)>velcheck)then
                         sumvel=sumvel+velvol(igate,jazim,kv)
                         sumvel2=sumvel2+                                  &
                           velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                         mvelok=mvelok+1
                         veldiff=velvol(igate,jazim,kv) - prevvel
                         sumtvel=sumtvel+veldiff*veldiff
                         prevvel=velvol(igate,jazim,kv)
                       endif
                     enddo ! igate
                   enddo ! jazim
                   do jazim = 1,jazmax
                     prevvel = velvol(irngmin,jazim,kv)
                     do igate = irngmin,irngmax
                       if(velvol(igate,jazim,kv)>velcheck)then
                         sumvel=sumvel+velvol(igate,jazim,kv)
                         sumvel2=sumvel2+                                   &
                           velvol(igate,jazim,kv)*velvol(igate,jazim,kv)
                         mvelok=mvelok+1
                         veldiff= velvol(igate,jazim,kv) - prevvel
                         sumtvel=sumtvel+veldiff*veldiff
                         prevvel = velvol(igate,jazim,kv)
                       endif
                     enddo ! igate
                   enddo ! jazim
                 endif
                 if(mvelok>1)then
                   stdvel = (sumvel2-sumvel*sumvel/mvelok)/(mvelok-1)
                   if(stdvel .lt. 0.0) stdvel = 0.0
                   stdvel = sqrt(stdvel)
                   tvel= sumtvel/mvelok
                   mvel= abs(sumvel/mvelok)
                 else
                   stdvel=999.
                   tvel=999.
                   mvel=999.
                 endif
               else  ! if ivgate is too far away from iigate
                 stdvel=999.
                 tvel=999.
                 mvel=999.
               endif
   
   !
   !  If calculated values match two or more characteristics,
   !  set the AP marker, tmp=1.  The actual resetting of refvol
   !  is deferred until all points have been checked on this level,
   !  so that the calculation of spin at subsequent gates is not affected.
   !
   
               flags=0
               if(spinchange >= spinthr) flags=flags+1
               if(deltdbz <= ddbzthr) flags=flags+1
               if(deltdbz <= ddbzthr2) flags=flags+1
               if(abs(mvel) < mvelthr ) flags=flags+1
   
               if(flags > 1 ) tmp(iigate,jjazim)=1.0
   
             endif  ! a valid reflectivity at iigate,jjazim,kk
           enddo  !iigate
         enddo  !jjazim
   
   !
   !   Where the AP marker has been set, set for reflectivity to apflag,
   !   and set corresponding velocities to apflag
   !
   
         do jjazim=1,kntrazim(kk)
   !
           jazmvel=1
           azmdiff = 999.0
           do jazim = 1,kntvazim(kv)
             if(abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))<azmdiff)then
               azmdiff = abs(azmvvol(jazim,kv)-azmrvol(jjazim,kk))
               jazmvel = jazim
             endif
           enddo
   !
           do iigate=1,kntrgate(jjazim,kk)
             if(tmp(iigate,jjazim) > 0. ) then
               kntapref=kntapref+1
               refvol(iigate,jjazim,kk) = apflag
               do igate=2,kntvgate(jazmvel,kv)-1
                 if(rngvvol(igate,kv) >= rngrvol(iigate,kk)) EXIT
               enddo
   
               if(abs(rngvvol(igate,kv)-rngrvol(iigate,kk)) <              &
                  abs(rngvvol(igate-1,kv)-rngrvol(iigate,kk)) ) then
                 ivgate=igate
               else
                 ivgate=igate-1
               endif
   
               if(abs(rngvvol(ivgate,kv)-rngrvol(iigate,kk))<              &
                  4*max(irefgatsp,ivelgatsp) .and.                         &
                  abs(azmvvol(jazmvel,kv)-azmrvol(jjazim,kk))<             &
                  float(iwinszazm))then
               irngmin = ivgate - int(irefgatsp/float(2*ivelgatsp))
               irngmin = max(irngmin,1)
               irngmax = ivgate + int(irefgatsp/float(2*ivelgatsp))
               irngmax = min(irngmax,kntvgate(jazmvel,kv))
               do igate=irngmin,irngmax
                 kntapvel=kntapvel+1
                 velvol(igate,jazmvel,kv)=apflag
               enddo
               ENDif
             endif
           enddo  ! iigate loop
         enddo  ! jjazim loop
   !
   !  For gcopt=2
   !     Zhao Kun Add Residual AP clutter Remove 2008.12.9
   !
         drnglim=0.5*irefgatsp
         if(gcopt == 2) then
           do ku=1,kntrelev
             do jj=1,kntrazim(ku)
               avgelvr2= elvrvol(jj,kk)
               do ii = 1,kntrgate(jj,ku)
   !
   !             If the reflectivity is suspected to be clutter or
   !             AP clutter then contintue check velocity
   !
                 if(    refvol(ii,jj,ku)  < dbzclutterL .OR.               &
                    abs(refvol(ii,jj,ku)) > dbzclutterH) then
   !               print*,'zhaokun',refvol(ii,jj,ku)
                   refvol(ii,jj,ku)=apflag
                   do jjj = 1,kntvazim(ku)
                     if(abs(azmrvol(jj,ku)-azmvvol(jjj,ku)) < 1.0 .AND.    &
                        abs(elvrvol(jj,ku)-elvvvol(jjj,ku)) < 0.3 ) then
                       do iii = 1,kntvgate(jjj,ku)
                         if(abs(rngrvol(ii,ku)-rngvvol(iii,ku)) <=         &
                                                           drnglim) then
                           if(abs(velvol(iii,jjj,ku))< velclutter)         &
                             velvol(iii,jjj,ku)= apflag
                         endif
                       enddo
                     endif
                   enddo
   
                 endif
               enddo
             enddo
           enddo
         endif   ! gcopt == 2
   
         if ( kntcheck > 0 ) then
           appct=100.*float(kntapref)/float(kntcheck)
         else
           appct=0.
         endif
         write(*,'(a,i6,a,f6.2,/a,i8,/a,i8,f9.2,a,/a,i8)')                 &
          ' AP detect completed for level ',kk,' elev=',avgelvr,           &
          '   Reflectivity gates checked:',kntcheck,                        &
          '               AP flagged ref:',kntapref,appct,' percent',      &
          '               AP flagged vel:',kntapvel
   
   !
   !   Apply despekl to the edited data
   !   This should help catch AP residue.
   !   Note: here force no median filter.
   !
         !CALL despekl(nrefgates,maxazim,nrefgates,maxazim,refcheck,0,refvol(1,1,kk))
         call despeckle(refvol,refcheck)
                      
       enddo  ! kk
   
   ! open(21,file='ap.dat',form='unformatted',status='unknown')
   ! write(21) nrefgates
   ! write(21) maxazim
   ! write(21) maxelev
   ! write(21) kntrelev
   ! write(21) kntrazim
   ! write(21) kntrgate
   ! write(21) rngrvol
   ! write(21) azmvvol
   ! write(21) elvrvol
   ! write(21) refvol
   ! close(21)
       else
         write(*,'(1x,a,f9.2,a,f9.2,/3x,a//)') 'Highest elev not above ', &
            elvmincmp,' deg: ',avgelvr2,'Skipping AP Detect'
       endif
   
     else
       write(*,'(a)') ' Need at least two reflectivity levels for AP detection'
       write(*,'(a)') ' Skipping AP detect'
     endif
   
     if(gcopt > 0) then
       do kk=1,kntvelev
         elvmax=-99.
         do jazim=1,kntvazim(kk)
           elvmax=max(elvmax,elvvvol(jazim,kk))
         enddo
   
         tmp = 0.0
         do jazim=1,kntvazim(kk)
           do igate=1, kntvgate(jazim,kk)
             tmp(igate,jazim) = velvol(igate,jazim,kk)
           enddo
         enddo
   
         do jazim=2,kntvazim(kk)-1
           do igate=2, kntvgate(jazim,kk)-1
             if( abs(tmp(igate+1,jazim)) < 1e-1 .AND.        &
                 abs(tmp(igate-1,jazim)) < 1e-1 .AND.        &
                 abs(tmp(igate,jazim+1)) < 1e-1 .AND.        &
                 abs(tmp(igate,jazim-1)) < 1e-1 )    then
                   velvol(igate,jazim,kk) = gcflag
             endif
           enddo
         enddo
   
         if(elvmax < elvmincmp) then
           do jazim=1,kntvazim(kk)
             if(elvvvol(jazim,kk) < elvmincmp) then
               do igate=1,kntvgate(jazim,kk)
                 if(abs(velvol(igate,jazim,kk)) < gcvrlim) then
                   velvol(igate,jazim,kk) = gcflag
                 endif
               enddo
             endif
           enddo
         endif
       enddo
     endif   ! gcchkopt
   
     RETURN
   end subroutine apdetect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine quad_refill_unfold_radar(rdat)
   implicit none
   type(t_radar_data), intent(inout) :: rdat


   integer, parameter :: n = 6

   real, parameter :: eps = 1.0E-25

   real :: avar(n,n)
   real :: rhsvar(n)
   real :: avel(n,n)
   real :: rhsvel(n)
   real :: sol(n)
   real :: work(n,n+1)
   real :: work1d(n+1)

   real :: array(4,4)
   real :: rhsv(4)
   real :: solv(4)

   integer :: ii,jj,kk,i,j,k,knt,kinbox,kntall,kntdaz
   integer :: kok,isort,jsort,mid
   integer :: kbgn,kend
   integer :: igate,jazim,kelev,jazmin,jmirror,jend,jknt
   integer :: iigate,jjazim,kkelev
   integer :: istatal,istatwrt
   integer :: nbeam,nrang
   integer :: kntchk,kntfold,kntflag
   integer :: maxgate
   integer :: istatus

   ! jsrchmn is the minimum number of radials to search to find data
   ! that are within the distance limits of the grid center.
   !
   integer, parameter :: jsrchmn=2, nsort=601, iorder=2
   real,    parameter :: dxthr0 = 1000.    ! minimum dx processing radius
   real :: dazlim=30.,vmedlim=15., dsort=0.5, sortmin=-100.
   integer :: kntbin(nsort)

   real :: deg2rad,rad2deg
   real :: twonyq,inv2nyq
   real :: sortscale
   real :: delx,dely,delz,azdiff
   real :: ddx,ddxy,ddx2,ddy,ddy2
   real :: cosaz,sinaz,sfcr,zagl
   real :: sum,sum2,sdev,thresh,slrange,elijk,azimijk,time
   real :: varmax,varmin,varavg,varmean,varmed
   real :: daz,sumdaz,azspc,dxthr
   real :: xs,ys,zps
   real :: radarx, radary, radalt
   real :: rngmin, rngmax

   real :: origvel,tstdev,fitvel,tstvel,tstdif,vardif

   integer, dimension(3) :: dims
   character(len=200) :: strFormat

   real, allocatable :: rngvol(:,:), rxvol(:,:,:), ryvol(:,:,:), rzvol(:,:,:)

   write(*,*) 'Quadric Fill Unfold Rejected Vel...'

   maxgate=maxvgate

   radarx=0.
   radary=0.
   radalt=rdat%altitude

   allocate(rngvol (maxgate,          maxelev))
   allocate(rxvol  (maxgate, maxazim, maxelev))
   allocate(ryvol  (maxgate, maxazim, maxelev))
   allocate(rzvol  (maxgate, maxazim, maxelev))
   rngvol=rngvvol

   rngmin=10000.

   sortscale=1./dsort

   time=0.

   rngmax=rngmin
   rxvol=rxvvol
   ryvol=ryvvol
   rzvol=rzvvol
   do kelev=1,rdat%ntilt
      if(.NOT.rdat%ifvel(kelev))cycle
      rngmax=rngmaxv(kelev)
   enddo

   do kkelev = 1, rdat%ntilt
      if(.NOT.rdat%ifvel(kkelev))cycle
      kntdaz=0
      sumdaz=0.
      do jjazim = 2,rdat%nazim(kkelev)
         daz=abs(rdat%razim(jjazim,kkelev)-rdat%razim(jjazim-1,kkelev))
         daz=min(daz,(360.-daz))
         kntdaz=kntdaz+1
         sumdaz=sumdaz+daz
      enddo
      if( kntdaz > 0 ) then
         azspc=sumdaz/float(kntdaz)
         write(*,*) ' Azimuth spacing: ',azspc,' degrees'
         azspc=deg2rad*azspc
      else
         azspc=1.0
      endif
      kntchk=0
      kntfold=0
      kntflag=0
      do jjazim = 1,rdat%nazim(kkelev)
         twonyq=2.0*rdat%vmax(kkelev)
         inv2nyq=1./twonyq
         do iigate= 1,rdat%nvgate(jjazim,kkelev)
            if(rdat%vel(iigate,jjazim,kkelev) < -1500.0) then
               kntchk=kntchk+1
               origvel=rdat%vel(iigate,jjazim,kkelev)+2000.

               xs  = rxvol(iigate,jjazim,kkelev)
               ys  = ryvol(iigate,jjazim,kkelev)
               zps = rzvol(iigate,jjazim,kkelev)

               kok=0
               sum=0.
               sum2=0.
               sdev=0.
               varavg=999999.
               kntbin=0
               delx=xs -radarx
               dely=ys -radary
               delz=zps-radalt

               slrange=rngvol    (iigate,kkelev)
               elijk  =rdat%rtilt(jjazim,kkelev)
               azimijk=rdat%razim(jjazim,kkelev)

               dxthr=max(dxthr0,(2.1*azspc*rngvol(iigate,kkelev)))
               write(*,*) ' dxthr =',(0.001*dxthr),' km'

               varmax=-999.
               varmin=999.
               do jj=1,n
                  do ii=1,n
                     avar(ii,jj)=0.
                  enddo
               enddo

               do ii=1,n
                  rhsvar(ii)=0.
               enddo

               kbgn=kkelev
               do k=kkelev-1,1,-1
                  if(.NOT.rdat%ifvel(k))cycle
                  if((elvmnvol(kkelev)-elvmnvol(k)) > 0.1) then
                     kbgn=k
                     EXIT
                  endif
               enddo
                                  
               kend=kkelev
               do k=kkelev+1,rdat%ntilt
                  if(.NOT.rdat%ifvel(k))cycle
                  if((elvmnvol(k)-elvmnvol(kkelev)) > 0.1) then
                     kend=k
                     EXIT
                  endif
               enddo
               write(*,*) ' Using levels',kbgn,' to ',kend,' at klevel:',kkelev
               write(*,*) '    elvmnvol(kkelev):',elvmnvol(kkelev)
          
               !  First pass, find min,max,mean,median.

               do kk=kbgn,kend
                  if(.NOT.rdat%ifvel(kk))cycle

                  !  Find nearest azimuth at this level

                  if(kk == kkelev) then
                     jazmin=jjazim
                  else
                     azdiff=181.
                     jazmin=1
                     do jazim=1,rdat%nazim(kk)
                        daz=rdat%razim(jazim,kk)-azimijk
                        if(daz >  180.) daz=daz-360.
                        if(daz < -180.) daz=daz+360.
                        daz=abs(daz)
                        if(daz < azdiff) then
                           azdiff=daz
                           jazmin=jazim
                        endif
                     enddo
                  endif

                  jmirror=jazmin+(rdat%nazim(kk)/2)
                  if(jmirror > rdat%nazim(kk)) jmirror=jmirror-rdat%nazim(kk)

                  !  Loop forward from jazmin

                  jend=rdat%nazim(kk)
                  if(jmirror > jazmin) jend=jmirror-1
                  jknt=0
                  do jazim=jazmin,jend
                     kinbox=0
                     jknt=jknt+1
                     daz=rdat%razim(jazim,kk)-azimijk
                     if(daz >  180.) daz=daz-360.
                     if(daz < -180.) daz=daz+360.
                     if(abs(daz) > dazlim) EXIT
                     do igate=1,rdat%nvgate(jazim,kk)
                        ddx=rxvol(igate,jazim,kk)-xs
                        ddy=ryvol(igate,jazim,kk)-ys

                        if( rngvol(igate,kk) > rngmin .AND.                    &
                            rngvol(igate,kk) < rngmax .AND.                    &
                            abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then
                           kinbox=kinbox+1
                           if(rdat%vel(igate,jazim,kk) > velcheck ) then
                              isort=1+NINT((rdat%vel(igate,jazim,kk)-sortmin)*sortscale)
                              isort=max(min(isort,nsort),1)
                              kntbin(isort)=kntbin(isort)+1
                              sum=sum+rdat%vel(igate,jazim,kk)
                              sum2=sum2+(rdat%vel(igate,jazim,kk)*rdat%vel(igate,jazim,kk))
                              kok=kok+1
                           endif   ! data ok
                        endif  ! inside box
                     enddo ! igate
                     if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                  enddo ! jazim

                 !  if kinbox > 0 continue from jazim=1

                  if((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==rdat%nazim(kk)) then
                     do jazim=1,jmirror-1
                        kinbox=0
                        jknt=jknt+1
                        daz=rdat%razim(jazim,kk)-azimijk
                        if(daz >  180.) daz=daz-360.
                        if(daz < -180.) daz=daz+360.
                        if(abs(daz) > dazlim) EXIT
                        do igate=1,rdat%nvgate(jazim,kk)

                           ddx=rxvol(igate,jazim,kk)-xs
                           ddy=ryvol(igate,jazim,kk)-ys

                           if( rngvol(igate,kk) > rngmin .AND.                    &
                               rngvol(igate,kk) < rngmax .AND.                    &
                               abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then
                              kinbox=kinbox+1
                              if(rdat%vel(igate,jazim,kk) > velcheck ) then
                                 isort=1+NINT((rdat%vel(igate,jazim,kk)-sortmin)*sortscale)
                                 isort=max(min(isort,nsort),1)
                                 kntbin(isort)=kntbin(isort)+1
                                 sum=sum+rdat%vel(igate,jazim,kk)
                                 sum2=sum2+(rdat%vel(igate,jazim,kk)*rdat%vel(igate,jazim,kk))
                                 kok=kok+1
                              endif
                           endif
                        enddo
                        if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                     enddo
                  endif

                  ! Loop backward from jazmin

                  jend= 1
                  if(jmirror < jazmin) jend=jmirror
                  jknt=0
                  do jazim=jazmin-1,jend,-1
                     kinbox=0
                     jknt=jknt+1
                     daz=rdat%razim(jazim,kk)-azimijk
                     if(daz >  180.) daz=daz-360.
                     if(daz < -180.) daz=daz+360.
                     if(abs(daz) > dazlim) EXIT
                     do igate=1,rdat%nvgate(jazim,kk)

                        ddx=rxvol(igate,jazim,kk)-xs
                        ddy=ryvol(igate,jazim,kk)-ys

                        if( rngvol(igate,kk) > rngmin .AND.                    &
                            rngvol(igate,kk) < rngmax .AND.                    &
                            abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then

                           kinbox=kinbox+1

                           if(rdat%vel(igate,jazim,kk) > velcheck ) then
                              isort=1+NINT((rdat%vel(igate,jazim,kk)-sortmin)*sortscale)
                              isort=max(min(isort,nsort),1)
                              kntbin(isort)=kntbin(isort)+1
                              sum=sum+rdat%vel(igate,jazim,kk)
                              sum2=sum2+(rdat%vel(igate,jazim,kk)*rdat%vel(igate,jazim,kk))
                              kok=kok+1
                           endif
                        endif
                     enddo
                     if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                  enddo

                  ! If not yet outside box, continue from last radial.

                  if((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) then
                     do jazim=rdat%nazim(kk),jmirror,-1
                        kinbox=0
                        jknt=jknt+1
                        daz=rdat%razim(jazim,kk)-azimijk
                        if(daz >  180.) daz=daz-360.
                        if(daz < -180.) daz=daz+360.
                        if(abs(daz) > dazlim) EXIT
                        do igate=1,rdat%nvgate(jazim,kk)

                           ddx=rxvol(igate,jazim,kk)-xs
                           ddy=ryvol(igate,jazim,kk)-ys

                           if( rngvol(igate,kk) > rngmin .AND.                    &
                               rngvol(igate,kk) < rngmax .AND.                    &
                               abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then
                              kinbox=kinbox+1
                              if(rdat%vel(igate,jazim,kk) > velcheck ) then
                                 isort=1+NINT((rdat%vel(igate,jazim,kk)-sortmin)*sortscale)
                                 isort=max(min(isort,nsort),1)
                                 kntbin(isort)=kntbin(isort)+1
                                 sum=sum+rdat%vel(igate,jazim,kk)
                                 sum2=sum2+(rdat%vel(igate,jazim,kk)*rdat%vel(igate,jazim,kk))
                                 kok=kok+1
                              endif
                           endif
                        enddo  ! igate
                        if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                     enddo ! jazim
                  endif
               enddo ! kk

               if( kok == 0 ) then
                  kntflag=kntflag+1
                  write(*,'(A,3I7,A,F10.2)') ' No neighbors near',           &
                         iigate,jjazim,kkelev,' Setting missing. Orig vel=', &
                         origvel
                  rdat%vel(iigate,jjazim,kkelev)=-777.
                  CYCLE
               else
                  varavg=sum/float(kok)
                  mid=(kok/2)+1
                  kntall=0
                  do isort=1,nsort-1
                     kntall=kntall+kntbin(isort)
                     if(kntall >= mid) EXIT
                  enddo
                  varmed=sortmin+((isort-1)*dsort)
                  if ( kok > 1 ) then
                     sdev=sqrt((sum2-(sum*sum/float(kok)))/float(kok-1))
                     thresh=max((2.*sdev),vmedlim)
                  else
                     thresh=vmedlim
                  endif
                  !write(*,*) ' Velocity difference threshold:',thresh
               endif

               !  Process data for local quadratic fit

               do kk=kbgn,kend 

                  !  Find nearest azimuth at this level

                  if(kk == kkelev) then
                     jazmin=jjazim
                  else
                     azdiff=181.
                     jazmin=1
                     do jazim=1,rdat%nazim(kk)
                       daz=rdat%razim(jazim,kk)-azimijk
                       if(daz >  180.) daz=daz-360.
                       if(daz < -180.) daz=daz+360.
                       daz=abs(daz)
                       if(daz < azdiff) then
                          azdiff=daz
                          jazmin=jazim
                       endif
                    enddo
                  endif

                  jmirror=jazmin+(rdat%nazim(kk)/2)
                  if(jmirror > rdat%nazim(kk)) jmirror=jmirror-rdat%nazim(kk)

                  !  Loop forward from jazmin

                  jend=rdat%nazim(kk)
                  if(jmirror > jazmin) jend=jmirror-1
                  jknt=0
                  do jazim=jazmin,jend
                     kinbox=0
                     jknt=jknt+1
                     daz=rdat%razim(jazim,kk)-azimijk
                     if(daz >  180.) daz=daz-360.
                     if(daz < -180.) daz=daz+360.
                     if(abs(daz) > dazlim) EXIT
                     do igate=1,rdat%nvgate(jazim,kk)

                        ddx=rxvol(igate,jazim,kk)-xs
                        ddy=ryvol(igate,jazim,kk)-ys

                        if( rngvol(igate,kk) > rngmin .AND.                    &
                            rngvol(igate,kk) < rngmax .AND.                    &
                            abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then

                           kinbox=kinbox+1
                           ddxy=ddx*ddy
                           ddx2=ddx*ddx
                           ddy2=ddy*ddy

                           if(    rdat%vel(igate,jazim,kk) > velcheck .AND.  &
                              abs(rdat%vel(igate,jazim,kk)-varmed) < thresh ) then

                              varmax=max(varmax,rdat%vel(igate,jazim,kk))
                              varmin=min(varmin,rdat%vel(igate,jazim,kk))

                              rhsvar(1)=rhsvar(1)+rdat%vel(igate,jazim,kk)
                              rhsvar(2)=rhsvar(2)+rdat%vel(igate,jazim,kk)*ddx
                              rhsvar(3)=rhsvar(3)+rdat%vel(igate,jazim,kk)*ddy
                              rhsvar(4)=rhsvar(4)+rdat%vel(igate,jazim,kk)*ddxy
                              rhsvar(5)=rhsvar(5)+rdat%vel(igate,jazim,kk)*ddx2
                              rhsvar(6)=rhsvar(6)+rdat%vel(igate,jazim,kk)*ddy2

                              avar(1,1)=avar(1,1)+1.
                              avar(1,2)=avar(1,2)+ddx
                              avar(1,3)=avar(1,3)+ddy
                              avar(1,4)=avar(1,4)+ddxy
                              avar(1,5)=avar(1,5)+ddx2
                              avar(1,6)=avar(1,6)+ddy2

                              avar(2,1)=avar(2,1)+ddx
                              avar(2,2)=avar(2,2)+ddx2
                              avar(2,3)=avar(2,3)+ddx*ddy
                              avar(2,4)=avar(2,4)+ddx*ddxy
                              avar(2,5)=avar(2,5)+ddx*ddx2
                              avar(2,6)=avar(2,6)+ddx*ddy2

                              avar(3,1)=avar(3,1)+ddy
                              avar(3,2)=avar(3,2)+ddy*ddx
                              avar(3,3)=avar(3,3)+ddy2
                              avar(3,4)=avar(3,4)+ddy*ddx2
                              avar(3,5)=avar(3,5)+ddy*ddx2
                              avar(3,6)=avar(3,6)+ddy*ddy2

                              avar(4,1)=avar(4,1)+ddxy
                              avar(4,2)=avar(4,2)+ddxy*ddx
                              avar(4,3)=avar(4,3)+ddxy*ddy
                              avar(4,4)=avar(4,4)+ddxy*ddxy
                              avar(4,5)=avar(4,5)+ddxy*ddx2
                              avar(4,6)=avar(4,6)+ddxy*ddy2

                              avar(5,1)=avar(5,1)+ddx2
                              avar(5,2)=avar(5,2)+ddx2*ddx
                              avar(5,3)=avar(5,3)+ddx2*ddy
                              avar(5,4)=avar(5,4)+ddx2*ddxy
                              avar(5,5)=avar(5,5)+ddx2*ddx2
                              avar(5,6)=avar(5,6)+ddx2*ddy2

                              avar(6,1)=avar(6,1)+ddy2
                              avar(6,2)=avar(6,2)+ddy2*ddx
                              avar(6,3)=avar(6,3)+ddy2*ddy
                              avar(6,4)=avar(6,4)+ddy2*ddxy
                              avar(6,5)=avar(6,5)+ddy2*ddx2
                              avar(6,6)=avar(6,6)+ddy2*ddy2

                           endif

                        endif
                     enddo  ! igate
                     if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                  enddo ! jazim

                  !  if kinbox > 0 continue from jazim=1

                  if((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==rdat%nazim(kk)) then
                     do jazim=1,jmirror-1
                        kinbox=0
                        jknt=jknt+1
                        daz=rdat%razim(jazim,kk)-azimijk
                        if(daz >  180.) daz=daz-360.
                        if(daz < -180.) daz=daz+360.
                        if(abs(daz) > dazlim) EXIT
                        do igate=1,rdat%nvgate(jazim,kk)
!
                           ddx=rxvol(igate,jazim,kk)-xs
                           ddy=ryvol(igate,jazim,kk)-ys
!
                           if( rngvol(igate,kk) > rngmin .AND.                    &
                               rngvol(igate,kk) < rngmax .AND.                    &
                               abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then

                              kinbox=kinbox+1
                              ddxy=ddx*ddy
                              ddx2=ddx*ddx
                              ddy2=ddy*ddy

                              if(    rdat%vel(igate,jazim,kk) > velcheck .AND.             &
                                 abs(rdat%vel(igate,jazim,kk)-varmed) < thresh ) then
!
                                 varmax=max(varmax,rdat%vel(igate,jazim,kk))
                                 varmin=min(varmin,rdat%vel(igate,jazim,kk))
!
                                 rhsvar(1)=rhsvar(1)+rdat%vel(igate,jazim,kk)
                                 rhsvar(2)=rhsvar(2)+rdat%vel(igate,jazim,kk)*ddx
                                 rhsvar(3)=rhsvar(3)+rdat%vel(igate,jazim,kk)*ddy
                                 rhsvar(4)=rhsvar(4)+rdat%vel(igate,jazim,kk)*ddxy
                                 rhsvar(5)=rhsvar(5)+rdat%vel(igate,jazim,kk)*ddx2
                                 rhsvar(6)=rhsvar(6)+rdat%vel(igate,jazim,kk)*ddy2
!
                                 avar(1,1)=avar(1,1)+1.
                                 avar(1,2)=avar(1,2)+ddx
                                 avar(1,3)=avar(1,3)+ddy
                                 avar(1,4)=avar(1,4)+ddxy
                                 avar(1,5)=avar(1,5)+ddx2
                                 avar(1,6)=avar(1,6)+ddy2
!
                                 avar(2,1)=avar(2,1)+ddx
                                 avar(2,2)=avar(2,2)+ddx2
                                 avar(2,3)=avar(2,3)+ddx*ddy
                                 avar(2,4)=avar(2,4)+ddx*ddxy
                                 avar(2,5)=avar(2,5)+ddx*ddx2
                                 avar(2,6)=avar(2,6)+ddx*ddy2

                                 avar(3,1)=avar(3,1)+ddy
                                 avar(3,2)=avar(3,2)+ddy*ddx
                                 avar(3,3)=avar(3,3)+ddy2
                                 avar(3,4)=avar(3,4)+ddy*ddxy
                                 avar(3,5)=avar(3,5)+ddy*ddx2
                                 avar(3,6)=avar(3,6)+ddy*ddy2
!
                                 avar(4,1)=avar(4,1)+ddxy
                                 avar(4,2)=avar(4,2)+ddxy*ddx
                                 avar(4,3)=avar(4,3)+ddxy*ddy
                                 avar(4,4)=avar(4,4)+ddxy*ddxy
                                 avar(4,5)=avar(4,5)+ddxy*ddx2
                                 avar(4,6)=avar(4,6)+ddxy*ddy2

                                 avar(5,1)=avar(5,1)+ddx2
                                 avar(5,2)=avar(5,2)+ddx2*ddx
                                 avar(5,3)=avar(5,3)+ddx2*ddy
                                 avar(5,4)=avar(5,4)+ddx2*ddxy
                                 avar(5,5)=avar(5,5)+ddx2*ddx2
                                 avar(5,6)=avar(5,6)+ddx2*ddy2

                                 avar(6,1)=avar(6,1)+ddy2
                                 avar(6,2)=avar(6,2)+ddy2*ddx
                                 avar(6,3)=avar(6,3)+ddy2*ddy
                                 avar(6,4)=avar(6,4)+ddy2*ddxy
                                 avar(6,5)=avar(6,5)+ddy2*ddx2
                                 avar(6,6)=avar(6,6)+ddy2*ddy2
!
                              endif

                           endif
                        enddo  ! igate
                        if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                     enddo ! jazim
                  endif

                  ! Loop backward from jazmin

                  jend= 1
                  if(jmirror < jazmin) jend=jmirror
                  jknt=0
                  do jazim=jazmin-1,jend,-1
                     kinbox=0
                     jknt=jknt+1
                     daz=rdat%razim(jazim,kk)-azimijk
                     if(daz >  180.) daz=daz-360.
                     if(daz < -180.) daz=daz+360.
                     if(abs(daz) > dazlim) EXIT
                     do igate=1,rdat%nvgate(jazim,kk)
!
                        ddx=rxvol(igate,jazim,kk)-xs
                        ddy=ryvol(igate,jazim,kk)-ys
!
                        if( rngvol(igate,kk) > rngmin .AND.                    &
                            rngvol(igate,kk) < rngmax .AND.                    &
                            abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then

                           kinbox=kinbox+1
                           ddxy=ddx*ddy
                           ddx2=ddx*ddx
                           ddy2=ddy*ddy

                           if(    rdat%vel(igate,jazim,kk) > velcheck .AND.             &
                              abs(rdat%vel(igate,jazim,kk)-varmed) < thresh ) then
!
                              varmax=max(varmax,rdat%vel(igate,jazim,kk))
                              varmin=min(varmin,rdat%vel(igate,jazim,kk))
!
                              rhsvar(1)=rhsvar(1)+rdat%vel(igate,jazim,kk)
                              rhsvar(2)=rhsvar(2)+rdat%vel(igate,jazim,kk)*ddx
                              rhsvar(3)=rhsvar(3)+rdat%vel(igate,jazim,kk)*ddy
                              rhsvar(4)=rhsvar(4)+rdat%vel(igate,jazim,kk)*ddxy
                              rhsvar(5)=rhsvar(5)+rdat%vel(igate,jazim,kk)*ddx2
                              rhsvar(6)=rhsvar(6)+rdat%vel(igate,jazim,kk)*ddy2
!
                              avar(1,1)=avar(1,1)+1.
                              avar(1,2)=avar(1,2)+ddx
                              avar(1,3)=avar(1,3)+ddy
                              avar(1,4)=avar(1,4)+ddxy
                              avar(1,5)=avar(1,5)+ddx2
                              avar(1,6)=avar(1,6)+ddy2
!
                              avar(2,1)=avar(2,1)+ddx
                              avar(2,2)=avar(2,2)+ddx2
                              avar(2,3)=avar(2,3)+ddx*ddy
                              avar(2,4)=avar(2,4)+ddx*ddxy
                              avar(2,5)=avar(2,5)+ddx*ddx2
                              avar(2,6)=avar(2,6)+ddx*ddy2

                              avar(3,1)=avar(3,1)+ddy
                              avar(3,2)=avar(3,2)+ddy*ddx
                              avar(3,3)=avar(3,3)+ddy2
                              avar(3,4)=avar(3,4)+ddy*ddxy
                              avar(3,5)=avar(3,5)+ddy*ddx2
                              avar(3,6)=avar(3,6)+ddy*ddy2

                              avar(4,1)=avar(4,1)+ddxy
                              avar(4,2)=avar(4,2)+ddxy*ddx
                              avar(4,3)=avar(4,3)+ddxy*ddy
                              avar(4,4)=avar(4,4)+ddxy*ddxy
                              avar(4,5)=avar(4,5)+ddxy*ddx2
                              avar(4,6)=avar(4,6)+ddxy*ddy2
!
                              avar(5,1)=avar(5,1)+ddx2
                              avar(5,2)=avar(5,2)+ddx2*ddx
                              avar(5,3)=avar(5,3)+ddx2*ddy
                              avar(5,4)=avar(5,4)+ddx2*ddxy
                              avar(5,5)=avar(5,5)+ddx2*ddx2
                              avar(5,6)=avar(5,6)+ddx2*ddy2
!
                              avar(6,1)=avar(6,1)+ddy2
                              avar(6,2)=avar(6,2)+ddy2*ddx
                              avar(6,3)=avar(6,3)+ddy2*ddy
                              avar(6,4)=avar(6,4)+ddy2*ddxy
                              avar(6,5)=avar(6,5)+ddy2*ddx2
                              avar(6,6)=avar(6,6)+ddy2*ddy2
!
                           endif
!
                        endif
                     enddo  ! igate
                     if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                  enddo ! jazim

                  ! If not yet outside box, continue from last radial.

                  if((kinbox > 0 .OR. jknt <= jsrchmn) .and. jend==1 ) then
                     do jazim=rdat%nazim(kk),jmirror,-1
                        kinbox=0
                        jknt=jknt+1
                        daz=rdat%razim(jazim,kk)-azimijk
                        if(daz >  180.) daz=daz-360.
                        if(daz < -180.) daz=daz+360.
                        if(abs(daz) > dazlim) EXIT
                        do igate=1,rdat%nvgate(jazim,kk)
!
                           ddx=rxvol(igate,jazim,kk)-xs
                           ddy=ryvol(igate,jazim,kk)-ys
!
                           if( rngvol(igate,kk) > rngmin .AND.                    &
                               rngvol(igate,kk) < rngmax .AND.                    &
                               abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) then

                              kinbox=kinbox+1
                              ddxy=ddx*ddy
                              ddx2=ddx*ddx
                              ddy2=ddy*ddy

                              if(    rdat%vel(igate,jazim,kk) > velcheck .AND.             &
                                 abs(rdat%vel(igate,jazim,kk)-varmed) < thresh ) then
!
                                 varmax=max(varmax,rdat%vel(igate,jazim,kk))
                                 varmin=min(varmin,rdat%vel(igate,jazim,kk))
!
                                 rhsvar(1)=rhsvar(1)+rdat%vel(igate,jazim,kk)
                                 rhsvar(2)=rhsvar(2)+rdat%vel(igate,jazim,kk)*ddx
                                 rhsvar(3)=rhsvar(3)+rdat%vel(igate,jazim,kk)*ddy
                                 rhsvar(4)=rhsvar(4)+rdat%vel(igate,jazim,kk)*ddxy
                                 rhsvar(5)=rhsvar(5)+rdat%vel(igate,jazim,kk)*ddx2
                                 rhsvar(6)=rhsvar(6)+rdat%vel(igate,jazim,kk)*ddy2
!
                                 avar(1,1)=avar(1,1)+1.
                                 avar(1,2)=avar(1,2)+ddx
                                 avar(1,3)=avar(1,3)+ddy
                                 avar(1,4)=avar(1,4)+ddxy
                                 avar(1,5)=avar(1,5)+ddx2
                                 avar(1,6)=avar(1,6)+ddy2
!
                                 avar(2,1)=avar(2,1)+ddx
                                 avar(2,2)=avar(2,2)+ddx2
                                 avar(2,3)=avar(2,3)+ddx*ddy
                                 avar(2,4)=avar(2,4)+ddx*ddxy
                                 avar(2,5)=avar(2,5)+ddx*ddx2
                                 avar(2,6)=avar(2,6)+ddx*ddy2
!
                                 avar(3,1)=avar(3,1)+ddy
                                 avar(3,2)=avar(3,2)+ddy*ddx
                                 avar(3,3)=avar(3,3)+ddy2
                                 avar(3,4)=avar(3,4)+ddy*ddxy
                                 avar(3,5)=avar(3,5)+ddy*ddx2
                                 avar(3,6)=avar(3,6)+ddy*ddy2
!
                                 avar(4,1)=avar(4,1)+ddx2
                                 avar(4,2)=avar(4,2)+ddx2*ddx
                                 avar(4,3)=avar(4,3)+ddx2*ddy
                                 avar(4,4)=avar(4,4)+ddx2*ddxy
                                 avar(4,5)=avar(4,5)+ddx2*ddx2
                                 avar(4,6)=avar(4,6)+ddx2*ddy2
!
                                 avar(5,1)=avar(5,1)+ddx2
                                 avar(5,2)=avar(5,2)+ddx2*ddx
                                 avar(5,3)=avar(5,3)+ddx2*ddy
                                 avar(5,4)=avar(5,4)+ddx2*ddxy
                                 avar(5,5)=avar(5,5)+ddx2*ddx2
                                 avar(5,6)=avar(5,6)+ddx2*ddy2
!
                                 avar(6,1)=avar(6,1)+ddy2
                                 avar(6,2)=avar(6,2)+ddy2*ddx
                                 avar(6,3)=avar(6,3)+ddy2*ddy
                                 avar(6,4)=avar(6,4)+ddy2*ddxy
                                 avar(6,5)=avar(6,5)+ddy2*ddx2
                                 avar(6,6)=avar(6,6)+ddy2*ddy2
!
                              endif

                           endif
                        enddo  ! igate
                        if(kinbox == 0 .AND. jknt > jsrchmn) EXIT
                     enddo ! jazim
                  endif

               enddo ! kk

               !   Solve for variable at grid point

               knt=nint(avar(1,1))
               if ( iorder > 1 .and. knt > 6 ) then
                  varmean=rhsvar(1)/avar(1,1)
                  CALL GJELIM(n,avar,rhsvar,sol,work,work1d,eps,istatus)
                  fitvel=min(varmax,max(varmin,sol(1)))
                  write(*,'(3F7.1,A,I8,4F7.1)') rdat%razim(jjazim,kkelev),       &
                       (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),        &
                       ' Qf Analysis1:',knt,  &
                       varmin,varmean,varmax,fitvel
               elseif ( iorder > 0 .and. knt > 5 ) then
                  do jj=1,4
                     do ii=1,4
                        array(ii,jj)=avar(ii,jj)
                     enddo
                  enddo
                  do ii=1,4
                     rhsv(ii)=rhsvar(ii)
                  enddo
                  CALL GJELIM(4,array,rhsv,solv,work,work1d,eps,istatus)
                  fitvel=min(varmax,max(varmin,solv(1)))
                  write(*,'(3F7.1,A,I8,4F7.1)') rdat%razim(jjazim,kkelev),      &
                       (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),       &
                       ' Qf Analysis2:',knt,  &
                       varmin,varmean,varmax,fitvel
               elseif ( knt > 0 ) then
                  varmean=rhsvar(1)/avar(1,1)
                  fitvel=varmean
                  write(*,'(3F7.1,A,I8,4F7.1)') rdat%razim(jjazim,kkelev),      &
                       (0.001*rngvol(iigate,kkelev)),elvmnvol(kkelev),       &
                       ' Qf Analysis3:',knt,  &
                       varmin,varmean,varmax,fitvel
               endif
               tstdev=twonyq*NINT((fitvel-origvel)*inv2nyq)
               tstvel=origvel+tstdev
               vardif=abs(tstvel-fitvel)

               if(abs(vardif) < thresh) then
                  kntfold=kntfold+1
                  rdat%vel(iigate,jjazim,kkelev) = tstvel
                  write(*,'(A,F10.1,A,F10.1,A,F10.1)')            &
                          ' Qf Unfold: Meas=',origvel,' Qfitvel=',     &
                                               fitvel,'  New=',tstvel
               else
                  kntflag=kntflag+1
                  rdat%vel(iigate,jjazim,kkelev) = -777.0
                  write(*,'(A,F10.1,A,F10.1,A)')               &
                          ' Qf Marked bad: Meas=',origvel,' Qfitvel=', &
                                                   fitvel,'  New=-777.0'
               endif
            endif
         enddo
      enddo

      write(*,'(/A/,A/,3I12)')                                     &
              '       After Least-Squares Quad Fit Folding Check',   &
              '       Gates Checked   Unfolded   Flagged: ',         &
                      kntchk,kntfold,kntflag
   enddo  ! level loop

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine unfold_radar(rdat, bkgvel)
   implicit none
   type(t_radar_data), intent(inout) :: rdat
   real, dimension(:,:,:), optional :: bkgvel
   
   integer :: maxgate, ngate, nazim

!  real,    allocatable :: bkgvel(:,:), tmp1(:,:), unfvel(:,:)
   real,    allocatable :: tmp1(:,:), unfvel(:,:), tmp2(:,:)
   integer, allocatable :: bgate(:),    egate(:)

   integer :: i,iray,jray,kray,eray,ipass
   integer :: js1,js2,jstart
   integer :: igate,kgate,kkgate
   integer :: bkgate,ekgate,lgate,kgate1,kgate2
   integer :: k,kh,kclose,kstart,kend
   integer :: knt,knt1,knt2,kntall,kntgood,kntfold,kntrej
   integer :: istatwrt,rem_points,found,nhlfazim
   integer :: igmin,igmax,irayd
   real    :: bknt,fknt,tknt,tpoint
   real    :: twonyq,inv2nyq,thrpri,thrpr2
   real    :: dbkvel,sum,vavg,rvwgt,tstdev,tstvel,sum1,vavg1,range
   real    :: uvel,vvel,bkrvel
   real    :: hlfbeam,elvtop,elvbot,hgtagl,hgttop,hgtbot,hgtmsl
   real    :: refvel,elevavg,sfcr
   real    :: dtr,umean,vmean,veloc,xcomp,ycomp,dotp,dotmin
   real    :: diff,shdiff,shdif2,shdif3,azdiff,azd,azd1,azd2,azd3,azim2
   real    :: rngfrst,gatespc
   real    :: hgtmin,hgtmax
   real    :: vnyq, radalt, lat1, lon1
   integer :: irfirst,igatesp,kelv

   integer, dimension(2) :: js
   integer :: loop

   real,    parameter :: beamwid = 1.0
   real,    parameter :: bkgwgt  = 0.33
   integer, parameter :: lukbak  = 5
   integer, parameter :: lukbak2 = 10
   integer, parameter :: lukfwd2 = 5
   integer, parameter :: lukbakr = 5
   
   real               :: validratio, absvmean, wmax
   real, dimension(:), allocatable :: rvmean

   integer, dimension(3) :: dims

   write(*,*) "Unfold and QC Vel..."

   dims     = ubound(rdat%vel)
   maxazim  = dims(2)
   maxelev  = rdat%ntilt
   maxvgate = dims(1)

   maxgate=maxvgate


!  allocate(bkgvel(maxgate,maxazim))
   allocate(unfvel(maxgate,maxazim))
   allocate(tmp1  (maxgate,maxazim))
!   allocate(tmp2  (maxgate,maxazim))
   allocate(bgate (maxazim))
   allocate(egate (maxazim))
   allocate(rvmean(maxazim))

   if(present(bkgvel))then
      bkgopt=1
   endif

   hlfbeam  = 0.5*beamwid
   rvwgt    = 1.0-bkgwgt
   dtr      = atan(1.)/45.
   radalt   = rdat%altitude

   do kelv=1, rdat%ntilt
 
      if(.NOT.rdat%ifvel(kelv)) cycle
      nazim    = rdat%nazim    (kelv)
      irfirst  = rdat%vgatesp  (kelv)
      igatesp  = rdat%vgatesp  (kelv)
      ngate    = rdat%nvgate (1,kelv)

      kntall   = 0
      kntgood  = 0
      kntfold  = 0
      kntrej   = 0
      rngfrst  = float(irfirst)
      gatespc  = float(igatesp)
      nhlfazim = (nazim/2)+1

      ! Set up tmp1 to be a quality array. 1=good 0=bad/missing
      igmin=ngate/2
      igmax=igmin+1
   
      do iray=1,maxazim
         bgate(iray)=1
         egate(iray)=1
      enddo
      do iray=1,maxazim
         do igate=1,maxgate
            tmp1  (igate,iray)=0.
         enddo
      enddo
   
      do iray=1,nazim
         bgate(iray)=0
         do igate=1,ngate
            unfvel(igate,iray)=rdat%vel(igate,iray,kelv)
            if(abs(unfvel(igate,iray)) < abs(velcheck)) then
               tmp1 (igate,iray)=1.
               egate(iray)=igate
               if(bgate(iray) == 0) bgate(iray)=igate
            endif
         enddo
         bgate(iray)=max(bgate(iray),1)
         igmin=min(igmin,bgate(iray))
         igmax=max(igmax,egate(iray))
      enddo
!      tmp2=tmp1
!      !tmp1=0.
!      do iray=1,nazim
!         do igate=bgate(iray),egate(iray)
!            if(tmp1(igate,iray)<1.and.iray>1) then
!write(605,*) igate,iray,tmp1(igate,iray),tmp2(igate,iray),unfvel(igate-1,iray),unfvel(igate,iray)
!               unfvel(igate,iray)=unfvel(igate-1,iray)
!               rdat%vel(igate,iray,kelv)=unfvel(igate,iray)
!               tmp1 (igate,iray)=1.
!            endif
!         enddo
!      enddo
   
      ! Determine jstart
      if (bkgopt>0) then
         elevavg=0.
         do iray=1,nazim
            elevavg=elevavg+rdat%rtilt(iray,kelv)
         enddo
         elevavg=elevavg/float(nazim)
         write(*,'(a,f10.2)') ' Average Elevation Angle: ',elevavg

         ! Determine mean environmental wind in layer observed by radar.
         ! Used to determine starting radial for unfolding.

         range=rngfrst+(igmin-1)*gatespc
         CALL beamhgt(elevavg,range,hgtmin,sfcr)
                       
         hgtmin=hgtmin+radalt
   
         !call get_wrf_sound(rdat%latitude,rdat%longitude,range,nzsnd,zsnd,usnd,vsnd)
         range=rngfrst+(igmax-1)*gatespc
         CALL beamhgt(elevavg,range,hgtmax,sfcr)
                      
         hgtmax=hgtmax+radalt
   
         write(*,'(A,F12.2,F12.2)') ' Height range (MSL): ',hgtmin,hgtmax
   
         rvmean=0.
         do iray=1, nazim
            knt=0
            do igate=1, maxgate
!              write(*,*) igate,iray,kelv, size(bkgvel)
               if(bkgvel(igate,iray,kelv)>value_invalid)then
                  knt=knt+1
                  rvmean(iray)=rvmean(iray)+abs(bkgvel(igate,iray,kelv))
               endif
            enddo
            knt=max(knt,1)
!           write(*,*) iray, size(rvmean)
!stop
            rvmean(iray)=rvmean(iray)/float(knt)
!           write(*,'(A,F12.2,A,2I6)') ' Mean velocity, rvmean: ',rvmean(iray),' iray,kelv: ',iray,kelv
         enddo
         
   
         ! Find observed radial with direction most perpendicular to the mean wind.
         ! Smallest dot product between mean wind and azimuth vector.
   
         dotmin=99999.
         js1=1
         do iray=1,nazim
            dotp=abs(rvmean(iray))
            if( dotp < dotmin) then
               dotmin=dotp
               js1=iray
            endif
         enddo
   
         write(*,'(A,I6,F12.2)') ' Most perpendicular radial: ',js1,rdat%razim(js1,kelv)
          
     
         ! Next find the radial opposite of this (+/-180)
     
         azim2=rdat%razim(js1,kelv)+180.
         if(azim2 > 360.) azim2=azim2-360.
         azdiff=90.
         js2=js1
         dotmin=99999.
         do iray=1,nazim
            azd1=abs(rdat%razim(iray,kelv)-azim2)
            azd2=abs(rdat%razim(iray,kelv)-azim2+360.)
            azd3=abs(rdat%razim(iray,kelv)-azim2-360.)
            azd=min(azd1,azd2,azd3)
            dotp=abs(rvmean(iray))
            if( azd < azdiff .and. dotp < dotmin) then
               dotmin=dotp
               js2=iray
            endif
         enddo
   
         write(*,'(A,I6,F12.2)') ' Opposite radial: ',js2, rdat%razim(js2,kelv)
          
     
         ! Find which zone around js1 or js2 has the most valid data.
         ! This is the best place to start
         
         knt1=0
         knt2=0
         do iray=1,nazim
            irayd=abs(iray-js1)
            if ( irayd < 6 ) then
               do igate=bgate(iray),egate(iray)
                  if(tmp1(igate,iray) > 0.) knt1=knt1+1
               enddo
            endif
            irayd=abs(iray-js2)
            if ( irayd < 6 ) then
               do igate=bgate(iray),egate(iray)
                  if(tmp1(igate,iray) > 0.) knt2=knt2+1
               enddo
            endif
         enddo
   
         write(*,'(A,I6,A,I6)') ' Data count js1: ',knt1,'  Data count js2: ',knt2 
           
   
         !if( knt2 > knt1 ) then
         !   jstart=js2
         !   js(1)=js2
         !   js(2)=js1
         !else
            jstart=js1
            js(1)=js1
            js(2)=js2
         !endif
      else
         js1=1
         js2=180
         js(1)=js1
         js(2)=js2
         jstart=1
         wmax=0
         do iray=1,nazim
            validratio=0
            absvmean  =0
            do igate=bgate(iray), egate(iray)
               if(tmp1(igate,iray) > 0.)then
                  validratio=validratio+1
                  absvmean  =absvmean  +abs(unfvel(igate,iray))
               endif
            enddo
            if(validratio>0)then
               absvmean=absvmean/validratio/rdat%vmax(kelv)
            endif
            validratio=validratio/ngate
            if(wmax<(validratio-absvmean) )then
               wmax=validratio-absvmean
               jstart=iray
            endif
         enddo
         write(*,'(A,I6,3F12.2)') ' Shear Check start at radial: ',jstart,rdat%razim(jstart,kelv), validratio, absvmean
      endif
      ! End of Determine jstart

      ! Determine background radial velocity at each valid data point.
      !if (bkgopt>0)then
      !   do iray=1,nazim
      !      do igate=bgate(iray),egate(iray)
      !         if( tmp1(igate,iray) > 0. ) then
      !            range=rngfrst+(igate-1)*gatespc
      !            CALL beamhgt(rdat%rtilt(iray,kelv),range,hgtagl,sfcr)
      !                          
      !            !call gcircle(rdat%latitude,rdat%longitude,rdat%razim(iray,kelv),sfcr,lat1,lon1)
      !            !call get_wrf_sound(lat1,lon1,range,nzsnd,zsnd,usnd,vsnd)

      !            hgtmsl=hgtagl+radalt
      !            elvtop=rdat%rtilt(iray,kelv)+hlfbeam
      !            CALL beamhgt(elvtop,range,hgtagl,sfcr)
      !                          
      !            hgttop=hgtagl+radalt
      !            hgttop=max(hgttop,(hgtmsl+200.))
      !            elvbot=rdat%rtilt(iray,kelv)-hlfbeam
      !            CALL beamhgt(elvbot,range,hgtagl,sfcr)
      !                          
      !            hgtbot=hgtagl+radalt
      !            hgtbot=min(hgtbot,(hgtmsl-200.))
      !            !knt=0
      !            !uvel=0.
      !            !vvel=0.
      !            !do k=1,nzsnd
      !            !   if( zsnd(k) > hgtbot .AND. zsnd(k) < hgttop ) then
      !            !      knt=knt+1
      !            !      uvel=uvel+usnd(k)
      !            !      vvel=vvel+vsnd(k)
      !            !   endif
      !            !   if (zsnd(k) > hgttop) EXIT
      !            !enddo
      !            !if(knt > 0) then
      !            !   uvel=uvel/float(knt)
      !            !   vvel=vvel/float(knt)
      !            !   CALL uv2vr(range,rdat%rtilt(iray,kelv),rdat%razim(iray,kelv),                       &
      !            !              uvel,vvel,bkgvel(igate,iray))
      !            !else
      !            !   diff = 99999.0
      !            !   do k=1,nzsnd
      !            !      if( abs(zsnd(k)-hgtmsl) < diff ) then
      !            !         diff = abs(zsnd(k)-hgtmsl)
      !            !         kclose = k
      !            !      endif
      !            !   enddo
      !            !   CALL uv2vr(range,rdat%rtilt(iray,kelv),rdat%razim(iray,kelv),                       &
      !            !              usnd(kclose),vsnd(kclose),bkgvel(igate,iray))
      !            !endif
      !         endif
      !      enddo
      !   enddo
      !endif ! bkgopt > 0
   

      ! ! Checking of data vs background radial velocity

      if( bkgopt > 0 ) then
         write(*,'(A,I4)') ' Unfolding using background wind profile, bkgopt=',bkgopt
           
         do iray=1,nazim
           vnyq    = rdat%vmax(kelv)
           twonyq  = 2.0*vnyq
           inv2nyq = 1./twonyq
           thrpri  = min(max((0.5*vnyq),10.),vnyq)
           thrpr2  = min(vnyq,(1.5*thrpri))

           do igate=bgate(iray),egate(iray)
              if( tmp1(igate,iray) > 0. ) then
                 kntall=kntall+1
                 if(abs(bkgvel(igate,iray,kelv))<200)then
                    tstdev=twonyq*NINT((bkgvel(igate,iray,kelv)-rdat%vel(igate,iray,kelv))*inv2nyq)
                    tstvel=rdat%vel(igate,iray,kelv)+tstdev
                    if(abs(tstdev) > 0.) then
                       kntfold=kntfold+1
                       write(600,"(A,3F12.2)") "bkg:",rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),tstvel
                    else
                       kntgood = kntgood+1
                    endif
                    unfvel(igate,iray)=tstvel
                    if(abs(unfvel(igate,iray))>100)then
                       write(701,*) "unfold check1",tstdev,tstvel,rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),vnyq
                    endif
                 else
                    bkgvel(igate,iray,kelv)=rdat%vel(igate,iray,kelv)
                 endif
              endif
           enddo
         enddo
         write(*,'(/A/,A/,3i12)')                                            &
               '       Comparison to background wind profile:',              &
               '       Gates  Good-Gates     Folded',                        &
                      kntall,kntgood,kntfold
      endif  ! bkgopt > 0

do loop=1, 2
   jstart=js(loop)
   
   
      if (shropt > 0 ) then
         do ipass = 1,2
            write(*,'(A,I4,A,I3)')                                              &
              ' Unfolding using local shear check, shropt=',shropt,             &
              '  ipass=',ipass
            write(*,'(A,I6,F12.2)')                                             &
              ' Shear-based unfolding will start at: ',jstart,rdat%razim(jstart,kelv)
            
            !   First subpass, increasing azimuth index
            
            !   Apply gate-to-gate shear checking for velocity folding.
            !   Loop forward from jstart

            kntall = 0
            kntgood = 0
            kntfold = 0

            do i=1,nazim
               iray=(jstart+i)-1
               if(iray > nazim) iray=((jstart+i)-1)-nazim
               
               vnyq    = rdat%vmax(kelv)
               twonyq  = 2.0*vnyq
               inv2nyq = 1./twonyq
               thrpri  = min(max((0.5*vnyq),10.),vnyq)
               thrpr2  = min(vnyq,(1.5*thrpri))
               
               do igate=bgate(iray),egate(iray)
                  if( tmp1(igate,iray) > 0. ) then
                     dbkvel=unfvel(igate,iray)-bkgvel(igate,iray,kelv)
                     kntall=kntall+1
                     lgate = 0
                     if(igate > bgate(iray)) then
                        ekgate = max(bgate(iray),igate-lukbak)
                        do kgate = igate-1,ekgate,-1
                           if(tmp1(kgate,iray) >0. ) then
                              lgate = kgate
                              EXIT
                           endif
                        enddo
                     endif
                     if(lgate >0) then
                        refvel = unfvel(lgate,iray)-bkgvel(lgate,iray,kelv)
                     else
                        refvel = 999.
                     endif
                     shdiff=abs(dbkvel-refvel)
                     if(shdiff > thrpri) then
                        sum=0.
                        bknt=0.
                        if(igate > bgate(iray)) then
                           ekgate=min(max((igate-lukbak2),bgate(iray)),(igate-1))
                           do kgate=igate-1,ekgate,-1
                              bknt=bknt+tmp1(kgate,iray)
                              sum=sum+tmp1(kgate,iray)*                              &
                                  (unfvel(kgate,iray)-bkgvel(kgate,iray,kelv))
                              if (bknt > 2. ) EXIT
                           enddo
                        endif
                        fknt=0.
                        eray=max((iray-lukbakr),1)
                        do kray=iray-1,eray,-1
                           ekgate=min((igate+lukfwd2),egate(kray))
                           do kgate=igate,ekgate
                              fknt=fknt+tmp1(kgate,kray)
                              sum=sum+tmp1(kgate,kray)*                             &
                                  (unfvel(kgate,kray)-bkgvel(kgate,kray,kelv))
                              if ( fknt > 4. ) EXIT
                           enddo
                        enddo
                        tknt=bknt+fknt
                        if( tknt > 2. ) then
                           vavg=(rvwgt*sum/tknt)+bkgvel(igate,iray,kelv)
                           shdif2=abs(unfvel(igate,iray)-vavg)
                           if( shdif2 < thrpr2 ) then
                              kntgood=kntgood+1
                           else
                              tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                              tstvel=unfvel(igate,iray)+tstdev
                              shdif3=abs(tstvel-vavg)
                              if( shdif3 < thrpr2 .and.abs(shdif3)<abs(shdiff)) then
                                 if( abs(tstdev) > 0. ) then
                                    write(600,"(A,3I5,4F12.2)") "shr:",igate,iray,kelv,rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),tstvel,vavg
                                    kntfold=kntfold+1
                 if(abs(tstvel)>2*vnyq)then
                    write(701,*) "unfold check2",tstdev,tstvel,vavg,rdat%vel(igate,iray,kelv),unfvel(igate,iray),vnyq
                 endif
                                    unfvel(igate,iray)=tstvel

                                 else
                                    kntgood=kntgood+1
                                 endif
                              endif
                           endif               ! unfvel-vavg >thresh
                        endif
                     else
                       kntgood=kntgood+1
                     endif                !(unfvel(igate,iray) - refvel) > thrpri
                  endif                  ! tmp1(igate,iray) > 0.
               enddo
            enddo
            write(*,'(/a/,a/,3i12)')                                              &
                  '       After forward pass shear consistency check',            &
                  '       Gates  Passed Check   Folded',                          &
                          kntall,kntgood,kntfold
   
   
            !   Apply gate-to-gate shear checking for velocity folding.
            !   Loop forward from jstart. No background used

            kntall = 0
            kntgood = 0
            kntfold = 0

            do i=1,nazim
               iray=(jstart+i)-1
               if(iray > nazim) iray=((jstart+i)-1)-nazim
               
               vnyq    = rdat%vmax(kelv)
               twonyq  = 2.0*vnyq
               inv2nyq = 1./twonyq
               thrpri  = min(max((0.5*vnyq),10.),vnyq)
               thrpr2  = min(vnyq,(1.5*thrpri))
               
               do igate=bgate(iray),egate(iray)
                  if( tmp1(igate,iray) > 0. ) then
                     dbkvel=unfvel(igate,iray)
                     kntall=kntall+1
                     lgate = 0
                     if(igate > bgate(iray)) then
                        ekgate = max(bgate(iray),igate-lukbak)
                        do kgate = igate-1,ekgate,-1
                           if(tmp1(kgate,iray) >0. ) then
                              lgate = kgate
                              EXIT
                           endif
                        enddo
                     endif
                     if(lgate >0) then
                        refvel = unfvel(lgate,iray)
                     else
                        refvel = 999.
                     endif
                     shdiff=abs(dbkvel-refvel)
                     if(shdiff > thrpri) then
                        sum=0.
                        bknt=0.
                        if(igate > bgate(iray)) then
                           ekgate=min(max((igate-lukbak2),bgate(iray)),(igate-1))
                           do kgate=igate-1,ekgate,-1
                              bknt=bknt+tmp1(kgate,iray)
                              sum=sum+tmp1(kgate,iray)*                              &
                                  (unfvel(kgate,iray))
                              if (bknt > 2. ) EXIT
                           enddo
                        endif
                        fknt=0.
                        eray=max((iray-lukbakr),1)
                        do kray=iray-1,eray,-1
                           ekgate=min((igate+lukfwd2),egate(kray))
                           do kgate=igate,ekgate
                              fknt=fknt+tmp1(kgate,kray)
                              sum=sum+tmp1(kgate,kray)*                             &
                                  (unfvel(kgate,kray))
                              if ( fknt > 4. ) EXIT
                           enddo
                        enddo
                        tknt=bknt+fknt
                        if( tknt > 2. ) then
                           vavg=(sum/tknt)
                           shdif2=abs(unfvel(igate,iray)-vavg)
                           if( shdif2 < thrpr2 ) then
                              kntgood=kntgood+1
                           else
                              tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                              tstvel=unfvel(igate,iray)+tstdev
                              shdif3=abs(tstvel-vavg)
                              if( shdif3 < thrpr2 .and.abs(shdif3)<abs(shdiff)) then
                                 if( abs(tstdev) > 0. ) then
                                    write(600,"(A,3I5,3F12.2)") "nobkg shr:",igate,iray,kelv,rdat%vel(igate,iray,kelv),tstvel,vavg
                                    kntfold=kntfold+1
                 if(abs(tstvel)>2*vnyq)then
                    write(701,*) "unfold check3",tstdev,tstvel,vavg,rdat%vel(igate,iray,kelv),unfvel(igate,iray),vnyq
                 endif

                                    unfvel(igate,iray)=tstvel
                                 else
                                    kntgood=kntgood+1
                                 endif
                              endif
                           endif               ! unfvel-vavg >thresh
                        endif
                     else
                       kntgood=kntgood+1
                     endif                !(unfvel(igate,iray) - refvel) > thrpri
                  endif                  ! tmp1(igate,iray) > 0.
               enddo
            enddo
            write(*,'(/a/,a/,3i12)')                                              &
                  '       After forward pass shear consistency check(nobkg)',            &
                  '       Gates  Passed Check   Folded',                          &
                          kntall,kntgood,kntfold
   
            !   Apply gate-to-gate shear checking for velocity folding.
            !   Loop backward from jstart-1

            kntall = 0
            kntgood = 0
            kntfold = 0
            kntrej = 0
            do i=1,nazim
               iray=(jstart-i)+1
               if(iray <= 0) iray=(jstart-i)+nazim
               
               vnyq = rdat%vmax(kelv)
               twonyq = 2.0*vnyq
               inv2nyq = 1./twonyq
               thrpri = min(max((0.5*vnyq),10.),vnyq)
               thrpr2 = min(vnyq,(1.5*thrpri))
               
               do igate=bgate(iray),egate(iray)
                  if( tmp1(igate,iray) > 0. ) then
                     dbkvel=unfvel(igate,iray)-bkgvel(igate,iray,kelv)
                     kntall=kntall+1
                     lgate = 0
                     if(igate > bgate(iray)) then
                        ekgate = max((igate-lukbak),bgate(iray))
                        do kgate = igate-1,ekgate,-1
                           if(tmp1(kgate,iray) >0. ) then
                              lgate = kgate
                              EXIT
                           endif
                        enddo
                     endif
                     if(lgate > 0) then
                        refvel = unfvel(lgate,iray)-bkgvel(lgate,iray,kelv)
                     else
                        refvel = 999.
                     endif
                     shdiff=abs(dbkvel-refvel)
                     if( shdiff > thrpri) then
                        sum=0.
                        fknt=0.
                        bknt=0.
                        if(igate > bgate(iray)) then
                           ekgate=max((igate-lukbak2),bgate(iray))
                           do kgate=igate-1,ekgate,-1
                              bknt=bknt+tmp1(kgate,iray)
                              sum=sum+tmp1(kgate,iray)*                                &
                                  (unfvel(kgate,iray)-bkgvel(kgate,iray,kelv))
                              if (bknt > 2. ) EXIT
                           enddo
                        endif
                        eray=min((iray+lukbakr),nazim)
                        do kray=iray+1,eray
                           ekgate=min((igate+lukfwd2),egate(kray))
                           do kgate=igate,ekgate
                              fknt=fknt+tmp1(kgate,kray)
                              sum=sum+tmp1(kgate,kray)*                                &
                                  (unfvel(kgate,kray)-bkgvel(kgate,kray,kelv))
                              if ( fknt > 4. ) EXIT
                           enddo
                        enddo
   
                        tknt=bknt+fknt
                        if(tknt > 2. ) then
                           vavg=(rvwgt*sum/tknt)+bkgvel(igate,iray,kelv)
                           shdif2=abs(unfvel(igate,iray)-vavg)
                           if( shdif2 < thrpr2 ) then
                              kntgood=kntgood+1
                           else
                              tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                              tstvel=unfvel(igate,iray)+tstdev
                              shdif3=(tstvel-vavg)
                              if( shdif3 < thrpr2 .and.abs(shdif3)<abs(shdiff)) then
                                 if( abs(tstdev) > 0. ) then
                                    write(600,"(A,3I5,4F12.2)") "shr:",igate,iray,kelv,rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),tstvel,vavg
                                    !write(600,"(A,4F12.2)") "shr:",rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),tstvel,vavg
                                    kntfold=kntfold+1
                 if(abs(tstvel)>2*vnyq)then
                    write(701,*) "unfold check4",tstdev,tstvel,vavg,rdat%vel(igate,iray,kelv),unfvel(igate,iray),vnyq
                 endif
                                    unfvel(igate,iray)=tstvel
                                 else
                                    kntgood=kntgood+1
                                 endif
                              else
                                 kntrej=kntrej+1
                                 !write(*,'(a,i5,a,i5,a,f10.1,a,f10.1,a,f10.1)')  &
                                 ! 'igate=',igate,' iray=',iray,                  &
                                 ! ' Unresolved: Meas=',unfvel(igate,iray),      &
                                 ! ' Avgvel=',vavg,'  Try=',tstvel
                                 unfvel(igate,iray)=unfvel(igate,iray)-2000.
                                 tmp1(igate,iray)=0.
                              endif              ! tstvel-vavg >thresh
                           endif               ! unfvel-vavg >thresh
                        endif           ! tknt > 0
                     else
                        kntgood=kntgood+1
                     endif                !(unfvel(igate,iray) - refvel) > thrpri
                  endif                  ! tmp1(igate,iray) > 0.
               enddo
            enddo
   
            write(*,'(/a/,a/,4i12)')                                            &
                  '       After reverse pass shear consistency check',          &
                  '       Gates  Passed Check   Folded     Rejected',           &
                          kntall,kntgood,kntfold,kntrej

            !   Apply gate-to-gate shear checking for velocity folding.
            !   Loop backward from jstart-1 not use background vel

            kntall = 0
            kntgood = 0
            kntfold = 0
            kntrej = 0
            do i=1,nazim
               iray=(jstart-i)+1
               if(iray <= 0) iray=(jstart-i)+nazim
               
               vnyq = rdat%vmax(kelv)
               twonyq = 2.0*vnyq
               inv2nyq = 1./twonyq
               thrpri = min(max((0.5*vnyq),10.),vnyq)
               thrpr2 = min(vnyq,(1.5*thrpri))
               
               do igate=bgate(iray),egate(iray)
                  if( tmp1(igate,iray) > 0. ) then
                     dbkvel=unfvel(igate,iray)
                     kntall=kntall+1
                     lgate = 0
                     if(igate > bgate(iray)) then
                        ekgate = max((igate-lukbak),bgate(iray))
                        do kgate = igate-1,ekgate,-1
                           if(tmp1(kgate,iray) >0. ) then
                              lgate = kgate
                              EXIT
                           endif
                        enddo
                     endif
                     if(lgate > 0) then
                        refvel = unfvel(lgate,iray)
                     else
                        refvel = 999.
                     endif
                     shdiff=abs(dbkvel-refvel)
                     if( shdiff > thrpri) then
                        sum=0.
                        fknt=0.
                        bknt=0.
                        if(igate > bgate(iray)) then
                           ekgate=max((igate-lukbak2),bgate(iray))
                           do kgate=igate-1,ekgate,-1
                              bknt=bknt+tmp1(kgate,iray)
                              sum=sum+tmp1(kgate,iray)*                                &
                                  (unfvel(kgate,iray))
                              if (bknt > 2. ) EXIT
                           enddo
                        endif
                        eray=min((iray+lukbakr),nazim)
                        do kray=iray+1,eray
                           ekgate=min((igate+lukfwd2),egate(kray))
                           do kgate=igate,ekgate
                              fknt=fknt+tmp1(kgate,kray)
                              sum=sum+tmp1(kgate,kray)*                                &
                                  (unfvel(kgate,kray))
                              if ( fknt > 4. ) EXIT
                           enddo
                        enddo
   
                        tknt=bknt+fknt
                        if(tknt > 2. ) then
                           vavg=(sum/tknt)
                           shdif2=abs(unfvel(igate,iray)-vavg)
                           if( shdif2 < thrpr2 ) then
                              kntgood=kntgood+1
                           else
                              tstdev=twonyq*NINT((vavg-unfvel(igate,iray))*inv2nyq)
                              tstvel=unfvel(igate,iray)+tstdev
                              shdif3=(tstvel-vavg)
                              if( (shdif3 < thrpr2 ).and.abs(shdif3)<abs(shdiff)) then
                                 if( abs(tstdev) > 0. ) then
                                    write(600,"(A,3I5,3F12.2)") "nobkg shr:",igate,iray,kelv,rdat%vel(igate,iray,kelv),tstvel,vavg
                                    !write(600,"(A,4F12.2)") "shr:",rdat%vel(igate,iray,kelv),bkgvel(igate,iray,kelv),tstvel,vavg
                                    kntfold=kntfold+1
                 if(abs(tstvel)>2*vnyq)then
                    write(701,*) "unfold check5-1",tstdev,tstvel,vavg,rdat%vel(igate,iray,kelv),unfvel(igate,iray),vnyq, fknt,bknt,tknt,refvel,shdiff,shdif3
                    write(701,*) "unfold check5-2",igate,ekgate,iray,eray,kelv
                    do kray=iray,eray
                    write(701,*) (int(tmp1(kgate,kray)), unfvel(kgate,kray),kgate=igate,ekgate)
                    enddo
                 endif
                                    unfvel(igate,iray)=tstvel
                                 else
                                    kntgood=kntgood+1
                                 endif
                              else
                                 kntrej=kntrej+1
                                 !write(*,'(a,i5,a,i5,a,f10.1,a,f10.1,a,f10.1)')  &
                                 ! 'igate=',igate,' iray=',iray,                  &
                                 ! ' Unresolved: Meas=',unfvel(igate,iray),      &
                                 ! ' Avgvel=',vavg,'  Try=',tstvel
                                 unfvel(igate,iray)=unfvel(igate,iray)-2000.
                                 tmp1(igate,iray)=0.
                              endif              ! tstvel-vavg >thresh
                           endif               ! unfvel-vavg >thresh
                        endif           ! tknt > 0
                     else
                        kntgood=kntgood+1
                     endif                !(unfvel(igate,iray) - refvel) > thrpri
                  endif                  ! tmp1(igate,iray) > 0.
               enddo
            enddo
   
            write(*,'(/a/,a/,4i12)')                                            &
                  '       After reverse pass shear consistency check(nobkg)',          &
                  '       Gates  Passed Check   Folded     Rejected',           &
                          kntall,kntgood,kntfold,kntrej

            !endif
         enddo !ipass
   
      endif ! shropt > 0
enddo
      do iray=1,nazim
         do igate=1,ngate
            !if(tmp2(igate,iray)>0)then
               rdat%vel(igate,iray,kelv)=unfvel(igate,iray)
            !else
            !   rdat%vel(igate,iray,kelv)=value_invalid
            !endif
         enddo
      enddo
   enddo ! kelv
   
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine unfold_4dd(rvVolume, missingVal, filt, success, soundVolume, lastVolume)
   use radar_data, only: copy_radar_data
   implicit none
   !/*
   !**
   !**  UW Radial Velocity Dealiasing Algorithm
   !**  Four-Dimensional Dealiasing (4DD)
   !**
   !**  DESCRIPTION:
   !**     This algorithm unfolds a volume of single Doppler radial velocity data.
   !**  The algorithm uses a previously unfolded volume (or VAD if previous volume
   !**  is unavailable) and the previous elevation sweep to unfold some of gates 
   !**  in each sweep. Then, it spreads outvalward from the 'good' gates, completing
   !**  the unfolding using gate-to-gate continuity. Gates that still remain
   !**  unfolded are compared to an areal average of neighboring dealiased gates.
   !**  Isolated echoes that still remain uncorrected are dealiased against a VAD
   !**  (as a last resort).
   !**
   !**  DEVELOPER:
   !**	Curtis N. James     25 Jan 99
   !**
   !**
   !*/
   !#include "4DD.h"
   !
   !#include <stdio.h>
   !#include <math.h>
   !#include <stdlib.h>
   !#include "rsl.h" /* Sweep */
   !
   !/* This routvaline performs a preliminary unfold on the data using the bin in the
   !** next highest sweep (already unfolded), and the last unfolded volume. If the 
   !** last unfolded volume is not available, the VAD is used (or sounding if VAD
   !** not available). If this is successful, the bin is considered GOOD. Every 
   !** other bin is set to GOOD=0 or GOOD=-1 (if the bin is missing or bad).
   !**
   !** Then, the algorithm scans azimuthally and radially. If the majority of the
   !** GOOD=1 bins (up to 8) adjacent to a particular GOOD=0 bin are within 
   !** THRESH*Vnyq it is saved and considered GOOD as well. Otherwise, Nyquist 
   !** intervals are added/subtracted until the condition is met. If unsuccessful, 
   !** then GOOD is set to -2 for that bin.
   !**
   !** When all GOOD=0 bins that are next to GOOD=1 bins have been examined, then
   !** GOOD=0 and GOOD=-2 bins are unfolded against a window of GOOD=1 values. If
   !** there are too few GOOD=1 values inside the window, then the data are
   !** unfolded against the VAD, if available.
   !** 
   !** PARAMETERS:
   !** rvVolume: The radial velocity field to be unfolded.
   !** soundVolume: A first guess radial velocity field using VAD (or sounding) 
   !** lastVolume: The last radial velocity field unfolded (time of lastVolume 
   !**    should be within 10 minutes prior to rvVolume)
   !** missingVal: The value for missing radial velocity data.
   !** filt: A flag that specifies whether or not to use Bergen/Albers filter.
   !** success: flag indicating whether or not unfolding is possible.
   !*/
   !
   
         type(t_radar_data), intent(inout) :: rvVolume
         real, dimension(:,:,:),optional,intent(in) :: soundVolume
         type(t_radar_data),optional,intent(in) :: lastVolume
         real, intent(in) :: missingVal
         integer, intent(in) :: filt
         integer, intent(out) :: success 
        
        
         integer :: sweepIndex, currIndex, i, l, m, n, direction, numSweeps, numRays
         integer :: numBins, numdo, left, right, next, prev
         integer :: rayindex(8), binindex(8)
         integer :: countindex, numneg, numpos, in, outval, startray, endray, firstbin
         integer :: lastbin, step = -1, startindex, endindex, prevIndex, abIndex, loopcount
         integer :: countbins, abSweepIndex
         
         integer :: numtimes, dcase, flag=1, wsuccess
        
         integer, dimension(:,:), allocatable :: GOOD !(MAXBINS,MAXRAYS);
        
         real :: NyqVelocity, NyqInterval, val, diff, fraction, finalval, initval
         real :: valcheck, goodval, winval, vdiff, fraction2
         real :: diffs(8)
        
         real ::  prevval, abval, pfraction, backval, cval, soundval, std
        
         real, dimension(:,:), allocatable :: OUT
        
         real :: elev, azimuth, abazimuth,rayazimuth
        
         character(len=50) :: ofname
         integer, parameter :: fpo=114
         type(t_radar_data) :: VALS
        
         ! 
         real,    parameter :: COMPTHRESH=0.25  !/* The threshold for performing initial dealiasing 
                                                !** using a previously unfolded volume. */
         real,    parameter :: COMPTHRESH2=0.49 !/* The threshold for performing initial dealiasing 
                                                !** using sounding (or VAD). */
         real,    parameter :: THRESH=0.4       !/* The unfolding threshold for unfolding using horizontal
                                                ! ** continuity. */
         logical, parameter :: VERBOSE=.FALSE.  !/* Verbose=1 for detailed printout during execution */
        
         integer, parameter :: DELNUM=20        !/* The first DELNUM velocity bins will be deleted along each
                                                !**  ray (should be between 0 and 5). */
         integer, parameter :: MAXCOUNT=10      !/* Maximum number of folds. */
         real,    parameter :: CKVAL=1.0        !/* If absolute value of the radial velocity gate is less 
                                                !** than this value, it will not be used as a PRELIM gate. */
        
         integer, parameter :: MINGOOD=5        !/* Number of good values required within unfolding window
                                                !**  to unfold the current bin. */
         real,    parameter :: STDTHRESH=0.8    !/* Fraction of the Nyquist velocity to use as a standard
                                                !**  deviation threshold when windowing. */
         integer, parameter :: PROXIMITY=5      !/* For unfolding using windowing.*/
   
         integer, parameter :: RM=0             !/* If soundvolume is not available, remove cells left over after
                                                !**  first pass. */
         integer, parameter :: PASS2=1          !/* Flag specifying the use of a second pass using only the
                                                !**   sounding (or VAD).*/
   
         integer :: maxbins, maxrays, ierr
         character(len=800) :: message
   
         numSweeps = rvVolume%ntilt
         maxbins= ubound(rvVolume%vel,dim=1)
         maxrays= ubound(rvVolume%vel,dim=2)
         write(*,*) maxbins, maxrays
         allocate(GOOD(maxbins,maxrays))
         !allocate( OUT(maxbins,maxrays),stat=ierr,errmsg=message)
         !if(ierr/=0)then
         !   write(*,*) "ALLOCATE OUT",trim(message)
         !   stop
         !endif
        
         call copy_radar_data(rvVolume,VALS)
        
         if (COMPTHRESH>1.0 .or. COMPTHRESH<=0.0 )then
            fraction=0.25
         else
            fraction=COMPTHRESH
         endif
        
         if (COMPTHRESH2>1.0 .or. COMPTHRESH2<=0.0 )then
            fraction2=0.25
         else 
            fraction2=COMPTHRESH2
         endif
        
         if (THRESH>1.0 .or. THRESH<=0.0 )then
            pfraction=0.5
         else 
            pfraction=THRESH
         endif
        
         if (present(soundVolume) .or. present(lastVolume)) then
   
         write(*,*) "soundVolume lastVolume"
   
         !/*###################### SWEEP LOOP ########################*/
   
         do sweepIndex=numSweeps,1,-1
   
            if(.NOT.rvVolume%ifvel(sweepIndex)) cycle
            NyqVelocity = rvVolume%vmax(sweepIndex) 
            
            NyqInterval = 2.0 * NyqVelocity
            numRays = rvVolume%nazim(sweepIndex)
            numBins = maxval(rvVolume%nvgate(1:numRays,sweepIndex))
            allocate(OUT(numBins,numRays))
            
            !/* First, unfold bins where vertical and temporal continuity
            !** produces the same value (i.e., bins where the radial velocity
            !** agrees with both the previous sweep and previous volume within
            !** a COMPTHRESH of the Nyquist velocity). This produces a number
            !** of good data points from which to begin unfolding assuming spatial
            !** continuity. */
            
            elev= sum(rvVolume%rtilt(1:numRays,sweepIndex))/numRays
            
            write(*,*) "############### SWEEP##############: ", sweepIndex,elev,numRays,numBins
            
            if (VERBOSE) write(*,*) "NyqVelocity: ",NyqVelocity,", missingVal: ",missingVal
        
            flag=1
            
            write(*,*) "FILTERING INITIAL"
            
            do  currIndex=1,numRays
         
               if (present(lastVolume))then
                  prevIndex=findRay(rvVolume, lastVolume,sweepIndex, sweepIndex,currIndex)
               endif
                  
         
               if (sweepIndex<numSweeps)then
                  abIndex=findRay(rvVolume, rvVolume,sweepIndex, abSweepIndex, currIndex)
               endif
                  
         
               do i=1,DELNUM
                  !/* Initialize Output Sweep with missing values: */
         
                  !write(*,*) i, currIndex, ubound(OUT)
                  OUT(i,currIndex)=missingVal
         
                  initval=missingVal
         
                  rvVolume%vel(i,currIndex,sweepIndex)=initval
                    
               enddo
         
         
               do i=DELNUM+1,numBins
                  !/* Initialize Output Sweep with missing values: */
         
                  OUT(i,currIndex)=missingVal
         
                  initval=missingVal
         
                  rvVolume%vel(i,currIndex,sweepIndex)=initval
         
                  !/* Assign uncorrect velocity bins to the array VALS: */
         
                  !write(*,*) i,currIndex,sweepIndex, ubound(VALS%vel)
                  val= VALS%vel(i,currIndex,sweepIndex)
         
                  valcheck=val
         
                  if (val==missingVal.or.val==131071.00)then
                     GOOD(i,currIndex)=-1
                  else 
         
                  !/* ########################################################################## */
                  !/* ############################ STEP 2    ################################### */
                  !/* ############################ FILTERING (3X3) ############################# */
                  !/* ########################################################################## */
         
                     if (filt==1) then
                        !/*Perform a 3x3 filter, as proposed by Bergen & Albers 1988*/
                         
                        countindex=0
                        if (currIndex==1)then
                           left=numRays
                        else
                           left=currIndex-1
                        endif
                
                        if (currIndex==numRays)then
                           right=1
                        else 
                           right=currIndex+1
                        endif
                
                        next=i+1
                        prev=i-1
                
                        !/* Look at all bins adjacent to current bin in question: */
         
                        !/* Three bin inward */
         
                        !write(*,*) i, currIndex,DELNUM,prev,next,left,right,sweepIndex, ubound(GOOD)
                        if (i>(DELNUM+1)) then
                           if (VALS%vel(prev,left,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                
                           if (VALS%vel(prev,currIndex,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                
                           if (VALS%vel(prev,right,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                        endif
                
   
                        ! /*  Two bin same distance */
   
                        if (VALS%vel(i,left,sweepIndex)/=missingVal) then
                           countindex=countindex+1;
                        endif
   
                        if (VALS%vel(i,currIndex,sweepIndex)/=missingVal) then
                           countindex=countindex+1;
                        endif
                
                        if (VALS%vel(i,right,sweepIndex)/=missingVal) then
                           countindex=countindex+1;
                        endif
                
                        !/* Three bin outward */
         
                        if (i<numBins) then  
                           if (VALS%vel(next,left,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                         
                           if (VALS%vel(next,currIndex,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                         
                           if (VALS%vel(next,right,sweepIndex)/=missingVal) then
                              countindex=countindex+1;
                           endif
                        endif
   
                        if (((i==numBins .or. i==DELNUM) .and. countindex>=3) .or.(countindex>=5)) then
                           !/* Save the bin for dealiasing: */
                           GOOD(i,currIndex)=0
                        else
   
                           !/* ++++++++++++++++++++ FILTERING ++++++++++++++++++++*/
   
                           !/* Assign missing value to the current bin. */
                           GOOD(i,currIndex)=-1
             
                        endif
   
                     !/* End 3 x 3 filter */
                     else 
                     ! /* If no filter is being applied save bin for dealiasing: */
                        GOOD(i,currIndex)=0
                     endif ! filt==1
   
                     !/* ########################################################################## */
                     !/* ##########################  STEP 3 ####################################### */
                     !/* ######## INITIAL DEALIASING ( LASTvolume, ABOVE TILT, SOUNDING )########## */
                     !/* ########################################################################## */
   
                     if (GOOD(i,currIndex)==0) then
                        !/* Try to dealias the bin using vertical and temporal 
                        !**   continuity (this is initial dealiasing). */
                     
                        if (val/=missingVal .and. present(lastVolume)) then
                           prevval=lastVolume%vel(i,prevIndex,sweepIndex)
                        else
                           prevval=missingVal;
                        endif
                     
                        if (val/=missingVal .and. (.NOT.present(lastVolume)) .and. present(soundVolume))then
                          soundval=soundVolume(i,currIndex,sweepIndex)
                        else
                          soundval=missingVal
                        endif
                     
                        if (val/=missingVal .and. sweepIndex<numSweeps) then
                           abval=rvVolume%vel(i,abIndex,sweepIndex)
                        else
                          abval=missingVal;
                        endif
   
                        !/* ###### use sounding volume ##### */
   
                        if (val  /=missingVal .and. prevval  == missingval .and. &
                            abval==missingVal .and. soundval /= missingVal ) then
                                    
                           cval=soundval
                           dcase=1  
   
                        !/* ###### use above tilt ###### */
   
                        else if (val  /=missingVal .and. prevval  == missingval .and. &
                                 abval/=missingVal ) then
                           cval=abval
                           dcase=2 
   
                        !/* ###### use last volume ##### */
   
                        else if (val  /=missingVal .and. prevval /= missingval ) then
                           cval=prevval
                           dcase=3
   
                        else
                           cval=missingVal
                           dcase=0
                        endif
   
                        !/* If gate data at last volume, sounding, or above tilt is good, 
                        !    do dealiasing for this bin we trust last volume most */
   
                        if (dcase>0) then
                        
                           diff=cval-val
                        
                           if (diff<0.0) then
                              diff=-diff
                              direction=-1
                           else
                              direction=1
                           endif
                        
                           numtimes=0
                        
                           do while (diff>0.99999*NyqVelocity .and. numtimes<=MAXCOUNT) 
                        
                              val=val+NyqInterval*direction
                        
                              numtimes=numtimes+1
                              diff=cval-val
                        
                              if (diff<0.0) then
                                 diff=-diff
                                 direction=-1
                              else
                                 direction=1
                              endif
                           enddo 
                        
                           if (dcase==1 .and. diff<fraction*NyqVelocity .and. &
                               abs(valcheck)>CKVAL) then
                        
                               if (VERBOSE) write(*,*) "GOOD1: ", val
                        
                               finalval=val
                        
                               rvVolume%vel(i,currIndex,sweepIndex)=finalval
                        
                               GOOD(i,currIndex)=1
                               OUT(i,currIndex)=val
                        
                           else if (dcase==2 .and. diff<fraction*NyqVelocity .and. &
                                    abs(soundval-val)<fraction*NyqVelocity .and. abs(valcheck)>CKVAL) then
                        
                               if (VERBOSE) write(*,*) "GOOD2: ", val
                        
                               finalval=val
                        
                               rvVolume%vel(i,currIndex,sweepIndex)=finalval
                        
                               GOOD(i,currIndex)=1
                               OUT(i,currIndex)=val
                        
                           else if (dcase==3 .and. diff<fraction*NyqVelocity .and. &
                                    abs(abval-val)<fraction*NyqVelocity .and. abs(valcheck)>CKVAL) then 
                        
                        
                               if (VERBOSE) write(*,*) "GOOD3: ", val
                        
                               finalval=val
                        
                               rvVolume%vel(i,currIndex,sweepIndex)=finalval
                        
                               GOOD(i,currIndex)=1
                               OUT(i,currIndex)=val
                           endif 
                        endif 
                     endif  ! GOOD(i,currIndex)==0
                  endif !val==missingVal.or.val==131071.00
               enddo ! i=DELNUM+1,numBins
            enddo ! currIndex=1,numRays
   
            !/* ########################################################################## */
            !/* #################################  STEP 4  ############################### */
            !/* ###################### SPATIAL DEALIASING  continuity  ################### */
            !/* ########################################################################## */
   
            !/*continue;*/
   
            write(*,*) "SPATIAL DEALIASING continuity" 
   
            !/* Now, unfold GOOD=0 bins assuming spatial continuity: */
            loopcount=0;
            
            do while (flag==1) 
            
               loopcount=loopcount+1
               flag=0
            
               if (step==1)then
                 step=-1
                 startindex=numRays
                 endindex=1
               else 
                 step=1
                 startindex=1
                 endindex=numRays
               endif
   
               do i=DELNUM+1,numBins
   
                  do currIndex=startindex,endindex,step
                  !/* Bug fix 1/31/03. */
                     val=VALS%vel(i,currIndex,sweepIndex)
   
                     !/*          if (GOOD(i,currIndex)==0 && val!=missingVal) {  */
                     if (GOOD(i,currIndex)==0 .and. val/=missingVal .and. val/=131071.00) then 
   
                        countindex=0
                        countbins=0
                        if (VERBOSE) write(*,*) "Startval: ", val
                        if (currIndex==1)then
                           left=numRays
                        else
                           left=currIndex-1
                        endif
                        if (currIndex==numRays)then
                           right=1
                        else 
                           right=currIndex+1
                        endif
                        next=i+1;
                        prev=i-1;
                       
                        !/* Look at all bins adjacent to current bin in question: */
                        !write(*,*) i, currIndex,DELNUM,prev,next,left,right,sweepIndex, ubound(GOOD)
                        if (i>(DELNUM+1)) then
                       
                           if (GOOD(prev,left)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=prev
                              rayindex(countindex)=left
                              if (VERBOSE) write(*,*)"pl "
                           endif
                       
                           if (GOOD(prev,left)==0) then
                              countbins=countbins+1
                           endif
                       
                           if (GOOD(prev,currIndex)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=prev
                              rayindex(countindex)=currIndex
                              if (VERBOSE) write(*,*)"pc "
                           endif
                       
                           if (GOOD(prev,currIndex)==0) then
                              countbins=countbins+1
                           endif
                       
                           if (GOOD(prev,right)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=prev
                              rayindex(countindex)=right
                              if (VERBOSE) write(*,*)"pr "
                           endif
                       
                           if (GOOD(prev,right)==0) then
                             countbins=countbins+1
                           endif
                       
                        endif
                       
                        if (GOOD(i,left)==1) then
                           if (flag==0) flag=1
                           countindex=countindex+1
                           binindex(countindex)=i
                           rayindex(countindex)=left
                           if (VERBOSE) write(*,*)"il "
                        endif
                       
                        if (GOOD(i,left)==0) then
                           countbins=countbins+1
                        endif
                       
                        if (GOOD(i,right)==1) then
                           if (flag==0) flag=1
                           countindex=countindex+1
                           binindex(countindex)=i
                           rayindex(countindex)=right
                           if (VERBOSE) write(*,*)"ir "
                        endif
                       
                        if (GOOD(i,right)==0) then
                           countbins=countbins+1
                        endif
                       
                        if (i<numBins) then  
                       
                           if (GOOD(next,left)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=next
                              rayindex(countindex)=left
                              if (VERBOSE) write(*,*)"nl "
                           endif
                       
                           if (GOOD(next,left)==0) then
                              countbins=countbins+1
                           endif
                       
                           if (GOOD(next,currIndex)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=next
                              rayindex(countindex)=currIndex
                              if (VERBOSE) write(*,*)"nc "
                           endif
                       
                           if (GOOD(next,currIndex)==0) then
                              countbins=countbins+1
                           endif
                       
                           if (GOOD(next,right)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=next
                              rayindex(countindex)=right
                              if (VERBOSE) write(*,*)"nr "
                           endif
                       
                           if (GOOD(next,right)==0) then
                              countbins=countbins+1
                           endif
                       
                        endif
                       
                        !/* Perform last step of Bergen and Albers filter: */
                        if (loopcount==1 .and. countbins+countindex<1) GOOD(i,currIndex)=-1
                              
                        if (countindex>=1) then
                           !/* Unfold against all adjacent values where GOOD==1 */
                           numtimes=0;
                       
                           do while(val/=missingVal .and. GOOD(i,currIndex)==0)
                              numtimes=numtimes+1
                              in=0
                              outval=0
                              numneg=0
                              numpos=0
   
                              if (VERBOSE) write(*,*)"countindex: ", countindex
                              
                              do l=1,countindex
                              
                                 if (VERBOSE) write(*,*)"",binindex(l),rayindex(l),goodval
   
                                 goodval=rvVolume%vel(binindex(l),rayindex(l),sweepIndex)
                                 
                                 diffs(l)=goodval-val;
                                 
                                 if (abs(diffs(l))<pfraction*NyqVelocity)then
                                    in=in+1
                                 else
                                    outval=outval+1
                                    if (diffs(l)>NyqVelocity)then
                                       numpos=numpos+1
                                    else if (diffs(l)<-NyqVelocity) then 
                                       numneg=numneg+1
                                    endif
   
                                 endif 
                              enddo
   
                              if (in>0 .and. outval==0) then
                              
                                 finalval=val
                              
                                 rvVolume%vel(i,currIndex,sweepIndex)=finalval
                              
                                 GOOD(i,currIndex)=1
                                 OUT(i,currIndex)=val
                                 if (VERBOSE) write(*,*)"Finalval: ", val
                                         
                              else
                              
                                 if (numtimes<=MAXCOUNT) then
                                    if ((numpos+numneg)<(in+outval-(numpos+numneg))) then
   
                                       if (loopcount>2)then
                                          !/* Keep the value after two passes through
                                          !** data. */
                                          finalval=val
                                     
                                          rvVolume%vel(i,currIndex,sweepIndex)=finalval
                                     
                                          GOOD(i,currIndex)=1
                                          OUT(i,currIndex)=val
                                          if (VERBOSE) write(*,*)"Keeping after two passes: ",val
                                       endif
   
                                    else if (numpos>numneg) then
   
                                       do l=1,countindex
                                          goodval=rvVolume%vel(binindex(l),rayindex(l),sweepIndex)
                                          diffs(l)=goodval-val
                                       enddo
   
                                       val=val+NyqInterval
   
                                    else if (numneg>numpos) then
   
                                       do l=1,countindex
                                          goodval=rvVolume%vel(binindex(l),rayindex(l),sweepIndex)
                                          diffs(l)=goodval-val
                                       enddo
   
                                       val=val-NyqInterval
                             
                                    else 
                                       !/* Save bin for windowing if unsuccessful after
                                       !** four passes: */
                                       if (loopcount>4) GOOD(i,currIndex)=-2
                                       if (VERBOSE) write(*,*)"Saving for windowing"
                                    endif
                                     
                                 else
                                    !   /* Remove bin: */
                                    GOOD(i,currIndex)=-2
                                    if (VERBOSE) write(*,*)"Saving for windowing"
                                 endif
                                        
                              endif  !in>0 .and. outval==0
                           enddo !while(val/=missingVal .and. GOOD(i,currIndex)==0)
                        endif !countindex>=1
                     endif !GOOD(i,currIndex)==0 .and. val/=missingVal .and. val/=131071.00
                  enddo !currIndex=startindex,endindex,step
               enddo !i=DELNUM+1,numBins
   
   
               !/* ########################################################################## */
               !/* #################################  STEP 5  ############################### */
               !/* #######################  WINDOW DEALIASING   average   ################### */
               !/* ########################################################################## */
   
               write(*,*)"yli"
   
   
               !/*continue;*/
   
               write(*,*)"WINDOW DEALIASING  average "
   
               ! /* Unfold remote bins or those that were previously unsuccessful
               ! **   using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
               ! **   if still no luck delete data (or unfold against VAD if PASS2). */
               !
               do i=DELNUM+1,numBins
   
                  do currIndex=1,numRays
   
                     if (GOOD(i,currIndex)==0 .or. GOOD(i,currIndex)==-2)then
                        val=VALS%vel(i,currIndex,sweepIndex)
   
                        startray=currIndex-PROXIMITY;
                        endray=currIndex+PROXIMITY;
                        firstbin=i-PROXIMITY;
                        lastbin=i+PROXIMITY;
   
                        if (startray<1) startray=numRays+startray
                        if (endray>numRays) endray=endray-numRays
                        if (firstbin<1) firstbin=1
                        if (lastbin>numBins) lastbin=numBins
   
                        call window(rvVolume, sweepIndex, startray, endray, firstbin, lastbin, winval, std, wsuccess)
                        
   
                        if (winval==missingVal .and. wsuccess==1)then !/* Expand the window: */
   
                           startray=currIndex-2*PROXIMITY;
                           endray=currIndex+2*PROXIMITY;
                           firstbin=i-2*PROXIMITY;
                           lastbin=i+2*PROXIMITY;
   
                           if (startray<1) startray=numRays+startray
                           if (endray>numRays) endray=endray-numRays
                           if (firstbin<1) firstbin=1
                           if (lastbin>numBins) lastbin=numBins
         
                           call window(rvVolume, sweepIndex, startray, endray, firstbin, lastbin, winval, std, wsuccess)
                            
                        endif      
         
                        if (winval/=missingVal) then
                           diff=winval-val;
         
                           if (diff<0.0) then
                              diff=-diff
                              direction=-1
                           else
                              direction=1
                           endif
   
                           numtimes=0
                          
                           do while (diff>0.99999*NyqVelocity .and. numtimes<=MAXCOUNT) 
                          
                              val=val+NyqInterval*direction
                          
                              numtimes=numtimes+1
                          
                              diff=winval-val
                          
                              if (diff<0.0) then
                                 diff=-diff
                                 direction=-1
                              else
                                 direction=1
                              endif
                           enddo
                                                         
   
                           if (diff<pfraction*NyqVelocity) then
                              !/* Return the value. */
                         
                              finalval=val
                         
                              rvVolume%vel(i,currIndex,sweepIndex)=finalval
                         
                              GOOD(i,currIndex)=1
                              OUT(i,currIndex)=val
                         
                           else if (diff<(1.0 - (1.0 - pfraction)/2.0)*NyqVelocity) then
                              !/* If within relaxed threshold, then return value, but
                              !**   do not use to dealias other bins. */
   
                              finalval=val
                            
                              rvVolume%vel(i,currIndex,sweepIndex)=finalval
                              GOOD(i,currIndex)=-1
                            
                           else 
                              !/* Remove bin */
                              GOOD(i,currIndex)=-1
                           endif
                            
                        else
                           if (wsuccess==0) then
                              !/* Remove bin */
                              GOOD(i,currIndex)=-1
                           else if ((.NOT.present(soundVolume)) .or. (.NOT.present(lastVolume))) then
                           
                              if (GOOD(i,currIndex)==0 .and. RM/=1) then
                           
                                 !/* Leave bin untouched. */
                                 val=VALS%vel(i,currIndex,sweepIndex)
                                 finalval=val
                           
                                 rvVolume%vel(i,currIndex,sweepIndex)=finalval
                                 GOOD(i,currIndex)=-1 !/* Don't use to unfold other bins*/
   
                              else
                                 !/* Remove bin */
                                 GOOD(i,currIndex)=-1;
                              endif
   
                           else if (GOOD(i,currIndex)==0 .and. PASS2>0 .and. present(soundVolume) .and. present(lastVolume)) then
   
                               !/* Leave GOOD(i,currIndex)=0 bins for a second pass.
                               !** In the second pass, we repeat unfolding, except this
                               !** time we use soundVolume for comparison instead of
                               !** lastVolume. */
   
                           else
                              !/* Remove bin */
                              GOOD(i,currIndex)=-1
                           endif !wsuccess==0
                        endif !winval/=missingVal
                     endif !GOOD(i,currIndex)==0 .or. GOOD(i,currIndex)==-2
                  enddo !currIndex=1,numRays
               enddo  !i=DELNUM+1,numBins                       
            enddo
   
   
            !/* ########################################################################## */
            !/* #################################  STEP 6.1  ############################# */
            !/* ####################### AUXILIARY DEALIASING sounding #################### */ 
            !/* ########################################################################## */
   
            !/*continue;*/
   
            write(*,*)"AUXILIARY DEALIASING sounding"
   
            !/* Beginning second pass through the data, this time using 
            !   **  soundVolume only: */
   
            if (present(lastVolume) .or. present(soundVolume))then  
   
               flag=1
   
               do currIndex=1,numRays
   
                  if (sweepIndex<numSweeps) abIndex=findRay(rvVolume, rvVolume,sweepIndex, abSweepIndex, currIndex)
                      
                  do i=DELNUM+1,numBins
                     if (GOOD(i,currIndex)==0) then
   
                        val=VALS%vel(i,currIndex,sweepIndex)
                        valcheck=val
                        
                        soundval= soundVolume(i,currIndex,sweepIndex)
                                
                        if (soundval/=missingVal .and. val/=missingVal) then
   
                           diff=soundval-val
   
                           if (diff<0.0) then
                              diff=-diff
                              direction=-1
   
                           else
                              direction=1
                           endif
   
                           numtimes=0
   
                           do while (diff>0.99999*NyqVelocity .and. numtimes<=MAXCOUNT)
   
                              val=val+NyqInterval*direction;
   
                              numtimes=numtimes+1;
                              diff=soundval-val;
   
                              if (diff<0.0) then
                                 diff=-diff
                                 direction=-1
                              else 
                                 direction=1
                              endif
                           enddo
   
                           if (diff<fraction2*NyqVelocity .and. abs(valcheck)>CKVAL) then
   
                              finalval=val
                              rvVolume%vel(i,currIndex,sweepIndex)=finalval
   
                              GOOD(i,currIndex)=1
                              OUT(i,currIndex)=val
                           endif
                        endif
                     endif
                  enddo
               enddo
   
   
               !/* ########################################################################## */
               !/* #################################  STEP 6.2  ############################# */
               !/* ############## AUXILIARY DEALIASING spatial continuity  ################## */ 
               !/* ########################################################################## */
   
               !/*continue;*/
   
               write(*,*)"AUXILIARY DEALIASING spatial continuity"
   
               !/* Now, try to unfold the rest of the GOOD=0 bins, assuming spatial
               !   continuity: */
               
               loopcount=0
               
               do while (flag==1) 
                  loopcount=loopcount+1
                  flag=0;
               
                  if (step==1) then
                     step=-1
                     startindex=numRays
                     endindex=1
                  else
                     step=1
                     startindex=1
                     endindex=numRays
                  endif
               
               
                  do currIndex=startindex,endindex,step
               
                     do i=DELNUM+1,numBins
               
                        if (GOOD(i,currIndex)==0) then
                           countindex=0;
                           val=VALS%vel(i,currIndex,sweepIndex)
               
                           valcheck=val
                           if (currIndex==1)then
                              left=numRays
                           else
                              left=currIndex-1
                           endif
               
                           if (currIndex==numRays)then
                              right=1
                           else 
                              right=currIndex+1
                           endif
               
                           next=i+1;
                           prev=i-1;
               
                           !/* Look at all bins adjacent to current bin in
                           !   question: */
               
                           if (i>DELNUM+1) then
                              if (GOOD(prev,left)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=prev
                                 rayindex(countindex)=left
                                 if (VERBOSE) write(*,*)"pl "
                              endif
               
                              if (GOOD(prev,currIndex)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=prev
                                 rayindex(countindex)=currIndex
                                 if (VERBOSE) write(*,*)"pc "
                              endif
               
                              if (GOOD(prev,right)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=prev
                                 rayindex(countindex)=right
                                 if (VERBOSE) write(*,*)"pr "
                              endif
                           endif
               
                           if (GOOD(i,left)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=i
                              rayindex(countindex)=left
                              if (VERBOSE) write(*,*)"il "
                           endif
               
                           if (GOOD(i,right)==1) then
                              if (flag==0) flag=1
                              countindex=countindex+1
                              binindex(countindex)=i
                              rayindex(countindex)=right
                              if (VERBOSE) write(*,*)"ir "
                           endif
               
                           if (i<numBins) then
               
                              if (GOOD(next,left)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=next
                                 rayindex(countindex)=left
                                 if (VERBOSE) write(*,*)"nl "
                              endif
               
                              if (GOOD(next,currIndex)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=next;
                                 rayindex(countindex)=currIndex;
                                 if (VERBOSE) write(*,*)"nc "
                              endif
               
                              if (GOOD(next,right)==1) then
                                 if (flag==0) flag=1
                                 countindex=countindex+1
                                 binindex(countindex)=next
                                 rayindex(countindex)=right
                                 if (VERBOSE) write(*,*)"nr "
                              endif
                           endif
   
                           !/* Unfold against all adjacent values with GOOD==1*/
                           if (countindex>=1) then
                              numtimes=0;
   
                              do while(val/=missingVal .and. GOOD(i,currIndex)==0)
                                 numtimes=numtimes+1
                                 in=0
                                 outval=0
                                 numneg=0
                                 numpos=0
                                 if (VERBOSE) write(*,*)"countindex: ", countindex
   
                                 do l=1,countindex
                                    goodval=rvVolume%vel(binindex(l),rayindex(l),sweepIndex)
                                    
                                    diffs(l)=goodval-val;
                                    if (abs(diffs(l))<pfraction*NyqVelocity)then
                                       in=in+1
                                    else 
                                       outval=outval+1
                                    
                                       if (diffs(l)>NyqVelocity) then
                                         numpos=numpos+1
                                       else if (diffs(l)<-NyqVelocity)then
                                         numneg=numneg+1
                                       endif
                                    endif
                                 enddo
   
                                 if (in>outval) then
                                    finalval=val
                                    rvVolume%vel(i,currIndex,sweepIndex)=finalval
   
                                    GOOD(i,currIndex)=1
                                    OUT(i,currIndex)=val
                                    if (VERBOSE) write(*,*)"Val: ", val
                                 else
   
                                    if (numtimes<=MAXCOUNT) then
   
                                       if (numpos+numneg<(in+outval-(numpos+numneg))) then
   
                                          if (loopcount<=2) then 
                                             val=missingVal !/* Try later */
   
                                          else
                                             !/* Keep the value after two passes through
                                             !** data */
   
                                             finalval=val
   
                                             rvVolume%vel(i,currIndex,sweepIndex)=finalval
                                             GOOD(i,currIndex)=1
                                             OUT(i,currIndex)=val
                                          endif
                                       else if (numpos>numneg) then
   
                                           val=val+NyqInterval
   
                                       else if (numneg>numpos) then
   
                                          val=val-NyqInterval
   
                                       else
                                       !/* Remove bin after four passes through data: */
                                          if (loopcount>4) GOOD(i,currIndex)=-1
                                       endif
                                    else
                                       !/* Remove bin: */
                                       GOOD(i,currIndex)=-1
                                    endif !numtimes<=MAXCOUNT
                                 endif ! in>outval
                              enddo !while(val/=missingVal .and. GOOD(i,currIndex)==0)
                           endif !countindex>=1
                        endif !GOOD(i,currIndex)==0
                     enddo !i=DELNUM+1,numBins
                  enddo !currIndex=startindex,endindex,step
               enddo !while (flag==1)
            endif !if(present(lastVolume) .or. present(soundVolume)
            
            deallocate (OUT)
            if(rvVolume%ifvel(sweepIndex))then
               abSweepIndex=sweepIndex
            endif
         enddo !sweepIndex=numSweeps,1,-1
   
         success=1
      else
         write(*,*)"First guess not available."
         success=0 
      endif
   
   contains
         integer function findRay (rvVolume1, rvVolume2, sweepIndex1, & 
                                    sweepIndex2, currIndex ) 
         implicit none
     
         type(t_radar_data), intent(in) :: rvVolume1, rvVolume2 
         integer, intent(in) :: sweepIndex1, sweepIndex2, currIndex
   
         integer ::  numSweeps, numRays, numBins, gatesz0, gatesz1, rayIndex1, i1
         real ::  prevval, az0, az1, diffaz, range0, range1
         real ::  startrange0,startrange1,spacing
         integer ::  direction, lastdir
        
         numRays = rvVolume2%nazim(sweepIndex2)
   
         az0 = rvVolume1%razim(currIndex,sweepIndex1)
   
         if (currIndex<numRays) then
             rayIndex1=currIndex
         else
             rayIndex1=numRays
         endif
         az1 = rvVolume2%razim(rayIndex1,sweepIndex2)
   
         if (az0==az1)then
            findRay= rayIndex1
            return
         endif
   
         !/* Compute the difference in azimuth between the two rays: */
   
         diffaz=az0-az1
         if (diffaz>=180.0) then
            diffaz=diffaz-360.0
         else if (diffaz<-180.0)then
            diffaz=diffaz+360.0
         endif
         
         !/* Get close to the correct index: */
         !/* Since the beamwidth is not necessarily the spacing between rays: */
   
         spacing = abs(rvVolume2%razim(1,sweepIndex2)-rvVolume2%razim(51,sweepIndex2))
   
         if (spacing>180) spacing=360.0-spacing
         spacing=spacing/50.0
         rayIndex1=rayIndex1+int(diffaz/spacing)
   
   !     write(*,*) numrays,rayIndex1,sweepIndex2,ubound(rvVolume2%razim)
         if (rayIndex1>numRays) rayIndex1=rayIndex1-numRays
         if (rayIndex1<1)       rayIndex1=rayIndex1+numRays
   !     write(*,*) numrays,rayIndex1,sweepIndex2,ubound(rvVolume2%razim)
         az1=rvVolume2%razim(rayIndex1,sweepIndex2)
   
         !/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   
         diffaz=az0-az1;
         if (diffaz>=180.0)then
            diffaz=diffaz-360.0
         else if (diffaz<-180.0)then
            diffaz=diffaz+360.0
         endif
   
         !/* Now add or subtract indices until the nearest ray is found: */
         if (diffaz>=0)then
            lastdir=1
         else 
            lastdir=-1
         endif
   
         do while (abs(diffaz)>spacing/2.0) 
   
            if (diffaz>=0) then
              rayIndex1=rayIndex1+1
              direction=1;
           
            else 
              rayIndex1=rayIndex1-1
              direction=-1
            endif
           
            if (rayIndex1>numRays) rayIndex1=rayIndex1-numRays
            if (rayIndex1<1)       rayIndex1=rayIndex1+numRays
            az1=rvVolume2%razim(rayIndex1,sweepIndex2)
           
            diffaz=az0-az1
            if (diffaz>=180.0)then
               diffaz=diffaz-360.0
            else if (diffaz<-180.0)then
               diffaz=diffaz+360.0
            endif
           
            if (direction/=lastdir)then
               exit
            else 
               lastdir=direction;
            endif
         enddo
   
         findRay= rayIndex1
         end function
   
   
         subroutine window(rvVolume, sweepIndex, startray, endray, firstbin, lastbin, ave, std, success)
         implicit none
         type(t_radar_data), intent(in) :: rvVolume
         integer, intent(in) :: sweepIndex, startray, endray, firstbin, lastbin
         real, intent(out) :: ave, std
         integer, intent(out) :: success
   
         integer :: num, currIndex, rangeIndex, numRays, numBins
         real :: val, sum, sumsq, NyqVelocity
        
         success=0
         NyqVelocity = rvVolume%vmax(sweepIndex)
   
         numRays = rvVolume%nazim(sweepIndex)
         numBins = rvVolume%nvgate(1,sweepIndex)
   
         !/* Now, sum the data in the window region between startray, 
         !**  endray, firstbin, lastbin. */
         std=0.0
         ave=0.0
         num=0
         sum=0.0
         sumsq=0.0
           
         if (firstbin>numBins .or. lastbin>numBins .or. firstbin<1 .or. lastbin<1)then
            ave= missingVal
            return
         endif
         if (startray>endray)then
            do currIndex=startray,numRays
               do rangeIndex=firstbin,lastbin
                  val=rvVolume%vel(rangeIndex,currIndex,sweepIndex)
                  if (val/=missingVal) then
                    num=num+1
                    sum=sum+val
                    sumsq=sumsq+val*val
                  endif
               enddo
            enddo
            do currIndex=1,endray
               do rangeIndex=firstbin,lastbin
                  val=rvVolume%vel(rangeIndex,currIndex,sweepIndex)
                  if (val/=missingVal) then
                     num=num+1
                     sum=sum+val
                     sumsq=sumsq+val*val
                  endif
               enddo
            enddo
         else
            do currIndex=startray,endray
               do rangeIndex=firstbin, lastbin
                  val=rvVolume%vel(rangeIndex,currIndex,sweepIndex)
                  if (val/=missingVal) then
                     num=num+1
                     sum=sum+val
                     sumsq=sumsq+val*val
                  endif
               enddo
            enddo
         endif
         if (num>=MINGOOD) then
            ave=sum/num
            std=sqrt(abs((sumsq-(sum*sum)/num)/(num-1)))
            if (std<=STDTHRESH*NyqVelocity)then
               success=1
               !/* printf("ave=%0.2f, std=%0.2f, sum=%0.2f\n", ave, *std, sum); */
            endif
         else 
            ave=missingVal
            std=0.0
         endif
         end subroutine
   end subroutine
end module
