module radar_data
use date_pack, only: date
implicit none

   integer, parameter :: radar_unit = 901
   integer, parameter :: RGates_SA =  460, VGates_SA =  920, WGates_SA =  920, MaxCuts_SA = 20 !SA, SB
   integer, parameter :: RGates_CB =  800, VGates_CB = 1600, WGates_CB = 1600, MaxCuts_CB = 20 !CA, CB
   integer, parameter :: RGates_SC = 1000, VGates_SC = 1000, WGates_SC = 1000, MaxCuts_SC = 30 !SC, CD
   integer, parameter :: RGates_CC =  500, VGates_CC =  500, WGates_CC =  500, MaxCuts_CC = 30 !CC, CCJ
   
   ! 98D
   integer, parameter ::   CODE_INVALID = 0,  CODE_RANFOLD = 1
   integer, parameter :: RES_POINT_FIVE = 2, RES_ONE_POINT = 4
   integer, parameter :: VOL_BEG=3, VOL_END=4, ELV_BEG=0, ELV_END=2

   ! General Usage
   real   , parameter :: VALUE_INVALID=-999., VALUE_RANFOLD=999. 
   real   , parameter :: RADIAN=3.14159/180.
   
   integer :: MaxRads, MaxCuts, RGates, VGates, WGates
   type(date), parameter :: d1970=date(1970,1,1,0,0,0,0)

   type t_radar_data
      
      character(len=20)  :: radar_name
      character(len=6 )  :: radar_type, file_format, radar_id
   
      real           :: latitude
      real           :: longitude
      real           :: altitude
   
      real           :: calibConst
      integer        :: vcp
      integer        :: year, month, day, hour, minute, second
   
      ! Dimension Gate, Azim, Tilt
      integer                                           :: ntilt
      integer      , dimension(:    )     , allocatable :: nazim             ! azimuth number of each tilt
      logical      , dimension(:    )     , allocatable :: ifref, ifvel, ifzdr, ifcc, iffdp, ifkdp, ifsnr
      integer      , dimension(: , :)     , allocatable :: nrgate  , nvgate  ! Velocity Resolution
      real         , dimension(:    )     , allocatable :: rgatesp , vgatesp, atmosAttenFactor
      real         , dimension(:    )     , allocatable ::  rmax   , vmax
      real         , dimension(: , :)     , allocatable :: rtilt   , razim   ! elevation angle and azimuth angle of each ray
      real(kind=8) , dimension(: , :)     , allocatable :: stime   , etime   ! scan start and end time: seconds since 1970-1-1 00:00:00.00
      real         , dimension(: , :)     , allocatable :: vres              ! Radical velocity resolution , (metstar 98d)
      real         , dimension(: , : , :) , allocatable :: vel  , ref  , spw, zdr, kdp, cc, fdp, snr
      real         , dimension(: , : , :) , allocatable :: vlon , vlat , valt
      real         , dimension(: , : , :) , allocatable :: rlon , rlat , ralt
      logical                                           :: if_have_loc = .false.
   end type


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_radar_filename_info(filename, radar_id, radar_time, radar_type)
   implicit none
   character(len=*), intent(in) :: filename
   character(len=*), intent(out):: radar_id, radar_time, radar_type

   character(len=46) :: base_name
   integer :: is, ie
  
   is=index(filename, "/", back=.TRUE.)
   is=is+1
   base_name =filename(is:)
   !write(*,*) trim(filename), is, trim(base_name)
   radar_id  ="Z----"

   if((base_name(1:6)/="Z_RADR").and.(base_name(1:6)/="T_RADR"))then
      return
   endif

!Z_RADR_I_Z9477_20130816000600_O_DOR_CB_CAP.bin
!T_RADR_I_Z9421_20170815000103_O_DOR_CD_CAP_FMT.bin
   radar_id  =base_name(10:14)
   radar_time=base_name(16:29)
   radar_type=base_name(37:38)
   
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_radar_data(filename,radar_type,radar_data)
   implicit none

   character(len=*),   intent(in)    :: filename
   character(len=*),   intent(inout)    :: radar_type
   type(t_radar_data), intent(inout) :: radar_data

   character(len=112) :: head
   integer          :: i,j,k

   integer :: nrad, nrmax, nvmax, nrgate, nvgate
   real :: min_z, max_z, min_v, max_v, min_w, max_w


   radar_data%ntilt=0
   radar_data%radar_type=radar_type

   open(radar_unit,file=filename,recl=112,status="old",access="direct")
   read(radar_unit,rec=1) head
   close(radar_unit)

   if(head(1:2)=="RD")then
      radar_data%file_format="CC2.0"
   elseif(head(1:7)=="CINRADC")then
      radar_data%file_format="CC1.0"
   elseif(head(101:109)=="CINRAD/CD".or.head(101:109)=="CINRAD/SC")then
      radar_data%file_format="784"
   elseif(head(1:4)=="RSTM")then
      radar_data%file_format="RSTM"
   else
      radar_data%file_format="98D"
   endif   
   !  radar_data%file_format="RSTM"

   select case(trim(radar_data%file_format))
   case ("98D")
      call read_radar_98d(filename,radar_type,radar_data)
   case ("CC2.0")
      call read_radar_rd (filename,radar_type,radar_data)
   case ("CC1.0")
      call read_radar_38 (filename,radar_type,radar_data)
   case ("784")
      call read_radar_784(filename,radar_type,radar_data)
   case ("RSTM")
      call read_radar_rstm(filename,radar_type,radar_data)
   case default
      write(*,*) "in read_radar_data:Unknown Radar Type:", trim(radar_type)
      radar_data%ntilt=0
      return
      !stop
   end select

   write(*,"(A,I4.4,5I2.2 )") "Radar Date :", radar_data%year   , &
                                              radar_data%month  , &
                                              radar_data%day    , &
                                              radar_data%hour   , &
                                              radar_data%minute , &
                                              radar_data%second

   write(*,"(A,2F12.4,F8.1)") "Lat,Lon,Alt:", radar_data%latitude, radar_data%longitude, radar_data%altitude 
   write(*,"(1(A,I10),2A  )") "Total Tilts:", radar_data%ntilt, ". File Format:",trim(radar_data%file_format)
   write(*,"(1X,A9,6(',',1X,A3),3(',',2X,A8),4(',',1X,A11),2(',',1X,A5))") &
           "Elevation","Z","V","Zdr","CC","Kdp","SNR","NAzimuth","NRefGate","NVelGate","RefGateSize","VelGateSize",&
           "MaxRange(m)","MaxVel(m/s)","NRMax","NVMax"
  
   nrgate=ubound(radar_data%ref,dim=1)
   nvgate=ubound(radar_data%vel,dim=1)
   nrad  =ubound(radar_data%vel,dim=2)
   do k=1, radar_data%ntilt
      if(.not.radar_data%ifvel(k)) radar_data%vmax(k)=VALUE_INVALID
      do j=radar_data%nazim(k)+1,nrad
         radar_data%razim(j,k)=VALUE_INVALID
         radar_data%rtilt(j,k)=VALUE_INVALID
      enddo
      nrmax=-999
      if(radar_data%ifref(k))then
         nrmax=int(radar_data%rmax(k)/radar_data%rgatesp(k)+0.5)
      !  write(*,*) nrmax, nrgate
         if(nrmax<nrgate)then
            radar_data%ref(nrmax+1:nrgate,:,k)=VALUE_INVALID
         endif
      endif
      nvmax=-999
      if(radar_data%ifvel(k))then
         nvmax=int(radar_data%rmax(k)/radar_data%vgatesp(k)+0.5)
         if(nvmax<nvgate)then
            radar_data%vel(nvmax+1:nvgate,:,k)=VALUE_INVALID
            radar_data%spw(nvmax+1:nvgate,:,k)=VALUE_INVALID
         endif
      endif
      write(*,"(F10.4,6(',',L4),3(',',I10),4(',',F12.2),2(',',I6))") & 
                         radar_data%rtilt (1,k), &
                         radar_data%ifref   (k), &
                         radar_data%ifvel   (k), &
                         radar_data%ifzdr   (k), &
                         radar_data%ifcc    (k), &
                         radar_data%ifkdp   (k), &
                         radar_data%ifsnr   (k), &
                         radar_data%nazim   (k), &
                         radar_data%nrgate(1,k), &
                         radar_data%nvgate(1,k), &
                         radar_data%rgatesp (k), &
                         radar_data%vgatesp (k), &
                         radar_data%rmax    (k), &
                         radar_data%vmax    (k), nrmax, nvmax
   enddo

   ! write(*,"(1X,A9,6(',',1X,A9))") "Elevation", "Zmin","Zmax","Vmin","Vmax","Wmin","Wmax"
   ! do k=1, radar_data%ntilt
   !    !min_z= 999.
   !    !max_z=-999.
   !    !min_v= 999.
   !    !max_v=-999.
   !    !min_w= 999.
   !    !max_w=-999.
   !    !do j=1, radar_data%nazim   (k)
   !    !   do i=1, radar_data%nrgate(j,k)
   !    !      if(radar_data%ref(i,j,k)/=VALUE_INVALID.and.radar_data%ref(i,j,k)/=VALUE_RANFOLD)then
   !    !         if(radar_data%ref(i,j,k)<min_z)then
   !    !            write(701,*) "min_z",i,j,k,min_z,radar_data%ref(i,j,k)
   !    !            min_z=radar_data%ref(i,j,k)
   !    !         endif
   !    !         if(radar_data%ref(i,j,k)>max_z)then
   !    !            write(701,*) "max_z",i,j,k,max_z,radar_data%ref(i,j,k)
   !    !            max_z=radar_data%ref(i,j,k)
   !    !         endif
   !    !      endif
   !    !   enddo
   !    !   do i=1, radar_data%nvgate(j,k)
   !    !      if(radar_data%vel(i,j,k)/=VALUE_INVALID.and.radar_data%vel(i,j,k)/=VALUE_RANFOLD)then
   !    !         if(radar_data%vel(i,j,k)<min_v)then
   !    !            write(701,*) "min_v",i,j,k,min_v,radar_data%vel(i,j,k)
   !    !            min_v=radar_data%vel(i,j,k)
   !    !         endif
   !    !         if(radar_data%vel(i,j,k)>max_v)then
   !    !            write(701,*) "max_v",i,j,k,max_v,radar_data%vel(i,j,k)
   !    !            max_v=radar_data%vel(i,j,k)
   !    !         endif
   !    !      endif
   !    !      if(radar_data%spw(i,j,k)/=VALUE_INVALID.and.radar_data%spw(i,j,k)/=VALUE_RANFOLD)then
   !    !         if(radar_data%spw(i,j,k)<min_w)then
   !    !            write(701,*) "min_w",i,j,k,min_w,radar_data%spw(i,j,k)
   !    !            min_w=radar_data%spw(i,j,k)
   !    !         endif
   !    !         if(radar_data%spw(i,j,k)>max_w)then
   !    !            write(701,*) "max_w",i,j,k,max_w,radar_data%spw(i,j,k)
   !    !            max_w=radar_data%spw(i,j,k)
   !    !         endif
   !    !      endif
   !    !   enddo
   !    !enddo
   !    !write(*,"(F10.2,6(',',F10.2))") radar_data%rtilt (1,k), min_z, max_z, min_v, max_v, min_w, max_w
   !    call get_maxmin_2d(radar_data%ref(:,:,k),min_z,max_z)
   !    call get_maxmin_2d(radar_data%vel(:,:,k),min_v,max_v)
   !    call get_maxmin_2d(radar_data%spw(:,:,k),min_w,max_w)
   !    write(*,"(F10.2,6(',',F10.2))") radar_data%rtilt (1,k), min_z, max_z, min_v, max_v, min_w, max_w
   ! enddo
!  stop

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_maxmin_2d(dat,minv,maxv)
  implicit none
  real, dimension(:,:), intent(in) :: dat
  real, intent(out) :: minv,maxv

  integer, dimension(2) :: n
  integer :: i, j

   n=ubound(dat)

   minv=999.
   maxv=-999.
   do i=1, n(1)
      do j=1, n(2)
         if(dat(i,j)/=VALUE_INVALID.and.dat(i,j)/=VALUE_RANFOLD.and.dat(i,j)>-200.)then
            if(dat(i,j)<minv)then
               minv=dat(i,j)
            endif
            if(dat(i,j)>maxv)then
               maxv=dat(i,j)
            endif
         endif
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist_3d(dat,dist) 
  implicit none
  real, dimension(:,:,:), intent(in) :: dat 
  integer, dimension(201),intent(out) :: dist

  integer, dimension(3) :: n
  integer :: i, j, k, l
  real :: r1, r2, minr, maxr

   n=ubound(dat)

   minr=-100
   maxr= 100
   dist=0
   do i=1, n(1)
      do j=1, n(2)
         do k=1, n(3)
            if(dat(i,j,k)>=(minr-0.5).and.dat(i,j,k)<(maxr+0.5))then
            l=int(dat(i,j,k)-minr+0.5)+1
            dist(l)=dist(l)+1
            endif
         enddo
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist_1d_int(n1,dat,dist)
  implicit none
  integer,               intent(in) :: n1
  integer, dimension(:), intent(in) :: dat
  integer, dimension(201),intent(out) :: dist

  integer, dimension(1) :: n
  integer :: i, j, k, l
  real :: r1, r2, minr, maxr

   minr=-100
   maxr= 100
   dist=0
   do i=1, n1
      if(dat(i)>=(minr-0.5).and.dat(i)<(maxr+0.5))then
      l=int(dat(i)-minr+0.5)+1
      dist(l)=dist(l)+1
      endif
   enddo
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mean_1d(n1,dat,mdat)
  implicit none
  integer,            intent(in) :: n1
  real, dimension(:), intent(in) :: dat
  real,               intent(out) :: mdat

  integer, dimension(1) :: n
  integer :: i, j, m 

   m=0
   mdat=0
   do i=1, n1
      if(dat(i)/=VALUE_INVALID.and.dat(i)/=VALUE_RANFOLD.and.dat(i)>-200.)then
         m=m+1
         mdat=mdat+dat(i) 
      endif
   enddo
   if(m>0)then
      mdat=mdat/m
   else
      mdat=value_invalid
   endif
   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allocate_radar_data(radar_data,RGates,VGates,WGates,MaxRads,MaxCuts)
   implicit none

   type(t_radar_data), intent(inout) :: radar_data
   integer           , intent(in)    :: RGates,VGates,WGates,MaxRads,MaxCuts
 
   !write(*,*) "in allocate"
   if(allocated(radar_data%vel))then
      deallocate(radar_data%nazim           )
      deallocate(radar_data%ifref           )
      deallocate(radar_data%ifvel           )
      deallocate(radar_data%ifzdr           )
      deallocate(radar_data%ifcc            )
      deallocate(radar_data%ifkdp           )
      deallocate(radar_data%ifsnr           )
      deallocate(radar_data%vmax            )
      deallocate(radar_data%rmax            )
      deallocate(radar_data%rgatesp         )
      deallocate(radar_data%vgatesp         )
      deallocate(radar_data%rtilt           )
      deallocate(radar_data%nrgate          )
      deallocate(radar_data%nvgate          )
      deallocate(radar_data%stime           )
      deallocate(radar_data%etime           )
      deallocate(radar_data%vres            )
      deallocate(radar_data%razim           )
      deallocate(radar_data%ref             )
      deallocate(radar_data%vel             )
      deallocate(radar_data%spw             )
      deallocate(radar_data%rlon            )
      deallocate(radar_data%rlat            )
      deallocate(radar_data%ralt            )
      deallocate(radar_data%vlon            )
      deallocate(radar_data%vlat            )
      deallocate(radar_data%valt            )
      deallocate(radar_data%atmosAttenFactor)
      deallocate(radar_data%zdr             )
      deallocate(radar_data%cc              )
      deallocate(radar_data%fdp             )
      deallocate(radar_data%kdp             )
      deallocate(radar_data%snr             )
   endif
   !write(*,*) "deallocate complete!!!"
   allocate(radar_data%nazim                   (MaxCuts))
   allocate(radar_data%ifvel                   (MaxCuts))
   allocate(radar_data%ifref                   (MaxCuts))
   allocate(radar_data%ifzdr                   (MaxCuts))
   allocate(radar_data%ifcc                    (MaxCuts))
   allocate(radar_data%ifkdp                   (MaxCuts))
   allocate(radar_data%ifsnr                   (MaxCuts))
   allocate(radar_data%vmax                    (MaxCuts))
   allocate(radar_data%rmax                    (MaxCuts))
   allocate(radar_data%rgatesp                 (MaxCuts))
   allocate(radar_data%vgatesp                 (MaxCuts))
   allocate(radar_data%atmosAttenFactor        (MaxCuts))
   allocate(radar_data%rtilt          (MaxRads, MaxCuts))
   allocate(radar_data%nrgate         (MaxRads, MaxCuts))
   allocate(radar_data%nvgate         (MaxRads, MaxCuts))
   allocate(radar_data%stime          (MaxRads, MaxCuts))
   allocate(radar_data%etime          (MaxRads, MaxCuts))
   allocate(radar_data%vres           (MaxRads, MaxCuts))
   allocate(radar_data%razim          (MaxRads, MaxCuts))
   allocate(radar_data%ref    (RGates, MaxRads, MaxCuts))
   allocate(radar_data%vel    (VGates, MaxRads, MaxCuts))
   allocate(radar_data%spw    (WGates, MaxRads, MaxCuts))
   allocate(radar_data%rlon   (RGates, MaxRads, MaxCuts))
   allocate(radar_data%rlat   (RGates, MaxRads, MaxCuts))
   allocate(radar_data%ralt   (RGates, MaxRads, MaxCuts))
   allocate(radar_data%vlon   (VGates, MaxRads, MaxCuts))
   allocate(radar_data%vlat   (VGates, MaxRads, MaxCuts))
   allocate(radar_data%valt   (VGates, MaxRads, MaxCuts))
   allocate(radar_data%zdr    (RGates, MaxRads, MaxCuts))
   allocate(radar_data%cc     (RGates, MaxRads, MaxCuts))
   allocate(radar_data%fdp    (RGates, MaxRads, MaxCuts))
   allocate(radar_data%kdp    (RGates, MaxRads, MaxCuts))
   allocate(radar_data%snr    (RGates, MaxRads, MaxCuts))

   !write(*,*) "allocate complete!!!"
   radar_data%ref   = VALUE_INVALID
   radar_data%vel   = VALUE_INVALID
   radar_data%spw   = VALUE_INVALID
   radar_data%zdr   = VALUE_INVALID
   radar_data%cc    = VALUE_INVALID
   radar_data%fdp   = VALUE_INVALID
   radar_data%kdp   = VALUE_INVALID
   radar_data%snr   = VALUE_INVALID
   radar_data%rlon  = VALUE_INVALID
   radar_data%rlat  = VALUE_INVALID
   radar_data%ralt  = VALUE_INVALID
   radar_data%vlon  = VALUE_INVALID
   radar_data%vlat  = VALUE_INVALID
   radar_data%valt  = VALUE_INVALID

   radar_data%nazim            = 0
   radar_data%nvgate           = 0
   radar_data%nrgate           = 0
   radar_data%stime            = 0
   radar_data%etime            = 0
   radar_data%vres             = 1
   radar_data%atmosAttenFactor = 0
   radar_data%calibConst       = 0

   radar_data%if_have_loc=.false.

   radar_data%ifref=.false.
   radar_data%ifvel=.false.
   radar_data%ifzdr=.false.
   radar_data%ifcc =.false.
   radar_data%ifkdp=.false.
   radar_data%ifsnr=.false.
   !write(*,*) "out allocate"
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_radar_grads(filename, radar_data)
   implicit none
   character(len=*),   intent(in) :: filename
   type(t_radar_data), intent(in) :: radar_data

   character(len=800) :: ctlfile, datfile
   real, dimension(:,:,:), allocatable :: ref, vel, spw

   integer :: max_ngate, max_nazim, ntilt, i, j, k

   max_ngate=max(maxval(radar_data%nrgate),maxval(radar_data%nvgate))
   max_nazim=maxval(radar_data%nazim(1:radar_data%ntilt))
   ntilt=radar_data%ntilt

   allocate(ref(max_ngate,max_nazim,ntilt))
   allocate(vel(max_ngate,max_nazim,ntilt))
   allocate(spw(max_ngate,max_nazim,ntilt))

   ref=value_invalid
   vel=value_invalid
   spw=value_invalid


   do k=1, ntilt
      do j=1,radar_data%nazim(k)
         do i=1,radar_data%nvgate(j,k)
            vel(i,j,k)=radar_data%vel(i,j,k)
            spw(i,j,k)=radar_data%spw(i,j,k)
         enddo
         !write(*,*) radar_data%nrgate(j,k)
         do i=1,radar_data%nrgate(j,k)
            ref(i,j,k)=radar_data%ref(i,j,k)
         enddo
      enddo
   enddo

   ctlfile=trim(filename)//".ctl"
   datfile=trim(filename)//".dat"

   OPEN(radar_unit,FILE=ctlfile,STATUS='unknown')
   write(*,*) "Writing grads ctl file:", trim(ctlfile)
   write(radar_unit,"(A      )") "DSET ^"//trim(datfile)
   write(radar_unit,"(A      )") "options big_endian sequential"
   write(radar_unit,"(A,F10.1)") "UNDEF ", value_invalid
   write(radar_unit,"(A      )") "TITLE  Sample Model Output"
   write(radar_unit,"(A,I5,A,F8.3 )") "XDEF ",max_nazim," linear 0.0 ",360./max_nazim
   write(radar_unit,"(A,I5,A )") "YDEF ",max_ngate," linear -90.0 0.001"
   write(radar_unit,"(A,I5,A )") "ZDEF ",ntilt    ," linear 1.00  1.00"
   write(radar_unit,"(A      )") "TDEF 1    LINEAR 00z01jan1995   1hr"
   write(radar_unit,"(A      )") "VARS 3"
   write(radar_unit,"(A,I5,A )") "Z ",ntilt," 99 Reflectivity dBz"
   write(radar_unit,"(A,I5,A )") "V ",ntilt," 99 Radical velocity m/s"
   write(radar_unit,"(A,I5,A )") "W ",ntilt," 99 Spectra Width"
   write(radar_unit,"(A      )") "ENDVARS"
   close(radar_unit)

   OPEN(radar_unit,FILE=datfile,STATUS='unknown',form="unformatted")
   write(*,*) "Writing grads dat file:", trim(datfile)
   WRITE(radar_unit) (((ref(i,j,k),j=1,max_nazim),i=1,max_ngate),k=1,ntilt) 
   WRITE(radar_unit) (((vel(i,j,k),j=1,max_nazim),i=1,max_ngate),k=1,ntilt) 
   WRITE(radar_unit) (((spw(i,j,k),j=1,max_nazim),i=1,max_ngate),k=1,ntilt) 
   close(radar_unit)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_csv(filename, radar_data)
   implicit none
   character(len=*),   intent(in) :: filename
   type(t_radar_data), intent(in) :: radar_data

   character(len=800) :: fileout

   integer :: i, j, k, n, iref, ivel, jvstart, iv, kv, jv
   real    :: ivratio, jvratio
   real, dimension(:), allocatable :: vel

   character(len=200) :: cNGate, cNAzim, strFormat

   n=0
   do k=1, radar_data%ntilt

      if(radar_data%ifvel(k))then
         n=n+1
         if(n>2)exit
         write(fileout,"(A,I2.2,A)") trim(filename)//".vel-",n,".csv"
         open(radar_unit,file=fileout,status="unknown")
         write(*,*) "Writing vel csv file:", trim(fileout)
         write(strFormat,*) "(A,",2*radar_data%nvgate(1,k),"(',',A,I4))"
         write(radar_unit,FMT=strFormat) "#,Ang",("V",i,"W",i,i=1,radar_data%nvgate(1,k))
         do j=1,radar_data%nazim(k)
            write(strFormat,*) "(I6,',',F8.3,",2*radar_data%nvgate(j,k),"(',',F8.2))"
            write(radar_unit,fmt=strFormat) j,radar_data%razim(j,k),&
                 (radar_data%vel(i,j,k),radar_data%spw(i,j,k),i=1,radar_data%nvgate(j,k))
         enddo
         close(radar_unit)
      endif
   enddo

   n=0
   do k=1, radar_data%ntilt

      if(radar_data%ifref(k))then
         n=n+1
         if(n>2)exit
         if(.NOT.radar_data%ifvel(k))then 
            kv=k+1
            jvratio=radar_data%nazim(k)/radar_data%nazim(kv)
            jvstart=1
            do jv=2,radar_data%nazim(kv)
               if(((radar_data%razim(jv  ,kv)- radar_data%razim( 1,k ))*   &
                   (radar_data%razim(jv-1,kv)- radar_data%razim( 1,k ))<=0.).and. &
                ABS(radar_data%razim(jv-1,kv)- radar_data%razim(jv,kv))<2*360./radar_data%nazim(kv).and. &
                   (radar_data%razim(jv-1,kv)/=radar_data%razim( 1,k )))then
                  jvstart=jv
                  exit
               endif
            enddo
            !write(*,*) "k,kv,jvstart:",k,kv,jvstart,&
            ! radar_data%razim(jvstart-1,kv),radar_data%razim(jvstart,kv),radar_data%razim(1,k)
         else
            kv=k
            jvratio=1
            jvstart=1
         endif
         ivratio=radar_data%rgatesp(k)/radar_data%vgatesp(kv)
         allocate(vel(radar_data%nrgate(1,k)))
         write(fileout,"(A,I2.2,A)") trim(filename)//".ref-",n,".csv"
         open(radar_unit,file=fileout,status="unknown")
         write(*,*) "Writing ref csv file:", trim(fileout)
         write(strFormat,*) "(A,",2*radar_data%nrgate(1,k),"(',',A,I4))"
         write(radar_unit,FMT=strFormat) "#,Ang1,Ang2",("Z",i,"V",i,i=1,radar_data%nrgate(1,k))
         do j=1,radar_data%nazim(k)
            vel=value_invalid
            jv=j*jvratio+jvstart-1
            if(jv>radar_data%nazim(kv))jv=jv-radar_data%nazim(kv)
            do i=1, radar_data%nrgate(j,k)
               iv=i*ivratio
               if(iv<=radar_data%nrgate(jv,kv)) vel(i)=radar_data%vel(iv,jv,kv)
            enddo
            write(strFormat,*) "(I6,2(',',F8.3),",2*radar_data%nrgate(j,k),"(',',F8.2))"
            write(radar_unit,fmt=strFormat) j,radar_data%razim(j,k),radar_data%razim(jv,kv),&
                 (radar_data%ref(i,j,k),vel(i),i=1,radar_data%nrgate(j,k))
         enddo
         deallocate(vel)
         close(radar_unit)
      endif
   enddo

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grads_grid2d(filename, nr, na, dat, missv)
   implicit none
   character(len=*),   intent(in) :: filename
   integer, intent(in) :: nr, na
   real, dimension(nr,na) :: dat
   real, intent(in) :: missv

   character(len=800) :: fileout, ctlfile, datfile, mapfile

   integer :: i, j, k, n
   real    :: lat1, lon1, alt, lat0, lon0, alt0, range, height, elev, sfcrng, azim

   integer           :: NFLAG, NLEV
   character(len=8)  :: STID
   real              :: TIM, VEL, SPW


   write(ctlfile,"(A,I2.2,A)") trim(filename)//".ctl"
   write(datfile,"(A,I2.2,A)") trim(filename)//".dat"

   OPEN(radar_unit,FILE=ctlfile,STATUS='unknown')
   write(*,*) "Writing grads ctl file:", trim(ctlfile)
   write(radar_unit,"(A      )") "DSET ^"//trim(datfile)
   write(radar_unit,"(A      )") "options big_endian sequential"
   write(radar_unit,"(A,F10.1)") "UNDEF ", value_invalid
   write(radar_unit,"(A      )") "TITLE  Sample Model Output"
   write(radar_unit,"(A,I5,A )") "XDEF ",nr," linear 1 ", nr
   write(radar_unit,"(A,I5,A )") "YDEF ",na," linear 1 ", na
   write(radar_unit,"(A,I5,A )") "ZDEF ",1    ," linear 1.00  1.00"
   write(radar_unit,"(A      )") "TDEF 1 1LINEAR 00z01jan1995   1hr"
   write(radar_unit,"(A      )") "VARS 1"
   write(radar_unit,"(A,I5,A )") "dat 0 99 Reflectivity dBz"
   write(radar_unit,"(A      )") "ENDVARS"
   close(radar_unit)

   OPEN(radar_unit,FILE=datfile,STATUS='unknown',form="unformatted")
   write(*,*) "Writing grads dat file:", trim(datfile)
   WRITE(radar_unit) ((dat(i,j),i=1,nr),j=1,na)
   close(radar_unit)

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   subroutine write_grads_lld(filename, nr, na, lat, lon, dat, missv)
!   implicit none
!   character(len=*),   intent(in) :: filename
!   integer, intent(in) :: nr, na
!   real, dimension(nr,na) :: lat, lon, dat
!   real, intent(in) :: missv
!
!   character(len=800) :: fileout, ctlfile, datfile, mapfile
!
!   integer :: i, j, k, n
!   real    :: lat1, lon1, alt, lat0, lon0, alt0, range, height, elev, sfcrng, azim
!
!   integer           :: NFLAG, NLEV
!   character(len=8)  :: STID
!   real              :: TIM, VEL, SPW
!
!
!   write(ctlfile,"(A,I2.2,A)") trim(filename)//".ctl"
!   write(datfile,"(A,I2.2,A)") trim(filename)//".dat"
!   write(mapfile,"(A,I2.2,A)") trim(filename)//".map"
!   open(radar_unit,file=ctlfile,status="unknown")
!   write(radar_unit,"(A)") "DSET ^"//trim(datfile)
!   write(radar_unit,"(A)") "DTYPE  station"
!   write(radar_unit,"(A)") "OPTIONS sequential big_endian"
!   write(radar_unit,"(A)") "STNMAP "//trim(mapfile) 
!   write(radar_unit,"(A,F15.1)") "UNDEF  ",missv
!   write(radar_unit,"(A)") "TITLE  Station Data Sample"
!   write(radar_unit,"(A)") "TDEF   1 linear 00z01Jan2000 1hr"
!   write(radar_unit,"(A)") "VARS 1"
!   write(radar_unit,"(A)") "dat    0  99  dat "
!   write(radar_unit,"(A)") "ENDVARS"
!
!   open(radar_unit,file=datfile,form='unformatted',status="unknown")
!   write(*,"(2A)") "Writing file:", trim(datfile)
!   TIM = 0.0
!   NLEV = 1
!   NFLAG = 1
!   do i=1, nr
!   do j=1, na
!      write(STID,'(I4,I4)') i
!      write(radar_unit) STID,lat(i,j),lon(i,j),TIM,NLEV,NFLAG
!      write(radar_unit) dat(i,j)
!   enddo
!   enddo
!   NLEV = 0
!   WRITE(radar_unit) STID,lat(1,1),lon(1,1),TIM,NLEV,NFLAG
!   close(radar_unit)
!
!   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_radar_38(filename,radar_type,rdat)
   use bytes, only: get_uc_value, get_i_value, get_us_value, get_s_value
   use date_pack, only: init_date, get_new_date, get_diff_date, one_hour
   implicit none

   character(len=*),   intent(in)    :: filename, radar_type
   type(t_radar_data), intent(inout) :: rdat
 

   type t_VPPISCANPARAMETER
      character(len=2 ) :: us_MaxV
      character(len=2 ) :: us_MaxL
      character(len=2 ) :: us_BindWidth
      character(len=2 ) :: us_BinNumber
      character(len=2 ) :: us_RecordNumber
      character(len=2 ) :: us_Arotate
      character(len=2 ) :: us_Prf1
      character(len=2 ) :: us_Prf2
      character(len=2 ) :: us_SpulseW
      character(len=2 ) :: us_Angle
      character(len=1 ) :: uc_SweepStatus
      character(len=1 ) :: uc_Ambiguousp
   end type

   type t_header_38
      character(len=16)         ::  c_FileType
      character(len=30)         ::  c_Country
      character(len=20)         ::  c_Province
      character(len=40)         ::  c_Station
      character(len=10)         ::  c_StationNumber
      character(len=20)         ::  c_RadarType
      character(len=16)         ::  c_Longitude
      character(len=16)         ::  c_Latitude
      character(len=4 )         ::  l_LongitudeValue
      character(len=4 )         ::  l_LatitudeValue
      character(len=4 )         ::  l_Height
      character(len=2 )         ::  s_MaxAngle
      character(len=2 )         ::  s_OptAngle
      character(len=1 )         :: uc_SYear1
      character(len=1 )         :: uc_SYear2
      character(len=1 )         :: uc_SMonth
      character(len=1 )         :: uc_SDay
      character(len=1 )         :: uc_SHour
      character(len=1 )         :: uc_SMinute
      character(len=1 )         :: uc_SSecond
      character(len=1 )         :: uc_TimeFrom
      character(len=1 )         :: uc_EYear1
      character(len=1 )         :: uc_EYear2
      character(len=1 )         :: uc_EMonth
      character(len=1 )         :: uc_EDay
      character(len=1 )         :: uc_EHour
      character(len=1 )         :: uc_EMinute
      character(len=1 )         :: uc_ESecond
      character(len=1 )         :: uc_ScanMode
      character(len=4 )         :: ul_SmilliSecond
      character(len=2 )         :: us_RHIA
      character(len=2 )         ::  s_RHIL
      character(len=2 )         ::  s_RHIH
      character(len=2 )         :: us_EchoType
      character(len=2 )         :: us_ProdCode
      character(len=1 )         :: uc_Calibration
      character(len=3 )         :: remain1
      type(t_VPPISCANPARAMETER) :: LayerInfo(30)
      character(len=4 )         :: lAntennaG
      character(len=4 )         :: lPower
      character(len=4 )         :: lWavelength
      character(len=2 )         :: us_BeamH
      character(len=2 )         :: us_BeamL
      character(len=2 )         :: us_Polarization
      character(len=2 )         :: us_LogA
      character(len=2 )         :: us_LineA
      character(len=2 )         :: us_AGCP
      character(len=2 )         :: us_FreqMode
      character(len=2 )         :: us_FreqRepeat
      character(len=2 )         :: us_PPPPulse
      character(len=2 )         :: us_FFTPoint
      character(len=2 )         :: us_ProcessType
      character(len=1 )         :: uc_ClutterT
      character(len=1 )         :: c_Sidelobe
      character(len=1 )         :: uc_VelocityT
      character(len=1 )         :: uc_FilderP
      character(len=1 )         :: uc_NoiseT
      character(len=1 )         :: uc_SQIT
      character(len=1 )         :: uc_IntensityC
      character(len=1 )         :: uc_IntensityR
      character(len=1 )         :: uc_CalNoise
      character(len=1 )         :: uc_CalPower
      character(len=1 )         :: uc_CalPulseWidth
      character(len=1 )         :: uc_CalWorkFreq
      character(len=1 )         :: uc_CalLog
      character(len=92)         :: remain3
      character(len=4 )         :: ul_DataOffset
   end type ! 1021 bytes 

   type(t_header_38) :: header_38

   type t_data_record
       character(len=2), dimension(:), allocatable :: z, v, w
   end type
   type(t_data_record), dimension(:,:), allocatable :: dat

   character(len=1024) :: header

   integer      :: TotalRays, i, j, k, n, ierr, value, FileSize
   type(date)   :: date_bj, date_utc, sdate, edate
   real         :: dtime
   real(kind=8) :: btime
   
   ! read in header only
   open(radar_unit,file=filename,recl=1021,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr)  &
                   header_38% c_FileType, &
                   header_38% c_Country, &
                   header_38% c_Province, &
                   header_38% c_Station, &
                   header_38% c_StationNumber, &
                   header_38% c_RadarType, &
                   header_38% c_Longitude, &
                   header_38% c_Latitude, &
                   header_38% l_LongitudeValue, &
                   header_38% l_LatitudeValue, &
                   header_38% l_Height, &
                   header_38% s_MaxAngle, &
                   header_38% s_OptAngle, &
                   header_38%uc_SYear1, &
                   header_38%uc_SYear2, &
                   header_38%uc_SMonth, &
                   header_38%uc_SDay, &
                   header_38%uc_SHour, &
                   header_38%uc_SMinute, &
                   header_38%uc_SSecond, &
                   header_38%uc_TimeFrom, &
                   header_38%uc_EYear1, &
                   header_38%uc_EYear2, &
                   header_38%uc_EMonth, &
                   header_38%uc_EDay, &
                   header_38%uc_EHour, &
                   header_38%uc_EMinute, &
                   header_38%uc_ESecond, &
                   header_38%uc_ScanMode, &
                   header_38%ul_SmilliSecond, &
                   header_38%us_RHIA, &
                   header_38% s_RHIL, &
                   header_38% s_RHIH, &
                   header_38%us_EchoType, &
                   header_38%us_ProdCode, &
                   header_38%uc_Calibration, &
                   header_38%remain1, &
                   header_38%LayerInfo    , &
                   header_38%lAntennaG, &
                   header_38%lPower, &
                   header_38%lWavelength, &
                   header_38%us_BeamH, &
                   header_38%us_BeamL, &
                   header_38%us_Polarization, &
                   header_38%us_LogA, &
                   header_38%us_LineA, &
                   header_38%us_AGCP, &
                   header_38%us_FreqMode, &
                   header_38%us_FreqRepeat, &
                   header_38%us_PPPPulse, &
                   header_38%us_FFTPoint, &
                   header_38%us_ProcessType, &
                   header_38%uc_ClutterT, &
                   header_38%c_Sidelobe, &
                   header_38%uc_VelocityT, &
                   header_38%uc_FilderP, &
                   header_38%uc_NoiseT, &
                   header_38%uc_SQIT, &
                   header_38%uc_IntensityC, &
                   header_38%uc_IntensityR, &
                   header_38%uc_CalNoise, &
                   header_38%uc_CalPower, &
                   header_38%uc_CalPulseWidth, &
                   header_38%uc_CalWorkFreq, &
                   header_38%uc_CalLog, &
                   header_38%remain3, &
                   header_38%ul_DataOffset

   close(radar_unit)

   sdate = init_date(get_uc_value(header_38%uc_SYear1 )*100+get_uc_value(header_38%uc_SYear2), &
                     get_uc_value(header_38%uc_SMonth ), &
                     get_uc_value(header_38%uc_SDay   ), &
                     get_uc_value(header_38%uc_SHour  ), &
                     get_uc_value(header_38%uc_SMinute), &
                     get_uc_value(header_38%uc_SSecond))

   edate = init_date(get_uc_value(header_38%uc_EYear1 )*100+get_uc_value(header_38%uc_EYear2), &
                     get_uc_value(header_38%uc_EMonth ), &
                     get_uc_value(header_38%uc_EDay   ), &
                     get_uc_value(header_38%uc_EHour  ), &
                     get_uc_value(header_38%uc_EMinute), &
                     get_uc_value(header_38%uc_ESecond))
   date_bj=edate
   date_utc = get_new_date(date_bj,-8*one_hour)

   btime=get_diff_date(get_new_date(sdate,-8*one_hour),d1970)

   rdat%Year   = date_utc%Year
   rdat%Month  = date_utc%Month
   rdat%Day    = date_utc%Day
   rdat%Hour   = date_utc%Hour
   rdat%Minute = date_utc%Minute
   rdat%Second = date_utc%Second

   rdat%latitude  = get_i_value(header_38%l_LatitudeValue )/3600000.
   rdat%longitude = get_i_value(header_38%l_LongitudeValue)/3600000.
   rdat%altitude  = get_i_value(header_38%l_Height        )/1000.

   rdat%ntilt = get_uc_value(header_38%uc_ScanMode)-100
   if(rdat%ntilt<0)then
      rdat%ntilt=0
      return
   endif

   MaxCuts = rdat%ntilt
   MaxRads = get_us_value(header_38%LayerInfo(1)%us_RecordNumber)
   RGates  = get_us_value(header_38%LayerInfo(1)%us_BinNumber  )

   do k=1, rdat%ntilt
       value=get_us_value(header_38%LayerInfo(k)%us_RecordNumber)
       if(MaxRads<value)then
          MaxRads=value
       endif
       value=get_us_value(header_38%LayerInfo(k)%us_BinNumber  )
       if(RGates<value)then
          RGates=value
       endif
   enddo

   VGates = RGates
   WGates = RGates

   call allocate_radar_data(rdat,RGates,VGates,WGates,MaxRads,MaxCuts)
   
   do k=1, rdat%ntilt

      rdat%nazim  (k) = get_us_value(header_38%LayerInfo(k)%us_RecordNumber)
      rdat%vmax   (k) = get_us_value(header_38%LayerInfo(k)%us_MaxV        )/100.
      rdat%rmax   (k) = get_us_value(header_38%LayerInfo(k)%us_MaxL        )*10.
      rdat%rgatesp(k) = get_us_value(header_38%LayerInfo(k)%us_BindWidth   )*2.
      rdat%vgatesp(k) = get_us_value(header_38%LayerInfo(k)%us_BindWidth   )*2.

      do j=1, rdat%nazim  (k) 

         rdat%nrgate(j,k) = get_us_value(header_38%LayerInfo(k)%us_BinNumber)
         rdat%nvgate(j,k) = get_us_value(header_38%LayerInfo(k)%us_BinNumber)
         rdat%rtilt (j,k) = get_us_value(header_38%LayerInfo(k)%us_Angle    )/100.
         rdat%razim (j,k) = (j-1)*360./rdat%nazim(k)
      enddo
   enddo

   allocate(dat(MaxRads, MaxCuts))
   do k=1, MaxCuts
      do j=1, MaxRads
         allocate(dat(j,k)%z(RGates))
         allocate(dat(j,k)%v(VGates))
         allocate(dat(j,k)%w(WGates))
      enddo
   enddo

   ! read in all data
   TotalRays=SUM(rdat%nazim(1:rdat%ntilt))
   FileSize=1024+(RGates+VGates+WGates)*2*TotalRays
   open(radar_unit,file=filename,recl=FileSize,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr) header, &
                                   (((dat(j,k)%z(i),i=1,RGates),&
                                     (dat(j,k)%v(i),i=1,VGates),&
                                     (dat(j,k)%w(i),i=1,WGates),j=1,rdat%nazim(k)),k=1,rdat%ntilt)
   close(radar_unit)

   dtime=get_diff_date(edate,sdate)/(TotalRays-1)

   write(*,"(2(A,I10))") "Total Rays:", TotalRays,". File Size:", FileSize 

   rdat%vcp=0
   if(rdat%ntilt==14.or.rdat%ntilt==16) rdat%vcp=11
   if(rdat%ntilt== 9.or.rdat%ntilt==11) rdat%vcp=21
   if(rdat%ntilt== 5.or.rdat%ntilt== 7) rdat%vcp=31

   n=0
   do k=1, rdat%ntilt
      rdat%ifref(k)=.true.
      rdat%ifvel(k)=.true.
      do j=1, rdat%nazim(k)
         n=n+1
         rdat%stime(j,k)=(n-1)*dtime+btime
         rdat%etime(j,k)=(n-1)*dtime+btime
         if(rdat%vmax (k)<63.5)then
            rdat%vres (j,k)=0.5
         else
            rdat%vres (j,k)=1.0
         endif
         do i=1, rdat%nrgate(j,k)
           rdat%ref(i,j,k)=Decode38(dat(j,k)%z(i))
         enddo
         do i=1, rdat%nvgate(j,k)
           rdat%vel(i,j,k)=Decode38(dat(j,k)%v(i))
           rdat%spw(i,j,k)=Decode38(dat(j,k)%w(i))
         enddo
      enddo
   enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains
      real function Decode38(s_code)
      implicit none
      character(len=2), intent(in) :: s_code
   
      integer :: idata
      ! 38
      integer, parameter :: NO_DATA = -32768
      
   
      idata=get_s_value(s_code)
      if(idata==NO_DATA)then
         Decode38=value_invalid
      else
         Decode38=idata/10.
      endif
      end function
   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_radar_784(filename,radar_type,rdat)
   use bytes, only: get_uc_value, get_i_value, get_us_value, get_s_value
   use date_pack, only: init_date, get_diff_date
   implicit none

   character(len=*),   intent(in)    :: filename, radar_type
   type(t_radar_data), intent(inout) :: rdat
 

   type t_RADARSITE  ! 170 Bytes
      character(len=30) :: str_country
      character(len=20) :: str_province
      character(len=40) :: str_station
      character(len=10) :: str_stationnumber
      character(len=20) :: str_radartype
      character(len=16) :: str_longitude
      character(len=16) :: str_latitude
      character(len=4) :: l_longitudevalue
      character(len=4) :: l_latitudevalue
      character(len=4) :: l_height
      character(len=2) :: s_Maxangle
      character(len=2) :: s_Opangle
      character(len=2) :: s_MangFreq
   end type

   type t_RADARPERFORMANCEPARAM ! 31 Bytes
      character(len=4) :: l_AntennaG
      character(len=2) :: us_BeamH
      character(len=2) :: us_BeamL
      character(len=1) :: uc_polarizations
      character(len=1) :: c_sidelobe
      character(len=4) :: l_Power
      character(len=4) :: l_wavelength
      character(len=2) :: us_logA
      character(len=2) :: us_LineA
      character(len=2) :: us_AGCP
      character(len=1) :: uc_clutterT
      character(len=1) :: uc_VelocityP
      character(len=1) :: uc_filderP
      character(len=1) :: uc_noiseT
      character(len=1) :: uc_SQIT
      character(len=1) :: uc_intensityC
      character(len=1) :: uc_intensityR
   end type

   type t_LAYERPARAM ! 21 Bytes
      character(len=1) :: uc_ambiguousp
      character(len=2) :: us_Arotate
      character(len=2) :: us_Prf1
      character(len=2) :: us_Prf2
      character(len=2) :: us_spulseW
      character(len=2) :: us_MaxV
      character(len=2) :: us_MaxL
      character(len=2) :: us_binWidth
      character(len=2) :: us_binnumber
      character(len=2) :: us_recordnumber
      character(len=2) :: s_Swangles
   end type

   type t_RADAROBSERVATIONPARAM ! 660 Bytes
      character(len=1) :: uc_stype
      character(len=2) :: us_Syear
      character(len=1) :: uc_Smonth
      character(len=1) :: uc_Sday
      character(len=1) :: uc_Shour
      character(len=1) :: uc_Sminute
      character(len=1) :: uc_Ssecond
      character(len=1) :: uc_Timep
      character(len=4) :: ul_Smillisecond
      character(len=1) :: u_calibration
      character(len=1) :: u_intensityI
      character(len=1) :: u_VelocityP
      type(t_LAYERPARAM), dimension(30) :: LayerParam
      character(len=2) :: us_RHIA
      character(len=2) :: s_RHIL
      character(len=2) :: s_RHIH
      character(len=2) :: us_Eyear
      character(len=1) :: uc_Emonth
      character(len=1) :: uc_Eday
      character(len=1) :: uc_Ehour
      character(len=1) :: uc_Eminute
      character(len=1) :: uc_Esecond
      character(len=1) :: uc_Etenth
   end type

   type t_RADARDATAFILEHEADER ! 1024 Bytes
      type(t_RADARSITE            ) :: RadarSiteInfo
      type(t_RADARPERFORMANCEPARAM) :: RadarPerformanceInfo
      type(t_RADAROBSERVATIONPARAM) :: RadarObservationInfo
      character(len=163)            :: Reserved
   end type

   type t_DATA
     character(len=1) :: uc_dBz
     character(len=1) :: uc_V
     character(len=1) :: uc_dBt
     character(len=1) :: uc_W
   end type

   type t_DATARECORD
      character(len=2) :: us_startaz , us_startel
      character(len=2) :: us_endaz   , us_endel
      type(t_DATA), dimension(1000) :: RawData
   end type

   character(len=1024) :: header

   type(t_RADARDATAFILEHEADER) :: header_784

   type(t_DATARECORD) ,dimension(:,:), allocatable :: dat

   integer      :: len_record
   integer      :: TotalRays, i, j, k, n, ierr, value, FileSize
   type(date)   :: date_bj, date_utc, sdate, edate
   real         :: dtime
   real(kind=8) :: btime
   
   ! read in header only
   open(radar_unit,file=filename,recl=1024,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr)  &
                   header_784%RadarSiteInfo       %str_country       ,&
                   header_784%RadarSiteInfo       %str_province      ,&
                   header_784%RadarSiteInfo       %str_station       ,&
                   header_784%RadarSiteInfo       %str_stationnumber ,&
                   header_784%RadarSiteInfo       %str_radartype     ,&
                   header_784%RadarSiteInfo       %str_longitude     ,&
                   header_784%RadarSiteInfo       %str_latitude      ,&
                   header_784%RadarSiteInfo       %l_longitudevalue  ,&
                   header_784%RadarSiteInfo       %l_latitudevalue   ,&
                   header_784%RadarSiteInfo       %l_height          ,&
                   header_784%RadarSiteInfo       %s_Maxangle        ,&
                   header_784%RadarSiteInfo       %s_Opangle         ,&
                   header_784%RadarSiteInfo       %s_MangFreq        ,&
                   header_784%RadarPerformanceInfo%l_AntennaG        ,&   
                   header_784%RadarPerformanceInfo%us_BeamH          ,& 
                   header_784%RadarPerformanceInfo%us_BeamL          ,& 
                   header_784%RadarPerformanceInfo%uc_polarizations  ,&         
                   header_784%RadarPerformanceInfo%c_sidelobe        ,&   
                   header_784%RadarPerformanceInfo%l_Power           ,&
                   header_784%RadarPerformanceInfo%l_wavelength      ,&     
                   header_784%RadarPerformanceInfo%us_logA           ,&
                   header_784%RadarPerformanceInfo%us_LineA          ,& 
                   header_784%RadarPerformanceInfo%us_AGCP           ,&
                   header_784%RadarPerformanceInfo%uc_clutterT       ,&    
                   header_784%RadarPerformanceInfo%uc_VelocityP      ,&     
                   header_784%RadarPerformanceInfo%uc_filderP        ,&   
                   header_784%RadarPerformanceInfo%uc_noiseT         ,&  
                   header_784%RadarPerformanceInfo%uc_SQIT           ,&
                   header_784%RadarPerformanceInfo%uc_intensityC     ,&      
                   header_784%RadarPerformanceInfo%uc_intensityR     ,&      
                   header_784%RadarObservationInfo%uc_stype          ,& 
                   header_784%RadarObservationInfo%us_Syear          ,& 
                   header_784%RadarObservationInfo%uc_Smonth         ,&  
                   header_784%RadarObservationInfo%uc_Sday           ,&
                   header_784%RadarObservationInfo%uc_Shour          ,& 
                   header_784%RadarObservationInfo%uc_Sminute        ,&   
                   header_784%RadarObservationInfo%uc_Ssecond        ,&   
                   header_784%RadarObservationInfo%uc_Timep          ,& 
                   header_784%RadarObservationInfo%ul_Smillisecond   ,&        
                   header_784%RadarObservationInfo%u_calibration     ,&      
                   header_784%RadarObservationInfo%u_intensityI      ,&     
                   header_784%RadarObservationInfo%u_VelocityP       ,( &    
                   header_784%RadarObservationInfo%LayerParam( i)%uc_ambiguousp                 ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_Arotate                    ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_Prf1                       ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_Prf2                       ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_spulseW                    ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_MaxV                       ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_MaxL                       ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_binWidth                   ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_binnumber                  ,&
                   header_784%RadarObservationInfo%LayerParam( i)%us_recordnumber               ,&
                   header_784%RadarObservationInfo%LayerParam( i)%s_Swangles          , i=1,30) ,&
                   header_784%RadarObservationInfo%us_RHIA           ,& 
                   header_784%RadarObservationInfo%s_RHIL            ,&
                   header_784%RadarObservationInfo%s_RHIH            ,&
                   header_784%RadarObservationInfo%us_Eyear          ,&  
                   header_784%RadarObservationInfo%uc_Emonth         ,&   
                   header_784%RadarObservationInfo%uc_Eday           ,& 
                   header_784%RadarObservationInfo%uc_Ehour          ,&  
                   header_784%RadarObservationInfo%uc_Eminute        ,&    
                   header_784%RadarObservationInfo%uc_Esecond        ,&    
                   header_784%RadarObservationInfo%uc_Etenth         ,&   
                   header_784%Reserved           

   close(radar_unit)

   sdate = init_date(get_us_value(header_784%RadarObservationInfo%us_SYear  ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_SMonth ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_SDay   ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_SHour  ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_SMinute), &
                     get_uc_value(header_784%RadarObservationInfo%uc_SSecond))

   edate = init_date(get_us_value(header_784%RadarObservationInfo%us_EYear  ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_EMonth ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_EDay   ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_EHour  ), &
                     get_uc_value(header_784%RadarObservationInfo%uc_EMinute), &
                     get_uc_value(header_784%RadarObservationInfo%uc_ESecond))
   !date_bj=edate
   !date_utc = get_new_date(date_bj,-8*one_hour)
   date_utc=edate

   !btime=get_diff_date(get_new_date(sdate,-8*one_hour),d1970)
   btime=get_diff_date(sdate,d1970)

   !write(*,*) "Year:",get_us_value(header_784%RadarObservationInfo%us_SYear), get_uc_value(header_784%RadarObservationInfo%us_SYear(1:1)), get_uc_value(header_784%RadarObservationInfo%us_SYear(2:2))
   write(*,*) "WaveLength:",get_i_value(header_784%RadarPerformanceInfo%l_wavelength)
   write(*,*) "BeamWidthV:",get_us_value(header_784%RadarPerformanceInfo%us_BeamH)
   write(*,*) "BeamWidthH:",get_us_value(header_784%RadarPerformanceInfo%us_BeamL)
   write(*,*) "MangFreq  :",get_s_value(header_784%RadarSiteInfo%s_MangFreq)

   rdat%Year   = date_utc%Year
   rdat%Month  = date_utc%Month
   rdat%Day    = date_utc%Day
   rdat%Hour   = date_utc%Hour
   rdat%Minute = date_utc%Minute
   rdat%Second = date_utc%Second

   rdat%latitude  = get_i_value(header_784%RadarSiteInfo%l_LatitudeValue )/100.
   rdat%longitude = get_i_value(header_784%RadarSiteInfo%l_LongitudeValue)/100.
   rdat%altitude  = get_i_value(header_784%RadarSiteInfo%l_Height        )/1000.

   rdat%ntilt = get_uc_value(header_784%RadarObservationInfo%uc_stype)-100
   if(rdat%ntilt<0)then
      rdat%ntilt=0
      return
   endif

   do k=1, rdat%ntilt
      write(*,*) k, get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_Prf1    ), &
                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_Prf2    ), &
                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_spulseW ), &
                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_MaxV    ), &
                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_MaxL    )
 
   enddo
   MaxCuts = rdat%ntilt
   MaxRads = get_us_value(header_784%RadarObservationInfo%LayerParam(1)%us_RecordNumber)
   RGates  = get_us_value(header_784%RadarObservationInfo%LayerParam(1)%us_BinNumber  )

   do k=1, rdat%ntilt
       value=get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_RecordNumber)
       if(MaxRads<value)then
          MaxRads=value
       endif
       value=get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_BinNumber  )
       if(RGates<value)then
          RGates=value
       endif
   enddo

   VGates = RGates
   WGates = RGates

   call allocate_radar_data(rdat,RGates,VGates,WGates,MaxRads,MaxCuts)
   
   do k=1, rdat%ntilt

      rdat%nazim  (k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_RecordNumber)
      rdat%vmax   (k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_MaxV        )/100.
      rdat%rmax   (k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_MaxL        )*10.
      rdat%rgatesp(k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_BinWidth    )/10.
      rdat%vgatesp(k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_BinWidth    )/10.

      do j=1, rdat%nazim  (k) 
         len_record=get_us_value(header_784%RadarObservationInfo%LayerParam(k)%us_BinNumber) ! =DataRecord Length
         rdat%nrgate(j,k) = (len_record-8)/4
         rdat%nvgate(j,k) = (len_record-8)/4
         rdat%rtilt (j,k) = get_us_value(header_784%RadarObservationInfo%LayerParam(k)%s_SWAngles  )/100.
!        rdat%razim (j,k) = (j-1)*360./rdat%nazim(k)
      enddo
   enddo

   allocate(dat(MaxRads, MaxCuts))
!   do k=1, MaxCuts
!      do j=1, MaxRads
!         allocate(dat(j,k)%z(RGates))
!         allocate(dat(j,k)%v(VGates))
!         allocate(dat(j,k)%w(WGates))
!      enddo
!   enddo

   ! read in all data
   TotalRays=SUM(rdat%nazim(1:rdat%ntilt))
   FileSize=1024+len_record*TotalRays
   open(radar_unit,file=filename,recl=FileSize,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr) header, &
                                   (( dat(j,k)%us_startaz, dat(j,k)%us_startel,&
                                      dat(j,k)%us_endaz  , dat(j,k)%us_endel  ,&
                                     (dat(j,k)%RawData(i)%uc_dBz, &
                                      dat(j,k)%RawData(i)%uc_V, &
                                      dat(j,k)%RawData(i)%uc_W, &
                                      dat(j,k)%RawData(i)%uc_dBt,i=1,rdat%nrgate(j,k)),j=1,rdat%nazim(k)),k=1,rdat%ntilt)
   close(radar_unit)

   dtime=get_diff_date(edate,sdate)/(TotalRays-1)

   write(*,"(2(A,I10))") "Total Rays:", TotalRays,". File Size:", FileSize 

!   do k=1, rdat%ntilt
!      write(*,"(F10.4,2(',',L6),3(',',I10),4(',',F12.2))") & 
!                         rdat%rtilt (1,k), &
!                         rdat%ifref   (k), &
!                         rdat%ifvel   (k), &
!                         rdat%nazim   (k), &
!                         rdat%nrgate(1,k), &
!                         rdat%nvgate(1,k), &
!                         rdat%rgatesp (k), &
!                         rdat%vgatesp (k), &
!                         rdat%rmax    (k), &
!                         rdat%vmax    (k)
!   enddo
   rdat%vcp=0
   if(rdat%ntilt==14.or.rdat%ntilt==16) rdat%vcp=11
   if(rdat%ntilt== 9.or.rdat%ntilt==11) rdat%vcp=21
   if(rdat%ntilt== 5.or.rdat%ntilt== 7) rdat%vcp=31

   n=0
   do k=1, rdat%ntilt
      rdat%ifref(k)=.true.
      rdat%ifvel(k)=.true.
      do j=1, rdat%nazim(k)
         n=n+1
         rdat%razim (j,k) = get_us_value(dat(j,k)%us_startaz)*360./65536. 
         rdat%stime(j,k)=(n-1)*dtime+btime
         rdat%etime(j,k)=(n-1)*dtime+btime
         rdat%vres (j,k)=rdat%vmax(k)/127.
         do i=1, rdat%nrgate(j,k)
           if(get_uc_value(dat(j,k)%RawData(i)%uc_dBz)/=0) rdat%ref(i,j,k)=(get_uc_value(dat(j,k)%RawData(i)%uc_dBz)-64.)/2
         enddo
         do i=1, rdat%nvgate(j,k)
           if(get_uc_value(dat(j,k)%RawData(i)%uc_V)/=0) rdat%vel(i,j,k)=rdat%vmax(k)*(get_uc_value(dat(j,k)%RawData(i)%uc_V)-128.)/127.
           if(get_uc_value(dat(j,k)%RawData(i)%uc_W)/=0) rdat%spw(i,j,k)=rdat%vmax(k)*get_uc_value(dat(j,k)%RawData(i)%uc_W)/512.
         enddo
      enddo
   enddo

   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_radar_98d(filename,radar_type,rdat)
   use bytes, only: get_s_value, get_f_value, &
                    get_uc_value, get_us_value, get_ui_value
   use date_pack, only: init_date, get_new_date, one_day, one_hour
   implicit none

   character(len=*),   intent(in)    :: filename
   character(len=*),   intent(inout) :: radar_type
   type(t_radar_data), intent(inout) :: rdat
 
   type t_radar_98d
     character(len=14)                            :: us_temp1
     character(len=2 )                            :: us_RadarStatus                   ! Radar Status: 1 - Radar Data
     character(len=12)                            :: us_temp2
     character(len=4 )                            :: ui_mSeconds                      ! MiliSeconds Of Radial Data Time
     character(len=2 )                            :: us_JulianDate                    ! Days From 1970/1/1
     character(len=2 )                            :: us_URange                        ! Rmax
     character(len=2 )                            :: us_Az                            ! Azimuth Angle
     character(len=2 )                            :: us_RadialNumber                  ! Radial Number
     character(len=2 )                            :: us_RadialStatus                  ! Radial Status
     character(len=2 )                            :: us_El                            ! Elevation Angle
     character(len=2 )                            :: us_ElNumber                      ! Elevation Number
     character(len=2 )                            ::  s_RangeToFirstGateOfRef         ! Range To First Gate Of Reflectivity unit:m
     character(len=2 )                            ::  s_RangeToFirstGateOfDop         ! Range To First Gate Of Doppler      unit:m
     character(len=2 )                            :: us_GateSizeOfReflectivity        ! Gate Size Of Reflectivity unit:m
     character(len=2 )                            :: us_GateSizeOfDoppler             ! Gate Size Of Doppler      unit:m
     character(len=2 )                            :: us_GatesNumberOfReflectivity     ! Gates Number Of Reflectivity
     character(len=2 )                            :: us_GatesNumberOfDoppler          ! Gates Number Of Doppler 
     character(len=2 )                            :: us_CutSectorNumber               ! Cut Sector Number
     character(len=4 )                            :: ui_CalibrationConst              ! Calibration Constant
     character(len=2 )                            :: us_PtrOfReflectivity             ! Pointer Of Reflectivity
     character(len=2 )                            :: us_PtrOfVelocity                 ! Pointer Of Velocity
     character(len=2 )                            :: us_PtrOfSpectrumWidth            ! Pointer Of SpectrumWidth
     character(len=2 )                            :: us_ResolutionOfVelocity          ! Resolution Of Velocity
     character(len=2 )                            :: us_VcpNumber                     ! Vcp Number
     character(len=8 )                            :: us_temp4
     character(len=2 )                            :: us_PtrOfArcReflectivity          ! Pointer Of Arc Reflectivity
     character(len=2 )                            :: us_PtrOfArcVelocity              ! Pointer Of Arc Velocity  
     character(len=2 )                            :: us_PtrOfArcWidth                 ! Pointer Of Arc Width 
     character(len=2 )                            :: us_Nyquist                       ! Nyquist
     character(len=2 )                            :: us_AAF                           ! Atmospheric attenuation factor (0.001dB/Km)
     character(len=2 )                            :: us_temp47
     character(len=2 )                            :: us_temp48
     character(len=2 )                            :: us_CircleTotal                   ! Total Circle Number
     character(len=30)                            :: uc_temp5
     character(len=1 ), allocatable, dimension(:) :: uc_Echodata                      ! [RGates+VGates+WGates]
     character(len=4 )                            :: uc_temp
   end type ! 128+RGates+VGates+WGates+4 Bytes
   ! type t_radar_784 == t_radar_98d

   type(t_radar_98d) :: data_98d

   type(date) :: date_now, date_bj, date_utc
   integer    :: irec, ierr, k
   integer    :: RadialStatus, ElIndex, AzIndex, FstBin, BinNum, LstBin, BnIndex, ptrPos, TotalRays, FileSize
   real       :: CurEl, CurAz, Vmax
   logical    :: VolBeg, RFlag, VFlag, WFlag, VolEnd, NewElv, if_checked
   integer(kind=8) :: btime, ctime
   
   integer :: i

   ! check radar type

   if_checked=.false.
   inquire(file=filename,size=filesize)
   if(MOD(filesize,132+RGates_SA+VGates_SA+WGates_SA)==0)then ! 2432
      if(radar_type/="SA".and.radar_type/="SB")then
         radar_type="SA"
         if_checked=.true.
      endif 
   elseif(MOD(filesize,132+RGates_CB+VGates_CB+WGates_CB)==0)then ! 4132
      if(radar_type/="CA".and.radar_type/="CB")then
         radar_type="CB"
         if_checked=.true.
      endif
   elseif(MOD(filesize,132+RGates_SC+VGates_SC+WGates_SC)==0)then ! 3132
      if(radar_type/="SC".and.radar_type/="CD")then
         radar_type="SC"
         if_checked=.true.
      endif
   elseif(MOD(filesize,132+RGates_CC+VGates_CC+WGates_CC)==0)then ! 1632
      if(radar_type/="CC")then
         radar_type="CC"
         if_checked=.true.
      endif
   else
      write(*,*) "Can't tell radar type by filesize:", filesize
      radar_type="unknown"
      return
   endif
   if(if_checked)then
      write(*,*) "Checked radar type:", radar_type
   endif

   MaxRads = 512
   select case(trim(radar_type))
   case ("SA","SB")
      RGates = RGates_SA 
      VGates = VGates_SA 
      WGates = WGates_SA 
      MaxCuts=MaxCuts_SA
   case ("CA","CB")
      RGates = RGates_CB 
      VGates = VGates_CB 
      WGates = WGates_CB 
      MaxCuts=MaxCuts_CB
   case ("SC","CD")
      RGates = RGates_SC 
      VGates = VGates_SC 
      WGates = WGates_SC 
      MaxCuts=MaxCuts_SC
   case ("CC")
      RGates = RGates_CC 
      VGates = VGates_CC 
      WGates = WGates_CC 
      MaxCuts=MaxCuts_CC
   case default
      write(*,*) "in read_radar_98d:Unknown Radar Type:", trim(radar_type)
      stop
   end select

   if(allocated (data_98d%uc_Echodata))then
      deallocate(data_98d%uc_Echodata)
   endif
   allocate(data_98d%uc_Echodata(RGates+VGates+WGates))

   open(radar_unit,file=filename,recl=128+RGates+VGates+WGates+4,status="old",access="direct")

   call allocate_radar_data(rdat,RGates,VGates,WGates,MaxRads,MaxCuts)
   
   rdat%nazim  = MaxRads
   rdat%nrgate = 0
   rdat%nvgate = 0

   VolBeg = .false.
   VolEnd = .false.

   irec=0
   do while(.true.)
      irec=irec+1
      ! Initialize flags
      RFlag = .false.
      VFlag = .false.
      WFlag = .false.

      read(radar_unit,rec=irec,iostat=ierr) data_98d%us_temp1, &
                                            data_98d%us_RadarStatus, &
                                            data_98d%us_temp2, &
                                            data_98d%ui_mSeconds, &
                                            data_98d%us_JulianDate, &
                                            data_98d%us_URange, &
                                            data_98d%us_Az, &
                                            data_98d%us_RadialNumber, &
                                            data_98d%us_RadialStatus, &
                                            data_98d%us_El, &
                                            data_98d%us_ElNumber, &
                                            data_98d% s_RangeToFirstGateOfRef, &
                                            data_98d% s_RangeToFirstGateOfDop, &
                                            data_98d%us_GateSizeOfReflectivity, &
                                            data_98d%us_GateSizeOfDoppler, &
                                            data_98d%us_GatesNumberOfReflectivity, &
                                            data_98d%us_GatesNumberOfDoppler, &
                                            data_98d%us_CutSectorNumber, &
                                            data_98d%ui_CalibrationConst, &
                                            data_98d%us_PtrOfReflectivity, &
                                            data_98d%us_PtrOfVelocity, &
                                            data_98d%us_PtrOfSpectrumWidth, &
                                            data_98d%us_ResolutionOfVelocity, &
                                            data_98d%us_VcpNumber, &
                                            data_98d%us_temp4, &
                                            data_98d%us_PtrOfArcReflectivity, &
                                            data_98d%us_PtrOfArcVelocity, &
                                            data_98d%us_PtrOfArcWidth, &
                                            data_98d%us_Nyquist, &
                                            data_98d%us_AAF, &
                                            data_98d%us_temp47, &
                                            data_98d%us_temp48, &
                                            data_98d%us_CircleTotal, &
                                            data_98d%uc_temp5, &
                                            data_98d%uc_Echodata, &
                                            data_98d%uc_temp

      if(ierr/=0) exit
      ! Check data
      if (get_us_value(data_98d%us_RadialStatus)==0.and.get_us_value(data_98d%us_El)==0.and.get_us_value(data_98d%us_JulianDate)==0)then
         cycle
      endif
      RadialStatus=get_us_value(data_98d%us_RadialStatus)

      !write(*,*) 'data_98d%ui_CalibrationConst"',get_f_value(data_98d%ui_CalibrationConst),',data_98d%us_AAF:"',get_s_value(data_98d%us_AAF)
      ! Start a volume scan
      write(601,*) irec, "RadialStatus:", RadialStatus, "VOL_BEG=3, VOL_END=4, ELV_BEG=0, ELV_END=2"

      
      if(RadialStatus == VOL_BEG)then

         write(602,*) "us_RadarStatus              : ",  get_us_value(data_98d%us_RadarStatus)
         write(602,*) "ui_mSeconds                 : ",  get_ui_value(data_98d%ui_mSeconds)
         write(602,*) "us_JulianDate               : ",  get_us_value(data_98d%us_JulianDate)
         write(602,*) "us_URange                   : ",  get_us_value(data_98d%us_URange)
         write(602,*) "us_Az                       : ",  get_us_value(data_98d%us_Az)
         write(602,*) "us_RadialNumber             : ",  get_us_value(data_98d%us_RadialNumber)
         write(602,*) "us_RadialStatus             : ",  get_us_value(data_98d%us_RadialStatus)
         write(602,*) "us_El                       : ",  get_us_value(data_98d%us_El)
         write(602,*) "us_ElNumber                 : ",  get_us_value(data_98d%us_ElNumber)
         write(602,*) " s_RangeToFirstGateOfRef    : ",  get_s_value (data_98d% s_RangeToFirstGateOfRef)
         write(602,*) " s_RangeToFirstGateOfDop    : ",  get_s_value (data_98d% s_RangeToFirstGateOfDop)
         write(602,*) "us_GateSizeOfReflectivity   : ",  get_us_value(data_98d%us_GateSizeOfReflectivity)
         write(602,*) "us_GateSizeOfDoppler        : ",  get_us_value(data_98d%us_GateSizeOfDoppler)
         write(602,*) "us_GatesNumberOfReflectivity: ",  get_us_value(data_98d%us_GatesNumberOfReflectivity)
         write(602,*) "us_GatesNumberOfDoppler     : ",  get_us_value(data_98d%us_GatesNumberOfDoppler)
         write(602,*) "us_CutSectorNumber          : ",  get_us_value(data_98d%us_CutSectorNumber)
         write(602,*) "ui_CalibrationConst         : ",  get_us_value(data_98d%ui_CalibrationConst)
         write(602,*) "us_PtrOfReflectivity        : ",  get_us_value(data_98d%us_PtrOfReflectivity)
         write(602,*) "us_PtrOfVelocity            : ",  get_us_value(data_98d%us_PtrOfVelocity)
         write(602,*) "us_PtrOfSpectrumWidth       : ",  get_us_value(data_98d%us_PtrOfSpectrumWidth)
         write(602,*) "us_ResolutionOfVelocity     : ",  get_us_value(data_98d%us_ResolutionOfVelocity)
         write(602,*) "us_VcpNumber                : ",  get_us_value(data_98d%us_VcpNumber)
         write(602,*) "us_PtrOfArcReflectivity     : ",  get_us_value(data_98d%us_PtrOfArcReflectivity)
         write(602,*) "us_PtrOfArcVelocity         : ",  get_us_value(data_98d%us_PtrOfArcVelocity)
         write(602,*) "us_PtrOfArcWidth)           : ",  get_us_value(data_98d%us_PtrOfArcWidth)
         write(602,*) "us_Nyquist                  : ",  get_us_value(data_98d%us_Nyquist)
         write(602,*) "us_AAF                      : ",  get_us_value(data_98d%us_AAF)
         write(602,*) "us_CircleTotal              : ",  get_us_value(data_98d%us_CircleTotal)
         write(602,*) "us_temp1                    : ",  trim(data_98d%us_temp1 ),',', (get_us_value(data_98d%us_temp1 (2*i-1:2*i)), i=1, 7)
         write(602,*) "us_temp2                    : ",  trim(data_98d%us_temp2 ),',', (get_us_value(data_98d%us_temp2 (2*i-1:2*i)), i=1, 6)
         write(602,*) "us_temp4                    : ",  trim(data_98d%us_temp4 ),',', (get_us_value(data_98d%us_temp4 (2*i-1:2*i)), i=1, 4)
         write(602,*) "us_temp47                   : ",  trim(data_98d%us_temp47),',', (get_us_value(data_98d%us_temp47(2*i-1:2*i)), i=1, 1)
         write(602,*) "us_temp48                   : ",  trim(data_98d%us_temp48),',', (get_us_value(data_98d%us_temp48(2*i-1:2*i)), i=1, 1)
         write(602,*) "uc_temp5                    : ",  trim(data_98d%uc_temp5 ),',', (get_us_value(data_98d%uc_temp5 (2*i-1:2*i)), i=1, 15)
         write(602,*) "uc_temp                     : ",  trim(data_98d%uc_temp  ),',', (get_us_value(data_98d%uc_temp  (2*i-1:2*i)), i=1, 2)

         VolBeg = .true.
         NewElv = .true.

         rdat%calibConst       = get_f_value(data_98d%ui_CalibrationConst)

         rdat%vcp = get_us_value(data_98d%us_VcpNumber)
         date_now       = get_new_date(d1970,   (get_us_value(data_98d%us_JulianDate)-1)*one_day)
         !write(601,*) "1",date_now, (get_us_value(data_98d%us_JulianDate)-1)*one_day
         date_now       = get_new_date(date_now, get_ui_value(data_98d%ui_mSeconds)/1000.)
         !write(601,*) "2",date_now, get_ui_value(data_98d%ui_mSeconds)/1000., get_ui_value(data_98d%ui_mSeconds)
         
         btime=(get_us_value(data_98d%us_JulianDate)-1)*DBLE(one_day)+get_ui_value(data_98d%ui_mSeconds)/1000.

         if(radar_type == "SA" .or. radar_type == "SB" .or. radar_type == "CA".or. radar_type == "CB")then
            date_utc=date_now
         else
            date_bj  = date_now
            date_utc = get_new_date(date_bj,-8*one_hour)
         endif
         rdat%Year   = date_utc%Year
         rdat%Month  = date_utc%Month
         rdat%Day    = date_utc%Day
         rdat%Hour   = date_utc%Hour
         rdat%Minute = date_utc%Minute
         rdat%Second = date_utc%Second

         write(601,*) "Start_Vol:", btime,  rdat%Year, rdat%Month, rdat%Day, rdat%Hour, rdat%Minute, rdat%Second
      endif
      !write(601,*) "Status:",RadialStatus, (get_us_value(data_98d%us_JulianDate)-1)*int(one_day)+get_ui_value(data_98d%ui_mSeconds)/1000,( (get_us_value(data_98d%us_JulianDate)-1)*one_day+get_ui_value(data_98d%ui_mSeconds)/1000.)-(8*one_hour)
      if(.NOT.VolBeg) cycle

      ! Start an elevation
      if(RadialStatus == ELV_BEG)then
         NewElv=.true.
      endif

      if(NewElv)then
         ! Elevation Angle
         CurEl   = (get_us_value(data_98d%us_El)/8.)*(180./4096.)
         ElIndex = get_us_value(data_98d%us_ElNumber)

         rdat%atmosAttenFactor(ElIndex) = get_s_value(data_98d%us_AAF)/1000.

         rdat%vmax   (ElIndex) = get_us_value(data_98d%us_Nyquist)/100.
         rdat%rmax   (ElIndex) = get_us_value(data_98d%us_URange )*100. 
         if(radar_type=="SC".or.radar_type=="CD") &
            rdat%rmax(ElIndex) =rdat%rmax(ElIndex) *10.
         rdat%rgatesp(ElIndex) = get_us_value(data_98d%us_GateSizeOfReflectivity)
         rdat%vgatesp(ElIndex) = get_us_value(data_98d%us_GateSizeOfDoppler)

         NewElv=.false.
         write(601,*) "New Elv:", ElIndex, CurEl,get_s_value(data_98d%us_AAF),rdat%atmosAttenFactor(ElIndex)
      endif

      ! Calculate azimuth angle and Azimuth Index
      CurAz = (get_us_value(data_98d%us_Az)/8.)*(180./4096.)
      if(CurAz >= 360.) CurAz = CurAz-360.
      AzIndex = get_us_value(data_98d%us_RadialNumber) 
      write(601,*) "Az Index:", AzIndex, CurAz

      if(RadialStatus == VOL_END)then
         VolEnd=.true.
         rdat%ntilt=ElIndex
         rdat%nazim(ElIndex) = AzIndex
!        write(*,*) "NAzim(",ElIndex,")=",AzIndex
      endif

      rdat%rtilt(AzIndex, ElIndex) = CurEl

      ctime=(get_us_value(data_98d%us_JulianDate)-1)*DBLE(one_day)+get_ui_value(data_98d%ui_mSeconds)/1000.
      !write(*,*) AzIndex, ElIndex,ctime, btime
      if(radar_type == "SA" .or. radar_type == "SB" .or. radar_type == "CA" .or. radar_type == "CB")then
         rdat%stime (AzIndex, ElIndex) = ctime
         rdat%etime (AzIndex, ElIndex) = ctime
      else
         rdat%stime (AzIndex, ElIndex) = ctime -DBLE(8*one_hour)
         rdat%etime (AzIndex, ElIndex) = ctime -DBLE(8*one_hour)
      endif

      write(601,*) AzIndex, ElIndex, ctime, rdat%stime (AzIndex, ElIndex)
      if(RadialStatus == ELV_END)then
         rdat%nazim(ElIndex) = AzIndex
!        write(*,*) "NAzim(",ElIndex,")=",AzIndex
      endif

      ! what kind of data in this cut
      if(get_us_value(data_98d%us_PtrOfReflectivity ) /=0) RFlag=.true.
      if(get_us_value(data_98d%us_PtrOfVelocity     ) /=0) VFlag=.true.
      if(get_us_value(data_98d%us_PtrOfSpectrumWidth) /=0) WFlag=.true.

      if(RFlag)then
         rdat%nrgate(AzIndex,ElIndex)=get_us_value(data_98d%us_GatesNumberOfReflectivity)
         rdat%ifref(ElIndex)=.true.
      endif
      if(VFlag)then
         rdat%nvgate(AzIndex,ElIndex)=get_us_value(data_98d%us_GatesNumberOfDoppler     )
         rdat%ifvel(ElIndex)=.true.
      endif
      rdat%razim (AzIndex,ElIndex)=CurAz
      rdat%vres  (AzIndex,ElIndex)=get_us_value(data_98d%us_ResolutionOfVelocity)/4.

      !Save reflectivity data into the array
      if(RFlag)then
         !Get first bin, last bin, and number of bins
         if(radar_type=="SA".or.radar_type=="SB".or. radar_type == "CA".or.radar_type=="CB")then
            FstBin = int(get_s_value(data_98d% s_RangeToFirstGateOfRef)/get_us_value(data_98d%us_GateSizeOfReflectivity)+0.5)
         else
            FstBin=0
         endif
         BinNum = get_us_value(data_98d%us_GatesNumberOfReflectivity)
         if(FstBin<0)then
            BinNum = FstBin+BinNum
            FstBin = -1*FstBin
         elseif(FstBin>0)then
            FstBin = 0 
         endif
         LstBin = FstBin + BinNum
         ptrPos = get_us_value(data_98d%us_PtrOfReflectivity)
    
         !write(*,*) "Save ref:",FstBin,LstBin,ptrPos,BinNum
         !Save data
         do BnIndex=FstBin+1, LstBin
            !write(*,*) BnIndex,AzIndex,ElIndex,ubound(rdat%ref),size(data_98d%uc_Echodata),ptrPos-100+BnIndex
            rdat%ref(BnIndex,AzIndex,ElIndex) = DecodeRef(data_98d%uc_Echodata(ptrPos-100+BnIndex));

            if(rdat%ref(BnIndex,AzIndex,ElIndex)==VALUE_RANFOLD)then
               rdat%ref(BnIndex,AzIndex,ElIndex)= VALUE_INVALID
            endif
         enddo
         !write(*,*) "ref:",(BnIndex,":",rdat%ref(BnIndex,AzIndex,ElIndex),",",BnIndex=FstBin+1, LstBin)
      endif
    
      !Save velocity data into the array
      if(VFlag)then
         !Get first bin, last bin, and number of bins
         if(radar_type=="SA".or.radar_type=="SB".or. radar_type == "CA".or.radar_type=="CB")then
            FstBin = int(get_s_value(data_98d% s_RangeToFirstGateOfDop)/get_us_value(data_98d%us_GateSizeOfDoppler)+0.5)
         else
            FstBin=0
         endif
         BinNum = get_us_value(data_98d%us_GatesNumberOfDoppler)
         if(FstBin<0)then
            BinNum = FstBin+BinNum
            FstBin = -1*FstBin
         elseif(FstBin>0)then
            FstBin = 0 
         endif
         LstBin = FstBin + BinNum
         ptrPos = get_us_value(data_98d%us_PtrOfVelocity)

         !write(*,*) "Save vel:",FstBin,LstBin,ptrPos,BinNum
         !Save data
         do BnIndex=FstBin+1, LstBin
            rdat%vel(BnIndex,AzIndex,ElIndex)=& 
               DecodeVel(data_98d%uc_Echodata(ptrPos-100+BnIndex),data_98d%us_ResolutionOfVelocity)

            if(rdat%vel(BnIndex,AzIndex,ElIndex)==VALUE_RANFOLD)then
               rdat%vel(BnIndex,AzIndex,ElIndex)= VALUE_INVALID
            endif
         enddo
         !write(*,*) "vel:",(BnIndex,":",rdat%vel(BnIndex,AzIndex,ElIndex),",",BnIndex=FstBin+1, LstBin)
      endif
    
      !Save spectrum width data into the array
      if(WFlag)then
         !Get first bin, last bin, and number of bins
         if(radar_type=="SA".or.radar_type=="SB".or. radar_type == "CA".or.radar_type=="CB")then
            FstBin = int(get_s_value(data_98d% s_RangeToFirstGateOfDop)/get_us_value(data_98d%us_GateSizeOfDoppler)+0.5)
         else
            FstBin=0
         endif
         BinNum = get_us_value(data_98d%us_GatesNumberOfDoppler) 
         if(FstBin<0)then
            BinNum = FstBin+BinNum
            FstBin = -1*FstBin
         elseif(FstBin>0)then
            FstBin = 0 
         endif
         LstBin = FstBin + BinNum
         ptrPos = get_us_value(data_98d%us_PtrOfSpectrumWidth)

         !write(*,*) "Save spw:",FstBin,LstBin,ptrPos,BinNum
         !Save data
         do BnIndex=FstBin+1, LstBin
            rdat%spw(BnIndex,AzIndex,ElIndex)=DecodeSpw(data_98d%uc_Echodata(ptrPos-100+BnIndex))

            if(    rdat%spw(BnIndex,AzIndex,ElIndex)==VALUE_RANFOLD)then
               !.or. ABS(rdat%spw(BnIndex,AzIndex,ElIndex))>rdat%vmax(ElIndex)
               !rdat%vel(BnIndex,AzIndex,ElIndex)= VALUE_INVALID
               rdat%spw(BnIndex,AzIndex,ElIndex)= VALUE_INVALID
            endif
         enddo
         !write(*,*) "spw:",(BnIndex,":",rdat%spw(BnIndex,AzIndex,ElIndex),",",BnIndex=FstBin+1, LstBin)
      endif

   enddo
   close(radar_unit)

   if(.NOT.VolEnd)then
      write(*,*) RadialStatus, ElIndex, AzIndex
      VolEnd=.true.
      rdat%ntilt=ElIndex
      rdat%nazim(ElIndex) = AzIndex
      write(*,*) "NAzim(",ElIndex,")=",AzIndex
   endif

   TotalRays = SUM(rdat%nazim(1:rdat%ntilt))
   FileSize  = (128+(RGates+VGates+WGates)+4)*TotalRays

   write(*,"(2(A,I10))") "Total Rays:", TotalRays,". File Size:", FileSize 

   contains
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    integer(kind=1) function EncodeRef(Ref)
  !    implicit none
  ! 
  !    real, intent(in) :: Ref
  ! 
  !    integer :: code
  ! 
  !    if(Ref==VALUE_RANFOLD)then
  !       code=CODE_RANFOLD
  !    elseif(Ref==VALUE_INVALID)then
  !       code=CODE_INVALID
  !    else
  !       code=(Ref+32.5)*2+2
  !    endif
  ! 
  !    EncodeRef=code
  ! 
  !    end function
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function DecodeRef(uc_code)
      implicit none
   
      character(len=1), intent(in) :: uc_code
   
      integer :: code
   
      code=get_uc_value(uc_code)
   
      if(code==CODE_INVALID.or.code==2)then
         DecodeRef=VALUE_INVALID
      else if(code==CODE_RANFOLD)then
         DecodeRef=VALUE_RANFOLD
      else
         DecodeRef=((code-2.)/2.-32.5)
      endif
      end function
   
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    integer(kind=1) function EncodeVel(Vel, Res)
  !    implicit none
  ! 
  !    real, intent(in) :: Vel
  !    real, intent(in) :: Res
  ! 
  !    integer :: code
  ! 
  !    if(Vel==VALUE_RANFOLD)then
  !       code=CODE_RANFOLD
  !    elseif(Vel==VALUE_INVALID)then
  !       code=CODE_INVALID
  !    else
  ! 
  !      if(Res==0.5)then  ! res:0.5 m/s
  !         code=(Vel+63.5)*2+2
  !      else              ! res:1.0 m/s
  !         code=(Vel+127)+2
  !      endif
  !    endif
  !    EncodeVel=code
  !    end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    integer(kind=1) function EncodeScaleOffset(Var, V_scale, V_offset)
  !    implicit none
  ! 
  !    real, intent(in) :: Var
  !    real, intent(in) :: V_scale, V_offset
  ! 
  !    integer :: code
  ! 
  !    if(Var==VALUE_RANFOLD)then
  !       code=CODE_RANFOLD
  !    elseif(Var==VALUE_INVALID)then
  !       code=CODE_INVALID
  !    else
  !       code=Var*V_scale+V_offset
  !       if(code==CODE_RANFOLD.or.code==CODE_INVALID) code=2
  !    endif
  !    EncodeScaleOffset=code
  !    end function
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function DecodeVel(uc_code, s_ResType)
      implicit none
   
      character(len=1), intent(in) :: uc_code
      character(len=2), intent(in) :: s_ResType
   
      integer :: code, ResType
   
      code   =get_uc_value(uc_code)
      ResType=get_us_value(s_ResType)
   
      if(code==CODE_INVALID)then
         DecodeVel=VALUE_INVALID
      else if(code==CODE_RANFOLD)then
         DecodeVel=VALUE_RANFOLD
      else
        if(ResType==RES_POINT_FIVE)then  ! res:0.5 m/s
           DecodeVel=((code-2.)/2.-63.5)
        else                             ! res:1.0 m/s
           DecodeVel=((code-2)-127.)
        endif
      endif
      end function
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   integer(kind=1) function EncodeSpw(Spw)
   !   implicit none
   !   real, intent(in) :: Spw
   !
   !   integer :: code
   !
   !   if(Spw==VALUE_RANFOLD)then
   !      code=CODE_RANFOLD
   !   elseif(Spw==VALUE_INVALID)then
   !      code=CODE_INVALID
   !   else
   !      code=(Spw+63.5)*2+2
   !   endif
   !
   !   EncodeSpw=code
   !   end function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function DecodeSpw(uc_code)
      implicit none
   
      character(len=1), intent(in) :: uc_code
   
      integer :: code
   
      code=get_uc_value(uc_code)
   
      if(code==CODE_INVALID)then
         DecodeSpw=VALUE_INVALID
      else if(code==CODE_RANFOLD)then
         DecodeSpw=VALUE_RANFOLD
      else
         DecodeSpw=((code-2.)/2.-63.5)
      endif
      end function

   end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_radar_rd(filename,radar_type,rdat)
   use bytes, only: get_i_value, get_c_value, get_s_value, &
                    get_us_value, get_uc_value
   use date_pack, only: init_date, get_new_date, get_diff_date, one_hour
   implicit none

   character(len=*),   intent(in)    :: filename, radar_type
   type(t_radar_data), intent(inout) :: rdat
 
! CC RD

   type t_LAYERPARAM ! 35 bytes  
      character(len=1) :: uc_DataType
      character(len=1) :: uc_Ambiguousp
      character(len=2) :: us_Arotate
      character(len=2) :: us_PRF1
      character(len=2) :: us_PRF2
      character(len=2) :: us_PulseW
      character(len=2) :: us_MaxV
      character(len=2) :: us_MaxL
      character(len=2) :: us_ZBinWidth
      character(len=2) :: us_VBinWidth
      character(len=2) :: us_WBinWidth
      character(len=2) :: us_ZBinNumber
      character(len=2) :: us_VBinNumber
      character(len=2) :: us_WBinNumber
      character(len=2) :: us_RecordNumber
      character(len=2) :: s_SwpAngles
      character(len=1) :: c_DataForm
      character(len=4) :: ul_DBegin
   end type
   
   type t_BINPARAM ! 8 bytes  
      character(len=2) :: s_Code
      character(len=2) :: s_Begin
      character(len=2) :: s_End
      character(len=2) :: s_BinLength
   end type

   type t_rd_header
      ! RADARDATAHEAD 12 bytes
      character(len=4  ) :: c_FileID
      character(len=4  ) :: f_VersionNo
      character(len=4  ) :: l_FileHeaderLength
       ! RADARSITE 168 bytes
      character(len=30 ) :: c_Country
      character(len=20 ) :: c_Province
      character(len=40 ) :: c_Station
      character(len=10 ) :: c_StationNumber
      character(len=20 ) :: c_RadarType
      character(len=16 ) :: c_Longitude
      character(len=16 ) :: c_Latitude
      character(len=4  ) :: l_LongitudeValue
      character(len=4  ) :: l_LatitudeValue
      character(len=4  ) :: l_Height
      character(len=2  ) :: s_MaxAngle
      character(len=2  ) :: s_OptiAngle
      ! RADARPERFORMANCEPARAM 36 bytes
      character(len=4  ) :: l_AntennaG
      character(len=2  ) :: us_VerBeamW
      character(len=2  ) :: us_HorBeamW
      character(len=1  ) :: uc_Polarizations
      character(len=2  ) :: us_SideLobe
      character(len=4  ) :: l_Power
      character(len=4  ) :: l_WaveLength
      character(len=2  ) :: us_LogA
      character(len=2  ) :: us_LineA
      character(len=2  ) :: us_AGCP   
      character(len=2  ) :: us_LogMinPower
      character(len=2  ) :: us_LineMinPower
      character(len=1  ) :: uc_ClutterT
      character(len=1  ) :: uc_VelocityP
      character(len=1  ) :: uc_FilterP
      character(len=1  ) :: uc_NoiseT
      character(len=1  ) :: uc_SQIT
      character(len=1  ) :: uc_IntensityC
      character(len=1  ) :: uc_IntersityR
      ! RADAROBSERVATIONPARAM 1282 bytes
      character(len=1  ) :: uc_SType
      character(len=2  ) :: us_SYear
      character(len=1  ) :: uc_SMonth
      character(len=1  ) :: uc_SDay
      character(len=1  ) :: uc_SHour
      character(len=1  ) :: uc_SMinute
      character(len=1  ) :: uc_SSecond
      character(len=1  ) :: uc_TimP
      character(len=4  ) :: ul_SMillisecond
      character(len=1  ) :: uc_Calibration
      character(len=1  ) :: uc_IntensityI
      character(len=1  ) :: uc_VelocityPn
      character(len=2  ) :: us_ZStartBin
      character(len=2  ) :: us_VStartBin
      character(len=2  ) :: us_WStartBin
      type(t_LAYERPARAM) :: LayerInfo(32)
      character(len=2  ) :: us_RHIA
      character(len=2  ) :: s_RHIL
      character(len=2  ) :: s_RHIH
      character(len=2  ) :: us_EYear
      character(len=1  ) :: uc_EMonth
      character(len=1  ) :: uc_EDay
      character(len=1  ) :: uc_EHour
      character(len=1  ) :: uc_EMinute
      character(len=1  ) :: uc_ESecond
      character(len=1  ) :: uc_ETenth
      character(len=2  ) :: us_ZBinByte
      type(t_BINPARAM  ) :: BinRange1(5)
      character(len=2  ) :: us_VBinByte
      type(t_BINPARAM  ) :: BinRange2(5)
      character(len=2  ) :: us_WBinByte
      type(t_BINPARAM  ) :: BinRange3(5)
      ! RADAROTHERINFORMATION 562 bytes
      character(len=2  ) :: c_StationID
      character(len=1  ) :: c_JHType
      character(len=559) :: c_Spare
   end type

   type(t_rd_header) :: rd_header
   
   type t_data_header ! 11 bytes  
      character(len=2) :: s_Elev
      character(len=2) :: us_Az
      character(len=1) :: uc_Hh ! Hour
      character(len=1) :: uc_Mm ! Minute
      character(len=1) :: uc_Ss ! Second
      character(len=4) :: ul_Ms ! MiliSecond
   end type

   character(len=1), dimension(:,:), allocatable :: uc_CorZ
   character(len=1), dimension(:,:), allocatable :: uc_UnZ 
   character(len=1), dimension(:,:), allocatable ::  c_V   
   character(len=1), dimension(:,:), allocatable :: uc_W   
   character(len=2), dimension(:,:), allocatable ::  s_ZDR 
   character(len=2), dimension(:,:), allocatable ::  s_PHDP
   character(len=2), dimension(:,:), allocatable ::  s_KDP 
   character(len=2), dimension(:,:), allocatable ::  s_LDRH
   character(len=2), dimension(:,:), allocatable ::  s_ROHV

   type(t_data_header), dimension(:), allocatable :: data_header

   character(len=2060) :: header

   integer    :: TotalRays, i, j, k, n, ierr, value, nazim, FileSize
   type(date) :: date_bj, date_utc
   integer(kind=8) :: atime, btime, ctime

   select case(trim(radar_type))
   case ("CC")
   case default
      write(*,*) "in read_radar_rd:Unknown Radar Type:", trim(radar_type)
!     stop
   end select

   ! read in header only
   open(radar_unit,file=filename,recl=2048+12,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr)            &
                   rd_header%c_FileID,           &
                   rd_header%f_VersionNo,        &
                   rd_header%l_FileHeaderLength, &
                   rd_header%c_Country,          &
                   rd_header%c_Province,         &
                   rd_header%c_Station,          &
                   rd_header%c_StationNumber,    &
                   rd_header%c_RadarType,        &
                   rd_header%c_Longitude,        &
                   rd_header%c_Latitude,         &
                   rd_header%l_LongitudeValue,   &
                   rd_header%l_LatitudeValue,    &
                   rd_header%l_Height,           &
                   rd_header%s_MaxAngle,         &
                   rd_header%s_OptiAngle,        &
                   rd_header%l_AntennaG,         &
                   rd_header%us_VerBeamW,        &
                   rd_header%us_HorBeamW,        &
                   rd_header%uc_Polarizations,   &
                   rd_header%us_SideLobe,        &
                   rd_header%l_Power,            &
                   rd_header%l_WaveLength,       &
                   rd_header%us_LogA,            &
                   rd_header%us_LineA,           &
                   rd_header%us_AGCP   ,         &
                   rd_header%us_LogMinPower,     &
                   rd_header%us_LineMinPower,    &
                   rd_header%uc_ClutterT,        &
                   rd_header%uc_VelocityP,       &
                   rd_header%uc_FilterP,         &
                   rd_header%uc_NoiseT,          &
                   rd_header%uc_SQIT,            &
                   rd_header%uc_IntensityC,      &
                   rd_header%uc_IntersityR,      &
                   rd_header%uc_SType,           &
                   rd_header%us_SYear,           &
                   rd_header%uc_SMonth,          &
                   rd_header%uc_SDay,            &
                   rd_header%uc_SHour,           &
                   rd_header%uc_SMinute,         &
                   rd_header%uc_SSecond,         &
                   rd_header%uc_TimP,            &
                   rd_header%ul_SMillisecond,    &
                   rd_header%uc_Calibration,     &
                   rd_header%uc_IntensityI,      &
                   rd_header%uc_VelocityPn,      &
                   rd_header%us_ZStartBin,       &
                   rd_header%us_VStartBin,       &
                   rd_header%us_WStartBin,       &
                   rd_header%LayerInfo,          &
                   rd_header%us_RHIA,            &
                   rd_header%s_RHIL,             &
                   rd_header%s_RHIH,             &
                   rd_header%us_EYear,           &
                   rd_header%uc_EMonth,          &
                   rd_header%uc_EDay,            &
                   rd_header%uc_EHour,           &
                   rd_header%uc_EMinute,         &
                   rd_header%uc_ESecond,         &
                   rd_header%uc_ETenth,          &
                   rd_header%us_ZBinByte,        &
                   rd_header%BinRange1,          &
                   rd_header%us_VBinByte,        &
                   rd_header%BinRange2,          &
                   rd_header%us_WBinByte,        &
                   rd_header%BinRange3,          &
                   rd_header%c_StationID,        &
                   rd_header%c_JHType,           &
                   rd_header%c_Spare
   close(radar_unit)

   date_bj = init_date(get_us_value(rd_header%us_EYear  ), &
                       get_uc_value(rd_header%uc_EMonth ), &
                       get_uc_value(rd_header%uc_EDay   ), &
                       get_uc_value(rd_header%uc_EHour  ), &
                       get_uc_value(rd_header%uc_EMinute), &
                       get_uc_value(rd_header%uc_ESecond))

   date_utc=get_new_date(date_bj,-8*one_hour)
   btime=get_diff_date(date_utc,d1970)

   rdat%Year   = date_utc%Year
   rdat%Month  = date_utc%Month
   rdat%Day    = date_utc%Day
   rdat%Hour   = date_utc%Hour
   rdat%Minute = date_utc%Minute
   rdat%Second = date_utc%Second

   rdat%latitude  = get_i_value(rd_header%l_LatitudeValue )/1000.
   rdat%longitude = get_i_value(rd_header%l_LongitudeValue)/1000.
   rdat%altitude  = get_i_value(rd_header%l_Height )/1000.

   rdat%ntilt=get_uc_value(rd_header%uc_Stype)-100
   if(rdat%ntilt<0)then
      rdat%ntilt=0
      return
   endif

   MaxCuts = rdat%ntilt
   MaxRads = get_us_value(rd_header%LayerInfo(1)%us_RecordNumber)
   RGates  = get_us_value(rd_header%LayerInfo(1)%us_ZBinNumber  )
   VGates  = get_us_value(rd_header%LayerInfo(1)%us_VBinNumber  )
   WGates  = get_us_value(rd_header%LayerInfo(1)%us_WBinNumber  )

   do k=1, rdat%ntilt
       value=get_us_value(rd_header%LayerInfo(k)%us_RecordNumber)
       if(MaxRads<value)then
          MaxRads=value
       endif
       value=get_us_value(rd_header%LayerInfo(1)%us_ZBinNumber  )
       if(RGates<value)then
          RGates=value
       endif
       value=get_us_value(rd_header%LayerInfo(1)%us_VBinNumber  )
       if(VGates<value)then
          VGates=value
       endif
       value=get_us_value(rd_header%LayerInfo(1)%us_WBinNumber  )
       if(WGates<value)then
          WGates=value
       endif
       if(get_uc_value(rd_header%LayerInfo(k)%c_DataForm)/=24)then
          write(*,*) "Other c_DataForm:",get_uc_value(rd_header%LayerInfo(k)%c_DataForm)," not implement"
          rdat%ntilt =0
          return
!         stop 
       endif
   enddo

   call allocate_radar_data(rdat,RGates,VGates,WGates,MaxRads,MaxCuts)
   
   do k = 1, rdat%ntilt
      rdat%nazim   (k) = get_us_value(rd_header%LayerInfo(k)%us_RecordNumber)
      rdat%rgatesp (k) = get_us_value(rd_header%LayerInfo(k)%us_ZBinWidth   )/10.
      rdat%vgatesp (k) = get_us_value(rd_header%LayerInfo(k)%us_VBinWidth   )/10.
      rdat%vmax    (k) = get_us_value(rd_header%LayerInfo(k)%us_MaxV        )/100.
      rdat%rmax    (k) = get_us_value(rd_header%LayerInfo(k)%us_MaxL        )*10.
      rdat%rtilt (:,k) = get_s_value (rd_header%LayerInfo(k)%s_SwpAngles    )/100.
      rdat%nrgate(:,k) = get_us_value(rd_header%LayerInfo(k)%us_ZBinNumber  )
      rdat%nvgate(:,k) = get_us_value(rd_header%LayerInfo(k)%us_VBinNumber  )
      if(rdat%nrgate(1,k)>0) rdat%ifref(k)=.true.
      if(rdat%nvgate(1,k)>0) rdat%ifvel(k)=.true.
   enddo

   TotalRays=SUM(rdat%nazim(1:rdat%ntilt))
   allocate(uc_CorZ(RGates,TotalRays))
   allocate(uc_UnZ (RGates,TotalRays))
   allocate( c_V   (VGates,TotalRays))
   allocate(uc_W   (WGates,TotalRays))
   allocate(data_header   (TotalRays))

   ! read in all data
   FileSize=12+2048+(RGates*2+VGates+WGates+11)*TotalRays
   open(radar_unit,file=filename,recl=FileSize,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr) header,                     &
                                     (data_header(n),             &
                                     (uc_CorZ  (i,n),i=1,RGates), &
                                     (uc_UnZ   (i,n),i=1,RGates), &
                                     (c_V      (i,n),i=1,VGates), &
                                     (uc_W     (i,n),i=1,WGates), n=1,TotalRays)
   close(radar_unit)

   write(*,"(2(A,I10))") "Total Rays:", TotalRays,". File Size:", FileSize 

   rdat%vcp=0
   if(rdat%ntilt==14.or.rdat%ntilt==16) rdat%vcp=11
   if(rdat%ntilt== 9.or.rdat%ntilt==11) rdat%vcp=21
   if(rdat%ntilt== 5.or.rdat%ntilt== 7) rdat%vcp=31

   n=0
   do k=1, rdat%ntilt
      do j=1, rdat%nazim(k)
         n=n+1
         if(n==1)atime=(get_uc_value(data_header(n)%uc_Hh)*60+get_uc_value(data_header(n)%uc_Mm))*60+get_uc_value(data_header(n)%uc_Ss)
         ctime=(get_uc_value(data_header(n)%uc_Hh)*60+get_uc_value(data_header(n)%uc_Mm))*60+get_uc_value(data_header(n)%uc_Ss)
         if(ctime<atime)then
            ctime=ctime+86400
         endif
         rdat%stime (j,k)=(ctime-atime)+btime
         rdat%etime (j,k)=(ctime-atime)+btime
         if(rdat%vmax (k)<63.5)then
            rdat%vres (j,k)=0.5
         else
            rdat%vres (j,k)=1.0
         endif
         rdat%razim(j,k)=get_us_value(data_header(n)%us_Az)/100.
         !write(*,*) "Elev,Az:", get_us_value(data_header(n)%us_Az), &
         !                       get_uc_value(data_header(n)%uc_Hh), &
         !                       get_uc_value(data_header(n)%uc_Mm), &
         !                       get_uc_value(data_header(n)%uc_Ss)
         do i=1, rdat%nrgate(j,k)
           value=get_uc_value(uc_CorZ(i,n))
           if(value/=0)then
              rdat%ref(i,j,k)=(value-64.)/2
           endif
           !write(99,*) i,j,k,get_uc_value(uc_CorZ(i,n))
         enddo
         do i=1, rdat%nvgate(j,k)
           value=get_c_value(c_V(i,n))
           if(value/=-128)then
              rdat%vel(i,j,k)=value*rdat%vmax(k)/127
           endif
           value=get_uc_value(uc_W(i,n))
           if(value/=0)then
              rdat%spw(i,j,k)=value*rdat%vmax(k)/512
           endif
         enddo
      enddo
   enddo

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_radar_rstm(filename,radar_type,rdat)
   use bytes, only: get_i_value, get_s_value, get_f_value, &
                    get_us_value, get_uc_value
   use date_pack, only: get_new_date
   implicit none

   character(len=*),   intent(in)    :: filename, radar_type
   type(t_radar_data), intent(inout) :: rdat
 

   type t_generic_header  ! 32 Bytes Table 2-2
      character(len=4 ) :: str_magic_number !"RSTM"
      character(len=2 ) :: s_major_version
      character(len=2 ) :: s_minor_version
      character(len=4 ) :: i_generic_type !1 base, 2 product
      character(len=4 ) :: i_product_type
      character(len=16) :: reserved
   end type

   type t_site_config  ! 128 Bytes Table 2-3
      character(len=8)  :: str_site_code
      character(len=32) :: str_site_name
      character(len=4)  :: f_latitude  ! Degree
      character(len=4)  :: f_longitude !Degree
      character(len=4)  :: i_antenna_height !m
      character(len=4)  :: i_ground_height !m
      character(len=4)  :: f_frequency ! MHz
      character(len=4)  :: f_beam_width_hori ! Degree
      character(len=4)  :: f_beam_width_vert ! Degree
      character(len=4)  :: i_rda_version ! 
      character(len=2)  :: s_radar_type  ! 1-SA,2-SB,3-SC,33-CA,34-CB,35-CC,36-CCJ,37-CD,65-XA 
      character(len=54) :: reserved ! 
   end type

   type t_task_config  ! 256 Bytes Table 2-4
      character(len=32)  :: str_task_name ! VCP21 ...
      character(len=128) :: str_task_description ! 
      character(len=4)   :: i_polarization_type ! 1-hori,2-vert,3-hori+vert,4-hori/vert
      character(len=4)   :: i_scan_type ! 0-vol,1-ppi,2-rhi,3-sector,4-sector vol,5-multi rhi,6-manual 
      character(len=4)   :: i_pulse_width !nanosecond
      character(len=4)   :: i_scan_start_time !second from 1970/1/1 0:0:0
      character(len=4)   :: i_cut_number !
      character(len=4)   :: f_horizontal_noise !dBm
      character(len=4)   :: f_vertical_noise !dBm
      character(len=4)   :: f_horizontal_calibration !dB
      character(len=4)   :: f_vertical_calibration !dB
      character(len=4)   :: f_horizontal_noise_temperature !K
      character(len=4)   :: f_vertical_noise_temperature !K
      character(len=4)   :: f_zdr_calibration !dB
      character(len=4)   :: f_phidp_calibration !Degree
      character(len=4)   :: f_ldr_calibration !dB
      character(len=40)  :: reserved 
   end type

   type t_cut_config  ! 256 Bytes. Table 2-5
      character(len=4 ) :: i_process_mode                      ! 1-PPP,2-FFT
      character(len=4 ) :: i_wave_form                         ! 0-CS,1-CD,2-CDX,3-Rx Test,4-Batch,5-Dual PRF,6-Staggered PRT
      character(len=4 ) :: f_prf1                              ! Hz
      character(len=4 ) :: f_prf2                              ! Hz
      character(len=4 ) :: i_dealiasing_mode                   ! 1-1PRF,2-2PRF3:2,3-2PRF4:3,4-2PRF5:4
      character(len=4 ) :: f_azimuth                           ! azimuth of RHI mode
      character(len=4 ) :: f_elevation                         ! elevation of PPI mode
      character(len=4 ) :: f_start_angle                       ! start azimuth of PPI, or high elevation of RHI
      character(len=4 ) :: f_end_angle                         ! end azimuth of PPI, or low elevation of RHI
      character(len=4 ) :: f_angular_resolution                ! Degree, PPI Only
      character(len=4 ) :: f_scan_speed                        ! Degree/Second, PPI azimuth or RHI elevation
      character(len=4 ) :: i_log_resolution                    ! m, distance resolution of power 
      character(len=4 ) :: i_doppler_resolution                ! m, distance resolution of velocity 
      character(len=4 ) :: i_maximum_range1                    ! m, max distance of PRF1
      character(len=4 ) :: i_maximum_range2                    ! m, max distance of PRF2
      character(len=4 ) :: i_start_range                       ! m
      character(len=4 ) :: i_sample1                           ! sample number of PRF1
      character(len=4 ) :: i_sample2                           ! sample number of PRF2
      character(len=4 ) :: i_phase_mode                        ! 1-fix,2-random,3-SZ code 
      character(len=4 ) :: f_atmosphere_loss                   ! dB/km 
      character(len=4 ) :: f_nyquist_speed                     ! m/s
      character(len=8 ) :: l_moments_mask                      ! 0-deny,1-allow, see table 2-6
      character(len=8 ) :: l_moments_size_mask                 ! 0-1Byte,1-2Bytes,
      character(len=4 ) :: i_misc_filter_mask                  ! 0-Not Apllied,1-Applied, see table 2-7
      character(len=4 ) :: f_sqi_threshold                     !
      character(len=4 ) :: f_sig_threshold                     !
      character(len=4 ) :: f_csr_threshold                     !
      character(len=4 ) :: f_log_threshold                     !
      character(len=4 ) :: f_cpa_threshold                     !
      character(len=4 ) :: f_pmi_threshold                     !
      character(len=4 ) :: thresholds_reserved                 !
      character(len=4 ) :: i_dbt_mask                          ! 0-Not Apllied,1-Applied, see table 2-8
      character(len=4 ) :: i_dbz_mask                          ! 0-Not Apllied,1-Applied, see table 2-8
      character(len=4 ) :: i_velocity_mask                     ! 0-Not Apllied,1-Applied, see table 2-8
      character(len=4 ) :: i_spectrum_width_mask               ! 0-Not Apllied,1-Applied, see table 2-8
      character(len=4 ) :: i_dp_width_mask                     ! 0-Not Apllied,1-Applied, see table 2-8
      character(len=12) :: mask_reserved                       ! 
      character(len=4 ) :: i_scan_sync                         !  reserved for multi radar sync
      character(len=4 ) :: i_direction                         ! 1-clock-wise,2-anti-clock
      character(len=2 ) :: s_ground_clutter_classifier_type    ! 1-no filter,2-all filter,3-dynamic filter,4-static filter
      character(len=2 ) :: s_ground_clutter_filter_type        ! 1-no filter,1-FDAF,2-Fix Band,3-Var Band,4-Var Std,5-IIR
      character(len=2 ) :: s_ground_clutter_filter_notch_width ! 0.1m/s
      character(len=2 ) :: s_ground_clutter_filter_window      ! 0-rect,1-Hanming,2-Blackman,3-Adaptive,4-No
      character(len=76) :: reserved                            !
   end type

   !Table 2-6 Moments Mask
   integer, parameter :: flag_moment_dbt = 1 ! Total Reflectivity (Before Filter)
   integer, parameter :: flag_moment_dbz = 2 ! Reflectivity (After Filter)
   integer, parameter :: flag_moment_v   = 3 ! Doppler Velocity
   integer, parameter :: flag_moment_w   = 4 ! Spectrum Width
   integer, parameter :: flag_moment_sqi = 5 ! Signal Quality Index
   integer, parameter :: flag_moment_cpa = 6 ! Clutter Phase Alignment
   integer, parameter :: flag_moment_zdr = 7 ! Differential Reflectivity
   integer, parameter :: flag_moment_ldr = 8 ! Liner Differential Ratio
   integer, parameter :: flag_moment_cc  = 9 ! Cross Correlation Coefficient
   integer, parameter :: flag_moment_fdp =10 ! Differential Phase
   integer, parameter :: flag_moment_kdp =11 ! Specific Differential Phase
   integer, parameter :: flag_moment_cp  =12 ! Clutter Probability
   integer, parameter :: flag_moment_hcl =14 ! Hydro Classification
   integer, parameter :: flag_moment_cf  =15 ! Clutter Flag
   integer, parameter :: flag_moment_snr =16 ! Signal Noise Ratio
   integer, parameter :: flag_moment_zc  =32 ! Corrected Reflectivity
   integer, parameter :: flag_moment_vc  =33 ! Corrected Doppler Velocity
   integer, parameter :: flag_moment_wc  =34 ! Corrected Spectrum Width
   integer, parameter :: flag_moment_zdrc=35 ! Corrected Differential Reflectivity
   
   type t_radial_header ! 64 Bytes
      character(len=4 ) :: i_radial_state     ! 0-elev start,1-medium,2-elev end,3-vol start,4-vol end,5-RHI start,6-RHI end
      character(len=4 ) :: i_spot_blank       ! 0-normal,1-hide
      character(len=4 ) :: i_sequence_number  ! radial number of vol scan
      character(len=4 ) :: i_radial_number    ! radial number of elevation
      character(len=4 ) :: i_elevation_number ! elevation number of vol scan 
      character(len=4 ) :: f_azimuth          ! Degree
      character(len=4 ) :: f_elevation        ! Degree
      character(len=4 ) :: i_seconds          ! From 1970/1/1 0:0:0
      character(len=4 ) :: i_microseconds     ! 
      character(len=4 ) :: i_length_of_data   ! Bytes 
      character(len=4 ) :: i_moment_number    ! 
      character(len=20) :: reserved           ! 
   end type

   type t_moment_header ! 32 Bytes
      character(len=4) :: i_data_type  ! see table 2-6
      character(len=4) :: i_scale      ! 
      character(len=4) :: i_offset     ! 
      character(len=2) :: s_bin_length ! Bytes 1 or 2
      character(len=2) :: s_flags      ! not used 
      character(len=4) :: i_length     ! Bytes (header not included) 
      character(len=12) :: reserved    !
   end type

   ! Var = ( Dat-offset)/Scale 
   ! 0-less than Threshold, 1- RangeFold, 2-NotScan, 3-Unknown, 4-Reserved
   type t_moment_data !
      character(len=1), dimension(:), allocatable :: dat
   end type

   type(t_generic_header) :: generic_header
   type(t_site_config   ) :: site_config
   type(t_task_config   ) :: task_config
   type(t_cut_config   ), dimension(:), allocatable :: cut_config

   type(t_radial_header), dimension(:  ), allocatable :: radial_header
   type(t_moment_header), dimension(:,:), allocatable :: moment_header
   type(t_moment_data  ), dimension(:,:), allocatable :: moment_data

   character(len=416) :: header

   integer, parameter                   :: max_radial=15000, max_moment=32
   integer                              :: n_radial, n_cut
   integer, dimension(max_radial)       :: n_moment, n_radial_length
   integer, dimension(max_radial,max_moment) :: n_moment_length

   integer, dimension(:), allocatable :: n_cut_radial, n_rgate, n_vgate

   integer      :: len_record, ival, ipos, jpos, ioffset, iscale, nbyte, itype
   integer      :: TotalRays, i, j, k, m, n, ierr, value, FileSize
   type(date)   :: date_bj, date_utc, sdate, edate
   real         :: dtime
   real(kind=8) :: btime

   character(len=1), dimension(:), allocatable :: str_file
   real, dimension(:), allocatable :: dat
   
   logical, dimension(:), allocatable :: ifvel, ifref, ifzdr, ifcc, ifkdp, ifsnr

   inquire(file=filename,size=filesize)
   allocate(str_file(filesize))
   open(radar_unit,file=filename,recl=filesize,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr) str_file
   close(radar_unit)

   !generic type:9-12
   if(get_i_value(str_file(9:12))/=1)then
      write(*,*) "Not Radar Base Data"
      return
   endif

   !scan type:325-328
   ival=get_i_value(str_file(325:328))
   if(ival/=1.and.ival/=0)then
      write(*,*) "Not Vol or PPI Scan"
       return
   endif

   !cut number:337-340
   n_cut=get_i_value(str_file(337:340))
   if(n_cut<0)then
      n_cut=0
      return
   endif
   
   allocate(n_cut_radial(n_cut))
   !32+128+256+n_cut*256
   ival=32+128+256+n_cut*256
   i=1
   do while(ival<filesize)
      ipos=ival+41
      jpos=ival+44
!     write(*,*) "i,ival,ipos:",i, ival, ipos
      !if(i>10) exit
      n_moment(i)=get_i_value(str_file(ipos:jpos))
      ipos=ival+64
      do j=1, n_moment(i)
         jpos=ipos+17
         n_moment_length(i,j)=get_i_value(str_file(jpos:jpos+3))
         ipos=ipos+n_moment_length(i,j)+32
!        write(*,*) j, n_moment_length(i,j)
      enddo
      ipos=ival+37
      jpos=ival+40
      n_radial_length(i)=get_i_value(str_file(ipos:jpos))
      ival=ival+n_radial_length(i)+64
!     write(*,*) i, n_radial_length(i)
      i=i+1
   enddo
   n_radial=i-1

   ival=maxval(n_moment(1:n_radial))
   write(601,*) "n_radial, n_moment:", n_radial, ival
   allocate(cut_config(n_cut))
   allocate(radial_header(n_radial))
   allocate(moment_header(n_radial,ival))
   allocate(moment_data  (n_radial,ival))

   allocate(ifref(n_cut))
   allocate(ifvel(n_cut))
   allocate(ifzdr(n_cut))
   allocate(ifcc (n_cut))
   allocate(ifkdp(n_cut))
   allocate(ifsnr(n_cut))

   do i=1, n_radial
!     write(*,*) "n_moment,length:",n_moment(i),n_radial_length(i)
      do j=1, n_moment(i)
!        write(*,*) "length:",n_moment_length(i,j)
         allocate(moment_data(i,j)%dat(n_moment_length(i,j)))
      enddo
   enddo

   ! read in all data
   open(radar_unit,file=filename,recl=filesize,status="old",access="direct")
   read(radar_unit,rec=1,iostat=ierr)  &
            generic_header%str_magic_number                   , &
            generic_header%s_major_version                    , &
            generic_header%s_minor_version                    , &
            generic_header%i_generic_type                     , &
            generic_header%i_product_type                     , &
            generic_header%reserved                           , &
            site_config%str_site_code                         , &
            site_config%str_site_name                         , &
            site_config%f_latitude                            , &
            site_config%f_longitude                           , &
            site_config%i_antenna_height                      , &
            site_config%i_ground_height                       , &
            site_config%f_frequency                           , &
            site_config%f_beam_width_hori                     , &
            site_config%f_beam_width_vert                     , &
            site_config%i_rda_version                         , &
            site_config%s_radar_type                          , &
            site_config%reserved                              , &
            task_config%str_task_name                         , &
            task_config%str_task_description                  , &
            task_config%i_polarization_type                   , &
            task_config%i_scan_type                           , &
            task_config%i_pulse_width                         , &
            task_config%i_scan_start_time                     , &
            task_config%i_cut_number                          , &
            task_config%f_horizontal_noise                    , &
            task_config%f_vertical_noise                      , &
            task_config%f_horizontal_calibration              , &
            task_config%f_vertical_calibration                , &
            task_config%f_horizontal_noise_temperature        , &
            task_config%f_vertical_noise_temperature          , &
            task_config%f_zdr_calibration                     , &
            task_config%f_phidp_calibration                   , &
            task_config%f_ldr_calibration                     , &
            task_config%reserved                              , &
           (cut_config(k)%i_process_mode                      , &
            cut_config(k)%i_wave_form                         , &
            cut_config(k)%f_prf1                              , &
            cut_config(k)%f_prf2                              , &
            cut_config(k)%i_dealiasing_mode                   , &
            cut_config(k)%f_azimuth                           , &
            cut_config(k)%f_elevation                         , &
            cut_config(k)%f_start_angle                       , &
            cut_config(k)%f_end_angle                         , &
            cut_config(k)%f_angular_resolution                , &
            cut_config(k)%f_scan_speed                        , &
            cut_config(k)%i_log_resolution                    , &
            cut_config(k)%i_doppler_resolution                , &
            cut_config(k)%i_maximum_range1                    , &
            cut_config(k)%i_maximum_range2                    , &
            cut_config(k)%i_start_range                       , &
            cut_config(k)%i_sample1                           , &
            cut_config(k)%i_sample2                           , &
            cut_config(k)%i_phase_mode                        , &
            cut_config(k)%f_atmosphere_loss                   , &
            cut_config(k)%f_nyquist_speed                     , &
            cut_config(k)%l_moments_mask                      , &
            cut_config(k)%l_moments_size_mask                 , &
            cut_config(k)%i_misc_filter_mask                  , &
            cut_config(k)%f_sqi_threshold                     , &
            cut_config(k)%f_sig_threshold                     , &
            cut_config(k)%f_csr_threshold                     , &
            cut_config(k)%f_log_threshold                     , &
            cut_config(k)%f_cpa_threshold                     , &
            cut_config(k)%f_pmi_threshold                     , &
            cut_config(k)%thresholds_reserved                 , &
            cut_config(k)%i_dbt_mask                          , &
            cut_config(k)%i_dbz_mask                          , &
            cut_config(k)%i_velocity_mask                     , &
            cut_config(k)%i_spectrum_width_mask               , &
            cut_config(k)%i_dp_width_mask                     , &
            cut_config(k)%mask_reserved                       , &
            cut_config(k)%i_scan_sync                         , &
            cut_config(k)%i_direction                         , &
            cut_config(k)%s_ground_clutter_classifier_type    , &
            cut_config(k)%s_ground_clutter_filter_type        , &
            cut_config(k)%s_ground_clutter_filter_notch_width , &
            cut_config(k)%s_ground_clutter_filter_window      , &
            cut_config(k)%reserved , k=1, n_cut)              , &
           ( radial_header(i)%i_radial_state                  , & 
             radial_header(i)%i_spot_blank                    , & 
             radial_header(i)%i_sequence_number               , & 
             radial_header(i)%i_radial_number                 , & 
             radial_header(i)%i_elevation_number              , & 
             radial_header(i)%f_azimuth                       , & 
             radial_header(i)%f_elevation                     , & 
             radial_header(i)%i_seconds                       , & 
             radial_header(i)%i_microseconds                  , & 
             radial_header(i)%i_length_of_data                , & 
             radial_header(i)%i_moment_number                 , & 
             radial_header(i)%reserved                        , & 
            (moment_header(i,j)%i_data_type                   , & 
             moment_header(i,j)%i_scale                       , & 
             moment_header(i,j)%i_offset                      , & 
             moment_header(i,j)%s_bin_length                  , & 
             moment_header(i,j)%s_flags                       , & 
             moment_header(i,j)%i_length                      , & 
             moment_header(i,j)%reserved                      , &
            (moment_data  (i,j)%dat(n),n=1,n_moment_length(i,j)) , j=1,n_moment(i)),i=1,n_radial)    !
   close(radar_unit)

   btime=get_i_value(radial_header(1)%i_seconds)
   sdate=get_new_date(d1970,DBLE(get_i_value(radial_header(1       )%i_seconds)))
   edate=get_new_date(d1970,DBLE(get_i_value(radial_header(n_radial)%i_seconds)))
   date_utc=edate
!
   rdat%Year   = date_utc%Year
   rdat%Month  = date_utc%Month
   rdat%Day    = date_utc%Day
   rdat%Hour   = date_utc%Hour
   rdat%Minute = date_utc%Minute
   rdat%Second = date_utc%Second
!
   rdat%latitude  = get_f_value(site_config%f_latitude)
   rdat%longitude = get_f_value(site_config%f_longitude)
   rdat%altitude  = get_i_value(site_config%i_antenna_height)
!
   rdat%ntilt = n_cut 
!
!   do k=1, rdat%ntilt
!      write(*,*) k, get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_Prf1    ), &
!                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_Prf2    ), &
!                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_spulseW ), &
!                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_MaxV    ), &
!                    get_us_value(header_784%RadarObservationInfo%LayerParam( k)%us_MaxL    )
! 
!   enddo
   ifref=.false.
   ifvel=.false.
   ifzdr=.false.
   ifcc =.false.
   ifkdp=.false.
   ifsnr=.false.
   MaxCuts = rdat%ntilt
   MaxRads = 0 
   allocate(n_rgate(n_cut))
   allocate(n_vgate(n_cut))
   n_rgate=0
   n_vgate=0
   do i=1, n_radial
      k=get_i_value(radial_header(i)%i_elevation_number)
      j=get_i_value(radial_header(i)%i_radial_number)
  !   write(*,*) "j,k:",j,k
      MaxRads= max(MaxRads,j)
      do j=1, n_moment(i)
         ival=get_i_value(moment_header(i,j)%i_data_type)
         if(ival==flag_moment_dbz)then
            n_rgate(k)=max(n_rgate(k),get_i_value(moment_header(i,j)%i_length)/get_s_value(moment_header(i,j)%s_bin_length))
            ifref(k)=.true.
         endif
         if(ival==flag_moment_v)then
            n_vgate(k)=max(n_vgate(k),get_i_value(moment_header(i,j)%i_length)/get_s_value(moment_header(i,j)%s_bin_length))
            ifvel(k)=.true.
         endif
         if(ival==flag_moment_zdr)then
            ifzdr(k)=.true.
         endif
         if(ival==flag_moment_cc)then
            ifcc(k)=.true.
         endif
         if(ival==flag_moment_kdp)then
            ifkdp(k)=.true.
         endif
         if(ival==flag_moment_snr)then
            ifsnr(k)=.true.
         endif
      enddo
   enddo

   RGates  = maxval(n_rgate) 
   VGates  = maxval(n_vgate) 
   WGates  = maxval(n_vgate) 

   write(601,*) "RGates,VGates,WGates,MaxRads,MaxCuts:",RGates,VGates,WGates,MaxRads,MaxCuts
   call allocate_radar_data(rdat,RGates,VGates,WGates,MaxRads,MaxCuts)
   rdat%ifref=ifref
   rdat%ifvel=ifvel
   rdat%ifzdr=ifzdr
   rdat%ifcc =ifcc 
   rdat%ifkdp=ifkdp
   rdat%ifsnr=ifsnr
 
!   
   write(601,"(9A10)") "k","range1","range2","rgatesp","vgatesp","vmax","Sangle","Eangle","ang_res"
   do k=1, rdat%ntilt
      rdat%vmax   (k) = get_f_value(cut_config(k)%f_nyquist_speed) 
      write(601,"(5I10,4F10.2)") k ,&
                 get_i_value(cut_config(k)%i_maximum_range1),&
                 get_i_value(cut_config(k)%i_maximum_range2),&
                 get_i_value(cut_config(k)%i_log_resolution), &
                 get_i_value(cut_config(k)%i_doppler_resolution), &
                 get_f_value(cut_config(k)%f_nyquist_speed),&
                 get_f_value(cut_config(k)%f_start_angle),& 
                 get_f_value(cut_config(k)%f_end_angle),&
                 get_f_value(cut_config(k)%f_angular_resolution)
      rdat%rmax   (k) = min(get_i_value(cut_config(k)%i_maximum_range1),get_i_value(cut_config(k)%i_maximum_range2)) 
      rdat%rgatesp(k) = get_i_value(cut_config(k)%i_log_resolution)
      rdat%vgatesp(k) = get_i_value(cut_config(k)%i_doppler_resolution)
   enddo

   filesize=32+128+256+256*n_cut
   rdat%nazim = 0 
   do i=1, n_radial 
      filesize=filesize+64
      j=get_i_value(radial_header(i)%i_radial_number)
      k=get_i_value(radial_header(i)%i_elevation_number)
      rdat%nazim  (k) = max(rdat%nazim(k),j) 
      do n=1, n_moment(i)
         len_record=get_i_value(moment_header(i,n)%i_length)
         filesize=filesize+32+len_record
         ival=get_i_value(moment_header(i,n)%i_data_type)
         if(ival==flag_moment_dbz)then
            rdat%nrgate(j,k) = get_i_value(moment_header(i,n)%i_length)/get_s_value(moment_header(i,n)%s_bin_length)
         endif
         if(ival==flag_moment_v)then
            rdat%nvgate(j,k) = get_i_value(moment_header(i,n)%i_length)/get_s_value(moment_header(i,n)%s_bin_length)
         endif
      enddo
      rdat%rtilt (j,k) = get_f_value(radial_header(i)%f_elevation)
      rdat%razim (j,k) = get_f_value(radial_header(i)%f_azimuth)
   enddo
!
   TotalRays=SUM(rdat%nazim(1:rdat%ntilt))
!
   write(*,"(2(A,I10))") "Total Rays:", TotalRays,". File Size:", FileSize 

   read(task_config%str_task_name(4:5),"(I2)") rdat%vcp
!
   n=max(Rgates, Vgates)
   allocate(dat(n))
   
   write(601,*) "n_radial:",n_radial
   write(601,"(5A8,2A10,A12,A8)") "m","state","spot","j","k","elev","azim","sec","ms" 
   write(603,"(3A5,8A10)") "I","J","K","Ref","Vel","Spw","Zdr","CC","Fdp","Kdp","SNR"
   do m=1, n_radial 
      j=get_i_value(radial_header(m)%i_radial_number)
      k=get_i_value(radial_header(m)%i_elevation_number)
      write(601,"(5I8,2F10.2,I12,I8)") m,get_i_value(radial_header(m)%i_radial_state),&
                                         get_i_value(radial_header(m)%i_spot_blank  ),&
                                     j,k,get_f_value(radial_header(m)%f_elevation   ),&
                                         get_f_value(radial_header(m)%f_azimuth     ),&
                                         get_i_value(radial_header(m)%i_seconds     ),&
                                         get_i_value(radial_header(m)%i_microseconds)
      do n=1, n_moment(m) 
            nbyte      = get_s_value(moment_header(m,n)%s_bin_length)
            len_record = get_i_value(moment_header(m,n)%i_length)/nbyte
            ioffset    = get_i_value(moment_header(m,n)%i_offset)
            iscale     = get_i_value(moment_header(m,n)%i_scale)
            ival       = get_i_value(moment_header(m,n)%i_data_type)
            write(602,*) "i_data_type:", m, n, ival,&
                                          "dbt = 1,"//&
                                          "dbz = 2,"//&
                                          "v   = 3,"//&
                                          "w   = 4,"//&
                                          "sqi = 5,"//&
                                          "cpa = 6,"//&
                                          "zdr = 7,"//&
                                          "ldr = 8,"//&
                                          "cc  = 9,"//&
                                          "dp  =10,"//&
                                          "kdp =11,"//&
                                          "cp  =12,"//&
                                          "hcl =14,"//&
                                          "cf  =15,"//&
                                          "snr =16,"//&
                                          "zc  =32,"//&
                                          "vc  =33,"//&
                                          "wc  =34,"//&
                                          "zdrc=35,"

            do i=1, len_record
               if(nbyte==1)then
                   dat(i)= get_rstm_value_1(moment_data(m,n)%dat(i*nbyte  :i*nbyte),iscale,ioffset)
               elseif(nbyte==2)then
                   dat(i)= get_rstm_value_2(moment_data(m,n)%dat(i*nbyte-1:i*nbyte),iscale,ioffset)
               else
                   write(*,*) "Moment nByte>2!!!",nbyte
               endif
            enddo
         if(ival==flag_moment_dbz)then
            rdat%ref(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_v)then
            rdat%vel(1:len_record,j,k)=dat(1:len_record)
            rdat%vres(j,k)=1./iscale
         endif
         if(ival==flag_moment_w)then
            rdat%spw(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_zdr)then
            rdat%zdr(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_cc)then
            rdat%cc(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_fdp)then
            rdat%fdp(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_kdp)then
            rdat%kdp(1:len_record,j,k)=dat(1:len_record)
         endif
         if(ival==flag_moment_snr)then
            rdat%snr(1:len_record,j,k)=dat(1:len_record)
         endif
         rdat%stime(j,k)=get_i_value(radial_header(m)%i_seconds)
         rdat%etime(j,k)=get_i_value(radial_header(m)%i_seconds)
         !do i=1, len_record
         !   if(rdat%ref(i,j,k)/=value_invalid.or. &
         !      rdat%vel(i,j,k)/=value_invalid.or. & 
         !      rdat%spw(i,j,k)/=value_invalid.or. & 
         !      rdat%zdr(i,j,k)/=value_invalid.or. &
         !      rdat%cc (i,j,k)/=value_invalid.or. & 
         !      rdat%fdp(i,j,k)/=value_invalid.or. & 
         !      rdat%kdp(i,j,k)/=value_invalid.or. & 
         !      rdat%snr(i,j,k)/=value_invalid)then
         !      write(603,"(3I5,8F10.2)") i, j, k, rdat%ref(i,j,k), rdat%vel(i,j,k), rdat%spw(i,j,k), rdat%zdr(i,j,k), rdat%cc(i,j,k), rdat%fdp(i,j,k), rdat%kdp(i,j,k), rdat%snr(i,j,k)
         !   endif
         !enddo
      enddo
   enddo

!  stop
   contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     real function get_rstm_value_1(code,iscale,ioffset)
     character(len=1), dimension(1), intent(in) :: code
     integer, intent(in) :: iscale, ioffset
   ! Var = ( Dat-offset)/Scale 
   ! 0-less than Threshold, 1- RangeFold,2-NotScan,3-Unknown,4-Reserved
     integer :: ival
   
     ival=get_uc_value(code(1))
     if(ival<=5)then
        get_rstm_value_1=VALUE_INVALID
        if(ival==1) get_rstm_value_1=VALUE_RANFOLD
        return
     else
        get_rstm_value_1=(DBLE(ival)-ioffset)/iscale
     endif
     end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     real function get_rstm_value_2(code,iscale,ioffset)
     character(len=1), dimension(2), intent(in) :: code
     integer, intent(in) :: iscale, ioffset
   ! Var = ( Dat-offset)/Scale 
   ! 0-less than Threshold, 1- RangeFold,2-NotScan,3-Unknown,4-Reserved
     integer :: ival
     character(len=2) :: scode

     scode(1:1)=code(1)
     scode(2:2)=code(2)
   
     ival=get_us_value(scode)
     if(ival<=5)then
        get_rstm_value_2=VALUE_INVALID
        if(ival==1) get_rstm_value_2=VALUE_RANFOLD
        return
     else
        get_rstm_value_2=(DBLE(ival)-ioffset)/iscale
     endif
     end function
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_radar_grads_station(filename, radar_data)
   use bytes, only: is_big_endian
   use date_pack, only: month_name
   use libradar, only: beamhgt, gcircle
   implicit none
   character(len=*),   intent(in) :: filename
   type(t_radar_data), intent(inout) :: radar_data

   character(len=800) :: fileout, ctlfile, datfile, mapfile

   integer :: i, j, k, n, iazim, i1
   real    :: lat1, lon1, alt, lat0, lon0, alt0, range, height, elev, sfcrng,azim, dazim
   real    :: lat2, lon2, lat3, lon3

   integer           :: NFLAG, NLEV
   character(len=8)  :: STID
   real              :: TIM, VEL, SPW, REF, CC, SNR, ZDR, FDP, KDP

   character(len=80) :: endian

   if(is_big_endian)then
      endian="big_endian"
   else
      endian="little_endian"
   endif

   call calculate_radar_loc(radar_data)
   lat0=radar_data%latitude
   lon0=radar_data%longitude
   alt0=radar_data%altitude

   n=0
   do k=1, radar_data%ntilt
      if(radar_data%ifvel(k))then
         n=n+1
         !if(n>2) exit
         write(ctlfile,"(A,I2.2,A)") trim(filename)//".vel.",n,".ctl"
         write(datfile,"(A,I2.2,A)") trim(filename)//".vel.",n,".dat"
         write(mapfile,"(A,I2.2,A)") trim(filename)//".vel.",n,".map"
         open(radar_unit,file=ctlfile,status="unknown")
         write(radar_unit,"(A)") "DSET ^"//trim(datfile)
         write(radar_unit,"(A)") "DTYPE  station"
         write(radar_unit,"(A)") "OPTIONS sequential "//trim(endian)
         write(radar_unit,"(A)") "STNMAP "//trim(mapfile)
         write(radar_unit,"(A,F8.1)") "UNDEF  ",value_invalid
         write(radar_unit,"(A)") "TITLE  Station Data Sample"
         write(radar_unit,"(A,I2.2,':',I2.2,'z',I2.2,A3,I4.4,A)") &
              "TDEF   1 linear ",radar_data%hour,radar_data%minute,&
              radar_data%day,month_name(radar_data%month),radar_data%year," 5mn"
         write(radar_unit,"(A)") "VARS 2"
         write(radar_unit,"(A)") "v    0  99  Radical Velocity "
         write(radar_unit,"(A)") "w    0  99  Spectrum Width "
         write(radar_unit,"(A)") "ENDVARS"

         open(radar_unit,file=datfile,form='unformatted',status="unknown")
         write(*,"(2A)") "Writing vel grads station file:", trim(datfile)
         TIM = 0.0
         NLEV = 1
         NFLAG = 1
         do j=1,radar_data%nazim(k)
            elev=radar_data%rtilt(j,k)
            azim=radar_data%razim(j,k)
            do i=1,radar_data%nvgate(j,k)
               VEL=radar_data%vel(i,j,k)
               SPW=radar_data%spw(i,j,k)
               if(ABS(VEL)>200)then
                  VEL=value_invalid
                  SPW=value_invalid
               endif
               if(VEL/=value_invalid)then
!.and.radar_data%spw(i,j,k)/=value_invalid
                  range=i*radar_data%vgatesp(k)
                  call beamhgt(elev,range,height,sfcrng)
                  !call gcircle(lat0,lon0,azim,sfcrng,lat1,lon1)
                  lat1=radar_data%vlat(i,j,k)
                  lon1=radar_data%vlon(i,j,k)
                  write(STID,'(I4,I4)') i,j
                  write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                  !if(ABS(radar_data%vel(i,j,k))<200.)then
                  !   vel=radar_data%vel(i,j,k)
                  !else
                  !   vel=value_invalid
                  !endif
                  !spw=radar_data%spw(i,j,k)
                  write(radar_unit) VEL, SPW
                  if(sfcrng*RADIAN>radar_data%vgatesp(k))then
                    iazim=ceiling(sfcrng*RADIAN/radar_data%vgatesp(k)/2)
                     do i1=iazim,0,-1
                        dazim=azim+(i1+0.0)*radar_data%vgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat1,lon1)
                        write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                        write(radar_unit) VEL, SPW
                        dazim=azim-(i1+0.0)*radar_data%vgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat1,lon1)
                        write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                        write(radar_unit) VEL, SPW
                     enddo
                  endif
               endif
            enddo
         enddo
         NLEV = 0
         WRITE(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
         close(radar_unit)
      endif
   enddo

   n=0
   do k=1, radar_data%ntilt
      if(radar_data%ifref(k))then
         n=n+1
         !if(n>1) exit
         write(ctlfile,"(A,I2.2,A)") trim(filename)//".ref.",n,".ctl"
         write(datfile,"(A,I2.2,A)") trim(filename)//".ref.",n,".dat"
         write(mapfile,"(A,I2.2,A)") trim(filename)//".ref.",n,".map"
         open(radar_unit,file=ctlfile,status="unknown")
         write(radar_unit,"(A)") "DSET ^"//trim(datfile)
         write(radar_unit,"(A)") "DTYPE  station"
         write(radar_unit,"(A)") "OPTIONS sequential "//trim(endian)
         write(radar_unit,"(A)") "STNMAP "//trim(mapfile)
         write(radar_unit,"(A,F8.1)") "UNDEF  ",value_invalid
         write(radar_unit,"(A)") "TITLE  Station Data Sample"
         write(radar_unit,"(A,I2.2,':',I2.2,'z',I2.2,A3,I4.4,A)") &
              "TDEF   1 linear ",radar_data%hour,radar_data%minute,&
              radar_data%day,month_name(radar_data%month),radar_data%year," 5mn"
         write(radar_unit,"(A)") "VARS 1"
         write(radar_unit,"(A)") "z    0  99  Reflectivity"
         write(radar_unit,"(A)") "ENDVARS"

         open(radar_unit,file=datfile,form='unformatted',status="unknown")
         write(*,"(2A)") "Writing ref grads station file:", trim(datfile)
         TIM = 0.0
         NLEV = 1
         NFLAG = 1
         do j=1,radar_data%nazim(k)
            elev=radar_data%rtilt(j,k)
            azim=radar_data%razim(j,k)
            do i=1,radar_data%nrgate(j,k)
               REF=radar_data%ref(i,j,k)
               if(ABS(REF)>200)then
                  REF=value_invalid
               endif
               if(REF/=value_invalid)then
                  range=i*radar_data%rgatesp(k)
                  call beamhgt(elev,range,height,sfcrng)
                  !call gcircle(lat0,lon0,azim,sfcrng,lat1,lon1)
                  lat1=radar_data%rlat(i,j,k)
                  lon1=radar_data%rlon(i,j,k)
                  write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                  write(radar_unit) REF 
!write(79,*) i,i,k,lat1,lon1,ref
                  if(sfcrng*RADIAN>radar_data%rgatesp(k))then
                    iazim=ceiling(sfcrng*RADIAN/radar_data%rgatesp(k)/2)
                     do i1=iazim,0,-1
                        dazim=azim+(i1+0.0)*radar_data%rgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat2,lon2)
                        write(radar_unit) STID,lat2,lon2,TIM,NLEV,NFLAG
                        write(radar_unit) REF 
                        dazim=azim-(i1+0.0)*radar_data%rgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat3,lon3)
                        write(radar_unit) STID,lat3,lon3,TIM,NLEV,NFLAG
                        write(radar_unit) REF 
!write(67,*) i,j,k,i1,iazim,lat1, lat2, lat3, lon1, lon2, lon3, radar_data%ref(i,j,k)
                     enddo
                  endif

               endif
            enddo
         enddo
         NLEV = 0
         WRITE(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
         close(radar_unit)
      endif

   enddo

   n=0
   do k=1, radar_data%ntilt
      if(radar_data%ifzdr(k))then
         n=n+1
         !if(n>1) exit
         write(ctlfile,"(A,I2.2,A)") trim(filename)//".zdr.",n,".ctl"
         write(datfile,"(A,I2.2,A)") trim(filename)//".zdr.",n,".dat"
         write(mapfile,"(A,I2.2,A)") trim(filename)//".zdr.",n,".map"
         open(radar_unit,file=ctlfile,status="unknown")
         write(radar_unit,"(A)") "DSET ^"//trim(datfile)
         write(radar_unit,"(A)") "DTYPE  station"
         write(radar_unit,"(A)") "OPTIONS sequential "//trim(endian)
         write(radar_unit,"(A)") "STNMAP "//trim(mapfile)
         write(radar_unit,"(A,F8.1)") "UNDEF  ",value_invalid
         write(radar_unit,"(A)") "TITLE  Station Data Sample"
         write(radar_unit,"(A,I2.2,':',I2.2,'z',I2.2,A3,I4.4,A)") &
              "TDEF   1 linear ",radar_data%hour,radar_data%minute,&
              radar_data%day,month_name(radar_data%month),radar_data%year," 5mn"
         write(radar_unit,"(A)") "VARS 5"
         write(radar_unit,"(A)") "zdr    0  99  Diferential Reflectivity"
         write(radar_unit,"(A)") "cc     0  99  Cross Correlation Coefficient"
         write(radar_unit,"(A)") "fdp    0  99  Differential Phase"
         write(radar_unit,"(A)") "kdp    0  99  Specific Differential Phase"
         write(radar_unit,"(A)") "snr    0  99  Signal Noise Ratio"
         write(radar_unit,"(A)") "ENDVARS"

         open(radar_unit,file=datfile,form='unformatted',status="unknown")
         write(*,"(2A)") "Writing zdr grads station file:", trim(datfile)
         TIM = 0.0
         NLEV = 1
         NFLAG = 1
         do j=1,radar_data%nazim(k)
            elev=radar_data%rtilt(j,k)
            azim=radar_data%razim(j,k)
            do i=1,radar_data%nrgate(j,k)
               ZDR=radar_data%zdr(i,j,k)
                CC=radar_data%cc (i,j,k)
               FDP=radar_data%fdp(i,j,k)
               KDP=radar_data%kdp(i,j,k)
               SNR=radar_data%snr(i,j,k)
               if(ABS(ZDR)>200)then
                  ZDR=value_invalid
                   CC=value_invalid
                  FDP=value_invalid
                  KDP=value_invalid
                  SNR=value_invalid
               endif
               if(ZDR/=value_invalid)then
                  range=i*radar_data%rgatesp(k)
                  call beamhgt(elev,range,height,sfcrng)
                  lat1=radar_data%rlat(i,j,k)
                  lon1=radar_data%rlon(i,j,k)
                  write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                  write(radar_unit) ZDR, CC, FDP, KDP, SNR
                  if(sfcrng*RADIAN>radar_data%rgatesp(k))then
                    iazim=ceiling(sfcrng*RADIAN/radar_data%rgatesp(k)/2)
                     do i1=iazim,0,-1
                        dazim=azim+(i1+0.0)*radar_data%rgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat1,lon1)
                        write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                        write(radar_unit) ZDR, CC, FDP, KDP, SNR
                        dazim=azim-(i1+0.0)*radar_data%rgatesp(k)/(sfcrng*RADIAN)
                        call gcircle(lat0,lon0,dazim,sfcrng,lat1,lon1)
                        write(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
                        write(radar_unit) ZDR, CC, FDP, KDP, SNR
                     enddo
                  endif
               endif
            enddo
         enddo
         NLEV = 0
         WRITE(radar_unit) STID,lat1,lon1,TIM,NLEV,NFLAG
         close(radar_unit)
      endif

   enddo

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine copy_radar_data(rdat_in,rdat_out)
   implicit none

   type(t_radar_data), intent(in ) :: rdat_in
   type(t_radar_data), intent(out) :: rdat_out

   integer :: RGates,VGates,WGates,MaxRads,MaxCuts
   integer :: i, j, k, nrad

   MaxCuts=rdat_in%ntilt
   MaxRads=maxval(rdat_in%nazim           (1:MaxCuts))
   RGates =maxval(rdat_in%nrgate(1:MaxRads,1:MaxCuts))
   VGates =maxval(rdat_in%nvgate(1:MaxRads,1:MaxCuts))
   WGates =maxval(rdat_in%nvgate(1:MaxRads,1:MaxCuts))

   call allocate_radar_data(rdat_out,RGates,VGates,WGates,MaxRads,MaxCuts)

   rdat_out%radar_name  = rdat_in%radar_name
   rdat_out%radar_type  = rdat_in%radar_type
   rdat_out%file_format = rdat_in%file_format
   rdat_out%radar_id    = rdat_in%radar_id
   rdat_out%latitude    = rdat_in%latitude
   rdat_out%longitude   = rdat_in%longitude
   rdat_out%altitude    = rdat_in%altitude
   rdat_out%vcp         = rdat_in%vcp
   rdat_out%year        = rdat_in%year
   rdat_out%month       = rdat_in%month
   rdat_out%day         = rdat_in%day
   rdat_out%hour        = rdat_in%hour
   rdat_out%minute      = rdat_in%minute
   rdat_out%second      = rdat_in%second
   rdat_out%ntilt       = rdat_in%ntilt
   rdat_out%calibConst  = rdat_in%calibConst
   rdat_out%if_have_loc = rdat_in%if_have_loc

   nrad  =ubound(rdat_out%vel,dim=2)
   do k=1, rdat_out%ntilt
      if(.not.rdat_out%ifvel(k))then
         rdat_out%vmax(k)=VALUE_INVALID
      endif
      do j=rdat_out%nazim(k)+1,nrad
         rdat_out%razim(j,k)=VALUE_INVALID
         rdat_out%rtilt(j,k)=VALUE_INVALID
      enddo
   enddo
   do k=1,rdat_in%ntilt

      rdat_out%nazim           (k)=rdat_in%nazim           (k)
      rdat_out%ifref           (k)=rdat_in%ifref           (k)
      rdat_out%ifvel           (k)=rdat_in%ifvel           (k)
      rdat_out%ifzdr           (k)=rdat_in%ifzdr           (k)
      rdat_out%ifcc            (k)=rdat_in%ifcc            (k)
      rdat_out%ifkdp           (k)=rdat_in%ifkdp           (k)
      rdat_out%ifsnr           (k)=rdat_in%ifsnr           (k)
      rdat_out%rmax            (k)=rdat_in%rmax            (k)
      rdat_out%vmax            (k)=rdat_in%vmax            (k)
      rdat_out%rgatesp         (k)=rdat_in%rgatesp         (k)
      rdat_out%vgatesp         (k)=rdat_in%vgatesp         (k)
      rdat_out%atmosAttenFactor(k)=rdat_in%atmosAttenFactor(k)

      do j=1, rdat_in%nazim(k)
         rdat_out%vres  (j,k)=rdat_in%vres  (j,k)
         rdat_out%nrgate(j,k)=rdat_in%nrgate(j,k)
         rdat_out%nvgate(j,k)=rdat_in%nvgate(j,k)
         rdat_out%stime (j,k)=rdat_in%stime (j,k)
         rdat_out%etime (j,k)=rdat_in%etime (j,k)
         rdat_out%rtilt (j,k)=rdat_in%rtilt (j,k)
         rdat_out%razim (j,k)=rdat_in%razim (j,k)

         if(rdat_out%ifref(k))then
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%ref(i,j,k)=rdat_in%ref(i,j,k)
            enddo
         endif
         if(rdat_out%ifvel(k))then
            do i=1, rdat_in%nvgate(j,k)
               rdat_out%vel(i,j,k)=rdat_in%vel(i,j,k)
               rdat_out%spw(i,j,k)=rdat_in%spw(i,j,k)
            enddo
         endif

         if(rdat_out%ifzdr(k))then
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%zdr(i,j,k)=rdat_in%zdr(i,j,k)
            enddo
         endif
         if(rdat_out%ifcc (k))then
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%cc (i,j,k)=rdat_in%cc (i,j,k)
            enddo
         endif
         if(rdat_out%ifkdp(k))then
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%kdp(i,j,k)=rdat_in%kdp(i,j,k)
               rdat_out%fdp(i,j,k)=rdat_in%fdp(i,j,k)
            enddo
         endif
         if(rdat_out%ifsnr(k))then
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%snr(i,j,k)=rdat_in%snr(i,j,k)
            enddo
         endif
      enddo
   enddo

   if(rdat_out%if_have_loc)then
      do k=1,rdat_in%ntilt
         do j=1, rdat_in%nazim(k)
            do i=1, rdat_in%nrgate(j,k)
               rdat_out%rlat(i,j,k)=rdat_in%rlat(i,j,k)
               rdat_out%rlon(i,j,k)=rdat_in%rlon(i,j,k)
               rdat_out%ralt(i,j,k)=rdat_in%ralt(i,j,k)
            enddo
            do i=1, rdat_in%nvgate(j,k)
               rdat_out%vlat(i,j,k)=rdat_in%vlat(i,j,k)
               rdat_out%vlon(i,j,k)=rdat_in%vlon(i,j,k)
               rdat_out%valt(i,j,k)=rdat_in%valt(i,j,k)
            enddo
         enddo
      enddo
   endif

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calculate_radar_loc(radar_data)
   use libradar, only: beamhgt, gcircle
   implicit none
   type(t_radar_data), intent(inout) :: radar_data

   integer :: i, j, k, n
   real    :: lat1, lon1, alt, lat0, lon0, alt0, range, height, elev, sfcrng,azim



   if(radar_data%if_have_loc)then
      return
   endif

   lat0=radar_data%latitude
   lon0=radar_data%longitude
   alt0=radar_data%altitude

   do k=1, radar_data%ntilt
      if(radar_data%ifvel(k))then
         do j=1,radar_data%nazim(k)
            elev=radar_data%rtilt(j,k)
            azim=radar_data%razim(j,k)
            do i=1,radar_data%nvgate(j,k)
               if(radar_data%vel(i,j,k)/=value_invalid)then
                  range=i*radar_data%vgatesp(k)
                  call beamhgt(elev,range,height,sfcrng)
                  call gcircle(lat0,lon0,azim,sfcrng,lat1,lon1)
                  radar_data%vlon(i,j,k)=lon1
                  radar_data%vlat(i,j,k)=lat1
                  radar_data%valt(i,j,k)=height+alt0
               endif
            enddo
         enddo
      endif
   enddo

   do k=1, radar_data%ntilt
      if(radar_data%ifref(k))then
         do j=1,radar_data%nazim(k)
            elev=radar_data%rtilt(j,k)
            azim=radar_data%razim(j,k)
            do i=1,radar_data%nrgate(j,k)
               if(radar_data%ref(i,j,k)/=value_invalid)then
                  range=i*radar_data%rgatesp(k)
                  call beamhgt(elev,range,height,sfcrng)
                  call gcircle(lat0,lon0,azim,sfcrng,lat1,lon1)
                  radar_data%rlon(i,j,k)=lon1
                  radar_data%rlat(i,j,k)=lat1
                  radar_data%ralt(i,j,k)=height+alt0
               endif
            enddo
         enddo
      endif
   enddo

   end subroutine

end module
