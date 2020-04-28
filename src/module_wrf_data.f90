MODULE wrf_data
USE map_utils, only: PROJ_INFO
IMPLICIT NONE

   SAVE

   INTEGER,                           PARAMETER   :: datestrlen = 19

   type t_wrf
      character(len=1024) :: filename

      TYPE(PROJ_INFO)     :: PROJ

      INTEGER :: nx, ny, nz, nt, domain, nsoil
      REAL    :: ds
      INTEGER :: nct   ! current time index

      CHARACTER(LEN=19) , DIMENSION(:) , ALLOCATABLE :: time_list
      REAL, DIMENSION(:,:  ), ALLOCATABLE :: u10, v10, psfc, q2, t2, ter
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p, h, u, v, w, t, th, q, ux, vx, qrn
   end type

   type(t_wrf) :: wrf

! Usage: 
! 1. call get_wrf_data (filename, istatus, [optional]get_time)

! 2 for each radar: call da_get_radar_rv

CONTAINS

   !--------------------------------------------------------------------
   SUBROUTINE get_wrf_data(filename, istatus, get_time)
   use netcdf, only: nf90_noerr
   IMPLICIT NONE
! 
   character(len=*), intent(in)  :: filename, get_time
   INTEGER,          INTENT(out) :: istatus
   
   optional :: get_time

   !
   INTEGER                       :: ncid_in
   integer :: i

   ! open wrf file
   wrf%filename=filename
   call open_wrf_file(filename, ncid_in, istatus)
   if(istatus/=nf90_noerr)then
      return
   endif

   ! read dimension and map
   call read_wrf_info(ncid_in, istatus)
   if(istatus/=nf90_noerr)then
      return
   endif

   ! allocate wrf data
   if(.not.allocated(wrf%p))then
      call allocate_wrf_data()
   endif

   ! check if want another time
   wrf%nct=1
   if(present(get_time))then
      do i=1, wrf%nt
         wrf%nct=i
         if(wrf%time_list(i)==get_time)then
            exit
         endif
      enddo
   endif

   ! read wrf data   
   call read_wrf_data_netcdf(ncid_in, istatus)
   if(istatus/=0) return

   ! close wrf file
   call close_wrf_file(ncid_in, istatus)

   END SUBROUTINE
   !--------------------------------------------------------------------
   subroutine open_wrf_file(filename, ncid_in, istatus)
   use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr
   implicit none

   CHARACTER(len=*), INTENT(in)  :: filename
   INTEGER,          intent(out) :: ncid_in
   INTEGER,          INTENT(out) :: istatus
   !
      istatus = nf90_open(filename, nf90_nowrite, ncid_in)

      IF(istatus /= nf90_noerr)THEN
         WRITE(*,*) "Error open wrf file: ", TRIM(filename)
      ENDIF
      WRITE(*,*) "Opening wrf file: ", TRIM(filename)

   end subroutine

   !--------------------------------------------------------------------
   subroutine  read_wrf_info(ncid_in, istatus)
   use netcdf, only: nf90_get_att, nf90_inq_dimid, nf90_Inquire_Dimension, &
                     nf90_inq_varid, nf90_get_var, NF90_GLOBAL
   use map_utils, only: PROJ_LC, PROJ_PS, PROJ_MERC, MAP_SET, IJ_TO_LATLON
   implicit none
   INTEGER, intent(in)   :: ncid_in
   INTEGER, INTENT(out)  :: istatus
   !
   INTEGER               :: dimid, varid, attid, nxid, nyid, nzid, ntid, nsoilid
   integer               :: i, j
   integer, dimension(4) :: start, count

   REAL                  :: LAT1, LON1, KNOWNI, KNOWNJ, DX, STDLON, TRUELAT1,TRUELAT2
   INTEGER               :: map_proj
   REAL                  :: CEN_LAT, CEN_LON, DY, X, Y
   INTEGER               :: CODE
   TYPE(PROJ_INFO)       :: PROJ1

      ! Attributes

      istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'GRID_ID',  wrf%domain )

      ! Dimensions
      istatus=nf90_inq_dimid(ncid_in,"west_east",nxid)
      istatus=nf90_Inquire_Dimension(ncid_in,nxid,len=wrf%nx)
      istatus=nf90_inq_dimid(ncid_in,"south_north",nyid)
      istatus=nf90_Inquire_Dimension(ncid_in,nyid,len=wrf%ny)
      istatus=nf90_inq_dimid(ncid_in,"bottom_top",nzid)
      istatus=nf90_Inquire_Dimension(ncid_in,nzid,len=wrf%nz)
      istatus=nf90_inq_dimid(ncid_in,"Time",ntid)
      istatus=nf90_Inquire_Dimension(ncid_in,ntid,len=wrf%nt)

      istatus=nf90_inq_dimid(ncid_in,"soil_layers_stag",nsoilid)
      istatus=nf90_Inquire_Dimension(ncid_in,nsoilid,len=wrf%nsoil)

      WRITE(*,"(5(A,I5))") "Dimensions: nx=",wrf%nx,", ny=",wrf%ny,", nz=",wrf%nz,", nt=",wrf%nt,", nsoil=",wrf%nsoil

      ! Setup Variables
      ! map
      istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'MAP_PROJ',  map_proj)
      if ( map_proj /= 0 ) then
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'DX'                        , DX      )
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'DY'                        , DY      )
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'TRUELAT1'                  , TRUELAT1)
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'TRUELAT2'                  , TRUELAT2)
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'STAND_LON'                 , STDLON  )
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'CEN_LAT'                   , CEN_LAT )
         istatus = nf90_get_att(ncid_in,NF90_GLOBAL,'CEN_LON'                   , CEN_LON )
      endif

      SELECT CASE(map_proj)
      CASE(1)
         CODE=PROJ_LC
      CASE(2)
         CODE=PROJ_PS
      CASE(3)
         CODE=PROJ_MERC
      CASE DEFAULT
         STOP "UNKNOWN PROJECT"
      END SELECT

      LAT1=CEN_LAT
      LON1=CEN_LON
      KNOWNI=(wrf%NX+1)/2.
      KNOWNJ=(wrf%NY+1)/2.
      CALL MAP_SET(CODE,LAT1,LON1,KNOWNI,KNOWNJ,DX,wrf%nx,wrf%ny,STDLON,TRUELAT1,TRUELAT2,PROJ1)

      X=2-KNOWNI
      Y=2-KNOWNJ
      CALL IJ_TO_LATLON(PROJ1, X, Y, LAT1, LON1)

      CALL MAP_SET(CODE,LAT1,LON1,1.,1.,DX,wrf%nx,wrf%ny,STDLON,TRUELAT1,TRUELAT2,wrf%PROJ)

      if(allocated (wrf%time_list))then
         deallocate(wrf%time_list)
      endif
      allocate(wrf%time_list(wrf%nt))

      write(*,*) "read Times..."
      ! Times wrf_date_list, wrf_time_list, wrf_start_date, init_time, run_length, output_interval, output_perhour
      start=1
      count=(/datestrlen,wrf%nt,1,1/)
      istatus = nf90_inq_varid(ncid_in, "Times", varid)
      istatus = nf90_get_var(ncid_in, varid, wrf%time_list, start=start,count=count)

      wrf%ds      = dx
   end subroutine

   !--------------------------------------------------------------------
   subroutine close_wrf_file(ncid_in, istatus)
   use netcdf, only: nf90_close
   implicit none

      integer, intent(in)  :: ncid_in
      integer, intent(out) :: istatus
      istatus=nf90_close(ncid_in)
   end subroutine
   !--------------------------------------------------------------------

   subroutine allocate_wrf_data()
   implicit none

      allocate (wrf% p   ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% h   ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% t   ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% th  ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% q   ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% qrn ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% u   ( wrf%nx+1 , wrf%ny   , wrf%nz   )  )
      allocate (wrf% v   ( wrf%nx   , wrf%ny+1 , wrf%nz   )  )
      allocate (wrf% w   ( wrf%nx   , wrf%ny   , wrf%nz+1 )  )
      allocate (wrf% ux  ( wrf%nx   , wrf%ny   , wrf%nz   )  )
      allocate (wrf% vx  ( wrf%nx   , wrf%ny   , wrf%nz   )  )
                                                     
      allocate (wrf% u10 ( wrf%nx   , wrf%ny              )  )
      allocate (wrf% v10 ( wrf%nx   , wrf%ny              )  )
      allocate (wrf% psfc( wrf%nx   , wrf%ny              )  )
      allocate (wrf% t2  ( wrf%nx   , wrf%ny              )  )
      allocate (wrf% q2  ( wrf%nx   , wrf%ny              )  )
      allocate (wrf% ter ( wrf%nx   , wrf%ny              )  )

   end subroutine

   !--------------------------------------------------------------------

   SUBROUTINE read_wrf_data_netcdf(ncid_in, istatus)
   use netcdf, only: nf90_inq_varid, nf90_get_var
   IMPLICIT NONE
! 
   INTEGER,           INTENT(in)       :: ncid_in
   INTEGER,           INTENT(out)      :: istatus

   INTEGER, DIMENSION(4)               :: start, count
   REAL, ALLOCATABLE, DIMENSION(:,:)   :: data2d1, data2d2
   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: data3d1, data3d2
   INTEGER                             :: varid

   REAL , PARAMETER :: r_d     = 287.
   REAL , PARAMETER :: cp      = 7.*r_d/2.
   REAL , PARAMETER :: rcp     = r_d/cp
   REAL , PARAMETER :: deg2rad = 3.14159265/180.
   REAL , PARAMETER :: c2k     = 273.15

   !  H
   allocate( data3d1(wrf%nx,wrf%ny,wrf%nz+1) )
   allocate( data3d2(wrf%nx,wrf%ny,wrf%nz+1) )

   ! get base and perturbation height at full eta levels:
   start    = 1
   count    = (/wrf%nx,wrf%ny,wrf%nz+1,1/)
   start(4) = wrf%nct
   istatus  = nf90_inq_varid(ncid_in, "PH", varid)
   istatus  = nf90_get_var(ncid_in, varid, data3d1, start = start, count = count)
   istatus  = nf90_inq_varid(ncid_in, "PHB", varid)
   istatus  = nf90_get_var(ncid_in, varid, data3d2, start = start, count = count)

   ! compute Z at half levels:
   data3d1 = (data3d1+data3d2)/9.81
   wrf%h       = 0.5*(data3d1(:,:,1:wrf%nz)+data3d1(:,:,2:wrf%nz+1))

   deallocate(data3d1, data3d2)

   ! P
   allocate( data3d1(wrf%nx,wrf%ny,wrf%nz) )
   allocate( data3d2(wrf%nx,wrf%ny,wrf%nz) )

   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx,wrf%ny,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "P", varid)
   istatus  = nf90_get_var(ncid_in, varid, data3d1, start=start, count=count)
   istatus  = nf90_inq_varid(ncid_in, "PB", varid)
   istatus  = nf90_get_var(ncid_in, varid, data3d2, start=start, count=count)

   wrf%p=data3d1+data3d2
   deallocate(data3d1, data3d2)

   ! T (perturbation potential temperature (theta-t0)
   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx,wrf%ny,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "T", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%th, start = start, count = count)
   wrf%th=wrf%th+300.
   wrf%t        = wrf%th*(wrf%p/100000.)**rcp

   ! U
   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx+1,wrf%ny,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "U", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%u, start = start, count = count)

   wrf%ux=.5*(wrf%u(1:wrf%nx,:,:)+wrf%u(2:wrf%nx+1,:,:))

   ! V
   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx,wrf%ny+1,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "V", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%v, start = start, count = count)

   wrf%vx=.5*(wrf%v(:,1:wrf%ny,:)+wrf%v(:,2:wrf%ny+1,:))

   ! W 
   start    = 1
   count    = (/wrf%nx,wrf%ny,wrf%nz+1,1/)
   start(4) = wrf%nct
   istatus  = nf90_inq_varid(ncid_in, "W", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%w, start = start, count = count)

   ! QVAPOR
   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx,wrf%ny,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "QVAPOR", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%q, start    = start, count = count)

   ! QRAIN
   start    = 1
   start(4) = wrf%nct
   count    = (/wrf%nx,wrf%ny,wrf%nz,1/)
   istatus  = nf90_inq_varid(ncid_in, "QRAIN", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%qrn, start    = start, count = count)

   ! Q2
   start    = 1
   count    = (/wrf%nx,wrf%ny,1,1/)
   start(3) = wrf%nct
   istatus  = nf90_inq_varid(ncid_in, "Q2", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%q2, start = start, count = count)

   ! PSFC
   start    = 1
   start(3) = wrf%nct
   count    = (/wrf%nx,wrf%ny,1,1/)
   istatus  = nf90_inq_varid(ncid_in, "PSFC", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%psfc, start = start, count = count)

   ! T2
   start    = 1
   start(3) = wrf%nct
   count    = (/wrf%nx,wrf%ny,1,1/)
   istatus  = nf90_inq_varid(ncid_in, "T2", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%t2, start = start, count = count)

   ! U10
   start    = 1
   start(3) = wrf%nct
   count    = (/wrf%nx,wrf%ny,1,1/)
   istatus  = nf90_inq_varid(ncid_in, "U10", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%u10, start = start, count = count)

   ! V10
   start    = 1
   start(3) = wrf%nct
   count    = (/wrf%nx,wrf%ny,1,1/)
   istatus  = nf90_inq_varid(ncid_in, "V10", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%v10, start = start, count = count)

   ! HGT (Terrain Height)
   start    = 1
   count    = (/wrf%nx,wrf%ny,1,1/)
   start(3) = wrf%nct
   istatus  = nf90_inq_varid(ncid_in, "HGT", varid)
   istatus  = nf90_get_var(ncid_in, varid, wrf%ter, start = start, count = count)

   end subroutine


   subroutine get_wrf_sound(lat,lon,range,nzs,hs,us,vs) !,ps,ts,qs
   use map_utils, only: rotate_uv, latlon_to_ij
   use libradar, only: interp2d
   implicit none

   real, intent(in) :: lat, lon, range
   integer, intent(out) :: nzs
  
   real,dimension(:), allocatable, intent(out) :: hs,us,vs !,ps,ts,qs

   real, dimension(1)    :: x, y
   integer :: i, j, k, ndx, ndy, count

   if(.NOT.wrf%proj%init) return
   ndx=range/wrf%ds
   ndy=range/wrf%ds
   call latlon_to_ij(wrf%proj,lat,lon,x(1),y(1))

   nzs=wrf%nz
   if(allocated(hs))then
      deallocate(hs)
      !deallocate(ps)
      deallocate(us)
      deallocate(vs)
      !deallocate(ts)
      !deallocate(qs)
   endif
   allocate(hs(nzs))
   allocate(us(nzs))
   allocate(vs(nzs))
   do k=1, wrf%nz
      ! hs(k)=0.
      ! us(k)=0.
      ! vs(k)=0.
      ! count=0
      ! do i=int(max((x-ndx),1.)), int(min((x+ndx),float(nx)))
      !    do j=int(max((y-ndy),1.)), int(min((y+ndy),float(ny)))
      !       hs(k)=hs(k)+h(i,j,k)
      !       us(k)=us(k)+u(i,j,k)
      !       vs(k)=vs(k)+v(i,j,k)
      !       count=count+1
      !    enddo
      ! enddo
      ! hs(k)=hs(k)/count
      ! us(k)=us(k)/count
      ! vs(k)=vs(k)/count
      call INTERP2D(wrf%NX, wrf%NY, wrf% h(:,:,k), 1, X, Y, hs(k), -999.)
      call INTERP2D(wrf%NX, wrf%NY, wrf%ux(:,:,k), 1, X, Y, us(k), -999.)
      call INTERP2D(wrf%NX, wrf%NY, wrf%vx(:,:,k), 1, X, Y, vs(k), -999.)
      call rotate_uv(wrf%proj,lon,us(k),vs(k))
   enddo
   end subroutine

END MODULE
