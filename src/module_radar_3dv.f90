module radar_3dv
implicit none

integer :: maxrgate, maxazim, maxelev, maxvgate
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_radar_3dv(filename,rdat)
   use config,     only: min_valid_nyquist
   use radar_data, only: t_radar_data, radar_unit, get_maxmin_2d
   use radar_qc,   only: rngrvol, rngvvol, azmvol, elvvol, &
                         rxvvol, rxrvol, ryvvol, ryrvol, gjelim, &
                         iordunf, nzsnd, zsnd, usnd, vsnd, rfrsnd, &
                         refcheck,rngmin,velcheck,velmedl,rngmaxr, rngmaxv
   use libradar, only: beamrng, xy2rd, gcircle
   implicit none
   character(len=*),   intent(in) :: filename
   type(t_radar_data), intent(in) :: rdat

   real, parameter :: dgrid =3000. ! 3km grid
   real, parameter :: miss = -888888.

   real,    dimension(:),     allocatable :: xs, ys, elvavg
   real,    dimension(:,:,:), allocatable :: refg, velg, hgtg, spwg, rngg
   integer, dimension(:,:),   allocatable :: iwrite, ks, ke
   real,    dimension(:,:),   allocatable :: latg, long

   integer, dimension(:,:,:), allocatable :: refg_qc, velg_qc  ! -1:ground, -2:ref noise, -3:spw noise
   integer :: nx, ny, nz, nhalf, istatus, ntotal
   integer :: i, j, k, igate, jazim, kelev, kref, kvel, jstart, kr, kv
   integer :: ngate, nazim, nelev

   real :: azim, sfcrng, range, height, min_var, max_var
   real, dimension(:), allocatable :: vn

   character(len=80), parameter :: loc_fmt_char="(A5,2X,A12,2(F9.4,1X),F8.1,2X,A19,2I6)"
   character(len=80), parameter :: col_fmt_char="(A12,3X,A19,2X,2(F12.4,2X),F8.1,2X,I6)"
   !character(len=80), parameter :: dat_fmt_char="( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )"
   character(len=80), parameter :: dat_fmt_char="( 3X, F12.1, 3(F12.3,I4,F12.3,2X) )"

   character(len=19) :: chtime
   integer, dimension(3) :: dims

   integer :: qc_flag=2
   real    :: vel_err=4., ref_err=7.5
                      
   write(chtime,"(I4.4,'-',I2.2,'-',I2.2,'_',I2.2,':',I2.2,':',I2.2)")  rdat%year,rdat%month ,rdat%day   , &
                                 rdat%hour,rdat%minute,rdat%second

   dims     = ubound(rdat%ref)
   maxrgate = dims(1)
   maxazim  = dims(2)
   maxelev  = rdat%ntilt
   dims     = ubound(rdat%vel)
   maxvgate = dims(1)

   nz=0
   do k=1, maxelev
      if(rdat%ifvel(k))then
         nz=nz+1
      endif
   enddo 

   nx=maxval(rdat%rmax(1:maxelev))/dgrid*2+1
   ny=nx
   
   allocate(    xs(nx      ))
   allocate(    ys(   ny   ))
   allocate(  latg(nx,ny   ))
   allocate(  long(nx,ny   ))
   allocate(iwrite(nx,ny   ))
   allocate(    ks(nx,ny   ))
   allocate(    ke(nx,ny   ))
   allocate(  refg(nx,ny,nz))
   allocate(  velg(nx,ny,nz))
   allocate(  hgtg(nx,ny,nz))
   allocate(  spwg(nx,ny,nz))
   allocate(velg_qc(nx,ny,nz))
   allocate(refg_qc(nx,ny,nz))
   allocate(  rngg(nx,ny,nz))
   allocate(vn(nz))

   allocate(elvavg(maxelev))

   do k=1, maxelev
      elvavg(k)=sum(rdat%rtilt(1:rdat%nazim(k),k))/rdat%nazim(k)
   enddo

   nhalf=nx/2+1
   do i=1, nx
      xs(i)=(i-nhalf)*dgrid
   enddo
   do j=1, ny
      ys(j)=(j-nhalf)*dgrid
   enddo
   iwrite=0
  
   write(*,*) " Remap radar ref and vel..."
   refg_qc=-88
   velg_qc=-88

   write(*,*) "maxrgate,maxvgate,maxazim:", maxrgate, maxvgate, maxazim
   kref=0
   kvel=0
   do kelev=1, maxelev
      write(*,*) " Remaping Elevation:", kelev
      if(rdat%ifref(kelev))then
         kref=kref+1
         !write(701,*) "remap ref",kref,kelev
         !min_var=999.
         !max_var=-999.
         !do i=1, maxrgate
         !   do j=1, maxazim
         !      if(rdat%ref(i,j,kelev)>refcheck)then
         !         if(rdat%ref(i,j,kelev)<min_var)then
         !            !write(701,*) "min_var",i,j,kelev,min_var,rdat%ref(i,j,kelev)
         !            min_var=rdat%ref(i,j,kelev)
         !         endif
         !         if(rdat%ref(i,j,kelev)>max_var)then
         !            write(701,*) "max_var",i,j,kelev,max_var,rdat%ref(i,j,kelev)
         !            max_var=rdat%ref(i,j,kelev)
         !         endif
         !      endif
         !   enddo
         !enddo
         !write(*,*) "ref min,max:",min_var, max_var
         call get_maxmin_2d(rdat%ref(:,:,kelev),min_var,max_var)
         write(*,*) "ref min,max:",min_var, max_var
!       stop 
         call remap2d(maxrgate,maxazim,1,nx,ny,nzsnd,                                           &
                      refcheck,miss,velmedl,5.,iordunf,                                         &
                      rdat%nrgate(1:maxazim,kelev),rdat%nazim(kelev),1,                        &
                      rdat%latitude,rdat%longitude,0.,0.,rdat%altitude,1.,                      &
                      rngmin,rngmaxr(kelev),                                                    &
                      rngrvol(:,kelev),azmvol(:,kelev),elvvol(:,kelev),                         &
                      rdat%ref(1:maxrgate,1:maxazim,kelev),rxrvol(:,:,kelev),ryrvol(:,:,kelev), &
                      xs,ys,zsnd,rfrsnd,refg(1:nx,1:ny,kref),istatus)
         do j=1, ny
            do i=1, nx
               sfcrng=sqrt(xs(i)*xs(i)+ys(j)*ys(j))
               call beamrng(elvavg(kelev),sfcrng,height,range)
!write(88,*) range, (-80.+20.*log(range)/log(10.)), refg(i,j,kref)
               !if(refg(i,j,kref)<(-80.+20.*log(range)/log(10.))) refg(i,j,kref)=miss 
               hgtg(i,j,kref)=height+rdat%altitude
               rngg(i,j,kref)=range
            enddo
         enddo
      endif
      if(rdat%ifvel(kelev))then
         kvel=kvel+1
         !write(701,*) "remap vel",kvel,kelev
         !min_var=999.
         !max_var=-999.
         !do i=1, maxvgate
         !   do j=1, maxazim
         !      if(rdat%vel(i,j,kelev)>velcheck)then
         !         if(rdat%vel(i,j,kelev)<min_var)then
         !            write(701,*) "min_var",i,j,kelev,min_var,rdat%vel(i,j,kelev)
         !            min_var=rdat%vel(i,j,kelev)
         !         endif
         !         if(rdat%vel(i,j,kelev)>max_var)then
         !            write(701,*) "max_var",i,j,kelev,max_var,rdat%vel(i,j,kelev)
         !            max_var=rdat%vel(i,j,kelev)
         !         endif
         !      endif
         !   enddo
         !enddo
         !write(*,*) "vel min,max:",min_var, max_var
         call get_maxmin_2d(rdat%vel(:,:,kelev),min_var,max_var)
         write(*,*) "vel min,max:",min_var, max_var
         call remap2d(maxvgate,maxazim,1,nx,ny,nzsnd,                                           &
                      velcheck,miss,velmedl,5.,iordunf,                                         &
                      rdat%nvgate(1:maxazim,kelev),rdat%nazim(kelev),1,                        &
                      rdat%latitude,rdat%longitude,0.,0.,rdat%altitude,1.,                      &
                      rngmin,rngmaxv(kelev),                                                    &
                      rngvvol(:,kelev),azmvol(:,kelev),elvvol(:,kelev),                         &
                      rdat%vel(1:maxvgate,1:maxazim,kelev),rxvvol(:,:,kelev),ryvvol(:,:,kelev), &
                      xs,ys,zsnd,rfrsnd,velg(:,:,kvel),istatus)
         !write(701,*) "remap spw",kvel,kelev
         call remap2d(maxvgate,maxazim,1,nx,ny,nzsnd,                                           &
                      velcheck,miss,velmedl,5.,iordunf,                                         &
                      rdat%nvgate(1:maxazim,kelev),rdat%nazim(kelev),1,                        &
                      rdat%latitude,rdat%longitude,0.,0.,rdat%altitude,1.,                      &
                      rngmin,rngmaxv(kelev),                                                    &
                      rngvvol(:,kelev),azmvol(:,kelev),elvvol(:,kelev),                         &
                      rdat%spw(1:maxvgate,1:maxazim,kelev),rxvvol(:,:,kelev),ryvvol(:,:,kelev), &
                      xs,ys,zsnd,rfrsnd,spwg(:,:,kvel),istatus)
         vn(kvel)=rdat%vmax(kelev)

      endif
  
   enddo

   ks=0
   ke=1 
   do j=1, ny
      do i=1, nx
         do k=1, min(nz,7) ! delete 2 high elevation
            if((abs(velg(i,j,k))<200.and.abs(spwg(i,j,k))<200).or.abs(refg(i,j,k))<200)then
               refg_qc(i,j,k)=0
               velg_qc(i,j,k)=0
            
               ! remove clearsky echo.
               ! refg <20-0.004*(hgtg-radar_alt) 
               height=hgtg(i,j,k)-rdat%altitude
               if(refg(i,j,k)<(20-0.004*height))then
                  refg_qc(i,j,k)=-1
                  velg_qc(i,j,k)=-1
               endif
               ! ref noise
               !if(refg(i,j,k)<(-80.+20.*log(rngg(i,j,k))/log(10.)))then
               !   refg_qc(i,j,k)=-2
               !   velg_qc(i,j,k)=-2
               !endif

               ! spw noise
               if(spwg(i,j,k)>8..or.spwg(i,j,k)<1.5)then
               !  refg_qc(i,j,k)=-2
                  velg_qc(i,j,k)=-2
               endif

               ! range check
               if(rngg(i,j,k)>rdat%rmax(k))then
                  refg_qc(i,j,k)=-3
                  velg_qc(i,j,k)=-3
               endif

               ! height check
               if(height>15000)then
                  refg_qc(i,j,k)=-4
                  velg_qc(i,j,k)=-4
               endif

               ! interpolate error
               if(abs(velg(i,j,k))>200.)then
                  velg_qc(i,j,k)=-6
               endif

               ! nyquist vel
               if(abs(vn(k))<min_valid_nyquist)then
                  velg_qc(i,j,k)=-7
               endif

               if(ks(i,j)==0) ks(i,j)=k
               ke(i,j)=k
            else
               refg_qc(i,j,k)=-5
               velg_qc(i,j,k)=-5
            endif
         enddo
         if(ks(i,j)/=0)then
            call xy2rd(xs(i),ys(j),sfcrng,azim)
            call gcircle(rdat%latitude,rdat%longitude,azim,sfcrng,latg(i,j),long(i,j))
            iwrite(i,j)=1
         endif

         do k=1, nz
            if(abs(velg(i,j,k))>200..and.velg(i,j,k)/=miss)then
               write(602,*) i,j,k,velg(i,j,k),latg(i,j),long(i,j),rngg(i,j,k), hgtg(i,j,k)
            endif
            if(velg(i,j,k)/=miss.and.spwg(i,j,k)/=miss.and.refg(i,j,k)/=miss)then
               write(91,*) latg(i,j),long(i,j), hgtg(i,j,k), velg(i,j,k), vn(k), refg(i,j,k), spwg(i,j,k)
            elseif(velg(i,j,k)/=miss.or.refg(i,j,k)/=miss)then
               write(92,*) latg(i,j),long(i,j), hgtg(i,j,k), velg(i,j,k), vn(k), refg(i,j,k), spwg(i,j,k)
            endif
         enddo

      enddo
   enddo

   ntotal=sum(iwrite)
   if(ntotal>0)then
      OPEN(radar_unit,FILE=filename,STATUS='unknown',FORM='formatted')

      write(*,*) "Writing radar_3dv file:", trim(filename)
      write(radar_unit,'(A12,I3)') "Total Radar:",1
      write(radar_unit,'(A80)')"================================================================================"
      write(radar_unit,'(A80)') loc_fmt_char
      write(radar_unit,fmt=loc_fmt_char) "RADAR",rdat%radar_id,rdat%longitude,rdat%latitude,rdat%altitude,chtime,ntotal,maxval(ke-ks+1)
      write(radar_unit,'(A80)') col_fmt_char
      write(radar_unit,'(A80)') dat_fmt_char
      DO i=1, nx 
         do j=1, ny
            if(iwrite(i,j)>0)then
               write(radar_unit,fmt=col_fmt_char) "FM-128      ",chtime,latg(i,j),long(i,j),hgtg(i,j,1),ke(i,j)-ks(i,j)+1
               DO k=ks(i,j),ke(i,j)
                  write(radar_unit,fmt=dat_fmt_char) hgtg(i,j,k),velg(i,j,k),velg_qc(i,j,k), vel_err, &
                                                                 refg(i,j,k),refg_qc(i,j,k), ref_err, &
                                                                 spwg(i,j,k),k, vn(k)

               ENDDO
            endif
         enddo
      enddo

      CLOSE(radar_unit)
   endif

   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE remap2d(maxgate,maxazim,maxelev,nx,ny,nzsnd,                 &
                      varchek,varmiss,vmedlim,dazlim,iorder,              &
                      kntgate,kntazim,kntelev,                            &
                      rdrlat,rdrlon,radarx,radary,rdralt,dazim,           &
                      rngmin,rngmax,                                      &
                      rngvol,azmvol,elvvol,                               &
                      varvol,rxvol,ryvol,                                 &
                      xs,ys,zsnd,rfrsnd,gridvar,istatus)
   use libradar, only: beamhgt
   use radar_qc, only: gjelim
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take data from the lowest tilt and remap it onto a plane for 
!  purposes of display in comparison to the ARPS k=2 reflectivity.
!  Uses a least squares a local quadratic fit to remap the data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS
!          June, 2002
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    varid    Radar variable ID for diagnostic file writing
!    varname  Name of radar variable for diagnostic file writing
!    varunit  Units of radar variable for diagnostic file writing
!
!    varchek  Threshold for checking data, good vs. flagged
!    varmiss  Value to assign to data for missing
!    vmedlim  Threshold limit for median check
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
!
!    rxvol    x-coordinate at radar data location
!    ryvol    y-coordinate at radar data location
!
!    gridvar  Radar variable remapped to scalar points
!
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
  INTEGER, INTENT(IN) :: maxgate
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nzsnd

  REAL, INTENT(IN)    :: varchek
  REAL, INTENT(IN)    :: varmiss
  REAL, INTENT(IN)    :: vmedlim
  REAL, INTENT(IN)    :: dazlim
  INTEGER, INTENT(IN) :: iorder

  INTEGER, INTENT(IN) :: kntgate(maxazim,maxelev)
  INTEGER, INTENT(IN) :: kntazim(maxelev)
  INTEGER, INTENT(IN) :: kntelev

  REAL, INTENT(IN)    :: rdrlat
  REAL, INTENT(IN)    :: rdrlon
  REAL, INTENT(IN)    :: radarx
  REAL, INTENT(IN)    :: radary
  REAL, INTENT(IN)    :: rdralt
  REAL, INTENT(IN)    :: dazim
  REAL, INTENT(IN)    :: rngmin
  REAL, INTENT(IN)    :: rngmax

  REAL, INTENT(IN)    :: rngvol(maxgate,maxelev)
  REAL, INTENT(IN)    :: azmvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: varvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)   :: rxvol(maxgate,maxazim,maxelev)
  REAL, INTENT(OUT)   :: ryvol(maxgate,maxazim,maxelev)

  REAL, INTENT(IN)    :: xs(nx)
  REAL, INTENT(IN)    :: ys(ny)
  REAL, INTENT(IN)    :: zsnd(nzsnd)
  REAL, INTENT(IN)    :: rfrsnd(nzsnd)
  REAL, INTENT(OUT)   :: gridvar(nx,ny)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
! Misc. Local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: n = 6
  REAL :: avar(n,n)
  REAL :: rhsvar(n)
  REAL :: avel(n,n)
  REAL :: rhsvel(n)
  REAL :: sol(n)
  REAL :: work(n,n+1)
  REAL :: work1d(n+1)

  REAL :: array(3,3)
  REAL :: rhsv(3)
  REAL :: solv(3)

  REAL, PARAMETER :: eps = 1.0E-25

  INTEGER :: ii,jj,kk,i,j,k,knt,kinbox
  INTEGER :: kok,isort,jsort,mid
  INTEGER :: kbgn,kend
  INTEGER :: igate,jazim,jazmin,jmirror,jend
  INTEGER :: istatal,istatwrt

  INTEGER, PARAMETER :: maxsort = 5000

  REAL :: deg2rad,rad2deg
  REAL :: delx,dely,delz,dazimr,daz,azdiff
  REAL :: ddx,ddxy,ddx2,ddy,ddy2,dxthr,dxthr0
  REAL :: azmrot,xcomp,ycomp,mapfct,sfcr,zagl
  REAL :: sum,sum2,sdev,thresh,slrange,elijk,azimijk,time
  REAL :: varmax,varmin,varavg,varmean,varmed

  REAL, ALLOCATABLE :: varsort(:)
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=80  ) :: runname  ! Name of this run
  INTEGER :: lfnkey
  CHARACTER (LEN=80 ) :: dirname  ! The name of output directory

  kk=1
  deg2rad = atan(1.)/45.
  rad2deg = 1./deg2rad
  dazimr = dazim*deg2rad
  dxthr0=0.6*max((xs(3)-xs(2)),(ys(3)-ys(2)))

  allocate(varsort(maxsort),stat=istatal)

  time=0.
  delz=0.
!
  DO j=1,ny
    DO i=1,nx
      gridvar(i,j)=varmiss
    END DO
  END DO
!
  write(*,*) kntazim(kk), kntgate(kntazim(kk),kk)
  DO jazim=1,kntazim(kk)
    xcomp=sin(deg2rad*azmvol(jazim,kk))
    ycomp=cos(deg2rad*azmvol(jazim,kk))

    DO igate=1,kntgate(jazim,kk)
      CALL beamhgt(elvvol(jazim,kk),rngvol(igate,kk),zagl,sfcr)
             
      rxvol(igate,jazim,kk)=radarx+xcomp*sfcr
      ryvol(igate,jazim,kk)=radary+ycomp*sfcr
    END DO
  END DO
!
  DO j=1,ny-1
    if(ys(j)> rngmax) cycle
    DO i=1,nx-1
      kok=0
      sum=0.
      sum2=0.
      sdev=0.
      varavg=999999.
      varsort=999999.
      delx=xs(i)-radarx
      dely=ys(j)-radary
      slrange=sqrt(delx*delx+dely*dely+delz*delz)
      dxthr=max(dxthr0,((slrange+dxthr0)*dazimr))
      if(slrange> rngmax) cycle
!
      IF( slrange > 0. ) THEN
!
!-----------------------------------------------------------------------
!
! Determine azimuth to this grid cell
!
!-----------------------------------------------------------------------
!
        IF(delx == 0.) THEN
          IF(dely >= 0.) THEN
            azimijk=0.
          ELSE
            azimijk=180.
          END IF
        ELSE
          azimijk=rad2deg*atan(delx/dely)
          IF(dely < 0.) azimijk=azimijk+180.
          IF(azimijk < 0.) azimijk=azimijk+360.
        END IF
!
        varmax=-999.
        varmin=999.
        DO jj=1,n
          DO ii=1,n
            avar(ii,jj)=0.
          END DO
        END DO 
!
        DO ii=1,n
          rhsvar(ii)=0.
        END DO
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
        azdiff=181.
        jazmin=1
        DO jazim=1,kntazim(kk)
          daz=azmvol(jazim,kk)-azimijk
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          daz=abs(daz)
          IF(daz < azdiff) THEN
            azdiff=daz
            jazmin=jazim
          END IF
        END DO

        jmirror=jazmin+(kntazim(kk)/2)
        IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  First pass, find median, avg, std dev.
!
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
        jend=kntazim(kk)
        IF(jmirror > jazmin) jend=jmirror-1
        
        DO jazim=jazmin,jend
          kinbox=0
          daz=azmvol(jazim,kk)-azimijk
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)
!
            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
              kinbox=kinbox+1
              !write(*,*) igate,jazim,kk
              IF(abs(varvol(igate,jazim,kk)) < abs(varchek) ) THEN
                DO isort=1,kok
                  IF(varvol(igate,jazim,kk) < varsort(isort)) EXIT
                END DO
                IF(kok < maxsort) THEN
                  DO jsort=kok,isort,-1
                    varsort(jsort+1)=varsort(jsort)
                  END DO
                  varsort(isort)=varvol(igate,jazim,kk)
                  sum=sum+varvol(igate,jazim,kk)
                  sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                  kok=kok+1
                ELSE
                  EXIT
                END IF
              END IF   ! data ok
            END IF  ! inside box
          END DO ! igate
          IF(kinbox == 0 .OR. kok == maxsort) EXIT
        END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
        IF(kinbox > 0 .and. jend==kntazim(kk)) THEN
          DO jazim=1,jmirror-1
            kinbox=0
            daz=azmvol(jazim,kk)-azimijk
            IF(daz > 180.) daz=daz-360.
            IF(daz < -180.) daz=daz+360.
            IF(abs(daz) > dazlim) EXIT
            DO igate=1,kntgate(jazim,kk)
!
              ddx=rxvol(igate,jazim,kk)-xs(i)
              ddy=ryvol(igate,jazim,kk)-ys(j)

              IF( rngvol(igate,kk) > rngmin .AND.                    &
                  rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                IF(abs(varvol(igate,jazim,kk)) < abs(varchek) ) THEN
                  IF(kok < maxsort) THEN
                    DO isort=1,kok
                      IF(varvol(igate,jazim,kk) < varsort(isort)) EXIT
                    END DO
                    DO jsort=kok,isort,-1
                      varsort(jsort+1)=varsort(jsort)
                    END DO
                    varsort(isort)=varvol(igate,jazim,kk)
                    sum=sum+varvol(igate,jazim,kk)
                    sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                    kok=kok+1
                  ELSE
                    EXIT
                  END IF
                END IF
              END IF
            END DO
            IF(kinbox == 0 .OR. kok == maxsort) EXIT
          END DO
        END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
        jend= 1
        IF(jmirror < jazmin) jend=jmirror
        DO jazim=jazmin-1,jend,-1
          kinbox=0
          daz=azmvol(jazim,kk)-azimijk
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
!
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)

            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

              kinbox=kinbox+1

              IF(abs(varvol(igate,jazim,kk)) < abs(varchek) ) THEN
                IF(kok < maxsort) THEN
                  DO isort=1,kok
                    IF(varvol(igate,jazim,kk) < varsort(isort)) EXIT
                  END DO
                  DO jsort=kok,isort,-1
                    varsort(jsort+1)=varsort(jsort)
                  END DO
                  varsort(isort)=varvol(igate,jazim,kk)
                  sum=sum+varvol(igate,jazim,kk)
                  sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                  kok=kok+1
                ELSE
                  EXIT
                END IF
              END IF
            END IF
          END DO
          IF(kinbox == 0 .OR. kok == maxsort) EXIT
        END DO
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
        IF(kinbox > 0 .and. jend==1 ) THEN
          DO jazim=kntazim(kk),jmirror,-1
            kinbox=0
            daz=azmvol(jazim,kk)-azimijk
            IF(daz > 180.) daz=daz-360.
            IF(daz < -180.) daz=daz+360.
            IF(abs(daz) > dazlim) EXIT
            DO igate=1,kntgate(jazim,kk)
!
              ddx=rxvol(igate,jazim,kk)-xs(i)
              ddy=ryvol(igate,jazim,kk)-ys(j)
! 
              IF( rngvol(igate,kk) > rngmin .AND.                    &
                  rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN
                kinbox=kinbox+1
                IF(abs(varvol(igate,jazim,kk)) < abs(varchek) ) THEN
                  IF(kok < maxsort ) THEN
                    DO isort=1,kok
                      IF(varvol(igate,jazim,kk) < varsort(isort)) EXIT
                    END DO
                    DO jsort=kok,isort,-1
                      varsort(jsort+1)=varsort(jsort)
                    END DO
                    varsort(isort)=varvol(igate,jazim,kk)
                    sum=sum+varvol(igate,jazim,kk)
                    sum2=sum2+(varvol(igate,jazim,kk)*varvol(igate,jazim,kk))
                    kok=kok+1
                  ELSE
                    EXIT
                  END IF
                END IF
              END IF
            END DO  ! igate
            IF(kinbox == 0 .OR. kok == maxsort) EXIT
          END DO ! jazim
        END IF
!
        mid=(kok/2)+1
        varmed=varsort(mid)
        IF(kok > 0) varavg=sum/float(kok)
        IF ( kok > 1 )                                               &
          sdev=sqrt((sum2-(sum*sum/float(kok)))/float(kok-1))
        thresh=max((2.*sdev),vmedlim)
!
!-----------------------------------------------------------------------
!
!  Process data for local quadratic fit
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Find nearest azimuth at this level
!
!-----------------------------------------------------------------------
!
        azdiff=181.
        jazmin=1
        DO jazim=1,kntazim(kk)
          daz=azmvol(jazim,kk)-azimijk
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          daz=abs(daz)
          IF(daz < azdiff) THEN
            azdiff=daz
            jazmin=jazim
          END IF
        END DO
             
        jmirror=jazmin+(kntazim(kk)/2)
        IF(jmirror > kntazim(kk)) jmirror=jmirror-kntazim(kk)
!
!-----------------------------------------------------------------------
!
!  Loop forward from jazmin
!
!-----------------------------------------------------------------------
!
        jend=kntazim(kk)
        IF(jmirror > jazmin) jend=jmirror-1
        DO jazim=jazmin,jend
          kinbox=0
          daz=azmvol(jazim,kk)-azimijk
          IF(daz > 180.) daz=daz-360.
          IF(daz < -180.) daz=daz+360.
          IF(abs(daz) > dazlim) EXIT
          DO igate=1,kntgate(jazim,kk)
!
            ddx=rxvol(igate,jazim,kk)-xs(i)
            ddy=ryvol(igate,jazim,kk)-ys(j)
! 
            IF( rngvol(igate,kk) > rngmin .AND.                    &
                rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

              kinbox=kinbox+1
              ddxy=ddx*ddy
              ddx2=ddx*ddx
              ddy2=ddy*ddy

              IF(abs(varvol(igate,jazim,kk)) < abs(varchek) .AND.  &
                 abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 if(varvol(igate,jazim,kk)>varmax)then
                    !write(701,*)"varmax",igate,jazim,kk,varvol(igate,jazim,kk),varmax
                 endif
                 varmax=max(varmax,varvol(igate,jazim,kk))
                 if(varvol(igate,jazim,kk)<varmin)then
                    !write(701,*)"varmin",igate,jazim,kk,varvol(igate,jazim,kk),varmin
                 endif
                 varmin=min(varmin,varvol(igate,jazim,kk))
                 
!
                 rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                 rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                 rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                 rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                 rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                 rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
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
                 avar(3,4)=avar(3,4)+ddy*ddx2
                 avar(3,5)=avar(3,5)+ddy*ddx2
                 avar(3,6)=avar(3,6)+ddy*ddy2
!
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
               END IF
!
             END IF
           END DO  ! igate
           IF(kinbox == 0) EXIT
         END DO ! jazim
!
!-----------------------------------------------------------------------
!
!  IF kinbox > 0 continue from jazim=1
!
!-----------------------------------------------------------------------
!
         IF(kinbox > 0 .and. jend==kntazim(kk)) THEN
           DO jazim=1,jmirror-1
             kinbox=0
             daz=azmvol(jazim,kk)-azimijk
             IF(daz > 180.) daz=daz-360.
             IF(daz < -180.) daz=daz+360.
             IF(abs(daz) > dazlim) EXIT
             DO igate=1,kntgate(jazim,kk)
!
               ddx=rxvol(igate,jazim,kk)-xs(i)
               ddy=ryvol(igate,jazim,kk)-ys(j)
!
               IF( rngvol(igate,kk) > rngmin .AND.                    &
                      rngvol(igate,kk) < rngmax .AND.                    &
                      abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

                 kinbox=kinbox+1
                 ddxy=ddx*ddy
                 ddx2=ddx*ddx
                 ddy2=ddy*ddy

                 IF(abs(varvol(igate,jazim,kk)) < abs(varchek) .AND.             &
                    abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 if(varvol(igate,jazim,kk)>varmax)then
                    !write(701,*)"varmax",igate,jazim,kk,varvol(igate,jazim,kk),varmax
                 endif
                   varmax=max(varmax,varvol(igate,jazim,kk))
                 if(varvol(igate,jazim,kk)<varmin)then
                    !write(701,*)"varmin",igate,jazim,kk,varvol(igate,jazim,kk),varmin
                 endif
                   varmin=min(varmin,varvol(igate,jazim,kk))
!
                   rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                   rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                   rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                   rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                   rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                   rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
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
                   avar(5,1)=avar(5,1)+ddxy
                   avar(5,2)=avar(5,2)+ddxy*ddx
                   avar(5,3)=avar(5,3)+ddxy*ddy
                   avar(5,4)=avar(5,4)+ddxy*ddxy
                   avar(5,5)=avar(5,5)+ddxy*ddx2
                   avar(5,6)=avar(5,6)+ddxy*ddy2

                   avar(6,1)=avar(6,1)+ddx2
                   avar(6,2)=avar(6,2)+ddx2*ddx
                   avar(6,3)=avar(6,3)+ddx2*ddy
                   avar(6,4)=avar(6,4)+ddx2*ddxy
                   avar(6,5)=avar(6,5)+ddx2*ddx2
                   avar(6,6)=avar(6,6)+ddx2*ddy2
!
                 END IF
!
               END IF
             END DO  ! igate
             IF(kinbox == 0) EXIT
           END DO ! jazim
         END IF
!
!-----------------------------------------------------------------------
!
! Loop backward from jazmin
!
!-----------------------------------------------------------------------
!
         jend= 1
         IF(jmirror < jazmin) jend=jmirror
         DO jazim=jazmin-1,jend,-1
           kinbox=0
           daz=azmvol(jazim,kk)-azimijk
           IF(daz > 180.) daz=daz-360.
           IF(daz < -180.) daz=daz+360.
           IF(abs(daz) > dazlim) EXIT
           DO igate=1,kntgate(jazim,kk)
!
             ddx=rxvol(igate,jazim,kk)-xs(i)
             ddy=ryvol(igate,jazim,kk)-ys(j)
!
             IF( rngvol(igate,kk) > rngmin .AND.                    &
                 rngvol(igate,kk) < rngmax .AND.                    &
                 abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

               kinbox=kinbox+1
               ddxy=ddx*ddy
               ddx2=ddx*ddx
               ddy2=ddy*ddy

               IF(abs(varvol(igate,jazim,kk)) < abs(varchek) .AND.             &
                  abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 if(varvol(igate,jazim,kk)>varmax)then
                    !write(701,*)"varmax",igate,jazim,kk,varvol(igate,jazim,kk),varmax
                 endif
                 varmax=max(varmax,varvol(igate,jazim,kk))
                 if(varvol(igate,jazim,kk)<varmin)then
                    !write(701,*)"varmin",igate,jazim,kk,varvol(igate,jazim,kk),varmin
                 endif
                 varmin=min(varmin,varvol(igate,jazim,kk))
!
                 rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                 rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                 rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                 rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                 rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                 rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
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
               END IF
!
             END IF
           END DO  ! igate
           IF(kinbox == 0) EXIT
         END DO ! jazim
!
!-----------------------------------------------------------------------
!
! If not yet outside box, continue from last radial.
!
!-----------------------------------------------------------------------
!
         IF(kinbox > 0 .and. jend==1 ) THEN
         DO jazim=kntazim(kk),jmirror,-1
           kinbox=0
           daz=azmvol(jazim,kk)-azimijk
           IF(daz > 180.) daz=daz-360.
           IF(daz < -180.) daz=daz+360.
           IF(abs(daz) > dazlim) EXIT
           DO igate=1,kntgate(jazim,kk)
!
             ddx=rxvol(igate,jazim,kk)-xs(i)
             ddy=ryvol(igate,jazim,kk)-ys(j)
!
             IF( rngvol(igate,kk) > rngmin .AND.                    &
                 rngvol(igate,kk) < rngmax .AND.                    &
                 abs(ddx) < dxthr .AND. abs(ddy) < dxthr ) THEN

               kinbox=kinbox+1
               ddxy=ddx*ddy
               ddx2=ddx*ddx
               ddy2=ddy*ddy

               IF(abs(varvol(igate,jazim,kk)) < abs(varchek) .AND.             &
                  abs(varvol(igate,jazim,kk)-varmed) < thresh ) THEN
!
                 if(varvol(igate,jazim,kk)>varmax)then
                    !write(701,*)"varmax",igate,jazim,kk,varvol(igate,jazim,kk),varmax
                 endif
                  varmax=max(varmax,varvol(igate,jazim,kk))
                 if(varvol(igate,jazim,kk)<varmin)then
                    !write(701,*)"varmin",igate,jazim,kk,varvol(igate,jazim,kk),varmin
                 endif
                  varmin=min(varmin,varvol(igate,jazim,kk))
!
                  rhsvar(1)=rhsvar(1)+varvol(igate,jazim,kk)
                  rhsvar(2)=rhsvar(2)+varvol(igate,jazim,kk)*ddx
                  rhsvar(3)=rhsvar(3)+varvol(igate,jazim,kk)*ddy
                  rhsvar(4)=rhsvar(4)+varvol(igate,jazim,kk)*ddxy
                  rhsvar(5)=rhsvar(5)+varvol(igate,jazim,kk)*ddx2
                  rhsvar(6)=rhsvar(6)+varvol(igate,jazim,kk)*ddy2
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
                  avar(5,1)=avar(5,1)+ddx2
                  avar(5,2)=avar(5,2)+ddx2*ddx
                  avar(5,3)=avar(5,3)+ddx2*ddy
                  avar(5,4)=avar(5,4)+ddx2*ddxy
                  avar(5,5)=avar(5,5)+ddx2*ddx2
                  avar(5,6)=avar(5,6)+ddx2*ddy2
!
                  avar(6,1)=avar(6,1)+ddx2
                  avar(6,2)=avar(6,2)+ddx2*ddx
                  avar(6,3)=avar(6,3)+ddx2*ddy
                  avar(6,4)=avar(6,4)+ddx2*ddxy
                  avar(6,5)=avar(6,5)+ddx2*ddx2
                  avar(6,6)=avar(6,6)+ddx2*ddy2
!
                END IF
 
              END IF
            END DO  ! igate
            IF(kinbox == 0) EXIT
          END DO ! jazim
        END IF
!
!-----------------------------------------------------------------------
!
!   Solve for variable at grid point
!
!-----------------------------------------------------------------------
!
        knt=nint(avar(1,1))
        if(varmin>=999.)then
            gridvar(i,j)=varmiss
        else 
        IF ( iorder > 1 .and. knt > 7 ) THEN          
          varmean=rhsvar(1)/avar(1,1)
          CALL GJELIM(n,avar,rhsvar,sol,work,work1d,eps,istatus)
          gridvar(i,j)=min(varmax,max(varmin,sol(1)))
          !if(abs(gridvar(i,j))>100)then
             !write(701,*) "remap2d",avar,rhsvar,varmin,varmax,solv(1)
             !stop
          !endif
        ELSE IF ( iorder > 0 .and. knt > 5 ) THEN
          DO jj=1,3
            DO ii=1,3
              array(ii,jj)=avar(ii,jj)
            END DO
          END DO
          DO ii=1,3
            rhsv(ii)=rhsvar(ii)
          END DO
          CALL GJELIM(3,array,rhsv,solv,work,work1d,eps,istatus)
          gridvar(i,j)=min(varmax,max(varmin,solv(1)))
          !if(abs(gridvar(i,j))>100)then
             !write(701,*) "remap2d",array,rhsv,varmin,varmax,solv(1)
             !stop
          !endif
        ELSE IF ( knt > 0 ) THEN
          !varmean=rhsvar(1)/avar(1,1)
          !gridvar(i,j)=varmean
          gridvar(i,j)=varmiss
        END IF
        endif

      END IF

    END DO   ! i loop
  END DO   ! j loop

  deallocate(varsort,stat=istatal)

  RETURN
  END SUBROUTINE remap2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module 
