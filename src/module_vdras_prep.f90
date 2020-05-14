module vdras_prep
implicit none

contains

   subroutine fill2d(dz,extrap,igrid, &
                     noct,rmnpts,mx,my,mz, &
                     nx,ny,nz,bad,ifill,pex)
      
   !  THIS ROUTINE PERFORMS 2-DIMENSIONAL DATA FILLING OF A CARTESIAN GRID
   !  USING A CONSTRAINED LEAST-SQUARES DATA FILL
   !
   !  DZ  -INPUT 3-DIMENSIONAL DATA SET TO BE FILLED.  ON OUTPUT DATA
   !         CONTAINS THE FILLED DATA.
   !  IGRID-  MAXIMUM NUMBER OF GRID POINTS TO SEARCH OUTWARD FROM A
   !         MISSING DATA LOCATION TO DETERMINE IF IT IS BOUNDED AND IF
   !         RMNPTS HAS BEEN SATISFIED.
   !  RMNPTS-MINIMUM NUMBER OF SURROUNDING POINTS REQUIRED
   !  NOCT-  MINIMUM NUMBER OF QUADRANTS OCCUPIED IN ORDER TO SATISFY THE SEARCH.
   !  NX-    NUMBER OF X GRID POINTS
   !  NY-    NUMBER OF Y GRID POINTS
   !  NZ-    NUMBER OF Z GRID POINTS
   !  BAD-   MISSING DATA FLAG
      
      implicit none
   
      integer, intent(in) :: igrid,noct,ifill, &
                             mx,my,mz, &
                             nx,ny,nz
      real, dimension(mx,my,mz), intent(inout) :: dz
      real, dimension(mx,my,mz), intent(out)   :: extrap
   
      real, intent(in) :: rmnpts, bad, pex
   
      real, dimension(mx,my) :: temp
   !yzm   integer, dimension(8) :: ioct(8)
      integer, dimension(8) :: ioct
      real, parameter :: eps=0.0001
   
      integer, dimension(mz) :: nfill, nexfl
   
      integer :: i, j, k, i0, j0, i1, j1, i2, j2, n, ix, iy, k0, k1, kq
      real    :: rnpts,sx,sx2,sy,sy2,t1,t2,t3,denom,rnum,sxy,sdx,sdy,sd
   
      print*
      print *,'PERFORM 2D FILLING'
      
      nfill(:)=0
      nexfl(:)=0
      
   !!$omp parallel do default(shared) &
   !!$omp private(k,i,j)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if((dz(i,j,k) - bad) < 1.0) then
            extrap(i,j,k)=0.
         else
            extrap(i,j,k)=1.
         end if
      enddo
      enddo
      enddo
   
   !!$omp parallel do default(shared) &
   !!$omp private(k,temp) &
   !!$omp private(i, j, i0, j0, i1, j1, i2, j2, ix, iy, nfill, nexfl, k0, k1, kq) &
   !!$omp private(rnpts,sx,sx2,sy,sy2,t1,t2,t3,denom,rnum,sxy,sdx,sdy,sd)
      k_loop: do k=1,nz
         do j=1,ny
         do i=1,nx
            temp(i,j)=dz(i,j,k)
         enddo
         enddo
      
      
      
   ! LOOP OVER GRID POINTS. IF GRID POINT VALUE IS MISSING FILL IT,
   ! OTHERWISE DO NOTHING.
      
        j0_loop: do j0=1,ny
          i0_loop: do i0=1,nx
            if(((temp(i0,j0) - bad) > 1.0).and.ifill.eq.0) cycle i0_loop
   
            ioct(1:4) = 0
   
   ! START SEARCHING OUTWARD FROM GRID POINT UNTIL THE NOCT AND RMNPTS
   ! CONDITIONS ARE MET.
            n_loop: do n=1,igrid
              j1=max0(1,j0-n)
              j2=min0(ny,j0+n)
              i1=max0(1,i0-n)
              i2=min0(nx,i0+n)
   
   ! INITIALIZE SUMMATION TERMS USED IN THE LEAST-SQUARES FIT
   
              rnpts=0.
              sx=0.
              sy=0.
              sx2=0.
              sy2=0.
              sxy=0.
              sdx=0.
              sdy=0.
              sd=0.
   
              j_loop: do j=j1,j2
                 iy=j-j0
                 i_loop: do i=i1,i2
                    if((temp(i,j) - bad) < 1.0) cycle i_loop
   
                    ix=i-i0
                    if(ix.ge.0.and.iy.gt.0.)ioct(1)=1
                    if(ix.gt.0.and.iy.le.0.)ioct(2)=1
                    if(ix.le.0.and.iy.lt.0.)ioct(3)=1
                    if(ix.lt.0.and.iy.ge.0.)ioct(4)=1
                    rnpts=rnpts+1.0
                    sx=sx+ix
                    sy=sy+iy
                    sx2=sx2+ix*ix
                    sy2=sy2+iy*iy
                    sxy=sxy+ix*iy
                    sd=sd+temp(i,j)
                    sdx=sdx+ix*temp(i,j)
                    sdy=sdy+iy*temp(i,j)
                 enddo i_loop
              enddo j_loop
   
              kq=sum(ioct(1:4))
   
              if(kq.lt.noct) cycle n_loop
              if(rnpts.lt.rmnpts) cycle n_loop
   
              nfill(k)=nfill(k)+1
              t1 = sx2*sy2-sxy*sxy
              t2 = sx*sy2-sxy*sy
              t3 = sx*sxy-sx2*sy
              denom=rnpts*t1-sx*t2+sy*t3
   
              if(denom.le.eps)then
                 write(unit=*, fmt='(1X,A,3I6,2X,E12.4)') &
                       'DENOM LT EPS-I,J,K,DENOM=',i0,j0,k0,denom
   !yzm          dz(i0,j0,k)=bad
                 cycle n_loop
              end if
   
              rnum=sd*t1-sdx*t2+sdy*t3
              dz(i0,j0,k)=rnum/denom
              if(kq.le.2) then
                nexfl(k)=nexfl(k)+1
                extrap(i0,j0,k)=pex
              else
                extrap(i0,j0,k)=1.
              end if
              exit n_loop
            enddo n_loop
   
         enddo i0_loop
         enddo j0_loop
      
      enddo k_loop
   
      print *,'NFILL=',sum(nfill(:)),' NEXFL=',sum(nexfl(:))
      
   end subroutine fill2d


!
! This routine has been modified a little - radar_dir1 is now
! not actually a directory, it is now the URL where data are to
! be read from. Niles Oien November 2001.
!
!
   subroutine readradarn(radar_dir1,yyyy,mm,dd,hh,mn,ss,unixtime, &
                         dt0,assimilation_window,rlond,rlatd,badpt, &
                         igridv,iedit,iproj_type,ecurv,pex,irealtime,vrthr,dtheta, &
                         uvad1,vvad1,num_rng,num_ang,ivad,iut, &
                         nrad, nxrad, nyrad, nzrad, nzd, nze, nze_top, nobs, max_times, &
                         ngrdp,nquad,rmnpts,intpbl,mpt, &
                         bdytop,nelev_vad, imdv_nx, imdv_ny, imdv_nz,irad,ivolfound, &
                         raddataz, raddatavr, raddatasd, xk,yl,zm,rlevel,rlonr,rlatr,raltr, &
                         zppi, pu, pz, pr, urs, cs, rs, std, &
                         zcart, urscart, cscart, pucart, pzcart, &
                         nyds, nyde, nxds, nxde, nzds, nzde, &
                         nyms, nyme, nxms, nxme, nzms, nzme, &
                         nyts, nyte, nxts, nxte, nzts, nzte, &
                         dx, dy, dz, x00, y00, z00, iunfold,uu,vv,topo0)
   !yzm                  dx, dy, dz, x00, y00, z00)
   
!  use dmp_util_module

   implicit none

   integer, intent(inout) :: iproj_type

   integer, intent(in) :: irealtime, igridv, iedit, &
                          nrad, nxrad, nyrad, nzrad, nze, nobs, max_times, &
                          num_rng,num_ang,ivad, ngrdp,nquad,intpbl,mpt, &
                          nelev_vad, irad
   integer, intent(inout) :: imdv_nx, imdv_ny, imdv_nz, nze_top

   integer, intent(in) :: nyds, nyde, nxds, nxde, nzds, nzde, &
                          nyms, nyme, nxms, nxme, nzms, nzme, &
                          nyts, nyte, nxts, nxte, nzts, nzte, &
                          nzd

   integer, intent(in) :: unixtime,yyyy,mm,dd,hh,mn,ss,assimilation_window,ecurv

   integer, dimension(nobs,nrad), intent(inout) :: iut
   integer, dimension(nrad),      intent(inout) :: ivolfound
   real, dimension(nze,nrad),     intent(out) :: rlevel           !yzm 

   real, dimension(nrad), intent(inout) :: xk, yl, zm

   real, dimension(nxrad, nyrad, nzrad), intent(inout) :: raddataz, raddatavr, raddatasd

!yzm   real, dimension(nzds:nzde,  nrad), intent(out) :: uvad1, vvad1
   real, dimension(num_rng,  nrad), intent(out) :: uvad1, vvad1

   real, dimension(nyms:nyme, nxms:nxme, nzms:nzme, nobs, nrad), intent(inout) :: urscart, cscart, pucart, pzcart
   real, dimension(nzms:nzme), intent(inout) :: zcart

   real, dimension(nyms:nyme, nxms:nxme, 1:nze, 1:nobs, 1:nrad), intent(inout) ::  &
                                               urs, pu, pz, pr, cs, rs, std
   real, dimension(nyms:nyme, nxms:nxme, 1:nze, 1:3, 1:nrad), intent(inout) :: zppi

   real, intent(in)  :: dt0, pex, vrthr, dtheta, bdytop, rmnpts
   real, intent(in)  :: rlond,rlatd,badpt

   real, intent(in)  :: dx, dy, dz, x00, y00, z00

   integer :: icheck_radar
   integer :: nfields
   integer :: have_dbz, have_vel
   character*8 :: dbz_field, vel_field
   
   real, dimension(nyms:nyme, nxms:nxme), intent(in) :: topo0

   real, intent(inout)  :: rlonr,rlatr,raltr

   character*64  proc_message
   character*128 radar_dir1,file1,filename1,file0
!
!     MAX_TIMES is the maximum number of radar scans to look at in the
!     assimilation window - Niles.
!
   integer*4 data_times
   dimension data_times(max_times)
   integer*4 current_data_time
   
   integer iyyd,immd,iddd,ihhd,imnd,issd
   character*128 filename

   real, dimension(nxms:nxme, nyms:nyme, nzms:nzme) :: workm
   real, dimension(nxms:nxme, nyms:nyme, 1:nze)     :: workd

   real, dimension(nrad) :: xloc, yloc, zloc
   
   real, dimension(nzrad) :: vlevel

   real    :: xminr,yminr,zminr
   real    :: dxkm,dykm,dzkm,dxkmr,dykmr,dzkmr
   real    :: el

   integer :: vlevelnz, iutcent, iucount

   integer :: i, j, k, l, m, il, ivol, ifound

   integer, intent(in) :: iunfold
   integer :: nscan, nscan_old                 !unfold
   real, dimension(nxds:nxde, nyds:nyde, nzds:nzde), intent(in) :: uu, vv

!----------------------------------------------------------
   integer, dimension(nxrad, nyrad) :: gmask        !yzm
   integer :: m0, imask
   real    :: badpt1              !yzm
!------------------------------------------

   integer :: nxp1, nyp1, nxp2, nyp2, nz

   integer :: nytsm1, nytem1, nxtsm1, nxtem1, &
              nytsp1, nytep1, nxtsp1, nxtep1, &
              nytsm2, nytem2, nxtsm2, nxtem2, &
              nytsp2, nytep2, nxtsp2, nxtep2

!   call get_ext_index(nytsm1, nytem1, nxtsm1, nxtem1, &
!                      nytsp1, nytep1, nxtsp1, nxtep1, &
!                      nytsm2, nytem2, nxtsm2, nxtem2, &
!                      nytsp2, nytep2, nxtsp2, nxtep2)

! Function to find end of strings.
!
   ifound = 0
   data_times = 0            !yzm
   print*,'in readradar, unixtime=',unixtime,max_times, ifound, data_times
   print*,'radar_dir1',radar_dir1
   print*,'unixtime',unixtime,assimilation_window,max_times, ifound

!   if(on_monitor) then
!      call find_mdv_in_interval(radar_dir1, unixtime, &
!        assimilation_window, max_times, ifound, data_times)
!   endif
!   
!   #ifdef DM_PARALLEL
!   call vdras_broadcast(ifound)
!   call vdras_broadcast(data_times, max_times)
!   #endif

   print*,'readradar data_times',(data_times(i),i=1,max_times)
!  file0 = radar_dir1(20:32)
!  print*,'file0',file0
!  open(unit=35,file=file0, form='formatted')
!   read(35,*)unixtime,assimilation_window,max_times, ifound
!   read(35,*)(data_times(i),i=1,max_times)
!  close(35)
   
   ivolfound(irad)=ifound
   print*,'Radar volumes found: ',ifound
   print*
   
   if (ifound .eq. 0) then
   
      print*, ' no data within ',assimilation_window
      print*, ' seconds of latest radar data file.'
      print*, ' at mdv url : '
      print*, 'radar_dir1=', trim(radar_dir1)
      print*
   
   else
   
!
! Otherwise we found something, sing about it ....
!
   print*, 'URL : radar_dir1=', trim(radar_dir1)
   print*
   
   do il=1, ifound
   
      iut(il,irad) = data_times(il)
   
!      call convert_from_unixtime(iyyd,immd,iddd,ihhd, &
!           imnd,issd,data_times(il))
   
      write (filename,45) iyyd,immd,iddd,ihhd,imnd,issd
   
      print*,' FILE : ',trim(filename)
   
   end do
   print*
   
 45   FORMAT(I4.4,I2.2,I2.2,'/',I2.2,I2.2,I2.2,'.mdv')
   
!
! --  Read MDV data
!
   
   do ivol = 1, ifound
   
   CURRENT_DATA_TIME = DATA_TIMES( IVOL )
   
   do m=1,nzrad
      do k=1,nyrad
         do l=1,nxrad
            raddataz(l,k,m)=badpt
            raddatavr(l,k,m)=badpt
         end do
      end do
   end do
   
!  proc_message  = 'Reading MDV file'//char(0)
!  call procmap_force_register( proc_message )
   
!   call convert_from_unixtime(iyyd,immd,iddd,ihhd, &
!        imnd,issd,current_data_time)
   
   write (filename,45) iyyd,immd,iddd,ihhd,imnd,issd
   
   print*,'READING FILE ',trim(filename)
   print*
   
!if(on_monitor) then

!  icheck_radar = 0 
   icheck_radar = 1 

   have_dbz = 0
   have_vel = 0
   nfields  = 0

   if(icheck_radar == 1) then

   dbz_field = 'DBZ'
   vel_field = 'VEL'

!   call check_mdv_fields(radar_dir1,current_data_time,0,dbz_field,vel_field,have_dbz,have_vel)

   if (have_dbz == 0) then
      print*, 'DBZ field included in this data.'
   else
      print*, 'DBZ field missing in this data.'
   endif

   if (have_vel == 0) then
      print*, 'VEL field included in this data.'
   else
      print*, 'VEL field missing in this data.'
   endif
   print*

   endif      !icheck_radar


   if (have_dbz .EQ. 0) then
!     call read_mdv(radar_dir1, current_data_time, 0, 0, &
!      call read_mdv(radar_dir1, current_data_time, 0, nfields, &
!           raddataz,imdv_nx,imdv_ny,imdv_nz, &
!           rlonr,rlatr,raltr,xminr,yminr,zminr,dxkmr,dykmr,dzkmr,badpt, &
!           iutcent,vlevel,vlevelnz,nscan)
!!yzm       iutcent,vlevel,vlevelnz)

       nfields = nfields + 1

   endif

!
!     Write out KIWX reflectivity in ascii format.
!
   
!$$$      if (irad.eq.1.and.IVOL.EQ.IFOUND) then
!$$$
!$$$         WRITE(FILENAME,FMT='(a,I10)') 'KIWX_REF-',IUT(IVOL,IRAD)
!$$$         open(unit=32,file=FILENAME,form='formatted',status='new')
!$$$         write(32,FMT='(f7.1)') (((RADDATAZ(k,l,m),k=1,nxrad),
!$$$     $        l=1,nyrad),m=1,nzrad)
!$$$         close(32)
!$$$
!$$$      endif
   
!
! Niles : Added check on array dimensions. Only for this
! case of read_mdv - one check is enough.
!
!
   if (imdv_nx .ne. nxrad .or. &
       imdv_ny .ne. nyrad) then
      print*,'NXRAD = ',nxrad,' NYRAD = ',nyrad
      print*,'imdv_nx = ',imdv_nx,' imdv_ny = ',imdv_ny
      print*, 'ERROR : Array dimensions do not match MDV input'
      stop
   endif
   
   if (imdv_nz .ne. nzrad) then
      print*
      write (6,fmt='(a,I3,a,I3)') ' WARNING: NZRAD = ',nzrad, &
           ' imdv_nz = ',imdv_nz
      print*
   endif
   
   if (have_vel .EQ. 0) then
!!     call read_mdv(radar_dir1, current_data_time, 0, 1, &
!      call read_mdv(radar_dir1, current_data_time, 0, nfields, &
!           raddatavr,imdv_nx,imdv_ny,imdv_nz, &
!           rlonr,rlatr,raltr,xminr,yminr,zminr,dxkmr,dykmr,dzkmr,badpt, &
!           iutcent,vlevel,vlevelnz,nscan)
!!yzm       iutcent,vlevel,vlevelnz)
   endif

!yzm---add mask-----under terrain is missing -8888.88 ---

   imask = 0
!  imask = 1

if(imask == 1) then

   badpt1 = -8888.88

! open(unit=73,file='/data/zhuming/data/taiwan/Meiyu/mask/mask.dat',form='formatted',status='old')
  open(unit=73, file='./mask.dat',form='unformatted',status='old')

    read(73) gmask
    write(6,'(i4,45i2)') (gmask(l,98), l=1,460)
    write(6,'(//)')
    write(6,'(i4,45i2)') (gmask(l,460), l=1,460)

!  do k = 1, nyrad
!   read(73,'(460i2)',end=92) (gmask(l,k), l=1,nxrad)
!   write(6,'(i4,45i2)') k,(gmask(l,k), l=1,45 )
!  enddo
!92 continue

  close(73)

   do k=1,nyrad
   do l=1,nxrad

      m0 = gmask(l,k)
      if(m0 >= nzrad+1) m0 = nzrad+1

      if(m0 .gt. 1) then

       do m = 1, m0-1
         raddatavr(l,k,m) = badpt1
         raddataz(l,k,m) = badpt1
       enddo
      endif

   enddo
   enddo

endif
!yzm----------------------------------------

!endif     !on_monitor


!   if (ifound .gt. 0 ) then
!      file1  = radar_dir1(25:32)
!      print*,'RADAR_DIR1',radar_dir1

!      print*, 'FILENAME', filename
!          write(filename1,fmt='(a8,a,a19)') file1,'/',filename
!          print*,'FILENAME1',filename1
!   open(unit=32,file=filename1,form='formatted')

!      read(32,fmt='5(i20)') imdv_nx,imdv_ny,imdv_nz, &
!                          iutcent, vlevelnz
!       
!      read(32,fmt='14(f8.2)')(vlevel(k),k=1,vlevelnz)
!         read(32,fmt='9(f12.5)') &
!        rlonr,rlatr,raltr,xminr,yminr,zminr,dxkmr,dykmr,dzkmr
!      read(32,fmt='10(f12.5)') (((raddataz(k,l,m),k=1,nxrad), &
!          l=1,nyrad),m=1,nzrad)
!         read(32,fmt='10(f12.5)') (((raddatavr(k,l,m),k=1,nxrad), &
!           l=1,nyrad),m=1,nzrad)
!       read(32,fmt='10(f12.5)') (((raddatasd(k,l,m),k=1,nxrad), &
!         l=1,nyrad),m=1,nzrad)
!      read(32,fmt='10(f12.5)') (((raddatasdz(k,l,m),k=1,nxrad), &
!            l=1,nyrad),m=1,nzrad)
!      close(32)
!   endif
!  
!-------------------------------------------------------------------
   IF (iunfold == 1) THEN
     if (nscan /= 11 .and. nscan /= 12 .and. vlevelnz == 14) then
        nscan_old = nscan
        if (vlevel(2) < 1.0) then
          nscan = 12
        else
          nscan = 11
        endif
        write(*,'(a,5i5)')'nscan is reset 1:',     &
               nscan_old,nscan,vlevelnz,ivol,irad
      elseif (nscan /= 31 .and. nscan /= 32 .and. vlevelnz == 5) then
        nscan_old = nscan
        nscan = 32
        write(*,'(a,5i5)')'nscan is reset 2:',   &
               nscan_old,nscan,vlevelnz,ivol,irad
      elseif (nscan /= 21 .and. vlevelnz == 9) then
        nscan_old = nscan
        nscan = 21
        write(*,'(a,5i5)')'nscan is reset 3:',   &
               nscan_old,nscan,vlevelnz,ivol,irad
      endif
      if (nscan /=  31 .and. nscan /= 32 .and. nscan /=  21 .and.   &
          nscan /=  11 .and. nscan /=  12) then
        write(*,'(a,4i5)') 'There is mismatch between scan_type and vlevelnz',  &
            nscan,vlevelnz,ivol,irad
        stop
      endif
    ENDIF          !iunfold
!---------------------------------------------------------------------
   
   xloc(irad) = (rlonr - rlond)*111.17*1000.0* &
           cos(3.14159*rlatd/180.)
   yloc(irad) = (rlatr - rlatd)*111.17*1000.0
   zloc(irad) = 0.0
   
   xminr=xminr+xloc(irad)/1000.
   yminr=yminr+yloc(irad)/1000.
   
   print*,'elevation angles:'
   do i=1,vlevelnz
    rlevel(i,irad) = vlevel(i)
     write (6,fmt='(I2,3X,f5.2)') i,vlevel(i)
   enddo
   print*
   
   print*,'radar altitude(km):','raltr=',raltr
   
!     FIND GRID INDICES OF RADAR POSITIONS
   

   xk(irad)=(xloc(irad)-x00)/dx
   yl(irad)=(yloc(irad)-y00)/dy
   zm(irad)=(zloc(irad)-z00)/dz
   
   if(ivol.eq.1) then
   print*
   write(6,114) irad,xloc(irad),yloc(irad),zloc(irad)
   write(6,119) irad,xk(irad),yl(irad),zm(irad)
   
   write(6,211) 'input data nx,ny,nz(km):',nxrad,nyrad,nzrad
   write(6,212) 'loc. of data grid (1,1,1(asl)):',xminr,yminr,zminr
   write(6,212) 'radar lat., lon., alt.:',rlatr,rlonr,raltr
   write(6,212) 'grid dx,dy,dz(km):',dxkmr,dykmr,dzkmr
   print*
   endif
   
 114  format('    radar ',I1,' loc. x,y,z:',3f12.6)
 119  FORMAT(' GRID INDEX OF RADAR ',I1,':',3F12.6)
 211  format(1x,a32,3i5)
 212  format(1x,a32,3f12.5)

!
! --- Fill data void, interpolate data to model grid and edit
!yzm-- move up 
!!$omp parallel do default(shared) & 
!!$omp private(k,l,m)
   do m=1,nzrad
   do k=1,nyrad
   do l=1,nxrad
      if((raddatavr(l,k,m)).gt.50.) then
         raddatavr(l,k,m)=badpt
      end if
      if((raddatavr(l,k,m)) .eq. badpt ) then
         raddataz(l,k,m)=badpt
      endif
!     if((raddataz(l,k,m)).gt.100) then
      if((raddataz(l,k,m)).gt.80.) then
         raddataz(l,k,m)=badpt
      end if
      
   enddo
   enddo
   enddo


   call pre_analyze(dxkmr,dykmr,dzkmr,xminr,yminr, &
                    dx,dy,dz,vrthr,dtheta, &
                    vlevel,yl,xk,zm,raltr, &
                    igridv,iedit,ivol, &
                    vlevelnz,iproj_type,ecurv,pex,irealtime, &
                    raddataz, raddatavr,  &
                    zppi, pu, pz, pr, cs, urs, rs, std, &
                    zcart, cscart, urscart, pucart, pzcart, &
                    nyds, nyde, nxds, nxde, nzds, nzde, &
                    nyms, nyme, nxms, nxme, nzms, nzme, &
                    nyts, nyte, nxts, nxte, nzts, nzte, &
                    nytsm1, nytep1, nxtsm1, nxtep1, &
                    nxrad,nyrad,nzrad, &
                    nrad,nobs,nzd,nze,nze_top, &
                    ngrdp,nquad,rmnpts,intpbl,mpt,irad, &
                    x00,y00,z00,badpt,iunfold,nscan,uu,vv, &
                    nelev_vad,num_ang,num_rng)        !yzm for unfolding
!
! --- vad analysis
!
   if(ivol.eq.1) then
      print*
      print*,'*****************************'
      print*,'** PERFORMING VAD ANALYSIS **'
      print*,'*****************************'
      print*
   
      if(iproj_type.eq.1)then
         el = vlevel(nelev_vad)
      endif
   
!     proc_message  = 'Perform VAD analysis'//char(0)
!     call procmap_auto_register( proc_message )
   
      dxkm=dx*0.001
      dykm=dy*0.001
      dzkm=dz*0.001
   
      if(ivad.eq.3) then
!!$omp parallel do default(shared) & 
!!$omp private(k,l,m)
         do m=1,nze
         do l=nytsm1,nytep1
         do k=nxtsm1,nxtep1
!           workd(k+1,l+1,m)=urs(l,k,m,ivol,irad)
            workd(k,l,m)=urs(l,k,m,ivol,irad)
         enddo
         enddo
         enddo
   
         call vad_anal(el,num_rng,num_ang, &
              uvad1,vvad1,nelev_vad,workd, &
              nrad, nxms, nxme, nyms, nyme, 1, nze, &
              dxkm,dykm,dzkm, &
              xk(irad)+1.0,yl(irad)+1.0,zm, &
              dz,z00,badpt,iproj_type,ecurv,irad)

      else if(ivad.eq.2) then
!!$omp parallel do default(shared) &
!!$omp private(k,l,m)
         do m=nzts, nzte
         do l=nytsm1,nytep1
         do k=nxtsm1,nxtep1
!           workm(k+1,l+1,m)=urscart(l,k,m,ivol,irad)
            workm(k,l,m)=urscart(l,k,m,ivol,irad)
         enddo
         enddo
         enddo
   
         iproj_type=2
         call vad_anal(el,num_rng,num_ang, &
              uvad1,vvad1,nelev_vad,workm, &
              nrad, nxms, nxme, nyms, nyme, nzms, nzme, &
              dxkm,dykm,dzkm, &
              xk(irad)+1.0,yl(irad)+1.0,zm, &
              dz,z00,badpt,iproj_type,ecurv,irad)
         iproj_type=1
   
      else
   
         call vad_anal(el,num_rng,num_ang, &
              uvad1,vvad1,nelev_vad,raddatavr, &
              nrad, 1, nxrad, 1, nyrad, 1, nzrad, &
              dxkmr,dykmr,dzkmr, &
              float(nxrad+1)/2.0,float(nyrad+1)/2.0,zm, &
              dz,z00,badpt,iproj_type,ecurv,irad)
   
      end if

   
   end if                    ! IF(IVOL.EQ.1)
   

!     REMOVE DATA ABOVE BDYTOP AND
!     CALCULATE THE NUMBER OF REMAINING GOOD DATA POINTS.
!
   
   iucount = 0
!!$omp parallel do default(shared) & 
!!$omp private(k,l,m)
      do m = 1, nze
      do k = nxms, nxme
      do l = nyms, nyme
   
         if(zppi(l,k,m,2,irad).gt.bdytop) then
!yzm        urs(l,k,m,ivol,irad)=badpt
!           CS(L,K,M,IVOL,IRAD)=BADPT                       ! leh
!           PU(L,K,M,IVOL,IRAD)=0.                          ! leh
!           PZ(L,K,M,IVOL,IRAD)=0.                          ! leh
         end if
!        if(urs(l,k,m,ivol,irad).ne.badpt) then
         if((urs(l,k,m,ivol,irad)-badpt) > 1.0) then
            iucount=iucount+1
         endif
   
      enddo
      enddo
      enddo


!     call vdras_sum(iucount)
   
      write(6,fmt='(a,I3,a,I3,a,I5)') &
           ' NUMBER OF GOOD VELOCITIES FROM RADAR ',irad, &
           ' VOLUME ',ivol,': ',iucount
      print*
   
!yzm added for setting value of badpt1=-888.88 on the points is badpt=-999.99 over sea 
!!$omp parallel do default(shared) & 
!!$omp private(k,l,m)
      do m = 1, nze
      do k = nxts, nxte
      do l = nyts, nyte
!     if((urs(l,k,m,ivol,irad)-badpt).lt.1.0 .and. topo0(l,k)<300.) then
!        urs(l,k,m,ivol,irad)=badpt1    
!     end if
!     if((cs(l,k,m,ivol,irad)-badpt)<1.0 .and. abs(topo0(l,k)-badpt)>1.0 .and. topo0(l,k)<300.) then
!     if((cs(l,k,m,ivol,irad)-badpt)<1.0 .and. abs(topo0(l,k)-badpt)>1.0 .and. topo0(l,k)<1.) then
!     if((cs(l,k,m,ivol,irad)-badpt)<1.0 .and. abs(topo0(l,k)-badpt)>1.0 .and. topo0(l,k)<10.) then
!        cs(l,k,m,ivol,irad)=badpt1      
!     end if
   enddo
   enddo
   enddo


   enddo
   
   endif                     
   

   end subroutine readradarn
   
   subroutine unfold_by_gradient_vr(RADDATAVR, &
                                    NX,NY,NZ,BADPT,vel_nyq)

!--------------------------------------------------------------
! Velocity unfolding program based on the horizontal gradient of vr
! Actually, it removes the folded area
!
! Eunha Lim                     April 2007
!--------------------------------------------------------------

   implicit none

   integer, intent(IN)         :: nx,ny,nz
   real,    intent(IN)         :: badpt

   real, dimension(nx,ny,nz), intent(INOUT)   :: raddatavr
   real, dimension(nz),       intent(IN)   :: vel_nyq

   integer, dimension(nx,ny)   :: mask

   real, parameter             :: eps = 1.0
   integer   :: i,j,k,kfold,count
   real      :: fold,grd



   do k = 1, nz
     count = 0
    mask = 0
   do j = 2, ny
loop1:   do i = 2, nx
    if (raddatavr(i,j,k) <= badpt .OR. &
       ABS(raddatavr(i,j,k)) <= vel_nyq(k)*0.9) CYCLE loop1

!--- Check the horizontal gradient of vr
     grd = ABS(raddatavr(i,j,k)+raddatavr(i,j,k-1))

    if (grd <= eps) mask(i,j) = 1
     if (mask(i-1,j) == 1 .AND.  mask(i,j) == 0) mask(i,j) = 1
     if (mask(i-1,j) == 1 .AND.  mask(i,j) == 1) mask(i,j) = 2

   enddo loop1
   enddo 

   do j = 1,ny
   do i = 1,nx
    if (mask(i,j) == 0) CYCLE
    raddatavr(i,j,k) = badpt
    count = count + 1
   enddo
   enddo

    write(*,'(a,2i10)')'Number of deleted data by unfolding',k,count
   enddo

   end subroutine unfold_by_gradient_vr

   subroutine unfold_preserve_hgradient(RADDATAVR, &
         NX,NY,NZ,BADPT,vel_nyq,mask)

!--------------------------------------------------------------
! Preserve small scale outflow by storms
!
! Eunha Lim                     April 2008
!--------------------------------------------------------------

   implicit none

   integer, intent(IN)         :: nx,ny,nz
   real,    intent(IN)         :: badpt

   real, dimension(nx,ny,nz), intent(IN)   :: raddatavr
   real, dimension(nz),       intent(IN)   :: vel_nyq

   integer, dimension(nx,ny,nz),intent(INOUT)   :: mask

   real, parameter             :: factor = 0.3
   integer, parameter          :: xgrid = 20, &       ! number of grid points
                          ygrid = 15

   integer   :: i,j,k,n1,n2,n



   mask = 0

   do 100 k = 1, nz

!--------------------------------------------------------------------
! Step1: Check the change of the sign between two adjecent grids
!--------------------------------------------------------------------

     do 1 j = 1, ny
       do i = 1, nx-1
        if (raddatavr(i,j,k) == badpt .OR. raddatavr(i+1,j,k) == badpt) CYCLE
!       if (k == 5) print*,'i,j,k=',i,j,k,raddatavr(i,j,k),raddatavr(i+1,j,k)
         if (raddatavr(i,j,k)*raddatavr(i+1,j,k) <= 0.)  mask(i,j,k) = 1
       enddo

       do i = nx, 2
        if (raddatavr(i,j,k) == badpt .OR. raddatavr(i-1,j,k) == badpt) CYCLE
         if (raddatavr(i,j,k)*raddatavr(i-1,j,k) <= 0.) mask(i,j,k) = 1
       enddo
1    continue

     do 2 i = 1, nx
       do j = 1, ny-1 
        if (raddatavr(i,j,k) == badpt .OR. raddatavr(i,j+1,k) == badpt) CYCLE
         if (raddatavr(i,j,k)*raddatavr(i,j+1,k) <= 0.) mask(i,j,k) = 1
       enddo

       do j = ny, 2 
        if (raddatavr(i,j,k) == badpt .OR. raddatavr(i,j-1,k) == badpt) CYCLE
         if (raddatavr(i,j,k)*raddatavr(i,j-1,k) <= 0.) mask(i,j,k) = 1
       enddo
2    continue


!--------------------------------------------------------------------
! Step2: Discriminate boundary of the out/inflow by storms 
!        from aliased boundary
!--------------------------------------------------------------------

     do 11 j = 1, ny
     do 11 i = 1, nx
      if (mask(i,j,k) == 1 .AND. ABS(raddatavr(i,j,k)) < vel_nyq(k)*factor) &
       mask(i,j,k) = mask(i,j,k) + 1
11   continue

!     do 12 j = 1, ny
!      n1 = 0
!     do 12 i = 1, nx
!      if (n1 == 0 .AND. mask(i,j,k) == 2) n1 = MIN(i+1,nx)
!      if (n1 /= 0 .AND. mask(i,j,k) == 0) then
!       if(k==5) print*,'j,n1,i',j,n1,i
!       mask(n1:i,j,k) = 0
!       n1 = 0
!      endif
!12   continue
!
!     do 13 i = 1, nx
!     do 13 j = 1, ny
!      n1 = MAX(j,1)
!      n2 = MIN(j+1,nx)
!      if (mask(i,n1,k)*mask(i,n2,k) >= 4) mask(i,j,k) = 3
!13   continue


!--------------------------------------------------------------------
! Step3: Check whether the strom boundary area is confined
!--------------------------------------------------------------------

! X- direction
     do 21 j = 1, ny
      n1 = 0
      n2 = 0
      do 21 i = 1, nx
         if (raddatavr(i,j,k) == badpt) CYCLE
       if (mask(i,j,k) == 2 .AND. n1 == 0) then
         n1 = MAX(i,1)
       elseif (mask(i,j,k) == 2 .AND. n1 /= 0) then
         n2 = MIN(i,nx)
       endif

         if (n1 == 0) CYCLE

       if (n2 > n1 .AND. (n2-n1) <= xgrid) then
!         if (k==5) print*,'j,n1,n2',j,n1,n2
         do n = n1, n2
           mask(n,j,k) = 3
           enddo
         n1 = n2
         n2 = 0
       elseif (n2 > n1 .AND. (n2-n1) > xgrid) then
          n1 = 0
          n2 = 0
       endif
21   continue

! Y-direction
     do 22 i = 1, nx
      n1 = 0
      n2 = 0
      do 22 j = 1, ny
         if (raddatavr(i,j,k) == badpt) CYCLE
       if (mask(i,j,k) == 3 .AND. n1 == 0) then
         n1 = MAX(j,1)
       elseif (mask(i,j,k) == 3 .AND. n1 /= 0) then
         n2 = MIN(j,ny)
       endif

         if (n1 == 0) CYCLE

       if (n2 > n1 .AND. (n2-n1) <= ygrid) then
         do n = n1, n2
           mask(i,n,k) = 3
           enddo
         n1 = n2
         n2 = 0
       elseif (n2 > n1 .AND. (n2-n1) > ygrid) then
          n1 = 0
          n2 = 0
       endif
22   continue

100 continue

!    write(31,*)mask
!   close(31)

!    k = 5
!   do j = 1, ny
!   do i = 1, nx
!     if (mask(i,j,k) /= 0) print*,'i,j,k,mask',i,j,k,mask(i,j,k)
!   enddo
!   enddo
       
   end subroutine unfold_preserve_hgradient

   subroutine unfold_using_specific_elev_data(RADDATAVR, &
         NX,NY,NZ,BADPT,vn,lev)

!--------------------------------------------------------------
! Velocity unfolding program based on the lower elevation data
!
! Eunha Lim                     April 2007
!--------------------------------------------------------------

   implicit none

   integer, intent(IN)         :: nx,ny,nz
   integer, intent(IN)         :: lev
   real,    intent(IN)         :: badpt

   real, dimension(nx,ny,nz), intent(INOUT)   :: raddatavr
   real, dimension(nz),       intent(IN)   :: vn

   real, parameter             :: gradient_factor = 1.0  !.8

   integer   :: i,j,k,kfold,count,klev,k1,k2,k3,iter
   real      :: fold



   do 1 iter = 1, 2
      k1 = lev
    if (iter == 1) then
      k2 = nz
      k3 = 1
    else
      k2 = 1
      k3 = -1
    endif
   do k = k1, k2, k3
     if (k == lev) CYCLE
    klev = k - k3
     count = 0
   do j = 1, ny
   do i = 1, nx
    if (raddatavr(i,j,klev) == badpt .OR. &
       raddatavr(i,j,k) == badpt) CYCLE

!--- Check the foldness

     fold = (raddatavr(i,j,klev) - raddatavr(i,j,k))*gradient_factor/vn(k)

    if (fold >= 0) then
      kfold = INT(fold + 0.5)
    else
      kfold = INT(fold - 0.5)
    endif

    if (kfold == 0) CYCLE 

    raddatavr(i,j,k) = raddatavr(i,j,k) + kfold*vn(k)
    count = count + 1

   enddo 
   enddo
    write(*,'(a,3i10)')'Number of unfolded data by lower elev.',k,klev,count
   enddo
1  continue

   end subroutine unfold_using_specific_elev_data

   subroutine unfold_velocity_bg3d(raddatavr,xk,yl,zm,idx_grdpt,scan_type,   &
             nelev_vad,num_ang,num_rng,nrad,irad,raltr,ecurv,iproj_type,    &
             u,v,dx,dy,dz,z00,badpt, &
             nxrad,nyrad,nzrad, &
             nyds, nyde, nxds, nxde, nzds, nzde)

!-------------------------------------------------------------------
! Unfold radial velocity using model background(3d) as a referenc velocity
!
! Usually forecast velocity is smaller than the observation,
! so sometimes unfolding is not correctly done.
! To avoid the situation, mutiply bg_amp_factor to the forecast vr.
!
! nyradquist velocity is hard coded for NEXRAD
!
!
! Eunha Lim         09/2006
!                               12/2006
!                               01/2007
!-------------------------------------------------------------------

   implicit none

! Radar data

   integer, intent(IN)        :: scan_type
   integer, intent(IN)        :: nxrad,nyrad,nzrad
   real, intent(IN)           :: badpt   
   real, intent(IN)           :: xk(nrad),yl(nrad),zm(nrad)
   real, intent(INOUT)        :: raddatavr(nxrad,nyrad,nzrad)
   integer, intent(IN)        :: idx_grdpt(nxrad,nyrad,nzrad,3)


!  Background wind  

   integer, intent(in) :: nyds, nyde, nxds, nxde, nzds, nzde

   real,intent(IN)            :: dx,dy,dz              ! m
   real,intent(IN)            :: z00                   ! m
   real,dimension(nxds:nxde, nyds:nyde, nzds:nzde),intent(IN) :: u, v

!  Nyquist velocity
   real,parameter              :: bg_amp_factor = 1.0  ! reference wind amplication factor
   character(len=80)           :: fname_vcp
   real,dimension(nzrad)       :: elev,vel_nyq,vn  

! VAD analysis related         
   integer,intent(IN)            :: num_ang,nrad,nelev_vad,irad
   integer,intent(IN)            :: num_rng,ecurv,iproj_type
   real,intent(IN)               :: raltr
   real,dimension(num_rng,nrad)  :: uvad,vvad
   real,dimension(num_rng)       :: z0

! Preserve horizontal gradient
   integer,dimension(nxrad,nyrad,nzrad)  :: mask

! Micelleneous
   integer   :: i,j,k,ii,jj,kk,count,kfold,ngrdp,nquad,VLEVELNZ
   real      :: xdist,ydist,disth,vr,vr_folded,fold,rmnpts


! Logical variables to control the methods to unfold
   logical, parameter           :: unfold_global = .TRUE., &
                                   unfold_local = .TRUE., &
                                   unfold_elev = .FALSE., &
                                   unfold_hgrd = .FALSE., &
                                   unfold_vad = .FALSE.,  &    
                                   unfold_preserve_hgrd = .TRUE.


!====================================================================
! Start Main program
!====================================================================

!--- Read Nyquist velocity 
   write(fname_vcp(:),'(a6,i2)')'./vcp.',scan_type

   print*,'scan_type=',scan_type
   print*,'VCP_profile =',TRIM(fname_vcp)
   open(2,file=TRIM(fname_vcp),form='formatted',status='old')

   read(2,*); read(2,*)
   do k = 1, nzrad
     read(2,*,end=10) elev(k),vel_nyq(k)
     vn(k) = 2.*vel_nyq(k)
     write(*,'(a,i5,2f8.1)')'vel_nyq= ',k,elev(k),vel_nyq(k)
   enddo
10 close(2)
   VLEVELNZ = k-1
!   write(*,'(a,i5)')'Number of vertical levels of Nyq_vel:',k-1


!------------------------------------------------------------------
!--- [0.0] Preserve horizontal gradient area
!------------------------------------------------------------------

    if (unfold_preserve_hgrd) then
      call unfold_preserve_hgradient(RADDATAVR, &
         NXRAD,NYRAD,NZRAD,BADPT,vel_nyq,mask)
   endif

!------------------------------------------------------------------
!--- [1.0] Global unfolding
!---       Unfolding using VDRAS first guess
!------------------------------------------------------------------

    if (unfold_global) then
!--- Check the foldness

   write(*,'(a,2f8.1)')'Max/Min before unfolding=', &
         MAXVAL(raddatavr,raddatavr/=badpt),MINVAL(raddatavr,raddatavr/=badpt)

   count = 0 

   do i = 1, nxrad
   do j = 1, nyrad

INNER :   do k = 1, nzrad

         if (raddatavr(i,j,k) == badpt) CYCLE INNER
         if (mask(i,j,k) == 3) CYCLE INNER

         ii = idx_grdpt(i,j,k,1)
         jj = idx_grdpt(i,j,k,2)
         kk = idx_grdpt(i,j,k,3)
         if (ii == INT(badpt) .or. jj == INT(badpt) .or. kk == INT(badpt)) CYCLE INNER

         xdist = (FLOAT(ii) - xk(irad)) * dx
         ydist = (FLOAT(jj) - yl(irad)) * dy
         disth = SQRT(xdist*xdist + ydist*ydist)

    if (disth == 0.) CYCLE INNER

     if (u(ii,jj,kk) == badpt .or. v(ii,jj,kk) == badpt) CYCLE INNER

    vr = (u(ii,jj,kk)*xdist + v(ii,jj,kk)*ydist)/disth * bg_amp_factor     


!--- Check the foldness

     fold = (vr - raddatavr(i,j,k))/vn(k)

    if (fold >= 0) then
      kfold = INT(fold + 0.5)
    else
      kfold = INT(fold - 0.5)
    endif

    if (kfold == 0) CYCLE INNER

         vr_folded = raddatavr(i,j,k)
    raddatavr(i,j,k) = raddatavr(i,j,k) + kfold*vn(k)
    count = count + 1

   enddo INNER
   enddo
   enddo

   write(*,'(a,2f8.1)')'Max/Min after unfolding =', &
         MAXVAL(raddatavr,raddatavr/=badpt),MINVAL(raddatavr,raddatavr/=badpt)
   write(*,'(a,i10)')'Total number of unfolded data =',count
   print*
   endif

!----------------------------------------
!--- [2] Additional Unfolding process
!----------------------------------------
!---------------------------
!--- [2.1] Local unfolding
!---------------------------
    if (unfold_local) then
      ngrdp = 8
      nquad = 4

!      do k = 1,5
         RMNPTS = (ngrdp*ngrdp)*0.8
         CALL local_unfold(RADDATAVR,NGRDP,NQUAD,RMNPTS, &
                      NXRAD,NYRAD,NZRAD,BADPT,vn)
!      enddo

!     ngrdp = 15
!      RMNPTS = (ngrdp*ngrdp)*0.8
!         CALL local_unfold(RADDATAVR,NGRDP,NQUAD,RMNPTS, &
!      NXRAD,NYRAD,NZRAD,BADPT,vn)
    endif

!------------------------------------------------------------------
!--- [2.2] Unfolding using specific elevation data
!          It is subjective, but it is worth for final brush it up
!------------------------------------------------------------------

    if (unfold_elev) then
      CALL unfold_using_specific_elev_data(RADDATAVR, &
         NXRAD,NYRAD,NZRAD,BADPT,vn,VLEVELNZ)
!--- Local unfolding
!      ngrdp = 15
!     RMNPTS = (ngrdp*ngrdp)*0.8
!      CALL local_unfold(RADDATAVR,NGRDP,NQUAD,RMNPTS, &
!         NXRAD,NYRAD,NZRAD,BADPT,vn)
    endif


!----------------------------------------------------------------------
!--- [2.3] Forced unfolding based on the horizontal gradient of the vr
!          It does not work very well
!----------------------------------------------------------------------

    if (unfold_hgrd) then
      call unfold_by_gradient_vr(RADDATAVR, &
         NXRAD,NYRAD,NZRAD,BADPT,vel_nyq)
    endif


!----------------------------------------------------------------------
!--- [2.4] Unfolding process using VAD analysis
!----------------------------------------------------------------------

    if (unfold_vad) then
      call vad_anal(elev(nelev_vad),num_rng,num_ang, &
             uvad,vvad,nelev_vad,raddatavr, &
             nrad,1,nxrad,1,nyrad,1,nzrad,  &
             dx*1.E-3,dy*1.E-3,dz*1.E-3,    &
             float(nxrad+1)/2.0,float(nyrad+1)/2.0,zm, &
             dz,z00,badpt,iproj_type,ecurv,irad)

      do k = 1, num_rng
        z0(k) = (Z00 + FLOAT(k)*DZ)*.001
        if(k /= 1 .AND. uvad(k,irad) == badpt) uvad(k,irad) = uvad(k-1,irad)
        if(k /= 1 .AND. vvad(k,irad) == badpt) vvad(k,irad) = vvad(k-1,irad)
        write(*,'(a,i5,3f8.2)')'k,uvad(k,irad),vvad(k,irad)',k,z0(k),uvad(k,irad),vvad(k,irad)
      enddo

   call unfold_velocity_using_vad(DX*1.E-3,DY*1.E-3,DZ*1.E-3, &
               raltr,NXRAD,NYRAD,NZRAD, &
               elev(1:VLEVELNZ),VLEVELNZ,IPROJ_TYPE,ECURV, &
               RADDATAVR,scan_type,badpt, &
               z0,uvad(:,irad),vvad(:,irad),num_rng)

!--- Local unfolding
      ngrdp = 15
     RMNPTS = (ngrdp*ngrdp)*0.8
      CALL local_unfold(RADDATAVR,NGRDP,NQUAD,RMNPTS, &
         NXRAD,NYRAD,NZRAD,BADPT,vn)
    endif

   end subroutine unfold_velocity_bg3d

   subroutine unfold_velocity_using_prev_vr(nx,ny,nz, &
             prev_vr,vr,badpt,scan_type)

!-------------------------------------------------------------------
! Unfold radial velocity using previously dealiased radial velocities
!
! Eunha Lim         08/2007
!
! NOTE: vr and prev_vr has different oder in argument
!-------------------------------------------------------------------

   implicit none

! Radar data
   integer,intent(IN)         :: nx,ny,nz
   integer,intent(IN)         :: scan_type
   real,intent(IN)            :: badpt   
   real,intent(IN),dimension(nx,ny,nz)::  prev_vr
   real,intent(INOUT),dimension(0:ny+1,0:nx+1,nz)::  vr

!  Nyquist velocity
   real,parameter              :: bg_amp_factor = 1.0  ! reference wind amplication factor
   character(len=80)           :: fname_vcp
   real,dimension(nz)       :: elev,vel_nyq,vn  

! Micelleneous
   integer   :: i,j,k,cnt,kfold,VLEVELNZ
   real      :: fold


!====================================================================
! Start Main program
!====================================================================

!--- Read Nyquist velocity 
   write(fname_vcp(:),'(a6,i2)')'./vcp.',scan_type

   print*,'scan_type=',scan_type
   print*,'VCP_profile =',TRIM(fname_vcp)
   open(2,file=TRIM(fname_vcp),form='formatted',status='old')

   read(2,*); read(2,*)
   do k = 1, nz
     read(2,*,end=10) elev(k),vel_nyq(k)
     vn(k) = 2.*vel_nyq(k)
     write(*,'(a,i5,2f8.1)')'vel_nyq= ',k,elev(k),vel_nyq(k)
   enddo
10 close(2)
   VLEVELNZ = k-1
!   write(*,'(a,i5)')'Number of vertical levels of Nyq_vel:',k-1
!   write(*,'(a,2f10.2)')'xk,yl',xk,yl
!   write(*,'(a,6i5)')'nx,ny,nz,nxyzrad',nx,ny,nz,nxrad,nyrad,nz


!------------------------------------------------------------------
!--- [1.0] Dealaising
!------------------------------------------------------------------

   write(*,'(a,2f8.1)')'Max/Min before unfolding=', &
         MAXVAL(vr,vr/=badpt),MINVAL(vr,vr/=badpt)

   cnt = 0 

   do i = 0, nx+1
   do j = 0, ny+1

INNER :   do k = 1, nz

         if (vr(j,i,k) == badpt .OR. prev_vr(i,j,k) == badpt) CYCLE INNER


!--- Check the foldness

         fold = (prev_vr(i,j,k) - vr(j,i,k))/vn(k)

        if (fold >= 0) then
          kfold = INT(fold + 0.5)
        else
          kfold = INT(fold - 0.5)
        endif

        vr(j,i,k) = vr(j,i,k) + kfold*vn(k)
        cnt = cnt + 1

   enddo INNER
   enddo
   enddo

   write(*,'(a,2f8.1)')'Max/Min after unfolding =', &
         MAXVAL(vr,vr/=badpt),MINVAL(vr,vr/=badpt)
   write(*,'(a,i10)')'Total number of unfolded data =',cnt
   print*

   end subroutine unfold_velocity_using_prev_vr

   subroutine unfold_velocity_using_vad(DXKMR,DYKMR,DZKMR, &
               raltkmr,NX,NY,NZ, &
               VLEVEL,VLEVELNZ,IPROJ_TYPE,ECURV, &
               RADDATAVR,scan_type,badpt, &
               hgt_bg,u_bg,v_bg,nz_bg)

!-------------------------------------------------------------------
! Unfold radial velocity using VAD as a referenc velocity
!
! Nyquist velocity is hard coded for NEXRAD
!
! Eunha Lim       May 2007   
!-------------------------------------------------------------------

   implicit none

! Radar data
   real,parameter         :: deg2rad = 0.017453292

   integer,intent(IN)     :: nx,ny,nz
   integer,intent(IN)     :: nz_bg
   integer,intent(IN)     :: ecurv,iproj_type,vlevelnz
   real,intent(IN)        :: dxkmr,dykmr,dzkmr
   real,intent(IN)        :: badpt    !vel_nyq
   real,intent(IN)        :: raltkmr
   real,dimension(nz_bg),intent(IN)       :: u_bg,v_bg,hgt_bg  ! height shoud be km

   integer,intent(IN)      :: scan_type
   real,dimension(vlevelnz),intent(IN)       :: vlevel
   real,dimension(nx,ny,nz),intent(INOUT)    :: raddatavr

!  Background wind  
   real,parameter                      :: bg_amp_factor = 1.0  ! reference wind amplication factor
   character(len=80)                   :: fname_bg,fname_vcp

! Micelleneous
   integer   :: i,j,k,l,k1,k2,flag,count,kfold
   real      :: lat,lon,elev
   real      :: xc,yc,xdist,ydist,disth,hgt_ppi,wu,wl,u,v,vr,vr_folded,fold,slope
   real, dimension(vlevelnz)                 :: vel_nyq,vn  


!====================================================================
! Start Main program
!====================================================================

!    print*,'DXKMR,DYKMR,DZKMR',DXKMR,DYKMR,DZKMR
!    print*,'raltkmr,NX,NY,NZ',raltkmr,NX,NY,NZ 
!    print*,'VLEVEL,VLEVELNZ,IPROJ_TYPE,ECURV',&
!          VLEVEL,VLEVELNZ,IPROJ_TYPE,ECURV
!    print*,'scan_type,badpt',scan_type,badpt
!    print*,'hgt_bg',hgt_bg
!    print*,'u_bg',u_bg
!    print*,'v_bg',v_bg

!--- For NEXRAD, temporary

     if (vlevelnz >= 15) then
      Print*,'# elevation is not match with VCP 11'
      stop
    endif


!--- Read Nyquist velocity 
   write(fname_vcp(:),'(a6,i2)')'./vcp.',scan_type
   open(2,file=TRIM(fname_vcp),form='formatted',status='old')
   print*,'VCP_profile =',TRIM(fname_vcp)

   read(2,*); read(2,*)
   do k = 1, vlevelnz
    read(2,*) elev,vel_nyq(k)
    vn(k) = 2.*vel_nyq(k)
    print*,'vel_nyq= ',k,vel_nyq(k)
   enddo
   close(2)

!--- Check the foldness

   xc = FLOAT(nx)/2. + 0.5
   yc = FLOAT(ny)/2. + 0.5
   count = 0

   do i = 1, nx
    xdist = (FLOAT(i)-xc)*dxkmr
   do j = 1, ny
    ydist = (FLOAT(j)-yc)*dykmr
INNER :   do k = 1, vlevelnz

    disth = SQRT(xdist*xdist + ydist*ydist)

    if (raddatavr(i,j,k) == badpt .OR. disth == 0.) CYCLE INNER

    if (ecurv == 0) then
      hgt_ppi = disth*TAN(vlevel(k)*deg2rad)                           !  + raltkmr
    else 
      hgt_ppi = disth*TAN(vlevel(k)*deg2rad) + 3./8.*disth**2/6370. !  + raltkmr
     endif


! Find the nearlest model height

     k2 = nz_bg
     do l = 1, nz_bg
      if (hgt_bg(l) >= hgt_ppi) then
       k2 = l
       exit
      endif
    enddo
    k1 = k2 - 1

    if( (k1 /= 0 .AND. u_bg(k1) == badpt) .OR. u_bg(k2) == badpt) CYCLE INNER

    if (k2 == 1) then    ! Extrapolation
!      slope = (u_bg(2)-u_bg(1))/(hgt_bg(2)-hgt_bg(1))
!      u = u_bg(1) + slope*(hgt_bg(1)-hgt_ppi)
!      v = v_bg(1) + slope*(hgt_bg(1)-hgt_ppi)
! assume the wind on the surface is 0
      wl = hgt_ppi/hgt_bg(k2)
      u = u_bg(k2)*wl
      v = v_bg(k2)*wl
!      print*,'ubg,vbg,u,v',u_bg(1),v_bg(1),u,v
    elseif (k2 == nz_bg) then  ! Put same value at the top
      u = u_bg(nz_bg) 
      v = v_bg(nz_bg) 
    else                 ! Interpolation
      wu = (hgt_bg(k2) - hgt_ppi)/(hgt_bg(k2) -hgt_bg(k1))
      wl = (hgt_ppi - hgt_bg(k1))/(hgt_bg(k2) -hgt_bg(k1))
      u = u_bg(k1)*wu + u_bg(k2)*wl
      v = v_bg(k1)*wu + v_bg(k2)*wl
    endif

    vr = (u*xdist + v*ydist)/disth * bg_amp_factor     
!    write(*,'(a,3i5,4f10.2)')'i,j,k',i,j,k,u,v,vr,raddatavr(i,j,k)
!    raddatavr(i,j,k) = vr                               ! For test base wind
!    CYCLE INNER


!--- Check the foldness

     fold = (vr - raddatavr(i,j,k))/vn(k)

    if (fold >= 0) then
      kfold = INT(fold + 0.5)
    else
      kfold = INT(fold - 0.5)
    endif

    if (kfold == 0) CYCLE INNER

     vr_folded = raddatavr(i,j,k)
    raddatavr(i,j,k) = raddatavr(i,j,k) + kfold*vn(k)
    count = count + 1

   enddo INNER
   enddo
   enddo

   Write(*,'(a,i10)')'Total number of unfolded data using VAD =',count

   end subroutine unfold_velocity_using_vad
   
   subroutine vad_anal(el,num_rng,num_ang, &
                    uvad,vvad,nelev_vad,vdradar, &
                    nrad, &
                    mx_bgn, mx_end, my_bgn, my_end, mz_bgn, mz_end, &
                    dx_cedkm,dy_cedkm,dz_cedkm, &
                    xradloc,yradloc,zm, &
                    dz,z00,rmissing,iproj_type,ecurv,irad)

   implicit none
   
   integer, intent(in) :: nrad
   integer, intent(in) :: mx_bgn, mx_end, my_bgn, my_end, mz_bgn, mz_end
   integer, intent(in) :: num_rng, num_ang, ecurv,irad, &
                          iproj_type,nelev_vad

   real, dimension(nrad), intent(in) :: zm
!yzm  real, dimension(num_rng,nrad), intent(in) :: uvad, vvad
   real, dimension(num_rng,nrad), intent(inout) :: uvad, vvad
   
   real, intent(in) :: el, dx_cedkm, dy_cedkm, dz_cedkm, xradloc,yradloc, &
                       dz, z00, rmissing

   real, dimension(mx_bgn:mx_end, my_bgn:my_end, mz_bgn:mz_end) :: vdradar
   
!     Local variables
   
   integer i,j,ielev,iloc,jloc
   
   real, dimension(num_rng,num_ang) :: vel
   real, dimension(num_rng) :: rng, z0, err, con
   real, dimension(num_ang) :: aza
   
   real :: tanel, delta_ang, &
           xloc, yloc, epx, epy
   
   real, parameter :: dtr = 0.017453293
   real, parameter :: rng_min = 25.0
   real, parameter :: cntmn = 20.0, &
                   tgap  = 120.0, &
                  rmsmx = 30.
!  real, parameter :: ALPHA=0.00005886970
   real, parameter :: alpha=5.88697e-5

!  print *, 'num_rng,num_ang,nrad,mx_bgn, mx_end, my_bgn, my_end, mz_bgn, mz_end,iproj_type,ecurv,irad=', &
!       num_rng,num_ang,nrad,mx_bgn, mx_end, my_bgn, my_end, mz_bgn, mz_end,iproj_type,ecurv,irad
   
!.. DATA IS INITIALLY IN CARTESIAN COORDINATES
!.. RESAMPLE RADIAL WIND ONTO CIRCLES OF RADIUS RNG_MIN
!.. AND HEIGHTS, Z=FLOAT(I)*DZ_CEDKM
   
   do i=1,num_rng
     if(iproj_type.eq.1)then
!.. IF PPI DATA, CALCULATE RANGE SO THAT VAD IS PERFORMED ON MODEL U,V SURFACES.
       tanel=tan(el*dtr)
       if(ecurv.eq.0) then
          rng(i)=(float(i)-zm(irad))*.001*dz/tanel
       else
          rng(i)=( -tanel + sqrt(tanel**2 + &
               4.*alpha*(float(i)-zm(irad))*.001*dz) )/(2.*alpha)
       end if
     elseif(iproj_type.eq.2)then
!.. IF CARTESIAN DATA, RANGE IS FIXED.
        rng(i)=rng_min
     endif
   enddo
   
   delta_ang = 360/num_ang
   do j=1,num_ang
      aza(j)=(j-1)*delta_ang
   enddo
   
   do i=1,num_rng
      if(iproj_type.eq.1)then
         ielev = nelev_vad
      elseif(iproj_type.eq.2)then
         ielev = i
      endif
      do j=1,num_ang
         xloc=xradloc+rng(i)*sin(aza(j)*dtr)/dx_cedkm
         yloc=yradloc+rng(i)*cos(aza(j)*dtr)/dy_cedkm
!yzm     iloc=int(xloc+0.0001)
!yzm     jloc=int(yloc+0.0001)
         iloc=int(xloc)
         jloc=int(yloc)
         if(iloc.lt.mx_bgn.or.iloc.ge.mx_end.or. &
            jloc.lt.my_bgn.or.jloc.ge.my_end)then
           vel(i,j) = rmissing
         else
           epx=xloc-iloc
           epy=yloc-jloc
           if(vdradar(iloc,jloc,ielev).eq.rmissing)then
             vel(i,j) = rmissing
           else
             vel(i,j)=vdradar(iloc,jloc,ielev)
           endif
!               if (i .eq. 1) then
!      print*,'iloc=,jloc=,vel',iloc,jloc,xloc,yloc,j,vel(i,j)
!               endif

!             IF(VDRADAR(ILOC,JLOC,IELEV).EQ.RMISSING.OR.
!    1           VDRADAR(ILOC+1,JLOC,IELEV).EQ.RMISSING.OR.
!    2           VDRADAR(ILOC,JLOC+1,IELEV).EQ.RMISSING.OR.
!    3           VDRADAR(ILOC+1,JLOC+1,IELEV).EQ.RMISSING)THEN
!                  VEL(I,J) = RMISSING
!             ELSE
!                  VEL(I,J)=(1.-EPY)*((1.-EPX)*VDRADAR(ILOC,JLOC,IELEV)
!    1                      +EPX*VDRADAR(ILOC+1,JLOC,IELEV))
!    2                      +EPY*((1.-EPX)*VDRADAR(ILOC,JLOC+1,IELEV)
!    3                      +EPX*VDRADAR(ILOC+1,JLOC+1,IELEV))
!             ENDIF
         endif
   enddo
   enddo
   
   write (6,fmt='(a,a,f6.1,a)') 'MAXIMUM ALLOWED AZIMUTH GAP FOR ', &
        'VAD ANALYSIS: ',tgap,' deg.'
   write (6,fmt='(a,a,I3)') 'MINIMUM NUMBER OF GOOD DATA POINTS ', &
        'FOR VAD ANALYSIS: ',int(cntmn)
   call vad(vel,aza,el,num_rng,num_ang,rng,cntmn, &
        tgap,rmsmx,uvad,vvad,err,con,rmissing,irad,nrad)
   
!..  set vad at lowest levels to missing
   
!      uvad(1,IRAD) = rmissing
!      vvad(1,IRAD) = rmissing
!      uvad(2,IRAD) = rmissing
!      vvad(2,IRAD) = rmissing
   
! EXTEND VAD PROFILE UP FROM LAST GOOD VALUE
   
!$$$      IGOODDATA=1
!$$$      DO 3000 I=NUM_RNG,1,-1
!$$$         IF(UVAD(I,IRAD).NE.RMISSING)THEN
!$$$            IGOODDATA=I
!$$$            GO TO 3010
!$$$         ENDIF
!$$$ 3000 CONTINUE
!$$$
!$$$ 3010 IF(IGOODDATA.EQ.1) GO TO 3030
!$$$
!$$$      DO 3020 I=IGOODDATA+1,NUM_RNG
!$$$         UVAD(I,IRAD)=UVAD(IGOODDATA,IRAD)
!$$$         VVAD(I,IRAD)=VVAD(IGOODDATA,IRAD)
!$$$ 3020 CONTINUE
!$$$
!$$$ 3030 CONTINUE
   
   do i=1,num_rng
      if(iproj_type.eq.1)then
         z0(i) = (z00 + float(i)*dz)*.001
      elseif(iproj_type.eq.2)then
         z0(i) = float(i-1)*dz_cedkm
      endif
      print 205,z0(i),rng(i),uvad(i,irad),vvad(i,irad),err(i),con(i)
   enddo
   
   print*
   print*,'***************************'
   print*,'** VAD ANALYSIS FINISHED **'
   print*,'***************************'
   print*
   
 205  FORMAT(1X,'Z= ',F6.3,2X,'RNG=',F6.1,2X,'UVAD,VVAD=',2F8.2, &
          '  ERR=',f8.2,'  CON=',f8.2)

   end subroutine vad_anal

!
!----------------------------------------------------------------------X
!
   subroutine vad(vel,aza,el,num_rng,num_ang,rng,cntmn,tgap,rmsmx, &
        uvad,vvad,err,con,rmissing,irad,nrad)
   implicit none
!
!  FUNCTION - COMPUTE MEAN RADIAL VELOCITY USING VAD ANALYSIS
!             F(I,J,OUT)=A0*A1*COS(A)*A2*COS(2*A)+B1*SIN(A)+B2*SIN(2*A)
!
!     VEL    - INPUT RADIAL VELOCITY FIELD. FIRST INDEX IS IN RANGE,
!                  SECOND IS FOR AZIMUTH ANGLE.
!     AZA    - ARRAY OF AZIMUTH ANGLES
!     EL     - ELEVATION ANGLE
!     NUM_RNG- NUMBER OF RANGE GATES
!     NUM_ANG- NUMBER OF AZIMUTH ANGLES
!     RNG    - RANGE (KM) TO EACH GATE OF DATA
!     CNTMN  - MINIMUM NUMBER OF GOOD DATA POINTS FOR VAD ANALYSIS
!     TGAP   - MAXIMUM ALLOWED AZIMUTH GAP         "   "      "   "
!     RMSMX  -    "       "    RMS DIFFERENCE BETWEEN INPUT AND VAD WINDS
   
!
!  VAD MEAN OUTPUT QUANTITIES
!
!     UVAD,VVAD  - HORIZONTAL WINDS       FOR THE ITH RANGE GATE
!     Z     - HEIGHT OF U,V POINT
!     CON    -     "      CONVERGENCE  "   "   "    "     "
!     ERR    - RMS DIFFERENCE BETWEEN MEASURED RADIAL VELOCITIES
!              AND VAD ANALYSIS
!
   
   integer irad, nrad, num_rng, num_ang
   
   real el, cntmn, tgap, rmsmx, rmissing
   
   real vel(num_rng,num_ang),aza(num_ang), &
      rng(num_rng),uvad(num_rng,nrad),vvad(num_rng,nrad), &
      err(num_rng),con(num_rng)
   
!     Local variables
   
   integer i, j
   
   real torad,todeg,sine,cose,z,a0,a1,a2,b1,b2,cnt,gapmx,gapmn,alft, &
        ang,angr,dat1,gap,sumsqdif,cntdif,vr_vad,vr_inp,rmsdif
   
   data torad,todeg/0.017453293,57.29577951/
   
   sine=sin(torad*el)
   cose=cos(torad*el)
   
   
!     DO 10 I=1,NUM_RNG-1
   do 10 i=1,num_rng
      uvad(i,irad) = rmissing
      vvad(i,irad) = rmissing
      con(i)= rmissing
      err(i)= rmissing
 10   CONTINUE
   
   
   
!
   
   
   do 100 i=1,num_rng
      if(rng(i).le.0.0)go to 100
      z=rng(i)*sine
      a0=0.0
      a1=0.0
      a2=0.0
      b1=0.0
      b2=0.0
      cnt=0.0
      gapmx=-999.0
      gapmn=999.0
      alft=rmissing
   
!        LOOP OVER ALL ANGLES TO CALCULATE THE FOURIER
!        COEFFICIENTS FOR ZEROTH, FIRST AND SECOND HARMONICS.
!
      do 90 j=1,num_ang
         ang=aza(j)
         angr=ang*torad
         dat1=vel(i,j)
         if(dat1.ne.rmissing)then
            a0=a0+dat1
            a1=a1+dat1*cos(angr)
            a2=a2+dat1*cos(angr*2.0)
            b1=b1+dat1*sin(angr)
!      print*,'b1=',j,b1,dat1,angr
            b2=b2+dat1*sin(angr*2.0)
            cnt=cnt+1.0
            if(alft.eq.rmissing)then
               alft=ang
            else
               gap=abs(ang-alft)
               alft=ang
               if(gap.gt.180.0)gap=abs(gap-360.0)
               if(gap.lt.gapmn)gapmn=gap
               if(gap.gt.gapmx)gapmx=gap
            end if
         end if
 90      CONTINUE
   
!     FROM FOURIER COEFFICIENTS:  CALCULATE THE VAD MEAN PARAMETERS
!     AND ANALYTIC WINDS FOR THIS RANGE, THEN GO ON TO NEXT RANGE
!
!     write (6,fmt='(a,I3,a,I4,a)') 'RANGE GATE: ',i, &
!          ' LARGEST GAP: ',int(gapmx),' deg.'
!     write (6,fmt='(a,a,I3)') 'NUMBER OF GOOD DATA POINTS ', &
!          'FOR VAD ANALYSIS: ',int(cnt)
   
      if(cnt.ge.cntmn.and.gapmx.le.tgap)then
         a0=    a0/cnt
         a1=2.0*a1/cnt
         a2=2.0*a2/cnt
         b1=2.0*b1/cnt
         b2=2.0*b2/cnt
!    print*,'a1=,b1=',a1,b1
         uvad(i,irad)=(b1/cose)
         vvad(i,irad)=(a1/cose)
         con(i)=-2.0*a0/(rng(i)*cose*cose)
         sumsqdif=0.0
         cntdif=0.0
         do 92 j=1,num_ang
            ang=aza(j)
   
            angr=ang*torad
            vr_vad=a0+a1*cos(angr)+a2*cos(angr*2.0) &
                            +b1*sin(angr)+b2*sin(angr*2.0)
            vr_inp=vel(i,j)
   
            if(vr_inp.ne.rmissing)then
               cntdif=cntdif+1.0
               sumsqdif=sumsqdif+(vr_inp-vr_vad)**2
            end if
 92         CONTINUE
   
!           CHECK RMS DIFFERENCE BETWEEN VAD AND INPUT WINDS.
   
!
         if(cntdif.gt.cntmn)then
            rmsdif=sqrt(sumsqdif/cntdif)
            err(i)=rmsdif
         end if
         if(cntdif.le.cntmn.or.rmsdif.gt.rmsmx)then
            uvad(i,irad) =rmissing
            vvad(i,irad) =rmissing
            con(i)=rmissing
         end if
      end if
 100  CONTINUE
   
   
   return
   end subroutine
   
   subroutine varia(data,var,igrid,mx,my,mz,bad)
   
   implicit none
   
   integer, intent(in) :: igrid,mx,my,mz
   
   real, dimension(mx,my,mz), intent(in)  :: data
   real, dimension(mx,my,mz), intent(out) :: var
   
   real, intent(in)    :: bad
   
!  Local variables
   
   integer :: l,k,m,k0,k1,l0,l1,ll,kk,ngpt
   real    :: rlocm
   
!
!     Calculate 2-d local variance
!     The size of the local domain is 2*IGRID
!

!!$omp parallel do default(shared) &
!!$omp private(k,l,m,k0,k1,l0,l1,kk,ll,ngpt,rlocm)
   m_loop: do m=1,mz
   k_loop: do k=1,my
   l_loop: do l=1,mx
!yzm  k_loop: do k=1,mx
!yzm  l_loop: do l=1,my
   
     k0=max0(1,k-igrid)
!yzm k1=min0(mx,k+igrid)
     k1=min0(my,k+igrid)
     l0=max0(1,l-igrid)
!yzm l1=min0(my,l+igrid)
     l1=min0(mx,l+igrid)
   
!     Calculate local mean
   
     ngpt=0
     rlocm=0.
     do kk=k0,k1
     do ll=l0,l1
        if((data(ll,kk,m) - bad) > 1.0) then
           rlocm=rlocm+data(ll,kk,m)
           ngpt=ngpt+1
        end if
     enddo
     enddo
   
     if(ngpt > (igrid+1)*(igrid+1)) then
        rlocm=rlocm/real(ngpt)
     else
        if(ngpt < 1) then
           var(l,k,m)=bad
        else
!----------Remove the data point if no data surrounding it.

           var(l,k,m)=999.
        endif

        cycle l_loop
     end if
   
!     Calculate local variance
   
     var(l,k,m)=0.
     do kk=k0,k1
     do ll=l0,l1
        if((data(ll,kk,m) - bad) > 1.0) then
           var(l,k,m)=var(l,k,m)+(data(ll,kk,m)-rlocm)**2
        end if
     enddo
     enddo
   
     var(l,k,m)=var(l,k,m)/real(ngpt)

   enddo l_loop
   enddo k_loop
   enddo m_loop
   
    print *,'var min/max:',minval(var,var>bad),maxval(var,var/=999.)

   end subroutine varia

   subroutine local_unfold (DZ,IGRID,NOCT,RMNPTS,   &
                              NX,NY,NZ,BAD,VN)

   implicit none

!  THIS ROUTINE PERFORMS LOCAL UNFOLDING 
!
!  DZ  -INPUT 3-DIMENSIONAL DATA SET TO BE UNFOLDED.
!  IGRID-  MAXIMUM NUMBER OF GRID POINTS TO SEARCH OUTWARD FROM A 
!         MISSING DATA LOCATION TO DETERMINE IF IT IS BOUNDED AND IF
!         RMNPTS HAS BEEN SATISFIED.
!  RMNPTS-MINIMUM NUMBER OF SURROUNDING POINTS REQUIRED
!  NOCT-  MINIMUM NUMBER OF QUADRANTS OCCUPIED IN ORDER TO SATISFY THE SEARCH.
!  NX-    NUMBER OF X GRID POINTS
!  NY-    NUMBER OF Y GRID POINTS
!  NZ-    NUMBER OF Z GRID POINTS
!  BAD-   MISSING DATA FLAG
!  VN-    2*(NYquist Velocity)

      integer :: igrid, nx, ny, nz, noct
      real :: rmnpts, bad
  
      real DZ(NX,NY,NZ),IOCT(8),VN(NZ)
      real, parameter :: EPS=0.0001

      integer :: i, j, k, i1, j1, i2, j2, l, ii, jj, ix, iy, kq, k1, n_unfold, kfold

      real :: denom, fold, RNPTS, RNUM, SD, SDX, SDY, SX, SX2, SXY, SY, SY2, &
              T1, T2, T3, vr

      PRINT*
      PRINT *,'PERFORM LOCAL UNFOLDING'


      DO 10 K=1,NZ
      n_unfold = 0
      DO 11 J=1,NY
      DO 11 I=1,NX

! LOOP OVER GRID POINTS. IF GRID POINT VALUE IS MISSING DO NOTHING, 
! OTHERWISE CHECK FOLDNESS.

       IF(DZ(I,J,K).EQ.BAD) CYCLE

       DO 55 L=1,4
55       IOCT(L)=0

! START SEARCHING OUTWARD FROM GRID POINT UNTIL THE NOCT AND RMNPTS
! CONDITIONS ARE MET.
         DO 60 L=1,IGRID
           J1=MAX0(1,J-L)
           J2=MIN0(NY,J+L)
           I1=MAX0(1,I-L)
           I2=MIN0(NX,I+L)
! INITIALIZE SUMMATION TERMS USED IN THE LEAST-SQUARES FIT
           RNPTS=0.
           SX=0.
           SY=0.
           SX2=0.
           SY2=0.
           SXY=0.
           SDX=0.
           SDY=0.
           SD=0.
           DO 90 JJ=J1,J2
             IY=J-JJ
           DO 90 II=I1,I2
             IX=I-II
             IF(DZ(II,JJ,K).EQ.BAD) CYCLE
             IF(IX.GE.0.AND.IY.GT.0.)IOCT(1)=1
             IF(IX.GT.0.AND.IY.LE.0.)IOCT(2)=1
             IF(IX.LE.0.AND.IY.LT.0.)IOCT(3)=1
             IF(IX.LT.0.AND.IY.GE.0.)IOCT(4)=1
             RNPTS=RNPTS+1.0
             SX=SX+IX
             SY=SY+IY
             SX2=SX2+IX*IX
             SY2=SY2+IY*IY
             SXY=SXY+IX*IY
             SD=SD+dz(ii,jj,k)
             SDX=SDX+IX*dz(ii,jj,k)
             SDY=SDY+IY*dz(ii,jj,k)
90         CONTINUE
           KQ=0
           DO 75 K1=1,4
75         KQ=KQ+IOCT(K1)

           IF(KQ.LT.NOCT) GO TO 60
           IF(RNPTS.LT.RMNPTS) GO TO 60
           T1 = SX2*SY2-SXY*SXY
           T2 = SX*SY2-SXY*SY
           T3 = SX*SXY-SX2*SY
           DENOM=RNPTS*T1-SX*T2+SY*T3
           IF(DENOM.LE.EPS) GO TO 60
           RNUM=SD*T1-SDX*T2+SDY*T3
           vr = RNUM/DENOM
           fold = (vr - dz(i,j,k))/vn(k)

           if (fold >= 0) then
             kfold = INT(fold + 0.5)
           else
             kfold = INT(fold - 0.5)
           endif

           if (kfold /= 0) then
                  dz(i,j,k) = dz(i,j,k) + kfold*vn(k)
                  n_unfold = n_unfold + 1
             GO TO 11
           endif
60      CONTINUE  
11    CONTINUE
      print*,'Number of unfolded data by local unfolding',k,n_unfold
10    CONTINUE

   end subroutine local_unfold

   subroutine pre_analyze(dxkmr,dykmr,dzkmr,xminr,yminr, &
                       dx,dy,dz,vrthr,dtheta, &
                       vlevel,yl,xk,zm,raltr, &
                       igridv,iedit,ivol, &
                       vlevelnz,iproj_type,ecurv,pex,irealtime, &
                       raddataz, raddatavr,  &
                       zppi, pu, pz, pr, cs, urs, rs, std,  &
                       zcart, cscart, urscart, pucart, pzcart, &
                       nyds, nyde, nxds, nxde, nzds, nzde, &
                       nyms, nyme, nxms, nxme, nzms, nzme, &
                       nyts, nyte, nxts, nxte, nzts, nzte, &
                       nytsm1, nytep1, nxtsm1, nxtep1, &
                       nxrad,nyrad,nzrad, &
                       nrad,nobs,nzd,nze,nze_top, &
                       ngrdp,nquad,rmnpts,intpbl,mpt,irad, &
                       x00,y00,z00,badpt,iunfold,nscan,uu,vv, &
                       nelev_vad,num_ang,num_rng)        !yzm for unfolding
!
!     Interpolate radar data to the model grid
!     in the horizontal (if iproj_type=1)
!     and in the vertical (if iproj_type=2)
!

!   use dmp_util_module
!  use module_define

   implicit none

   integer, intent(in) :: igridv,iedit,ivol, &
                          vlevelnz,iproj_type,ecurv,irealtime, &
                          nyds, nyde, nxds, nxde, nzds, nzde, &
                          nyms, nyme, nxms, nxme, nzms, nzme, &
                          nyts, nyte, nxts, nxte, nzts, nzte, &
                          nytsm1, nytep1, nxtsm1, nxtep1, &
                          nxrad,nyrad,nzrad, &
                          nrad,nobs,nzd,nze, &
                          ngrdp,nquad,intpbl,mpt,irad
   integer, intent(out) :: nze_top

   real,    intent(in) :: dxkmr,dykmr,dzkmr,xminr,yminr, &
                          dx,dy,dz,pex,vrthr,dtheta,raltr

   real,    intent(in) :: x00,y00,z00,badpt,rmnpts

   
   real, dimension(nzrad), intent(in)    :: vlevel 
   real, dimension(nrad),  intent(in)    ::  yl, xk, zm

   real, dimension(nxrad,nyrad,nzrad), intent(inout) :: raddataz,raddatavr
   real, dimension(nxrad,nyrad,nzrad) :: extrap1, extrap2, variac

   real, dimension(nyms:nyme, nxms:nxme, nze, nobs, nrad), &
         intent(inout)  :: pu, pz, pr, cs, urs, rs, std

   real, dimension(nyms:nyme, nxms:nxme, nze, 3, nrad), intent(inout) :: zppi

   real, dimension(nyms:nyme, nxms:nxme, nzms:nzme, nobs, nrad), intent(inout) :: cscart, urscart, pzcart, pucart

   real, dimension(nzms:nzme), intent(inout) :: zcart

   real, dimension(nyms:nyme, nxms:nxme, nzrad) :: radz, radvr, purad, pzrad, radsd

   integer :: l, k, m
   integer :: ipt, jpt, npt, icount, mcs, mppi
   integer :: lmindx, kmindx, lmaxdx, kmaxdx, nztop, ics, icsp, jcs, jcsp
   real    :: refd, xminrm, yminrm, dxr, dyr, dzr, xmindx, ymindx, zmindx
   real    :: xd, xd1, xd2, yd, yd1, yd2, zd, zd1, zd2, disth
   real    :: rad_x0, rad_x1, som, degrad, dzppi, ranged, stdmax, stdvr, stdqr
   real, dimension(nzrad) :: sum_vr

   integer, dimension(nzrad) :: ngpt
   character(len=256) :: msg
   

!----Velocity unfolding---------------------------
   integer, intent(in) :: iunfold,nscan
   real, dimension(nxds:nxde, nyds:nyde, nzds:nzde), intent(in) :: uu, vv
   integer :: idx_grdpt(nxrad,nyrad,nzrad,3)
   integer :: ii,jj,kk, i,j
   real,allocatable,dimension(:,:,:)  :: zppi2, local_zppi2
   integer :: nelev_vad,num_ang,num_rng

!----------------------------------------------------------
!  integer, dimension(nxrad, nyrad) :: gmask        !yzm
!  integer :: m0, imask 
   real    :: badpt1              !yzm
!------------------------------------------

!     Initialize some arrays
   
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nzrad
   do l=nytsm1,nytep1
   do k=nxtsm1,nxtep1
      radz(l,k,m)=badpt
      radvr(l,k,m)=badpt
      pzrad(l,k,m)=0.
      purad(l,k,m)=0.
   enddo
   enddo
   enddo
   
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nze
   do l=nytsm1,nytep1
   do k=nxtsm1,nxtep1
      pu(l,k,m,ivol,irad)=0.
      pr(l,k,m,ivol,irad)=0.
      pz(l,k,m,ivol,irad)=0.
   enddo
   enddo
   enddo

!   print *,'raddataz min/max:',minval(raddataz),maxval(raddataz)
!   print *,'raddatavr min/max:',minval(raddatavr),maxval(raddatavr)


!     Calculate the data grid index of the first model grid (0,0,1)
   
   xminrm=xminr*1000.
   yminrm=yminr*1000.
   dxr=dxkmr*1000.
   dyr=dykmr*1000.
   dzr=dzkmr*1000.
   ymindx=(y00-yminrm)/dyr+1.
   xmindx=(x00-xminrm)/dxr+1.
   zmindx=z00/dzr+1.
   lmindx=max(int(ymindx),1)
   kmindx=max(int(xmindx),1)
   lmaxdx=min(int(ymindx+(dy*(nyde-nyds+1)+dy)/dyr),nyrad)
   kmaxdx=min(int(xmindx+(dx*(nxde-nxds+1)+dx)/dxr),nxrad)

   print *, 'xminr, yminr=', xminr, yminr
   print *, 'dxkmr, dykmr, dzkmr=', dxkmr, dykmr, dzkmr
   print *, 'y00, x00, z00=', y00, x00, z00
   print *, 'xminrm, xminrm=', xminrm, xminrm
   print *, 'ymindx, xmindx, zmindx=', ymindx, xmindx, zmindx
   print *, 'kmindx,kmaxdx=', kmindx,kmaxdx
   print *, 'lmindx,lmaxdx=', lmindx,lmaxdx
   
   if(lmindx.gt.nyrad.or.kmindx.gt.nxrad.or. &
        lmaxdx.lt.1.or.kmaxdx.lt.1) then
      print*,'Radar # ',irad, &
           ' data are all outside of analysis domain.'
      print*,'Increase the analysis domain coverage or ', &
           'remove this radar.'
      stop
      !write(msg, fmt='(3a,i6)') 'Stopped in file: <', __FILE__, '>, at line:', __LINE__
      !call vdras_finalize(trim(msg))
   endif
   if(iproj_type.eq.2) then
      nztop=int(dz*(nzde-nzds+1)/dzr+2.)
   else
      nztop=vlevelnz
   end if

   nze_top = nztop

   if(nze_top > nze) nze_top = nze

   print *, 'nztop=', nztop

   if(nztop > nzrad) then
      print *, 'nzrad=', nzrad
      print *, "nztop > nzrad"
      stop
      !write(msg, fmt='(3a,i6)') 'Stopped in file: <', __FILE__, '>, at line:', __LINE__
      !call vdras_finalize(trim(msg))
   endif
   
   if (ivol.eq.1) then
      write(6,111) 'data grid index of model minx,maxx:' &
           ,kmindx,kmaxdx
      write(6,111) 'data grid index of model miny,maxy:' &
           ,lmindx,lmaxdx
      print*
   endif
 111  format(1x,a40,2i6)
   
      print *,'raddataz min/max:',minval(raddataz,raddataz/=badpt),maxval(raddataz,raddataz/=badpt)
      print *,'raddatavr min/max:',minval(raddatavr,raddatavr/=badpt),maxval(raddatavr,raddatavr/=badpt)

!----remove ground clutter---------
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
      do m=1,nztop
      do l=1,nyrad
      do k=1,nxrad
         if(abs(raddatavr(k,l,m)).le.vrthr) then
            raddataz(k,l,m)=badpt
            raddatavr(k,l,m)=badpt
         end if
      enddo
      enddo
      enddo
   
      print *,'aft vrthr raddataz min/max:',minval(raddataz,raddataz/=badpt),maxval(raddataz,raddataz/=badpt)
      print *,'aft vrthr raddatavr min/max:',minval(raddatavr,raddatavr/=badpt),maxval(raddatavr,raddatavr/=badpt)
!     Data filling and compute local variance. igrid can not exceed nztop-1
   
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nzrad
   do l=1,nyrad
   do k=1,nxrad
      extrap1(k,l,m)=0.
      extrap2(k,l,m)=0.
   enddo
   enddo
   enddo

   write(unit=*, fmt='(a,i6)') &
        'iproj_type=', iproj_type, &
       'irealtime =', irealtime
   
   if(iproj_type.eq.2.and.irealtime.ne.2) then
      print *, 'ngrdp,nquad,rmnpts=', ngrdp,nquad,rmnpts
      print *, 'nxrad,nyrad,nzrad=', nxrad,nyrad,nzrad
      print *, 'fill2d raddataz: iedit =', iedit
   
      call fill2d(raddataz,extrap2,ngrdp,nquad,rmnpts, &
                  nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,0,pex)
   
      print *,'aft fill2d raddataz min/max:',minval(raddataz,raddataz/=badpt),maxval(raddataz,raddataz/=badpt)
!
!     Remove data with high variance
!
      if(iedit.eq.1) then
   
         call varia(raddatavr,variac,igridv,nxrad,nyrad,nzrad,badpt)
         ngpt(:)=0
   
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
         do m=1,nztop
         do l=1,nyrad
         do k=1,nxrad
            if(variac(k,l,m).gt.60.) then
               raddatavr(k,l,m)=badpt
               ngpt(m)=ngpt(m)+1
            end if
         enddo
         enddo
         enddo
   
         print *,'nztop =', nztop
         print *,'Number of data removed due to high variance: ', sum(ngpt(:))
      end if
   
!     CALL FILL3D(raddatavr,TEMP,EXTRAP1,3,3,30.,NXRAD,NYRAD,NZRAD, &
!                 NXRAD,NYRAD,NZRAD,BADPT,1,PEX)

!  print *, 'fill2d raddatavr:'

      call fill2d(raddatavr,extrap1,ngrdp,nquad,rmnpts, &
                 nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,1,pex)
   
      print *,'aft fill2d raddatavr min/max:',minval(raddatavr,raddatavr/=badpt),maxval(raddatavr,raddatavr/=badpt)

   else if( iproj_type.eq.1.and.irealtime.ne.2 ) then
   
      print *, 'ngrdp,nquad,rmnpts=', ngrdp,nquad,rmnpts
      print *, 'nxrad,nyrad,nzrad=', nxrad,nyrad,nzrad
      print *, 'pex=', pex

   degrad=0.017453292
   do m=nzts,nzte
     zcart(m) = z00 + float(m)*dz
   enddo

!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nztop
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
     disth=sqrt(((float(k)-xk(irad))*dx)**2 &
                +((float(l)-yl(irad))*dy)**2)
     if(ecurv.eq.0) then
!yzm    dzppi=0.
        dzppi=zm(irad)*dz + z00
     else
!yzm    dzppi=3./8.*disth**2/6370000.
        dzppi=3./8.*disth**2/6370000. + zm(irad)*dz + z00
     end if
     zppi(l,k,m,2,irad)=disth*tan(vlevel(m)*degrad)+dzppi
     zppi(l,k,m,1,irad)=disth*tan((vlevel(m)-dtheta/2.) &
          *degrad)+dzppi
     if(zppi(l,k,m,1,irad).lt.0.0) then
        zppi(l,k,m,1,irad)=0.
     end if
     zppi(l,k,m,3,irad)=disth*tan((vlevel(m)+dtheta/2.) &
          *degrad)+dzppi
   enddo
   enddo
   enddo
   
    write(*,'(a,2f12.2)')'pre_analyze zppi: MIN/MAX',MINVAL(zppi),MAXVAL(zppi)


!yzm --- Velocity unfolding-------------
 
if (iunfold == 1) then


   allocate(local_zppi2(nxts:nxte, nyts:nyte, nztop))
   allocate(zppi2(nxds:nxde, nyds:nyde, nztop))

!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nztop
   do k=nxts,nxte
   do l=nyts,nyte
       local_zppi2(k,l,m)  = zppi(l,k,m,2,irad)
   enddo
   enddo
   enddo

!#ifdef DM_PARALLEL
! call vdras_collect(zppi2, local_zppi2, nztop)
!!call vdras_broadcast(zppi2, nxde, nyde, nztop)
!#else
zppi2 = local_zppi2
!#endif

  write(*,'(a,2f12.5)')'local_zppi2: MIN/MAX',MINVAL(local_zppi2,local_zppi2/=badpt),MAXVAL(local_zppi2)
  write(*,'(a,2f12.5)')'zppi2: MIN/MAX',MINVAL(zppi2,zppi2/=badpt),MAXVAL(zppi2)

! if(irad.eq.1) print*,'local_zppi2=', zppi2(nxts:nxte,nyts,2)

! if(on_monitor) then
!   if(irad.eq.1) print*,'zppi2=', zppi2(:,nyds,2)
! endif

    deallocate(local_zppi2)

!  if(on_monitor) then

!--- 1st step : Remove the noise based on high variance

        call varia(raddatavr,variac,igridv,nxrad,nyrad,nzrad,badpt)
         ngpt(:)=0

!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
         do m=1,nztop
         do l=1,nyrad
         do k=1,nxrad
            if(variac(k,l,m).gt.60.) then
               raddatavr(k,l,m)=badpt
               ngpt(m)=ngpt(m)+1
            end if
         enddo
         enddo
         enddo

         print *,'Number of data removed due to high variance: ', sum(ngpt(:))

!--- 2nd step : Determine the grid points of background with for each
!               radar data points

!!$omp parallel do default(shared) &
!!$omp private(i,j,k,ii,jj,kk)
        IDX_GRDPT(:,:,:,:) = INT(BADPT)
        DO  I = 1, NXRAD
          II = (FLOAT(I)*DXR + XMINRM - X00) / DX
!yzm      IF (II < nxds-1 .OR. II > nxde+1) CYCLE
          IF (II < nxds .OR. II > nxde) CYCLE
        DO  J = 1, NYRAD
          JJ = (FLOAT(J)*DYR + YMINRM - Y00) / DY
!yzm      IF (JJ < nyds-1 .OR. JJ > nyde+1) CYCLE
          IF (JJ < nyds .OR. JJ > nyde) CYCLE
!yzm    DO  K = 1, NZRAD
        DO  K = 1, nztop
!         kk = (zppi2(ii,jj,k) - Z00) / DZ + 1
          kk = INT((zppi2(ii,jj,k) - Z00) / DZ)
          kk = MAX(1,KK)
          IF (kk > nzde) CYCLE
          IDX_GRDPT(I,J,K,1) = II
          IDX_GRDPT(I,J,K,2) = JJ
          IDX_GRDPT(I,J,K,3) = KK
       enddo
       enddo
       enddo

    deallocate(zppi2)

!leh      print*,'max IDX_GRDPT(I,J,K,1)', MAXVAL(IDX_GRDPT(:,:,:,1))
!leh      print*,'max IDX_GRDPT(I,J,K,2)', MAXVAL(IDX_GRDPT(:,:,:,2))
!leh      print*,'max IDX_GRDPT(I,J,K,3)', MAXVAL(IDX_GRDPT(:,:,:,3))


!     call unfold_velocity_bg3d(raddatavr,xk(irad),yl(irad),idx_grdpt,nscan, &
!                               uu,vv,dx,dy,badpt, &
!                               nxrad,nyrad,nzrad, &
!                               nyds, nyde, nxds, nxde, nzds, nzde) 

      call unfold_velocity_bg3d(raddatavr,xk,yl,zm,idx_grdpt,nscan, &
                                nelev_vad,num_ang,num_rng,nrad,irad,raltr,ecurv,iproj_type, &
                                uu,vv,dx,dy,dz,z00,badpt, &
                                nxrad,nyrad,nzrad, &
                                nyds, nyde, nxds, nxde, nzds, nzde) 

!---------ifill = 0, Do not smooth, do only fill

      call fill2d(raddatavr,extrap1, 4, 2, 20., &
                  nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,0,pex)

!  endif

!#ifdef DM_PARALLEL
!call vdras_broadcast(raddatavr, nxrad, nyrad, nzrad)
!#endif

    print *,'aft unfold raddataz min/max:',minval(raddataz,raddataz>badpt),maxval(raddataz)
    print *,'aft unfold raddatavr min/max:',minval(raddatavr,raddatavr>badpt),maxval(raddatavr)

endif
!yzm--------------end unfold--------

   
!----remove ground clutter---------
!!$omp parallel do default(shared)  &
!!$omp private(l,k,m)
!      do m=1,nztop
!      do l=1,nyrad
!      do k=1,nxrad
!         if(abs(raddatavr(k,l,m)).le.vrthr) then
!            raddataz(k,l,m)=badpt
!            raddatavr(k,l,m)=badpt
!         end if
!      enddo
!      enddo
!      enddo
   
   
!     Remove data with high variance
!
      if(iedit.eq.1) then

         call varia(raddataz,variac,igridv,nxrad,nyrad,nzrad,badpt)

      print *,'variac for raddataz min/max:', minval(variac,variac/=badpt),maxval(variac)

         ngpt(:)=0
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
          do m=1,nztop
          do l=1,nyrad
          do k=1,nxrad
             if(variac(k,l,m).gt.150.) then
!            if(variac(k,l,m).gt.100.) then
                raddataz(k,l,m)=badpt
                ngpt(m)=ngpt(m)+1
             end if
          enddo
          enddo
          enddo
 
!       print *,'variac min/max:',minval(variac),maxval(variac)


          print *,'nztop =', nztop
          print *,'Number of ref data removed due to high variance: ',sum(ngpt(:))
       end if
 
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nzrad
   do l=1,nyrad
   do k=1,nxrad
      extrap2(k,l,m)=0.
   enddo
   enddo
   enddo

      call fill2d(raddataz,extrap2,ngrdp,nquad,rmnpts, &
                  nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,1,pex)

!
!     Remove data with high variance
!
      if(iedit.eq.1) then
   
         call varia(raddatavr,variac,igridv,nxrad,nyrad,nzrad,badpt)

         ngpt(:)=0
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
           do m=1,nztop
        do l=1,nyrad
        do k=1,nxrad
            if(variac(k,l,m).gt.60.) then
               raddatavr(k,l,m)=badpt
               ngpt(m)=ngpt(m)+1
            end if
         enddo
         enddo
         enddo

         print *,'nztop =', nztop
         print *,'Number of vr data removed due to high variance: ',sum(ngpt(:))
      end if
   
   do m=1,nzrad
   do l=1,nyrad
   do k=1,nxrad
      extrap1(k,l,m)=0.
   enddo
   enddo
   enddo
      call fill2d(variac,extrap1,ngrdp,nquad,rmnpts, &
           nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,1,pex)

   do m=1,nzrad
   do l=1,nyrad
   do k=1,nxrad
      extrap1(k,l,m)=0.
   enddo
   enddo
   enddo
      call fill2d(raddatavr,extrap1,ngrdp,nquad,rmnpts, &
           nxrad,nyrad,nzrad,nxrad,nyrad,nzrad,badpt,1,pex)

   end if

!   print *,'raddataz min/max:',minval(raddataz),maxval(raddataz)
!   print *,'raddatavr min/max:',minval(raddatavr),maxval(raddatavr)

      do m=1,nze_top
      sum_vr(m) = 0.0
      do l = 1,nyrad
      do k = 1,nxrad
!         if(raddatavr(k,l,m) .ne. badpt) then
         if(raddatavr(k,l,m) .gt. badpt) then
        sum_vr(m) = sum_vr(m) + raddatavr(k,l,m) 
         endif
           enddo
      enddo
     enddo
   
!yzm---add mask-----under terrain is missing -8888.88 ---

!  imask = 0
!  imask = 1

!if(imask == 1) then
!
!   badpt1 = -8888.88
!
!
! if(on_monitor) then
!
!! open(unit=73, file='/data/zhuming/data/taiwan/Meiyu/mask/mask.dat',form='formatted',status='old')
!  open(unit=73, file='./mask.dat',form='unformatted',status='old')
!
!    read(73) gmask
!    write(6,'(i4,45i2)') (gmask(l,98), l=1,460)
!    write(6,'(//)')
!    write(6,'(i4,45i2)') (gmask(l,460), l=1,460)
!
!!  do k = 1, nyrad
!!   read(73,'(460i2)',end=92) (gmask(l,k), l=1,nxrad)
!!   write(6,'(i4,45i2)') k,(gmask(l,k), l=1,45 )
!!  enddo
!
!!92 continue
!
!  close(73)
!
!!!$omp parallel do default(shared) &
!!!$omp private(k,l,m)
!   do k=1,nyrad
!   do l=1,nxrad
!
!      m0 = gmask(l,k)
!!     if(m0 >= nzrad+1) m0 = nzrad+1
!      if(m0 >= nztop+1) m0 = nztop+1
!!      if(k==98) print*,'m0=', l, gmask(l,k), m0
!
!      if(m0 .gt. 1) then
!      
!       do m = 1, m0-1
!         raddatavr(l,k,m) = badpt1
!         raddataz(l,k,m) = badpt1
!         variac(l,k,m) = badpt1
!!   if(raddataz(l,k,m).eq.badpt1) print*,'rad in badpt1= ',l,k,m,raddataz(l,k,m)
!       enddo
!      endif
!
!   enddo
!   enddo
!
! endif


!#ifdef DM_PARALLEL
!call vdras_broadcast(raddatavr, nxrad, nyrad, nzrad)
!call vdras_broadcast(raddataz, nxrad, nyrad, nzrad)
!call vdras_broadcast(variac, nxrad, nyrad, nzrad)
!#endif
!
!endif           

!yzm-------------------------------------
 
!     Horizontal interpolation of data to model grid;
!     intpbl=1, bilinear if more than 1 good data point, nearest point if
!     only 1 good data point. intpbl=0, average over the number of good
!     data point if the number of good grid point greater than mpt.
!
!     print*,'intpbl=',intpbl,mpt

!!$omp parallel do default(shared) &
!!$omp private(l,k,m,yd,jcs,jcsp)
   do m=1,nztop
   do l=nytsm1,nytep1
      yd=float(l)*dy/dyr+ymindx
      jcs=int(yd+0.00001)
      yd1=yd-float(jcs)
      yd2=1.0-yd1
      jcsp=jcs+1
   do k=nxtsm1,nxtep1
      xd=float(k)*dx/dxr+xmindx
      ics=int(xd+0.00001)
      xd1=xd-float(ics)
      xd2=1.0-xd1
      icsp=ics+1
   
      if (jcs.ge.1.and.jcsp.le.nyrad.and.ics.ge.1.and.icsp.le.nxrad) &
           then
   
      if(intpbl.eq.1) then
      if(((raddatavr(ics,jcs,m)-badpt) > 1.0).and.((raddatavr(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddatavr(ics,jcs,m)*xd2+raddatavr(icsp,jcs,m)*xd1
      else if(((raddatavr(ics,jcs,m)-badpt) < 1.0).and.((raddatavr(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddatavr(icsp,jcs,m)
      else if(((raddatavr(icsp,jcs,m)-badpt) < 1.0).and.((raddatavr(ics,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddatavr(ics,jcs,m)
      else if(((raddatavr(icsp,jcs,m)-badpt) < 1.0).and.((raddatavr(ics,jcs,m)-badpt) < 1.0)) then
         rad_x0=raddatavr(icsp,jcs,m)
      else
         rad_x0=badpt
      end if
      if(((raddatavr(ics,jcsp,m)-badpt) > 1.0).and.((raddatavr(icsp,jcsp,m) -badpt) > 1.0)) then
         rad_x1=raddatavr(ics,jcsp,m)*xd2+raddatavr(icsp,jcsp,m)*xd1
      else if(((raddatavr(ics,jcsp,m)-badpt) < 1.0).and.((raddatavr(icsp,jcsp,m)-badpt) > 1.0)) then
         rad_x1=raddatavr(icsp,jcsp,m)
      else if(((raddatavr(icsp,jcsp,m)-badpt) < 1.0).and.((raddatavr(ics,jcsp,m)-badpt) > 1.0)) then
         rad_x1=raddatavr(ics,jcsp,m)
      else if(((raddatavr(icsp,jcsp,m)-badpt) < 1.0).and.((raddatavr(ics,jcsp,m)-badpt) < 1.0)) then
         rad_x1=raddatavr(icsp,jcsp,m)
      else
         rad_x1=badpt
      end if
   
      if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) > 1.0)) then
         radvr(l,k,m)=rad_x0*yd2+rad_x1*yd1
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) > 1.0)) then
         radvr(l,k,m)=rad_x1
      else if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) < 1.0)) then
         radvr(l,k,m)=rad_x0
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) < 1.0)) then
         radvr(l,k,m)=rad_x0
      else
         radvr(l,k,m)=badpt
      end if

      if(((raddataz(ics,jcs,m)-badpt) > 1.0).and.((raddataz(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddataz(ics,jcs,m)*xd2+raddataz(icsp,jcs,m)*xd1
      else if(((raddataz(ics,jcs,m)-badpt) < 1.0).and.((raddataz(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddataz(icsp,jcs,m)
      else if(((raddataz(icsp,jcs,m)-badpt) < 1.0).and.((raddataz(ics,jcs,m)-badpt) > 1.0)) then
         rad_x0=raddataz(ics,jcs,m)
      else if(((raddataz(icsp,jcs,m)-badpt) < 1.0).and.((raddataz(ics,jcs,m)-badpt) < 1.0)) then
         rad_x0=raddataz(icsp,jcs,m)
      else
         rad_x0=badpt
      end if
      if(((raddataz(ics,jcsp,m)-badpt) > 1.0).and.((raddataz(icsp,jcsp,m)-badpt) > 1.0)) then
         rad_x1=raddataz(ics,jcsp,m)*xd2+raddataz(icsp,jcsp,m)*xd1
      else if(((raddataz(ics,jcsp,m)-badpt) < 1.0).and.((raddataz(icsp,jcsp,m)-badpt) > 1.0)) then
         rad_x1=raddataz(icsp,jcsp,m)
      else if(((raddataz(icsp,jcsp,m)-badpt) < 1.0).and.((raddataz(ics,jcsp,m)-badpt) > 1.0)) then
         rad_x1=raddataz(ics,jcsp,m)
      else if(((raddataz(icsp,jcsp,m)-badpt) < 1.0).and.((raddataz(ics,jcsp,m)-badpt) < 1.0)) then
         rad_x1=raddataz(icsp,jcsp,m)
      else
         rad_x1=badpt
      end if
   
      if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) > 1.0)) then
         radz(l,k,m)=rad_x0*yd2+rad_x1*yd1
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) > 1.0)) then
         radz(l,k,m)=rad_x1
      else if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) < 1.0)) then
         radz(l,k,m)=rad_x0
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) < 1.0)) then
         radz(l,k,m)=rad_x0
      else
         radz(l,k,m)=badpt
      end if


      if(((variac(ics,jcs,m)-badpt) > 1.0).and.((variac(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=variac(ics,jcs,m)*xd2+variac(icsp,jcs,m)*xd1
      else if(((variac(ics,jcs,m)-badpt) < 1.0).and.((variac(icsp,jcs,m)-badpt) > 1.0)) then
         rad_x0=variac(icsp,jcs,m)
      else if(((variac(icsp,jcs,m)-badpt) < 1.0).and.((variac(ics,jcs,m)-badpt) > 1.0)) then
         rad_x0=variac(ics,jcs,m)
      else if(((variac(icsp,jcs,m)-badpt) < 1.0).and.((variac(ics,jcs,m)-badpt) < 1.0)) then
         rad_x0=variac(icsp,jcs,m)
      else
         rad_x0=badpt
      end if

    if(((variac(ics,jcsp,m)-badpt) > 1.0).and.((variac(icsp,jcsp,m)-badpt) > 1.0)) then
         rad_x1=variac(ics,jcsp,m)*xd2+variac(icsp,jcsp,m)*xd1
      else if(((variac(ics,jcsp,m)-badpt) < 1.0).and.((variac(icsp,jcsp,m)-badpt) > 1.0)) then
         rad_x1=variac(icsp,jcsp,m)
      else if(((variac(icsp,jcsp,m)-badpt) < 1.0).and.((variac(ics,jcsp ,m)-badpt) > 1.0)) then
         rad_x1=variac(ics,jcsp,m)
      else if(((variac(icsp,jcsp,m)-badpt) < 1.0).and.((variac(ics,jcsp ,m)-badpt) < 1.0)) then
         rad_x1=variac(icsp,jcsp,m)
      else
         rad_x1=badpt
      end if

      if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) > 1.0)) then
         radsd(l,k,m)=rad_x0*yd2+rad_x1*yd1
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) > 1.0)) then
         radsd(l,k,m)=rad_x1
      else if(((rad_x0-badpt) > 1.0).and.((rad_x1-badpt) < 1.0)) then
         radsd(l,k,m)=rad_x0
      else if(((rad_x0-badpt) < 1.0).and.((rad_x1-badpt) < 1.0)) then
         radsd(l,k,m)=rad_x0
      else
         radsd(l,k,m)=badpt
      end if
   
      else
         npt=0
         som=0.
         do jpt=jcs,jcsp
            do ipt=ics,icsp
               if(((raddatavr(ipt,jpt,m)-badpt) > 1.0)) then
                  npt=npt+1
                  som=som+raddatavr(ipt,jpt,m)
               end if
            end do
         end do
         if(npt.ge.mpt) then
            radvr(l,k,m)=som/float(npt)
         end if
   
         npt=0
         som=0.
         do jpt=jcs,jcsp
            do ipt=ics,icsp
               if(((raddataz(ipt,jpt,m)-badpt) > 1.0)) then
                  npt=npt+1
                  som=som+raddataz(ipt,jpt,m)
               end if
            end do
         end do
         if(npt.ge.mpt) then
            radz(l,k,m)=som/float(npt)
         else
            radz(l,k,m)=badpt
         end if
   
         npt=0
         som=0.
         do jpt=jcs,jcsp
            do ipt=ics,icsp
               if(((variac(ipt,jpt,m)-badpt) > 1.0)) then
                  npt=npt+1
                  som=som+variac(ipt,jpt,m)
               end if
            end do
         end do
         if(npt.ge.mpt) then
            radsd(l,k,m)=som/float(npt)
         else
            radsd(l,k,m)=badpt
         end if
   
      end if
   
      purad(l,k,m)=extrap1(ics,jcs,m)*xd2*yd2 &
                  +extrap1(icsp,jcs,m)*xd1*yd2 &
                  +extrap1(ics,jcsp,m)*xd2*yd1 &
                  +extrap1(icsp,jcsp,m)*xd1*yd1
      pzrad(l,k,m)=extrap2(ics,jcs,m)*xd2*yd2 &
                  +extrap2(icsp,jcs,m)*xd1*yd2 &
                  +extrap2(ics,jcsp,m)*xd2*yd1 &
                  +extrap2(icsp,jcsp,m)*xd1*yd1
   
      endif
   
!   if(radz(l,k,m).eq.badpt1) print*,'radz in badpt1= ',l,k,m,radz(l,k,m)
      enddo
      enddo
      enddo
   
!  print*,' hori inter'
      do m=1,nze_top
         sum_vr(m) = 0.0
         do l = nytsm1,nytep1
         do k = nxtsm1,nxtep1
            if (radvr(l,k,m) .ne. badpt) then
          sum_vr(m) = sum_vr(m) + radvr(l,k,m) 
            endif
         enddo
         enddo
      enddo
!
!     Vertical linear interpolation
!
   if(iproj_type.eq.2) then
!!$omp parallel do default(shared) &
!!$omp private(l,k,m,zd,zd1,zd2)
     do m=nzts,nzte
       zd=real(m)*dz/dzr+zmindx
       mcs=int(zd+0.00001)
       zd1=zd-float(mcs)
       zd2=1.0-zd1

       if(mcs.eq.0) then
       do k=nxtsm1,nxtep1
       do l=nytsm1,nytep1
        if(((radvr(l,k,m)-badpt) > 1.0).and.((radvr(l,k,m+1)-badpt) > 1.0)) then
           urs(l,k,m,ivol,irad)=(2.*radvr(l,k,m)-radvr(l,k,m+1))*zd2+ &
                    radvr(l,k,m)*zd1
        elseif(((radvr(l,k,m)-badpt) > 1.0).and.((radvr(l,k,m+1)-badpt) < 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,m+1)
        elseif(((radvr(l,k,m)-badpt) < 1.0).and.((radvr(l,k,m+1)-badpt) > 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,m)
        elseif(((radvr(l,k,m)-badpt) < 1.0).and.((radvr(l,k,m+1)-badpt) < 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,m+1)
        else
           urs(l,k,m,ivol,irad)=badpt
        end if

        if(((radz(l,k,m)-badpt) > 1.0).and.((radz(l,k,m+1)-badpt) > 1.0)) then
           cs(l,k,m,ivol,irad)=(2.*radz(l,k,m)-radz(l,k,m+1))*zd2+ &
                    radz(l,k,m)*zd1
        elseif(((radz(l,k,m)-badpt) > 1.0).and.((radz(l,k,m+1)-badpt) < 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,m+1)
        elseif(((radz(l,k,m)-badpt) < 1.0).and.((radz(l,k,m+1)-badpt) > 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,m)
        elseif(((radz(l,k,m)-badpt) < 1.0).and.((radz(l,k,m+1)-badpt) < 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,m+1)
        else
           cs(l,k,m,ivol,irad)=badpt
        end if

!yzm    if(((radsd(l,k,m)-badpt) > 1.0).and.((radsd(l,k,m+1)-badpt) > 1.0)) then
!          std(l,k,m,ivol,irad)=(2.*radsd(l,k,m)-radsd(l,k,m+1))*zd2+ &
!                   radsd(l,k,m)*zd1
!       else
!          std(l,k,m,ivol,irad)=badpt
!       end if

        pu(l,k,m,ivol,irad)=(2.*purad(l,k,m)-purad(l,k,m+1))*zd2+ &
                    purad(l,k,m)*zd1
        pz(l,k,m,ivol,irad)=(2.*pzrad(l,k,m)-pzrad(l,k,m+1))*zd2+ &
                    pzrad(l,k,m)*zd1
      enddo
       enddo
       else
       do k=nxtsm1,nxtep1
       do l=nytsm1,nytep1
        if(((radvr(l,k,mcs)-badpt) > 1.0).and.((radvr(l,k,mcs+1)-badpt) > 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,mcs)*zd2+radvr(l,k,mcs+1)*zd1
        elseif(((radvr(l,k,mcs)-badpt) > 1.0).and.((radvr(l,k,mcs+1)-badpt) < 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,mcs+1)
        elseif(((radvr(l,k,mcs)-badpt) < 1.0).and.((radvr(l,k,mcs+1)-badpt) > 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,mcs)
        elseif(((radvr(l,k,mcs)-badpt) < 1.0).and.((radvr(l,k,mcs+1)-badpt) < 1.0)) then
           urs(l,k,m,ivol,irad)=radvr(l,k,mcs+1)
        else
           urs(l,k,m,ivol,irad)=badpt
        end if

        if(((radz(l,k,mcs)-badpt) > 1.0).and.((radz(l,k,mcs+1)-badpt) > 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,mcs)*zd2+radz(l,k,mcs+1)*zd1
        elseif(((radz(l,k,mcs)-badpt) > 1.0).and.((radz(l,k,mcs+1)-badpt) < 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,mcs+1)
        elseif(((radz(l,k,mcs)-badpt) < 1.0).and.((radz(l,k,mcs+1)-badpt) > 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,mcs)
        elseif(((radz(l,k,mcs)-badpt) < 1.0).and.((radz(l,k,mcs+1)-badpt) < 1.0)) then
           cs(l,k,m,ivol,irad)=radz(l,k,mcs+1)
        else
           cs(l,k,m,ivol,irad)=badpt
        end if

!       if(((radsd(l,k,mcs)-badpt) > 1.0).and.((radsd(l,k,mcs+1)-badpt) > 1.0)) then
!          std(l,k,m,ivol,irad)=radsd(l,k,mcs)*zd2+radsd(l,k,mcs+1)*zd1
!       else
!          std(l,k,m,ivol,irad)=badpt
!       endif

        pu(l,k,m,ivol,irad)=purad(l,k,mcs)*zd2+purad(l,k,mcs+1)*zd1
        pz(l,k,m,ivol,irad)=pzrad(l,k,mcs)*zd2+pzrad(l,k,mcs+1)*zd1
        enddo
        enddo
        end if

        do k=nxtsm1,nxtep1
        do l=nytsm1,nytep1
           urscart(l,k,m,ivol,irad)=urs(l,k,m,ivol,irad)
           cscart(l,k,m,ivol,irad)=cs(l,k,m,ivol,irad)

           pucart(l,k,m,ivol,irad)=pu(l,k,m,ivol,irad)        !yzm
           pzcart(l,k,m,ivol,irad)=pz(l,k,m,ivol,irad)        !yzm
        enddo
        enddo
      enddo
   
   else
   
!
!     The PPI data are stored in URS & CS. The interpolated Cartesian data
!     are stored in URSCART and CSCART
!
   
   if(iproj_type.eq.1) then
!!$omp parallel do default(shared) &
!!$omp private(l,k,m)
   do m=1,nze_top
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
      urs(l,k,m,ivol,irad)=radvr(l,k,m)
      cs(l,k,m,ivol,irad)=radz(l,k,m)
      std(l,k,m,ivol,irad)=radsd(l,k,m)
      pu(l,k,m,ivol,irad)=purad(l,k,m)
      pz(l,k,m,ivol,irad)=pzrad(l,k,m)
   enddo
   enddo
   enddo
   endif
!  do m=1,nze_top
!  print*,'pu,pz,IVOL,IRAD',ivol,irad, &
!        sum(pu(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad)), &
!        sum(pz(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad))
!  enddo

   do m=1,nze_top
      sum_vr(m) = 0.0
      do l = nytsm1,nytep1
      do k = nxtsm1,nxtep1
         if(urs(l,k,m,ivol,irad) .ne. badpt) then
            sum_vr(m) = sum_vr(m) + urs(l,k,m,ivol,irad) 
         endif
      enddo
      enddo
   enddo
   
!yzm---add data interpolation from PPI to cartesian for PU and PZ-   

   do m=nzts,nzte
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
      if(sqrt(((float(k)-xk(irad))*dx)**2 &
           +((float(l)-yl(irad))*dy)**2).lt.6000.) then
         urscart(l,k,m,ivol,irad)=badpt
         cscart(l,k,m,ivol,irad)=badpt
         pucart(l,k,m,ivol,irad)=0.                !yzm 
         pzcart(l,k,m,ivol,irad)=0.                !yzm 
      else
         if(zcart(m).lt.zppi(l,k,1,1,irad)) then
            urscart(l,k,m,ivol,irad)=radvr(l,k,1)
            cscart(l,k,m,ivol,irad)=radz(l,k,1)
            pucart(l,k,m,ivol,irad)=purad(l,k,1)        !yzm 
            pzcart(l,k,m,ivol,irad)=pzrad(l,k,1)        !yzm 
         elseif(zcart(m).ge.zppi(l,k,1,1,irad).and. &
                 zcart(m).lt.zppi(l,k,1,2,irad)) then
            if(((radvr(l,k,1)-badpt) > 1.0)) then
               urscart(l,k,m,ivol,irad)=radvr(l,k,1)
               pucart(l,k,m,ivol,irad)=purad(l,k,1)     !yzm 
            else
               urscart(l,k,m,ivol,irad)=badpt
               pucart(l,k,m,ivol,irad)=0.               !yzm 
            end if
            if(((radz(l,k,1)-badpt) > 1.0)) then
               cscart(l,k,m,ivol,irad)=radz(l,k,1)
               pzcart(l,k,m,ivol,irad)=pzrad(l,k,1)        !yzm 
            else
               cscart(l,k,m,ivol,irad)=badpt
               pzcart(l,k,m,ivol,irad)=0.               !yzm 
            end if
         else
            mppi_loop: do mppi=1,nztop-1
               if(zcart(m).ge.zppi(l,k,mppi,2,irad).and.zcart(m).lt. &
                    zppi(l,k,mppi+1,2,irad)) then
                  if(((radvr(l,k,mppi)-badpt) > 1.0).and.((radvr(l,k,mppi+1)-badpt) > 1.0)) then
                     urscart(l,k,m,ivol,irad)=radvr(l,k,mppi)+ &
                          (zcart(m)-zppi(l,k,mppi,2,irad))/ &
                          (zppi(l,k,mppi+1,2,irad)-zppi(l,k,mppi,2,irad)) &
                          *(radvr(l,k,mppi+1)-radvr(l,k,mppi))
                     pucart(l,k,m,ivol,irad)=purad(l,k,mppi)+ &
                          (zcart(m)-zppi(l,k,mppi,2,irad))/ &
                          (zppi(l,k,mppi+1,2,irad)-zppi(l,k,mppi,2,irad)) &
                          *(purad(l,k,mppi+1)-purad(l,k,mppi))                    !yzm
                  else
                     urscart(l,k,m,ivol,irad)=badpt
                     pucart(l,k,m,ivol,irad)=0.                     !yzm 
                  end if
                  if(((radz(l,k,mppi)-badpt) > 1.0).and.((radz(l,k,mppi+1)-badpt) > 1.0)) then
                     cscart(l,k,m,ivol,irad)=radz(l,k,mppi)+ &
                          (zcart(m)-zppi(l,k,mppi,2,irad))/ &
                          (zppi(l,k,mppi+1,2,irad)-zppi(l,k,mppi,2,irad)) &
                          *(radz(l,k,mppi+1)-radz(l,k,mppi))
                     pzcart(l,k,m,ivol,irad)=pzrad(l,k,mppi)+ &
                          (zcart(m)-zppi(l,k,mppi,2,irad))/ &
                          (zppi(l,k,mppi+1,2,irad)-zppi(l,k,mppi,2,irad)) &
                          *(pzrad(l,k,mppi+1)-pzrad(l,k,mppi))                  !yzm
                  else
                     cscart(l,k,m,ivol,irad)=badpt
                     pzcart(l,k,m,ivol,irad)=0.                     !yzm 
                  end if
                  exit mppi_loop
               end if
             enddo mppi_loop
         end if

         if(zcart(m).ge.zppi(l,k,nztop,2,irad).and. &
              zcart(m).le.zppi(l,k,nztop,3,irad)) then
            if(((radvr(l,k,nztop)-badpt) > 1.0)) then
               urscart(l,k,m,ivol,irad)=radvr(l,k,nztop)
                pucart(l,k,m,ivol,irad)=purad(l,k,nztop)      !yzm
            else
               urscart(l,k,m,ivol,irad)=badpt
                pucart(l,k,m,ivol,irad)=0.                     !yzm 
            end if
            if(((radz(l,k,nztop)-badpt) > 1.0)) then
               cscart(l,k,m,ivol,irad)=radz(l,k,nztop)
               pzcart(l,k,m,ivol,irad)=pzrad(l,k,nztop)        !yzm
            else
               cscart(l,k,m,ivol,irad)=badpt
               pzcart(l,k,m,ivol,irad)=0.            !yzm
            end if
         end if
      end if
   enddo
   enddo
   enddo
   

!  print*,'pu=',ivol,irad,pu(nyte/2,nxte/2,1:nze,ivol,irad) 
!  print*,'pz=',ivol,irad,pz(nyte/2,nxte/2,1:nze,ivol,irad) 
!  do m=1,nze_top
!  print*,'sum of pu,pz=',ivol,irad, &
!        sum(pu(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad)), &
!        sum(pz(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad))
!  enddo
!  do m=nzts,nzte
!  print*,'sum of pucart,pzcart=',ivol,irad, &
!        sum(pucart(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad)), &
!        sum(pzcart(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad))
!  enddo
!  print*,'pucart=',ivol,irad,pucart(nyte/2,nxte/2,nzts:nzte,ivol,irad) 
!  print*,'pzcart=',ivol,irad,pzcart(nyte/2,nxte/2,nzts:nzte,ivol,irad) 
    

   end if
   
   icount=0
   do m=1,nztop
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
      if(zppi(l,k,m,2,irad).le.zcart(nzde)) then
         if((urs(l,k,m,ivol,irad) - badpt) < 1.0) cycle
         if((urs(l,k,m,ivol,irad)).gt.30.) then
            icount = icount + 1
         end if
      end if
   enddo
   enddo
   enddo
   

!   call vdras_sum(icount)

   print*
   print*,'NUMBER OF URS GREATER THAN 30 M/S: ',icount
   
   icount=0
   do m=nzts,nzte
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
      if((urscart(l,k,m,ivol,irad)).gt.30..and. &
         ((urscart(l,k,m,ivol,irad)-badpt) > 1.0)) then
         icount = icount + 1
      end if
   enddo
   enddo
   enddo
   
!   call vdras_sum(icount)
   
   print*
   print*,'NUMBER OF URSCART GREATER THAN 30 M/S: ',icount
   print*
   
!  do m=1,nze_top
!  do k=nxtsm1,nxtep1
!  do l=nytsm1,nytep1
!      disth=sqrt(((float(k)-xk(irad))*dx)**2 &
!            +((float(l)-yl(irad))*dy)**2)
!.. IN REGIONS OF NO RADAR RETURN WITHIN THE RANGE OF THE RADAR,
!.. SET QR TO ZERO INSTEAD OF BAD.
!     if(((cs(l,k,m,ivol,irad)-badpt) < 1.0).and.disth.le.225000.) then
!     cs(l,k,m,ivol,irad)=-15.0
!     cs(l,k,m,ivol,irad)=0.0
!     pz(l,k,m,ivol,irad)=0.0
!     end if
!yzm   if(k.le.10) then
!       urs(l,k,m,ivol,irad)=badpt
!      endif
!  end do
!  end do
!  end do
   

!  stdmax=badpt
!  do m=1,nze_top
!  do k=nxtsm1,nxtep1
!  do l=nytsm1,nytep1
!      if(std(l,k,m,ivol,irad).gt.stdmax) then
!      stdmax=std(l,k,m,ivol,irad)
!      end if
!  end do
!  end do
!  end do

!  call vdras_maxval(stdmax)

!  ranged=20000.
   ranged=70000.

   do m=1,nze_top
   do k=nxtsm1,nxtep1
   do l=nytsm1,nytep1
       disth=sqrt(((float(k)-xk(irad))*dx)**2 &
             +((float(l)-yl(irad))*dy)**2)
       if(disth.le.ranged) then
       stdvr=(ranged-disth)/ranged*2.0+1.
       stdqr=(ranged-disth)/ranged*2.0+1.
       else
       stdvr=1.
       stdqr=1.
       end if
       pu(l,k,m,ivol,irad)=pu(l,k,m,ivol,irad)/stdvr**2
       pr(l,k,m,ivol,irad)=pz(l,k,m,ivol,irad)/stdqr**2
!      if(((std(l,k,m,ivol,irad)-badpt) > 1.0)) then
!      if(std(l,k,m,ivol,irad).gt.0.5) then
!      pu(l,k,m,ivol,irad)=pu(l,k,m,ivol,irad)/ &
!           std(l,k,m,ivol,irad)**2
!      else
!      pu(l,k,m,ivol,irad)=pu(l,k,m,ivol,irad)/0.25
!      end if
!      end if
!      pr(l,k,m,ivol,irad)=pz(l,k,m,ivol,irad)
   end do
   end do
   end do


!  do m=1,nze_top
!    print*,'pu,pr,ivol,irad',ivol,irad, &
!    sum(pu(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad)), &
!    sum(pr(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad)), &
!    sum(urs(nytsm1:nytep1,nxtsm1:nxtep1,m,ivol,irad))
!  enddo
   
!-----data thining
   
!      print *,'before thinning'
!      do l=0,nyp1
!      do k=0,nx
!      do m=1,nz,2
!          DISTH=SQRT(((FLOAT(K)-XK(IRAD))*DX)**2
!     *          +((FLOAT(L)-YL(IRAD))*DY)**2)
!         if(disth.le.30000.) then
!         dzppimin=dz*2
!         do mm=1,nze
!         dzppi=zcart(m)-zppi(l,k,mm,2,irad)
!         if(dzppi.lt.dzppimin.and.dzppi.ge.0.) then
!         urstmp=urs(l,k,mm,ivol,irad)
!         urs(l,k,mm,ivol,irad)=badpt
!         cstmp=cs(l,k,mm,ivol,irad)
!         cs(l,k,mm,ivol,irad)=badpt
!         dzppimin=dzppi
!         mz=mm
!         end if
!         end do
!         urs(l,k,mz,ivol,irad)=urstmp
!         cs(l,k,mz,ivol,irad)=cstmp
!         end if
!      end do
!      end do
!      end do
   
   end subroutine pre_analyze

end module
