module libradar
implicit none

real(kind=8), parameter :: PI         =3.141592653589793238 
real(kind=8), parameter :: EARTHRADIUS=6371000.
real(kind=8), parameter :: DEG2RAD    =(PI/180.)
real(kind=8), parameter :: RAD2DEG    =(180./PI)
real(kind=8), parameter :: FOURTHIRDE =(4.*EARTHRADIUS/3.)
real(kind=8), parameter :: EIGHTTHIRDE=(8.*EARTHRADIUS/3.)
real(kind=8), parameter :: FOURTHIRDSQ=(FOURTHIRDE*FOURTHIRDE)
real(kind=8), parameter :: ONE        =1.
real(kind=8), parameter :: MINUSONE   =-1.


contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE xy2rd(x, y, r, d )
      implicit none
      !Catesian Coordinate to Polar Coordinate
      real, intent(in ) :: x, y
      real, intent(out) :: r, d

      r=SQRT(x*x+y*y)

      IF( x == 0.) THEN
        IF( y >= 0.) THEN
          d=0.
        ELSE
          d=180.
        END IF
      ELSE
        if(y==0)then
           if(x>0)then
              d=90.
           else
              d=270.
           endif
        else
           d=rad2deg*atan(x/y)
        endif
        IF(y < 0.) d=d+180.
        IF(d < 0.) d=d+360.
      END IF

      END SUBROUTINE

      SUBROUTINE rd2xy( r, d, x, y)
      implicit none
      !Polar Coordinate to Catesian Coordinate
      real, intent(in ) :: r, d
      real, intent(out) :: x, y

      x=r*cos(deg2rad*(90-d)) 
      y=r*sin(deg2rad*(90-d)) 

      END SUBROUTINE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE dhdrange(elvang,range,dhdr)
   !
   !-----------------------------------------------------------------------
   !
   !  PURPOSE:
   !
   !  Calculate the local change in height of the radar
   !  beam with respect to a change in range.  Due to
   !  curvature of the beam and the earth's surface this is
   !  generally different what would be calculated from the
   !  elevation angle measured at the radar.  This derivative
   !  is needed for finding 3-d velocities from radial winds
   !  and accounting for terminal velocity of precipitation.
   !
   !  This formulation, consistent with subroutine beamhgt,
   !  assumes a 4/3 earth radius beam curvature.  This formula
   !  is obtained by differentiating Eq 2.28 of Doviak and
   !  Zrnic', Doppler Radar and Weather Observations, 1st Ed.
   !
   !-----------------------------------------------------------------------
   !
   !  INPUT:
   !
   !    elvang   Elevation angle (degrees) of radar beam
   !    range    Distance (meters) along radar beam from radar
   !
   !  OUTPUT:
   !    dhdr     Change in height per change in range (non-dimensional)
   !
   !-----------------------------------------------------------------------

   IMPLICIT NONE
   REAL, INTENT(IN) :: range
   REAL, INTENT(IN) :: elvang
   REAL, INTENT(OUT) :: dhdr
   !
   real(kind=8) :: sinelv,dhdrdb,drange
   !
   drange = DBLE(range)
   sinelv = SIN(DEG2RAD*DBLE(elvang))
   dhdrdb = (drange+FOURTHIRDE*sinelv)/                                     &
            SQRT(drange*drange + FOURTHIRDSQ + EIGHTTHIRDE*drange*sinelv)
   dhdr   = dhdrdb
   !
   END SUBROUTINE dhdrange

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine beamrng(elvang,sfcrng,height,range)
   implicit none
   real, intent(in)  :: elvang, sfcrng
   real, intent(out) :: height, range 

   real(kind=8) :: elvrad,hgtdb,rngdb,drange

   elvrad = DEG2RAD*elvang
   rngdb  = sfcrng/FOURTHIRDE
   hgtdb  = FOURTHIRDE/(cos(rngdb)-(sin(rngdb)*tan(elvrad)))
   drange = hgtdb*sin(rngdb)/cos(elvrad)
   height = hgtdb-FOURTHIRDE
   range  = drange

   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine beamhgt(elvang,range,height,sfcrng)
   implicit none
   real, intent(in)  :: elvang, range
   real, intent(out) :: height, sfcrng

   real(kind=8) :: elvrad,hgtdb,rngdb,drange
   
   elvrad = DEG2RAD*elvang
   drange = range
   hgtdb  = sqrt(drange*drange + FOURTHIRDSQ + EIGHTTHIRDE*drange*sin(elvrad)) - FOURTHIRDE
   height = hgtdb
   rngdb  = FOURTHIRDE * asin (drange*cos(elvrad)/(FOURTHIRDE + hgtdb) )
   sfcrng = rngdb
   
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine beamelv(height,sfcrng,elvang,range)
   implicit none
   real, intent(in)  :: height, sfcrng
   real, intent(out) :: elvang, range

   real(kind=8) :: elvrad,hgtdb,rngdb,drange

   if(sfcrng > 0.) then
      hgtdb  = FOURTHIRDE+height
      rngdb  = sfcrng/FOURTHIRDE

      elvrad = atan((hgtdb*cos(rngdb) -FOURTHIRDE )/(hgtdb * sin(rngdb)))
      drange = (hgtdb*sin(rngdb))/cos(elvrad)
      elvang = RAD2DEG*elvrad
      range  = drange
   else
      elvang = 90.
      range  = height
   endif
  
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine disthead(lat1,lon1,lat2,lon2,headng,dist)
   implicit none
   real, intent(in)  :: lat1, lon1, lat2, lon2
   real, intent(out) :: headng, dist

   real(kind=8) :: alat1,alat2,dlon,arcdst,cosdst,coshd,denom

   alat1  = DEG2RAD * lat1
   alat2  = DEG2RAD * lat2
   dlon   = DEG2RAD * (lon2-lon1)
   cosdst = sin(alat1) * sin(alat2) + cos(alat1) * cos(alat2) * cos(dlon)
   arcdst = acos(cosdst)
   dist   = EARTHRADIUS * arcdst
 
   denom  = cos(alat1)*sin(arcdst)
   headng=0.
   if(ABS(denom) > 1.0E-6) then
     coshd  = (sin(alat2) - sin(alat1)*cosdst) / denom
     coshd  = MAX(coshd,MINUSONE)
     coshd  = MIN(coshd,ONE)
     headng = RAD2DEG*acos(coshd)
     if ( sin(dlon) < 0 ) headng = 360.-headng
   elseif ( ABS(cos(alat1)) < 1.0E-6 .and. alat1 > 0.) then
     headng=180.
   endif 
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gcircle(lat1,lon1,head,dist,lat2,lon2)
   implicit none
   real, intent(in)  :: lat1, lon1, head, dist
   real, intent(out) :: lat2, lon2

   real(kind=8) :: alat1,alat2,dlon,arcdst,cosdst,coshd
   real(kind=8) :: denom,sinlat2,cosdlon
 
   alat1  = DEG2RAD*lat1
   arcdst = dist/EARTHRADIUS
   cosdst = cos(arcdst)
   coshd  = cos(DEG2RAD*head)
   sinlat2= coshd*cos(alat1)*sin(arcdst) + sin(alat1)*cosdst
   sinlat2= MAX(sinlat2,MINUSONE)
   sinlat2= MIN(sinlat2,ONE)
   alat2  = asin(sinlat2)
   lat2   = RAD2DEG*alat2
 
   denom=cos(alat1)*cos(alat2)
   if(denom /= 0.) then
     cosdlon = (cosdst - sin(alat1)*sinlat2)/(cos(alat1)*cos(alat2))
     cosdlon = MAX(cosdlon,MINUSONE)
     cosdlon = MIN(cosdlon,ONE)
     dlon    = RAD2DEG*acos(cosdlon)
     if(sin(DEG2RAD*head) < 0.) dlon=-dlon
     lon2=lon1+dlon
   else
     lon2=lon1
   endif
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine xy2rd(X,Y,radius,head)
   !implicit none
   !real, intent(in)  :: X, Y
   !real, intent(out) :: radius, head

   !real(kind=8) :: tana,alpha,degree

   !radius = sqrt((X*X)+(Y*Y))
   !tana   = abs(X/Y)
   !alpha  = atan(tana)
   !degree = RAD2DEG*alpha
   !if(X>0.and.Y>0)then
   !   degree=degree
   !elseif(X>0.and.Y<0)then
   !   degree=180.-degree
   !elseif(X<0.and.Y>0)then
   !   degree=360.-degree
   !elseif(X<0.and.Y<0)then
   !   degree=180.+degree
   !else
   !   degree=0.
   !endif 
   !head=degree
   !end subroutine


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function bilinear(v1,v2,v3,v4,x,y)
   implicit none
   real, intent(in) :: v1, v2, v3, v4, x, y

   real :: p1,p2

   p1=(v2-v1)*x+v1
   p2=(v4-v3)*x+v3
   bilinear=(p2-p1)*y+p1
   end function


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function trilinear(v1,v2,v3,v4,v5,v6,v7,v8,x,y,z)
   implicit none
   real, intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, x, y, z
   
   real :: p1,p2

   p1=bilinear(v1,v2,v3,v4,x,y)
   p2=bilinear(v5,v6,v7,v8,x,y)
   trilinear=(p2-p1)*z+p1
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE INTERP2D ( NX, NY, FIELD, NUM, XI, YJ, DATA_OUT, MISS )
   IMPLICIT NONE
!      
      INTEGER,                   INTENT(IN)  :: NUM, NX, NY
      REAL, DIMENSION( NX , NY ),INTENT(IN)  :: FIELD
      REAL,                      INTENT(IN)  :: MISS
      REAL, DIMENSION( NUM ),    INTENT(IN)  :: XI , YJ
      REAL, DIMENSION( NUM ),    INTENT(OUT) :: DATA_OUT
!
      REAL                                   :: X , Y
      INTEGER                                :: K , IS , JS , IE , JE


      !WRITE(*,*) FIELD
      DATA_OUT = MISS
      DO K =1, NUM
         IF( XI ( K ) > NX .OR. XI ( K ) <1 .OR. &
             YJ ( K ) > NY .OR. YJ ( K ) <1 )THEN
             CYCLE
         ENDIF
         IS =INT( XI ( K ) )
!        WRITE(*,*) XI(K), IS
         IF( IS >= NX ) IS = NX -1
         IF( IS < 1) IS =1
         X = XI ( K )-FLOAT( IS )
!        WRITE(*,*) XI(K), IS
         IE = IS +1

         JS =INT( YJ ( K ))
         IF( JS >= NY ) JS = NY -1
         IF( JS < 1) JS =1
         Y = YJ ( K )-FLOAT( JS )
         JE = JS +1

         DATA_OUT ( K )= (1.0- X )*(1.0- Y )* FIELD ( IS , JS )+ &
                         (1.0- X )*      Y  * FIELD ( IS , JE )+ &
                               X  *(1.0- Y) * FIELD ( IE , JS )+ &
                               X  *      Y  * FIELD ( IE , JE )
      ENDDO

   END SUBROUTINE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine uv2vr( range, elev , azimu , ugrid , vgrid , vr)
 !
 !-----------------------------------------------------------------------
 !
 !  PURPOSE:
 !   u,v to radial direction
 !
 !-----------------------------------------------------------------------
 !
 !  INPUT:
 !    range      distance to radar location
 !    elev       elevation angle
 !    azimu      azimuth angle
 !    ugrid      u
 !    vgrid      v
 !
 !  OUTPUT:
 !    vr         radial wind
 !
 !-----------------------------------------------------------------------
 !
   IMPLICIT NONE
 !
 !-----------------------------------------------------------------------
 !
 !  Variable Declarations:
 !
 !-----------------------------------------------------------------------
 !
   real, intent(in)  :: range, elev, azimu
   real, intent(in)  :: ugrid,vgrid
   real, intent(out) :: vr
 !
   real :: ange,dhdr,dsdr
   real :: uazmrad,vazmrad
   real :: DPI
   parameter ( DPI=3.1415926535898/180.0 )
 
   uazmrad   = -SIN( azimu*DPI )
   vazmrad   = -COS( azimu*DPI )
 
   CALL dhdrange(elev,range,dhdr)
   dsdr=SQRT(AMAX1(0.,(1.-dhdr*dhdr)))

   vr=(uazmrad*ugrid + vazmrad*vgrid) * dsdr
 !
   end subroutine uv2vr

! Example for using above function:
! 
! program example 
! use libradar
! implicit none
!  real :: elev,range
!  real :: sfcrng,height
!  real :: lon1,lon2,lat1,lat2
!  real :: azim,dist,head
! 
!  lat1=35.
!  lon1=-100.
!  azim=30.
! 
!  height=2000.
!  sfcrng=100000.
!  call beamelv(height,sfcrng,elev,range)
!  write(*,*) height,sfcrng,elev,range
!  call beamhgt(elev,range,height,sfcrng)
!  write(*,*) height,sfcrng,elev,range 
!  call gcircle(lat1,lon1,azim,sfcrng,lat2,lon2)
!  write(*,*) lat1,lon1,azim,sfcrng,lat2,lon2
!  call disthead(lat1,lon1,lat2,lon2,head,dist)
!  write(*,*) lat1,lon1,head,dist,lat2,lon2
!

!
recursive subroutine quick_sort(a,d,n,s,e) ! little - big
implicit none
  integer :: n    !
  real    :: a(n) !
  integer :: d(n) !
  integer :: s    !
  integer :: e    !
  integer :: l,r  !

  real    :: k    ! 
  real    :: temp ! 
  !
  l=s 
  r=e+1
  ! right > left
  if ( r<=l ) return
  k=a(s)  !
  do while(.true.)
    ! a(l)<k
    do while( .true. )
      l=l+1
      if ( (l>=e) ) exit
      if ( (a(l) > k) ) exit
    end do
    ! a(r)>k
    do while( .true. )
      r=r-1
      if ( (r<=s) ) exit
      if ( (a(r) < k) ) exit
    end do
    ! right left
    if ( r <= l ) exit
    ! a(l),a(r)
    temp=a(l)
    a(l)=a(r)
    a(r)=temp
    temp=d(l)
    d(l)=d(r)
    d(r)=temp
  end do
  ! a(s),a(r)
  temp=a(s)
  a(s)=a(r)
  a(r)=temp
  temp=d(s)
  d(s)=d(r)
  d(r)=temp
  ! r
  call quick_sort(a,d,n,s,r-1)
  ! r
  call quick_sort(a,d,n,r+1,e)
  return
end subroutine quick_sort
 
end module

