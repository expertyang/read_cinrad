module grd
implicit none

interface write_grd
   module procedure write_grd_real, write_grd_int, write_grd_int8, write_grd_real8
end interface

interface write_test_data
   module procedure write_test_data_real, write_test_data_int
end interface

character(len=20), private :: ctime =""
character(len=6),  private :: clevel=""
character(len=6),  private :: cid   =""

interface set_file
   module procedure set_file_char, set_file_int
end interface
contains

subroutine set_file_char(cname,cvar)
implicit none
character(len=*), intent(in) :: cname, cvar

if(trim(cname)=="time")then
   ctime=cvar
elseif(trim(cname)=="level")then
   clevel=cvar
elseif(trim(cname)=="id")then
   cid=cvar
endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_file_int(cname,ivar)
implicit none
character(len=*), intent(in) :: cname
integer, intent(in) :: ivar

if(trim(cname)=="level")then
   write(clevel,"(I2.2)") ivar
endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_filename(cname,filename,extname)
implicit none
character(len=*), intent(in ) :: cname
character(len=*), intent(out) :: filename
character(len=*), optional, intent(in ) :: extname

if(len_trim(ctime)/=0)then
   filename=trim(ctime)
endif
if(len_trim(cid)/=0)then
   if(len_trim(filename)/=0)then
      filename=trim(filename)//"_"
   endif
   filename=trim(filename)//trim(cid)
endif
if(len_trim(clevel)/=0)then
   if(len_trim(filename)/=0)then
      filename=trim(filename)//"_"
   endif
   filename=trim(filename)//trim(clevel)
endif
if(len_trim(cname)/=0)then
   if(len_trim(filename)/=0)then
      filename=trim(filename)//"_"
   endif
   filename=trim(filename)//trim(cname)
endif

if(present(extname)) filename=trim(filename)//"."//trim(extname)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_diff_data(nr,na,dat,rdiff,adiff,miss)
   implicit none
   integer, intent(in) :: nr, na
   real, dimension(nr,na), intent(in) :: dat
   real, dimension(nr,na), intent(out) :: rdiff, adiff
   real, intent(in) :: miss

   integer :: i, i1, j, j1, k

   rdiff=miss
   adiff=miss
   do j=1, na
      do i=1, nr
         i1=i+1
         if(i1<=nr) then
            if(dat(i,j)/=miss.and.dat(i1,j)/=miss)then
               rdiff(i,j)=dat(i1,j)-dat(i,j)
            endif 
         endif 
         j1=j+1
         if(j1>na)j1=1
         if(dat(i,j)/=miss.and.dat(i,j1)/=miss)then
            adiff(i,j)=dat(i,j1)-dat(i,j)
         endif 
      enddo
   enddo
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_diff_data0(nr,na,dat,minv,maxv,rdiff,adiff,miss)
   implicit none
   integer, intent(in) :: nr, na
   real, dimension(nr,na), intent(in) :: dat
   real, dimension(nr,na), intent(out) :: rdiff, adiff
   real, intent(in) :: miss, minv, maxv

   integer :: i, i1, j, j1, k

   rdiff=miss
   adiff=miss
   do j=1, na
      do i=1, nr
         if(dat(i,j)>=minv.and.dat(i,j)<maxv)then
         i1=i+1
         if(i1<=nr) then
            if(dat(i,j)/=miss.and.dat(i1,j)/=miss)then
               rdiff(i,j)=dat(i1,j)-dat(i,j)
            endif
         endif
         j1=j+1
         if(j1>na)j1=1
         if(dat(i,j)/=miss.and.dat(i,j1)/=miss)then
            adiff(i,j)=dat(i,j1)-dat(i,j)
         endif
         endif
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_diff_flag_data(nr,na,flag,dat,rdiff,adiff,miss)
   implicit none
   integer, intent(in) :: nr, na
   real, dimension(nr,na), intent(in) :: dat
   integer, dimension(nr,na), intent(in) :: flag 
   real, dimension(nr,na), intent(out) :: rdiff, adiff
   real, intent(in) :: miss

   integer :: i, i1, j, j1, k

   rdiff=miss
   adiff=miss
   do j=1, na
      do i=1, nr
         i1=i+1
         if(i1<=nr) then
            if(dat(i,j)/=miss.and.dat(i1,j)/=miss.and.flag(i,j)==1.and.flag(i1,j)==1)then
               rdiff(i,j)=dat(i1,j)-dat(i,j)
            endif
         endif
         j1=j+1
         if(j1>na)j1=1
         if(dat(i,j)/=miss.and.dat(i,j1)/=miss.and.flag(i,j)==1.and.flag(i,j1)==1)then
            adiff(i,j)=dat(i,j1)-dat(i,j)
         endif
      enddo
   enddo
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_dist_data(filein,nr,na,dat)
   implicit none
   character(len=*), intent(in) :: filein
   integer, intent(in) :: nr, na
   real, dimension(nr,na), intent(in) :: dat

   character(len=8000) :: filename
   integer :: i, j, idx
   integer, parameter :: nd=401
   real, parameter :: dd=0.5, sd=-100.
   integer, dimension(nd) :: dist
   character(len=20) :: c_format
  
   call get_filename(filein,filename)

   open(21,file=filename,form="formatted",status="unknown")
   write(*,*) "writing file:", trim(filename), nr, na

   dist=0
   do j=1, na
      do i=1, nr
         idx=floor((dat(i,j)-sd)/dd+0.5)+1
         if(idx>=1.and.idx<=nd)then
            dist(idx)=dist(idx)+1
         endif
      enddo
   enddo
   write(c_format,"(A,I3,A)") "(I10,",nd-1,"(',',I10))"
   write(21,FMT=c_format) dist
   close(21)
   end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine interp_radar_2d(nr1,na1,azim1,dis1,dat1,nr2,na2,azim2,dis2,dat2,miss,imeth)
   implicit none
   integer,                     intent(in) :: nr1, na1, nr2, na2, imeth
   real,                        intent(in) :: miss, dis1, dis2
   real,    dimension(nr1,na1), intent(in) :: dat1
   real,    dimension(nr2,na2), intent(out):: dat2
   real,    dimension(na1),     intent(in) :: azim1
   real,    dimension(na2),     intent(in) :: azim2

   integer :: i1, j1, i2, j2, i3, j3, ns, ne

   real,    dimension(na2) :: jy
   integer, dimension(na2) :: js, je
   integer, dimension(1) :: loc1, loc2

   real :: minv, maxv, ix
   real :: d1, d2, d3, d4, r0, r1
   real :: w1, w2, w3, w4 
! imeth=1,min; imeth=2, max; imeth=3, linear; imeth=4, nearest

   dat2=miss

   minv=minval(azim1)
   loc1=minloc(azim1)
   maxv=maxval(azim1)
   loc2=maxloc(azim1)
   ns=loc1(1)
   ne=loc2(1)
   if(abs(ns-ne)==1)then
      ns=1
      ne=na1
   else
      ns=min(loc1(1),loc2(1))
      ne=max(loc1(1),loc2(1))
   endif
   js=0
   !write(*,*) "azim1",azim1
   !write(*,*) "azim2",azim2
   write(*,*) loc1(1), minv ,loc2(1), maxv
   do j1=1, na2
      if(azim2(j1)<minv)then
         j2=loc2(1)
         j3=loc1(1)
         js(j1)=j2
         je(j1)=j3
         jy(j1)=(azim2(j1)+360-azim1(j2))/(azim1(j3)+360-azim1(j2))
!        write(*,*) "azim2(j1)<minv", azim2(j1), minv
      elseif(azim2(j1)>maxv)then
         j2=loc2(1)
         j3=loc1(1)
         js(j1)=j2
         je(j1)=j3
         jy(j1)=(azim2(j1)-azim1(j2))/(azim1(j3)+360-azim1(j2))
!        write(*,*) "azim2(j1)>maxv", azim2(j1), maxv
      else 
         do j2=ns, ne !1, na1
            j3=j2+1
            if(j3>na1) j3=1
            if(j2==loc2(1)) cycle
            if((azim2(j1)-azim1(j2))*(azim2(j1)-azim1(j3))<=0)then
!write(*,*) j1, j2, j3, azim2(j1), azim1(j2), azim1(j3)
               js(j1)=j2
               je(j1)=j3
               jy(j1)=(azim2(j1)-azim1(j2))/(azim1(j3)-azim1(j2))
               exit
            endif
         enddo
      endif
      !write(*,*) j1, minv, maxv, js(j1), je(j1), azim2(j1)
      !if(js(j1)==0.or.je(j1)==0)then
      !write(*,*) "js(j1)==0.or.je(j1)==0" 
      !stop
      !else
      !write(*,*) azim2(j1), js(j1), je(j1), azim1(js(j1)), azim1(je(j1)), jy(j1)
      !endif
      do i1=1, nr2
         ix=(i1-1)*dis2/dis1+1
         i2=floor(ix)
         i3=i2+1
         ix=(ix-i2)
         if(ix==0.and.i2==nr1)then
           i3=nr1
         endif
         if(i3>nr1) cycle
!write(*,*) ubound(dat1),i2,i3,js(j1),je(j1)
         d1=dat1(i2,js(j1))
         d2=dat1(i3,js(j1))
         d3=dat1(i2,je(j1))
         d4=dat1(i3,je(j1))
         w1=(1-ix)*(1-jy(j1))
         w2=(  ix)*(1-jy(j1))
         w3=(1-ix)*(  jy(j1))
         w4=(  ix)*(  jy(j1))
         if(imeth==1)then
            if(w1/=0.and.d1/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d1
               dat2(i1,j1)=min(dat2(i1,j1),d1)
            endif
            if(w2/=0.and.d2/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d2
               dat2(i1,j1)=min(dat2(i1,j1),d2)
            endif
            if(w3/=0.and.d3/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d3
               dat2(i1,j1)=min(dat2(i1,j1),d3)
            endif
            if(w4/=0.and.d4/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d4
               dat2(i1,j1)=min(dat2(i1,j1),d4)
            endif
         elseif(imeth==2)then
            if(w1/=0.and.d1/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d1
               dat2(i1,j1)=max(dat2(i1,j1),d1)
            endif
            if(w2/=0.and.d2/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d2
               dat2(i1,j1)=max(dat2(i1,j1),d2)
            endif
            if(w3/=0.and.d3/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d3
               dat2(i1,j1)=max(dat2(i1,j1),d3)
            endif
            if(w4/=0.and.d4/=miss)then
               if(dat2(i1,j1)==miss) dat2(i1,j1)=d4
               dat2(i1,j1)=max(dat2(i1,j1),d4)
            endif
         elseif(imeth==3)then
            dat2(i1,j1)=w1*d1+w2*d2+w3*d3+w4*d4
            if((d1==miss.and.w1/=0).or.&
               (d2==miss.and.w2/=0).or.&
               (d3==miss.and.w3/=0).or.&
               (d4==miss.and.w4/=0))then
              dat2(i1,j1)=miss
            endif
         elseif(imeth==4)then
            if(d1/=miss)then
               r1=sqrt(((ix  )**2)+((jy(j1)  )**2))
               if(dat2(i1,j1)==miss)then
                  dat2(i1,j1)=d1
                  r0=r1
               endif
               if(r1<r0) then
                  dat2(i1,j1)=d1
                  r0=r1
               endif
            endif
            if(d2/=miss)then
               r1=sqrt(((ix-1)**2)+((jy(j1)  )**2))
               if(dat2(i1,j1)==miss)then
                  dat2(i1,j1)=d2
                  r0=r1
               endif
               if(r1<r0) then
                  dat2(i1,j1)=d2
                  r0=r1
               endif
            endif
            if(d3/=miss)then
               r1=sqrt(((ix   )**2)+((jy(j1)-1)**2))
               if(dat2(i1,j1)==miss)then
                  dat2(i1,j1)=d3
                  r0=r1
               endif
               if(r1<r0) then
                  dat2(i1,j1)=d3
                  r0=r1
               endif
            endif
            if(d4/=miss)then
               r1=sqrt(((ix-1)**2)+((jy(j1)-1)**2))
               if(dat2(i1,j1)==miss)then
                  dat2(i1,j1)=d4
                  r0=r1
               endif
               if(r1<r0) then
                  dat2(i1,j1)=d4
                  r0=r1
               endif
            endif
         endif
      enddo
   enddo

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grid_csv(filein,nr,na,dr,azim,vel,ref,ref0,spw, flag,region,cluster,miss)
   implicit none
   character(len=*),          intent(in) :: filein
   integer,                   intent(in) :: nr, na
   real,                      intent(in) :: miss, dr
   real,    dimension(na),    intent(in) :: azim 
   real,    dimension(nr,na), intent(in) :: vel, ref, ref0, spw 
   integer, dimension(nr,na), intent(in) :: flag, region, cluster

   character(len=800)          :: filename
   integer :: i, j

   call get_filename(filein, filename)

   open(21,file=filename,form="formatted",status="unknown")
   write(*,*) "writing file:", trim(filename), miss

  write(21,"(A)") "POINT_NAME, Range,  Azim,    Vel,    Ref,   Ref0,    Spw,Flag, REGION,CLUSTER"
   do j=1, na
      do i=1, nr
         if(vel(i,j)/=miss)then
            write(21,"(A,',',F6.2,',',F6.1,4(',',F7.1),',',I4,',',A,',',A)") &
               point_name(i,j),i*dr/1000.,azim(j),vel(i,j),ref(i,j),ref0(i,j),spw(i,j),flag(i,j),region_name(region(i,j)),cluster_name(cluster(i,j))
         endif
      enddo
   enddo
   write(*,*) "write csv file end"
   close(21)
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   character(len=7) function cluster_name(i)
   implicit none
   integer, intent(in) :: i

   !write(*,*) "cluster_name", i
   write(cluster_name,"('C-',I5.5)") i

   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   character(len=7) function region_name(i)
   implicit none
   integer, intent(in) :: i

   write(region_name,"('R-',I5.5)") i

   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   character(len=10) function point_name(i,j)
   implicit none
   integer, intent(in) :: i, j

   write(point_name,"('P-',I4.4,'_',I3.3)") i,j

   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_region_csv(filein,num_region,cluster,num_neighbor,base_region,nyquist_number,miss)
   implicit none
   character(len=*),               intent(in) :: filein
   integer,                        intent(in) :: num_region
   real,                           intent(in) :: miss
   integer, dimension(num_region), intent(in) :: cluster,num_neighbor,base_region,nyquist_number

   integer :: i, j

   character(len=800) :: filename
   call get_filename(filein,filename)

   open(21,file=filename,form="formatted",status="unknown")
   write(*,*) "writing file:", trim(filename)

   write(21,"(A)") " REGION,CLUSTER,NEIGHBOR,BASE-REG,NYQUIST"
   do i=1, num_region
      write(21,"(A,', ',A,',',',I8,', ',A,',',',I7)") &
      region_name(i),cluster_name(cluster(i)),num_neighbor(i),region_name(base_region(i)),nyquist_number(i)
   enddo
   close(21)
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_test_data_real(filename,nr,na,azim,dis,dat,miss)
   implicit none
   character(len=*),       intent(in) :: filename
   integer,                intent(out) :: nr, na
   real,                   intent(out) :: miss, dis
   real, dimension(   :), allocatable, intent(out) :: azim
   real, dimension(: ,:), allocatable, intent(out) :: dat

   open(21,file=filename,form="unformatted",status="old")
   write(*,*) "reading file:", trim(filename)
   read(21) nr, na, miss, dis
   write(*,*) "nr,na:", nr, na, miss, dis
   allocate(azim(na))
   allocate(dat(nr,na))
   read(21) azim
   read(21) dat
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_test_data_real(filein,nr,na,azim,dis,dat,miss)
   implicit none
   character(len=*),       intent(in) :: filein
   integer,                intent(in) :: nr, na
   real,                   intent(in) :: miss, dis
   real, dimension(   na), intent(in) :: azim
   real, dimension(nr,na), intent(in) :: dat
   
   character(len=800) :: filename
   call get_filename(filein,filename)

   open(21,file=filename,form="unformatted",status="unknown")
   write(*,*) "writing file:", trim(filename), nr, na, miss, dis
   write(21) nr, na, miss, dis
   write(21) azim 
   write(21) dat 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_test_data_int(filein,nr,na,azim,dis,dat,miss)
   implicit none
   character(len=*),          intent(in) :: filein
   integer,                   intent(in) :: nr, na
   integer,                   intent(in) :: miss
   real,                      intent(in) :: dis
   real,    dimension(   na), intent(in) :: azim
   integer, dimension(nr,na), intent(in) :: dat

   character(len=800) :: filename
   call get_filename(filein,filename)

   open(21,file=filename,form="unformatted",status="unknown")
   write(*,*) "writing file:", trim(filename)
   write(21) nr, na, miss, dis
   write(21) azim 
   write(21) dat 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grd_real(filein,nx,ny,dat,miss)
   implicit none
   character(len=*), intent(in) :: filein
   integer, intent(in) :: nx, ny
   real, dimension(nx,ny), intent(in) :: dat
   real, intent(in) :: miss 

   character(len=80) :: fmt_str
   integer :: i, j
   real    :: minv, maxv

   character(len=800) :: filename
   call get_filename(filein,filename)

   minv=miss
   maxv=miss
   do i=1, nx
      do j=1, ny
         if(dat(i,j)/=miss)then
            if(minv==miss)then
               minv=dat(i,j)
            else 
               minv=min(minv,dat(i,j))
            endif
         
            if(maxv==miss)then
               maxv=dat(i,j)
            else
               maxv=max(maxv,dat(i,j))
            endif
         endif
      enddo
   enddo


   WRITE(*,*) " Writing grd file: "//TRIM(FILENAME)
   OPEN(21,FILE=FILENAME,STATUS='UNKNOWN')
   WRITE(21,"(A)") "DSAA"
   WRITE(21,"(2I15)") NX, NY
   WRITE(21,"(2I15)") 1,nx 
   WRITE(21,"(2I15)") 1,ny 
   WRITE(21,"(2ES15.6)") minv,maxv
   write(fmt_str,"(A,I5,A)") "(",nx,"ES15.6)"
   DO J=1,NY
      WRITE(21,FMT=fmt_str) (DAT(I,J),I=1,NX)
   ENDDO
   CLOSE(21)
!   write(*,*) "write end"

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grd_int(filein,nx,ny,dat,miss)
   implicit none
   character(len=*), intent(in) :: filein
   integer, intent(in) :: nx, ny
   integer, dimension(nx,ny), intent(in) :: dat
   integer, intent(in) :: miss

   character(len=80) :: fmt_str
   integer :: i, j
   integer :: minv, maxv

   character(len=800) :: filename
   call get_filename(filein,filename)

!   write(*,*) ubound(dat), nx ,ny
   minv=miss
   maxv=miss
   do i=1, nx
      do j=1, ny
         if(dat(i,j)/=miss)then
            if(minv==miss)then
               minv=dat(i,j)
            else
               minv=min(minv,dat(i,j))
            endif

            if(maxv==miss)then
               maxv=dat(i,j)
            else
               maxv=max(maxv,dat(i,j))
            endif
         endif
      enddo
   enddo


   WRITE(*,*) " Writing grd file: "//TRIM(FILENAME)
   OPEN(21,FILE=FILENAME,STATUS='UNKNOWN')
   WRITE(21,"(A)") "DSAA"
   WRITE(21,"(2I15)") NX, NY
   WRITE(21,"(2I15)") 1,nx
   WRITE(21,"(2I15)") 1,ny
   WRITE(21,"(2I15)") minv,maxv
   write(fmt_str,"(A,I5,A)") "(",nx,"I15)"
   DO J=1,NY
      WRITE(21,FMT=fmt_str) (DAT(I,J),I=1,NX)
   ENDDO
   CLOSE(21)
!   write(*,*) "write end"

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grd_int8(filein,nx,ny,xmin, xmax, ymin, ymax, dat)
   implicit none
   character(len=*), intent(in) :: filein
   integer, intent(in) :: nx, ny
   real,    intent(in) :: xmin, xmax, ymin, ymax 
   integer, dimension(nx,ny), intent(in) :: dat

   character(len=80) :: fmt_str
   integer :: i, j
   integer :: minv, maxv

   character(len=800) :: filename
   call get_filename(filein,filename)

!   write(*,*) ubound(dat), nx ,ny
   minv=minval(dat)
   maxv=maxval(dat)

   WRITE(*,*) " Writing grd file: "//TRIM(FILENAME)
   OPEN(21,FILE=FILENAME,STATUS='UNKNOWN')
   WRITE(21,"(A)") "DSAA"
   WRITE(21,"(2I15)") NX, NY
   WRITE(21,"(2F15.2)") xmin,xmax
   WRITE(21,"(2f15.2)") ymin,ymax
   WRITE(21,"(2I15)") minv,maxv
   write(fmt_str,"(A,I5,A)") "(",nx,"I15)"
   DO J=1,NY
      WRITE(21,FMT=fmt_str) (DAT(I,J),I=1,NX)
   ENDDO
   CLOSE(21)
!   write(*,*) "write end"

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grd_real8(filein,nx,ny,xmin, xmax, ymin, ymax, dat)
   implicit none
   character(len=*), intent(in) :: filein
   integer, intent(in) :: nx, ny
   real,    intent(in) :: xmin, xmax, ymin, ymax
   real, dimension(nx,ny), intent(in) :: dat

   character(len=80) :: fmt_str
   integer :: i, j
   real :: minv, maxv

   character(len=800) :: filename
   call get_filename(filein,filename)

!   write(*,*) ubound(dat), nx ,ny
   minv=minval(dat)
   maxv=maxval(dat)

   WRITE(*,*) " Writing grd file: "//TRIM(FILENAME)
   OPEN(21,FILE=FILENAME,STATUS='UNKNOWN')
   WRITE(21,"(A)") "DSAA"
   WRITE(21,"(2I15)") NX, NY
   WRITE(21,"(2F15.5)") xmin,xmax
   WRITE(21,"(2f15.5)") ymin,ymax
   WRITE(21,"(2F15.5)") minv,maxv
   write(fmt_str,"(A,I5,A)") "(",nx,"F15.5)"
   DO J=1,NY
      WRITE(21,FMT=fmt_str) (DAT(I,J),I=1,NX)
   ENDDO
   CLOSE(21)
!   write(*,*) "write end"

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_csv_real4(filein,nx,x, dat)
   implicit none
   character(len=*),    intent(in) :: filein
   integer,             intent(in) :: nx
   real, dimension(nx), intent(in) :: x
   real, dimension(nx), intent(in) :: dat

   character(len=80) :: fmt_str
   integer :: i, j
   real :: minv, maxv

   character(len=800) :: filename
   call get_filename(filein,filename)

   WRITE(*,*) " Writing csv file: "//TRIM(FILENAME)
   OPEN(21,FILE=FILENAME,STATUS='UNKNOWN')
   do i=1, nx
      WRITE(21,*) x(i),',',dat(i) 
   ENDDO
   CLOSE(21)
!   write(*,*) "write end"

   end subroutine

end module
