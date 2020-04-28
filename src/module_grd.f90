module grd

interface write_grd
   module procedure write_grd_real, write_grd_int
end interface

interface write_test_data
   module procedure write_test_data_real, write_test_data_int
end interface
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_test_data_real(filename,nr,na,azim,dis,dat,miss)
   implicit none
   character(len=*),       intent(in) :: filename
   integer,                intent(in) :: nr, na
   real,                   intent(in) :: miss, dis
   real, dimension(   na), intent(in) :: azim
   real, dimension(nr,na), intent(in) :: dat
   

   open(21,file=filename,status="unknown")
   write(*,*) "writing file:", trim(filename)
   write(21) nr, na, miss, dis
   write(21) azim 
   write(21) dat 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_test_data_int(filename,nr,na,azim,dis,dat,miss)
   implicit none
   character(len=*),          intent(in) :: filename
   integer,                   intent(in) :: nr, na
   integer,                   intent(in) :: miss
   real,                      intent(in) :: dis
   real,    dimension(   na), intent(in) :: azim
   integer, dimension(nr,na), intent(in) :: dat


   open(21,file=filename,status="unknown")
   write(*,*) "writing file:", trim(filename)
   write(21) nr, na, miss, dis
   write(21) azim 
   write(21) dat 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_grd_real(filename,nx,ny,dat,miss)
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in) :: nx, ny
   real, dimension(nx,ny), intent(in) :: dat
   real, intent(in) :: miss 

   character(len=80) :: fmt_str
   integer :: i, j
   real    :: minv, maxv

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
   subroutine write_grd_int(filename,nx,ny,dat,miss)
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in) :: nx, ny
   integer, dimension(nx,ny), intent(in) :: dat
   integer, intent(in) :: miss

   character(len=80) :: fmt_str
   integer :: i, j
   integer :: minv, maxv

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

end module
