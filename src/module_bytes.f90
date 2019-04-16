module bytes
! Input Code is Little Endian
   logical, private :: is_big_endian, if_checked=.false.

      INTERFACE get_i_value
         MODULE PROCEDURE get_i_value_1 , get_i_value_4 
      END INTERFACE
      INTERFACE get_f_value
         MODULE PROCEDURE get_f_value_1 , get_f_value_4 
      END INTERFACE

contains
   ! ui - unsigned int   : 4 bytes
   ! us - unsigned short : 2 bytes
   ! uc - unsigned char  : 1 bytes
   ! l  - long int       : 4 bytes
   ! i  - int            : 4 bytes
   ! s  - signed short   : 2 bytes
   ! c  - char           : 1 bytes
   ! f  - float          : 4 bytes   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   function char_to_int1(uc_code) result(rt_code)
!   implicit none
!   character(len=*), intent(in) :: uc_code
!   integer(kind=1),  dimension(:) :: rt_code
!
!   integer :: n, i
!
!   n=len(uc_code)
!   do i=1,n
!      rt_code(i)=get_uc_value(uc_code(i:i))
!   enddo
!   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_uc_value(uc_code)
   implicit none
   character(len=1), intent(in) :: uc_code

   get_uc_value=ichar(uc_code)
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_i_value_4(s_code)
   implicit none
   character(len=1), dimension(4),intent(in) :: s_code

   integer, dimension(4) :: code 
   integer(kind=8) :: i256=256

   code(4)=get_uc_value(s_code(1)(1:1))
   code(3)=get_uc_value(s_code(2)(1:1))
   code(2)=get_uc_value(s_code(3)(1:1))
   code(1)=get_uc_value(s_code(4)(1:1))

   if(code(1)>127)then ! negative number
      get_i_value_4=code(1)*256*256*256+code(2)*256*256+code(3)*256+code(4)-i256*256*256*256
   else
      get_i_value_4=code(1)*256*256*256+code(2)*256*256+code(3)*256+code(4)
   endif
   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_i_value_1(s_code)
   implicit none
   character(len=4), intent(in) :: s_code

   integer, dimension(4) :: code 
   integer(kind=8) :: i256=256

   code(4)=get_uc_value(s_code(1:1))
   code(3)=get_uc_value(s_code(2:2))
   code(2)=get_uc_value(s_code(3:3))
   code(1)=get_uc_value(s_code(4:4))

   if(code(1)>127)then ! negative number
      get_i_value_1=code(1)*256*256*256+code(2)*256*256+code(3)*256+code(4)-i256*256*256*256
   else
      get_i_value_1=code(1)*256*256*256+code(2)*256*256+code(3)*256+code(4)
   endif
   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function get_f_value_4(s_code)
   implicit none
   character(len=1), dimension(4), intent(in) :: s_code

   character(len=4) :: code
   real             :: value

   equivalence ( code, value)

   if(.not.if_checked)then
       is_big_endian= check_big_endian()
       if_checked=.true.
   endif
   if(is_big_endian)then
      code(4:4)=s_code(1)(1:1)
      code(3:3)=s_code(2)(1:1)
      code(2:2)=s_code(3)(1:1)
      code(1:1)=s_code(4)(1:1)
   else
      code(4:4)=s_code(4)(1:1)
      code(3:3)=s_code(3)(1:1)
      code(2:2)=s_code(2)(1:1)
      code(1:1)=s_code(1)(1:1)
   endif

   get_f_value_4=value
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real function get_f_value_1(s_code)
   implicit none
   character(len=4), intent(in) :: s_code

   character(len=4) :: code
   real             :: value
   
   equivalence ( code, value)

   if(.not.if_checked)then
       is_big_endian= check_big_endian()
       if_checked=.true.
   endif
   if(is_big_endian)then
      code(4:4)=s_code(1:1)
      code(3:3)=s_code(2:2)
      code(2:2)=s_code(3:3)
      code(1:1)=s_code(4:4)
   else
      code=s_code
   endif

   get_f_value_1=value
   end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_s_value(s_code)
   implicit none
   character(len=2), intent(in) :: s_code

   integer, dimension(2) :: code 

   code(2)=get_uc_value(s_code(1:1))
   code(1)=get_uc_value(s_code(2:2))

   if(code(1)>127)then ! negative number
      get_s_value=code(1)*256+code(2)-256*256
   else
      get_s_value=code(1)*256+code(2)
   endif
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_c_value(c_code)
   implicit none
   character(len=1), intent(in) :: c_code

   integer, dimension(1) :: code 

   code(1)=get_uc_value(c_code(1:1))

   if(code(1)>127)then ! negative number
      get_c_value=code(1)-256
   else
      get_c_value=code(1)
   endif
   end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_us_value(us_code)
   implicit none
   character(len=2), intent(in) :: us_code
  
   integer, dimension(2) :: code 

   code(2)=get_uc_value(us_code(1:1))
   code(1)=get_uc_value(us_code(2:2))
   get_us_value=code(1)*256+code(2)
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer(kind=8) function get_ui_value(ui_code)
   implicit none
   character(len=4), intent(in) :: ui_code
  
   integer(kind=8), dimension(2) :: code 

   code(2)=get_us_value(ui_code(1:2))
   code(1)=get_us_value(ui_code(3:4))
   get_ui_value=code(1)*256*256+code(2)
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   logical function check_big_endian()
   implicit none

   integer (kind=4) :: idata
   character(len=4) :: cdata
   equivalence(idata,cdata)

   idata=0
   cdata(4:4)=char(1)
   
   if(idata ==1)then
      check_big_endian=.true.
   else
      check_big_endian=.false.
   endif
   end function
end module
