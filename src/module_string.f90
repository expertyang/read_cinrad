MODULE string
IMPLICIT NONE

character(len=1), parameter :: endl=char(0)

CONTAINS

   !--------------------------------------------------------------------
   SUBROUTINE get_next_word(line,word,delim)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(INOUT)  :: line
   CHARACTER(LEN=*), INTENT(OUT)    :: word
   CHARACTER(LEN=*), INTENT(IN)     :: delim
!   
   INTEGER :: idx, len
!   
   line=ADJUSTL(line)
   len=LEN_TRIM(line)
   idx=SCAN(line, delim)
   IF(idx>0)THEN
      word=line(1:idx-1)
      line=line(idx+1:len)
   ELSE
      word=line
      line=""
   ENDIF
   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE count_word(line, delim, count)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)     :: line
   CHARACTER(LEN=*), INTENT(IN)     :: delim
   INTEGER,          INTENT(OUT)    :: count
   
   CHARACTER(LEN=LEN(line))   :: word
   CHARACTER(LEN=LEN(line))   :: tmp
   
   tmp=line
   count=0
   DO WHILE(LEN_TRIM(tmp)/=0)
      CALL get_next_word(tmp,word,delim)
!      write(*,*) trim(tmp)//trim(word)
      count=count+1
   ENDDO
   
   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE string2int(word,num,istat)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: word
   INTEGER,          INTENT(OUT) :: num, istat
!   
   READ(word,*,IOSTAT=istat) num
   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE string2real(word,num,istat)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN)  :: word
   REAL,             INTENT(OUT) :: num
   INTEGER,          INTENT(OUT) :: istat
!   
   READ(word,*,IOSTAT=istat) num
   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE to_upper(str)
   IMPLICIT NONE
   CHARACTER(*), INTENT(INOUT) :: str
!   
   INTEGER, PARAMETER :: ichar_a=ICHAR("a"), ichar_z=ICHAR("z"), ichar_ca=ICHAR("A")
   INTEGER :: code
   INTEGER :: i
!   
   DO i=1, LEN(str)
      code=ICHAR(str(i:i))
      IF(code>=ichar_a.AND.code<=ichar_z)THEN
         str(i:i)=CHAR(code+ichar_ca-ichar_a)
      ENDIF
   ENDDO

   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE to_lower(str)
   IMPLICIT NONE
   CHARACTER(*), INTENT(INOUT) :: str
!   
   INTEGER, PARAMETER :: ichar_a=ICHAR("a"), ichar_cz=ICHAR("Z"), ichar_ca=ICHAR("A")
   INTEGER :: code
   INTEGER :: i
!   
   DO i=1, LEN(str)
      code=ICHAR(str(i:i))
      IF(code>=ichar_ca.AND.code<=ichar_cz)THEN
         str(i:i)=CHAR(code+ichar_a-ichar_ca)
      ENDIF
   ENDDO

   END SUBROUTINE
   !--------------------------------------------------------------------
   SUBROUTINE name2mon(name, month)
   IMPLICIT NONE
   CHARACTER(*),  INTENT(INOUT)  :: name
   INTEGER,       INTENT(OUT)    :: month
   CALL to_upper(name)
   SELECT CASE(name)
   CASE("JAN")
      month=1
   CASE("FEB")
      month=2
   CASE("MAR")
      month=3
   CASE("APR")
      month=4
   CASE("MAY")
      month=5
   CASE("JUN")
      month=6
   CASE("JUL")
      month=7
   CASE("AUG")
      month=8
   CASE("SEP")
      month=9
   CASE("OCT")
      month=10
   CASE("NOV")
      month=11
   CASE("DEC")
      month=12
   END SELECT
   END SUBROUTINE
   !--------------------------------------------------------------------


   subroutine search_last_number(string, index)
   implicit none
   character(len=*), intent(in) :: string
   integer, intent(out) :: index

   integer :: i, j, n

   n=len_trim(string)
   index=0

   do i=1, n
      j=scan("0123456789", string(i:i))
      if(j==0)then
         index=i-1
         return
      endif
   enddo
   end subroutine

   SUBROUTINE count_line(iunit,nline)
   IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(out) :: nline
      
      INTEGER :: ierr
      
      nline=0
      DO WHILE(.TRUE.)
         READ(iunit,*, IOSTAT=ierr)
         IF(ierr/=0) EXIT
         nline=nline+1
      ENDDO
      REWIND(iunit)
   END SUBROUTINE

   subroutine search_next_blank(string, index)
   implicit none
   character(len=*), intent(in) :: string
   integer, intent(out) :: index

   integer :: i, j, n

   n=len_trim(string)
   index=0

   do i=1, n
      if(string(i:i)==" ")then
         index=i
         return
      endif
   enddo
   end subroutine

END MODULE
