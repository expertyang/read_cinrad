MODULE date_pack

      IMPLICIT NONE
      REAL ,      PARAMETER :: one_day = 86400. , one_hour = 3600. , &
                               one_minute = 60. , one_second = 1. , one_millisecond = 0.001
      
      CHARACTER(*) , PARAMETER :: date_char_fmt = '( I4.4,2I2.2,3I2.2,".",I3.3 )' 

      CHARACTER(*) , PARAMETER :: date_str_fmt = '( I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2,".",I3.3 )'
 
      INTEGER      , PARAMETER :: date_char_len = 18
                                  
      INTEGER      , PARAMETER :: date_str_len = 23

      INTEGER      , PARAMETER :: mday(12) = (/ 31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 , 31 /)
      character(len=3) , PARAMETER :: month_name (12) = (/ "Jan" , "Feb" , "Mar" , "Apr" , "May" , "Jun" , &
                                                           "Jul" , "Aug" , "Sep" , "Oct" , "Nov" , "Dec" /)
      PRIVATE mday
      
      
      
      TYPE date
      SEQUENCE
         INTEGER :: year         ,    &
                    month        ,    &
                    day          ,    &
                    hour         ,    &
                    minute       ,    &
                    second       ,    &
                    millisecond 
      END TYPE
      
      TYPE ( date ) , PARAMETER :: empty_date = DATE ( 0 , 0 , 0 , 0 , 0 , 0 , 0 )

      INTERFACE month_days 
         MODULE PROCEDURE month_days_date , month_days_ym
      END INTERFACE 
         
      INTERFACE init_date
         MODULE PROCEDURE ymdhms_to_date , ymd_hms_to_date , string_to_date
      END INTERFACE

      INTERFACE split_date
         MODULE PROCEDURE split_date_ymdhms , split_date_ymd_hms
      END INTERFACE
      
      INTERFACE compare
         MODULE PROCEDURE compare_integer , compare_real , compare_date
      END INTERFACE
      
      INTERFACE get_new_date
         MODULE PROCEDURE get_new_date_r4 , get_new_date_r8
      END INTERFACE


      PRIVATE get_diff_date_positive
CONTAINS

!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION date_to_string ( datein ) RESULT ( string )
      
      IMPLICIT NONE
      TYPE ( date ) ,     INTENT( IN ) :: datein
      CHARACTER( LEN = date_str_len ) :: string
      
      WRITE ( string , FMT = date_str_fmt ) datein
   
   END FUNCTION  
   
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION date_to_char( datein ) RESULT ( string )
      
      IMPLICIT NONE
      TYPE ( date ) ,     INTENT( IN ) :: datein
      CHARACTER( LEN = date_char_len ) :: string
      
      WRITE( string , FMT = date_char_fmt ) datein
      
   END FUNCTION  
   
!--------------------------------------------------------------------------------
   
   ELEMENTAL LOGICAL FUNCTION is_leap_year ( year ) 
   
      IMPLICIT NONE
      INTEGER , INTENT ( IN ) :: year
      
      is_leap_year = .FALSE.
      IF (MOD(year,4).eq.0) THEN  
         is_leap_year = .TRUE.
         IF (MOD(year,100).eq.0) THEN
            is_leap_year = .FALSE.
            IF (MOD(year,400).eq.0) THEN
               is_leap_year = .TRUE.
            END IF
         END IF
      END IF

   END FUNCTION  

!--------------------------------------------------------------------------------

   ELEMENTAL FUNCTION month_days_date ( datein ) RESULT (num_days)
   
      IMPLICIT NONE
   
      TYPE( date ) , INTENT(IN) :: datein
      INTEGER :: num_days
      
      INTEGER :: year , month

      year  = datein % year
      month = datein % month
      num_days = mday ( month )
      IF ( month .EQ. 2 ) THEN
         IF ( is_leap_year ( year ) ) THEN
                  num_days = mday(2) + 1  
         END IF
      END IF
      
   END FUNCTION 

!--------------------------------------------------------------------------------

   ELEMENTAL FUNCTION month_days_ym ( year , month ) RESULT (num_days)
   
      IMPLICIT NONE
   
      INTEGER , INTENT(IN) ::  year , month 
      INTEGER :: num_days

      num_days = mday ( month )
      IF ( month .EQ. 2 ) THEN
         IF ( is_leap_year ( year ) ) THEN
                  num_days = mday(2) + 1  
         END IF
      END IF

   END FUNCTION 

!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION increase_day ( datein ) RESULT ( dateout )
   
      IMPLICIT NONE
      TYPE ( date ) :: datein , dateout
      INTENT ( IN ) :: datein
      
      dateout = datein
      dateout % day = dateout % day + 1
      IF ( dateout % day > month_days ( dateout ) ) THEN
         dateout % day = 1
         dateout % month = dateout % month + 1
         IF ( dateout % month > 12 ) THEN
            dateout % month = 1
            dateout % year = dateout % year + 1
         ENDIF
      ENDIF
   END FUNCTION

!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION decrease_day ( datein ) RESULT ( dateout )
   
      IMPLICIT NONE
      TYPE ( date ) :: datein , dateout
      INTENT ( IN ) :: datein
      
      dateout = datein
      dateout % day = dateout % day - 1
      IF ( dateout % day < 1 ) THEN
         dateout % month = dateout % month - 1
         IF ( dateout % month < 1 ) THEN
            dateout % year  = dateout % year - 1
            dateout % month = 12
         ENDIF
         dateout % day   = month_days ( dateout )
      ENDIF
   END FUNCTION


!-------------------------------------------------------------------------------

   ELEMENTAL LOGICAL FUNCTION is_valid_date ( datein )
   
      IMPLICIT NONE
      TYPE ( date ) , INTENT ( IN ) :: datein
      
      is_valid_date = ( datein % month       >= 1 .AND. datein % month       <= 12   ) .AND. &
                      ( datein % day         >= 1 .AND. datein % day         <= 31   ) .AND. &
                      ( datein % hour        >= 0 .AND. datein % hour        <= 23   ) .AND. &
                      ( datein % minute      >= 0 .AND. datein % minute      <= 59   ) .AND. &
                      ( datein % second      >= 0 .AND. datein % second      <= 59   ) .AND. &
                      ( datein % millisecond >= 0 .AND. datein % millisecond <= 999  ) 
                      
      IF ( is_valid_date ) THEN            
         is_valid_date = ( datein % day      >= 1 .AND. datein % day         <= month_days ( datein ) )
      ENDIF
      
   END FUNCTION
      
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION get_new_date_r4 ( datein , ds ) RESULT ( dateout ) 

      IMPLICIT NONE
      TYPE ( date ) :: datein , dateout
      REAL :: ds     ! delta_seconds
              
      INTENT ( IN ) :: datein , ds 
      
      REAL :: ads  
      INTEGER :: nday , nhr , nmin , nsec , nms , i 
      
      dateout = datein
      
      ads = ABS ( ds)
      
      nday = ads / one_day
      nhr  = MOD ( ads , one_day ) / one_hour
      nmin = MOD ( ads , one_hour ) / one_minute
      nsec = MOD ( ads , one_minute ) / one_second
      nms  = INT ( MOD ( ads , one_second ) / one_millisecond + 0.5 ) 
 
      IF ( ds >= 0 ) THEN
      
         dateout % millisecond = dateout % millisecond + nms
         IF ( dateout % millisecond >= 1000 ) THEN
            dateout % millisecond = dateout % millisecond - 1000
            nsec = nsec + 1
         ENDIF

         dateout % second = dateout % second + nsec
         IF ( dateout % second >=60 )  THEN
            dateout % second = dateout % second - 60
            nmin = nmin + 1
         ENDIF
      
         dateout % minute = dateout % minute + nmin
         IF ( dateout % minute >= 60 )  THEN
            dateout % minute = dateout % minute - 60
            nhr = nhr + 1
         ENDIF

         dateout % hour = dateout % hour + nhr
         IF ( dateout % hour >= 24 )  THEN
            dateout % hour = dateout % hour - 24
            nday = nday + 1
         ENDIF

         DO i = 1 , nday
            dateout = increase_day ( dateout )
         ENDDO
         
      ELSE ! ds < 0
       
         dateout % millisecond = dateout % millisecond - nms
         IF ( dateout % millisecond < 0 )  THEN
            dateout % millisecond = dateout % millisecond + 1000
            nsec = nsec + 1
         ENDIF

         dateout % second = dateout % second - nsec
         IF ( dateout % second < 0 )  THEN
            dateout % second = dateout % second + 60
            nmin = nmin + 1
         ENDIF
      
         dateout % minute = dateout % minute - nmin
         IF ( dateout % minute < 0 )  THEN
            dateout % minute = dateout % minute + 60
            nhr = nhr + 1
         ENDIF

         dateout % hour = dateout % hour - nhr
         IF ( dateout % hour < 0 )  THEN
            dateout % hour = dateout % hour + 24
            nday = nday + 1
         ENDIF

         DO i = 1 , nday
            dateout = decrease_day ( dateout )
         ENDDO
         
      ENDIF        

   END FUNCTION
      
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION get_new_date_r8 ( datein , ds ) RESULT ( dateout ) 

      IMPLICIT NONE
      TYPE ( date ) :: datein , dateout
      REAL(kind=8) :: ds     ! delta_seconds
              
      INTENT ( IN ) :: datein , ds 
      
      REAL(KIND=8) :: ads  
      INTEGER :: nday , nhr , nmin , nsec , nms , i 
      
      dateout = datein
      
      ads = ABS ( ds)
      
      nday = ads / one_day
      nhr  = DMOD ( ads , DBLE(one_day   )) / one_hour
      nmin = DMOD ( ads , DBLE(one_hour  )) / one_minute
      nsec = DMOD ( ads , DBLE(one_minute)) / one_second
      nms  = INT ( DMOD ( ads , DBLE(one_second) ) / one_millisecond + 0.5 ) 
 
      IF ( ds >= 0 ) THEN
      
         dateout % millisecond = dateout % millisecond + nms
         IF ( dateout % millisecond >= 1000 ) THEN
            dateout % millisecond = dateout % millisecond - 1000
            nsec = nsec + 1
         ENDIF

         dateout % second = dateout % second + nsec
         IF ( dateout % second >=60 )  THEN
            dateout % second = dateout % second - 60
            nmin = nmin + 1
         ENDIF
      
         dateout % minute = dateout % minute + nmin
         IF ( dateout % minute >= 60 )  THEN
            dateout % minute = dateout % minute - 60
            nhr = nhr + 1
         ENDIF

         dateout % hour = dateout % hour + nhr
         IF ( dateout % hour >= 24 )  THEN
            dateout % hour = dateout % hour - 24
            nday = nday + 1
         ENDIF

         DO i = 1 , nday
            dateout = increase_day ( dateout )
         ENDDO
         
      ELSE ! ds < 0
       
         dateout % millisecond = dateout % millisecond - nms
         IF ( dateout % millisecond < 0 )  THEN
            dateout % millisecond = dateout % millisecond + 1000
            nsec = nsec + 1
         ENDIF

         dateout % second = dateout % second - nsec
         IF ( dateout % second < 0 )  THEN
            dateout % second = dateout % second + 60
            nmin = nmin + 1
         ENDIF
      
         dateout % minute = dateout % minute - nmin
         IF ( dateout % minute < 0 )  THEN
            dateout % minute = dateout % minute + 60
            nhr = nhr + 1
         ENDIF

         dateout % hour = dateout % hour - nhr
         IF ( dateout % hour < 0 )  THEN
            dateout % hour = dateout % hour + 24
            nday = nday + 1
         ENDIF

         DO i = 1 , nday
            dateout = decrease_day ( dateout )
         ENDDO
         
      ENDIF        

   END FUNCTION
      
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION get_diff_date ( new_date , old_date ) RESULT ( diff_date ) 
   
      IMPLICIT NONE
      TYPE ( date ) , INTENT( IN ) :: new_date , old_date
      REAL :: diff_date  ! in seconds

      IF ( compare_date ( new_date , old_date )  >= 0 ) THEN
         diff_date = get_diff_date_positive ( new_date , old_date )
      ELSE
         diff_date = - get_diff_date_positive ( old_date , new_date )
      ENDIF   
   END FUNCTION

!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION get_diff_date_positive ( new_date , old_date ) RESULT ( diff_date ) ! new_date - old_date
      
      IMPLICIT NONE
      TYPE ( date ) , INTENT( IN ) :: new_date , old_date
      REAL :: diff_date  ! in seconds
      
      INTEGER :: old_day , new_day , i
      
      old_day = 0
      IF ( old_date % month > 1 ) THEN
         DO i = 1 , old_date % month - 1
            old_day = old_day + month_days ( old_date % year , i )
         ENDDO
      ENDIF
      
      old_day = old_day + old_date % day - 1
      
      new_day = 0
      DO i = old_date % year , new_date % year - 1
         IF ( is_leap_year ( i ) ) THEN
            new_day = new_day + 366 
         ELSE
            new_day = new_day + 365 
         ENDIF
      ENDDO
      
      IF ( new_date % month > 1 ) THEN
         DO i = 1 , new_date % month - 1
            new_day = new_day + month_days ( new_date % year , i )
         ENDDO
      ENDIF
      
      new_day = new_day + new_date % day - 1
      
      diff_date = ( new_day - old_day ) * one_day
      diff_date = diff_date + ( new_date % hour        - old_date % hour        ) * one_hour
      diff_date = diff_date + ( new_date % minute      - old_date % minute      ) * one_minute
      diff_date = diff_date + ( new_date % second      - old_date % second      ) * one_second
      diff_date = diff_date + ( new_date % millisecond - old_date % millisecond ) * one_millisecond
      
      
   END FUNCTION

!-------------------------------------------------------------------------------

   SUBROUTINE split_date_ymdhms ( datein , year , month , day , hour , minute , second , millisecond ) 
     
      IMPLICIT NONE

      TYPE ( date ) , INTENT( IN ) :: datein 
   
      INTEGER , INTENT( OUT ) :: year , month , day , hour , minute , second , millisecond
      OPTIONAL hour , minute , second , millisecond
      
      year   = datein % year
      month  = datein % month
      day    = datein % day
      
      IF ( PRESENT( hour        ) ) hour         = datein % hour
      IF ( PRESENT( minute      ) ) minute       = datein % minute
      IF ( PRESENT( second      ) ) second       = datein % second
      IF ( PRESENT( millisecond ) ) millisecond  = datein % millisecond
      
   END SUBROUTINE 
   
!-------------------------------------------------------------------------------

   SUBROUTINE split_date_ymd_hms ( datein , ymd , hms ) 
     
      IMPLICIT NONE

      TYPE ( date ) , INTENT( IN ) :: datein 
      INTEGER        :: ymd
      REAL           :: hms
      INTENT( OUT )  :: ymd , hms
      OPTIONAL hms
      
      INTEGER :: year , month , day , hour , minute , second , millisecond
      
      ymd   = ( datein % year * 10000 ) + ( datein % month * 100 ) + datein % day
      
      IF ( PRESENT( hms        ) ) THEN
         hms         = datein % hour * 10000. + datein % minute * 100. + &
                       datein % second * 1.   + datein % millisecond * 0.001
      ENDIF      
   END SUBROUTINE 
   
!-------------------------------------------------------------------------------

   FUNCTION ymdhms_to_date ( year , month , day , hour , minute , second , millisecond ) RESULT ( dateout )
     
      IMPLICIT NONE
      
      INTEGER , INTENT( IN ) :: year , month , day , hour , minute , second , millisecond
      OPTIONAL hour , minute , second , millisecond
   
      TYPE ( date ) :: dateout
      
      dateout = empty_date

      dateout % year         = year
      dateout % month        = month
      dateout % day          = day
!      dateout % hour         = 0
!      dateout % minute       = 0
!      dateout % second       = 0
!      dateout % millisecond  = 0
      
      IF ( PRESENT( hour        ) ) dateout % hour         = hour
      IF ( PRESENT( minute      ) ) dateout % minute       = minute
      IF ( PRESENT( second      ) ) dateout % second       = second
      IF ( PRESENT( millisecond ) ) dateout % millisecond  = millisecond
      
!     WRITE(*,*) dateout
      IF ( .NOT. is_valid_date (dateout) ) THEN
         WRITE ( * , * ) 'invalid date : ' , date_to_string ( dateout )
         STOP
      ENDIF
   END FUNCTION

!--------------------------------------------------------------------------------


   FUNCTION ymd_hms_to_date ( ymd , hms ) RESULT ( dateout )
     
      IMPLICIT NONE
      
      INTEGER      :: ymd
      REAL         :: hms
      INTENT( IN ) :: ymd , hms
      OPTIONAL     :: hms
      
      TYPE ( date ) :: dateout
      
      INTEGER ::  year , month , day , ihms , hour , minute , second , millisecond
      
      dateout = empty_date

      day   = MOD ( ymd , 100 )
      month = ( MOD ( ymd , 10000 ) ) / 100
      year  = ( ymd ) / 10000 
      
      dateout % year         = year
      dateout % month        = month
      dateout % day          = day
!      dateout % hour         = 0
!      dateout % minute       = 0
!      dateout % second       = 0
!      dateout % millisecond  = 0
      
      IF ( PRESENT( hms ) ) THEN
         second      = INT ( MOD ( hms , 100. ) )
         minute      = INT ( MOD ( hms , 10000. ) / 100. ) 
         hour        = INT ( hms / 10000. ) 
         millisecond = INT ( ( MOD ( hms , 100. ) - second ) * 1000.)
         
         dateout % hour         = hour
         dateout % minute       = minute
         dateout % second       = second
         dateout % millisecond  = millisecond
      ENDIF
      
      IF ( .NOT. is_valid_date (dateout) ) THEN
         WRITE ( * , * ) 'invalid date : ' , date_to_string ( dateout )
         STOP
      ENDIF
   END FUNCTION

!--------------------------------------------------------------------------------

   FUNCTION string_to_date ( string ) RESULT ( dateout )
      
      IMPLICIT NONE
      CHARACTER( LEN = date_str_len ) , INTENT( IN ) :: string
      TYPE ( date ) :: dateout
      
      INTEGER :: ierr
      
      dateout = empty_date

      READ ( string , FMT = date_str_fmt , iostat = ierr) dateout
      
      IF ( ierr /= 0 ) THEN
         WRITE ( * , * ) 'error in string_to_date : ', string
         STOP
      ENDIF
   END FUNCTION  
   
!-------------------------------------------------------------------------------

   FUNCTION char_to_date ( string ) RESULT ( dateout )
      
      IMPLICIT NONE
      CHARACTER( LEN = date_char_len ) , INTENT( IN ) :: string
      TYPE ( date ) :: dateout
      
      INTEGER :: ierr
      
      dateout = empty_date

      READ ( string , FMT = date_char_fmt , iostat = ierr) dateout
      
      IF ( ierr /= 0 ) THEN
         WRITE ( * , * ) 'error in char_to_date : ', string
         STOP
      ENDIF
   END FUNCTION  
   
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION compare_date ( date_a , date_b ) RESULT ( value ) 
      
      IMPLICIT NONE
      TYPE ( date ) , INTENT ( IN ) :: date_a , date_b
      INTEGER :: value
      
      value = compare_integer ( date_a % year , date_b % year ) 
      
      IF ( value == 0 ) value = compare_integer ( date_a % month       , date_b % month       )
      IF ( value == 0 ) value = compare_integer ( date_a % day         , date_b % day         )
      IF ( value == 0 ) value = compare_integer ( date_a % hour        , date_b % hour        )
      IF ( value == 0 ) value = compare_integer ( date_a % minute      , date_b % minute      )
      IF ( value == 0 ) value = compare_integer ( date_a % second      , date_b % second      )
      IF ( value == 0 ) value = compare_integer ( date_a % millisecond , date_b % millisecond )
      
   END FUNCTION
   
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION compare_integer ( a , b ) RESULT ( value )
   
      IMPLICIT NONE
      INTEGER , INTENT ( IN ) :: a , b
      INTEGER :: value
      
      IF ( a == b )THEN
         value = 0   
      ELSE IF ( a < b ) THEN
         value = -1
      ELSE IF ( a > b  ) THEN
         value = 1
      ENDIF
      
   END FUNCTION
   
!-------------------------------------------------------------------------------

   ELEMENTAL FUNCTION compare_real ( a , b ) RESULT ( value )
   
      IMPLICIT NONE
      REAL , INTENT ( IN ) :: a , b
      INTEGER :: value
      
      IF ( a == b )THEN
         value = 0   
      ELSE IF ( a < b ) THEN
         value = -1
      ELSE IF ( a > b  ) THEN
         value = 1
      ENDIF
      
   END FUNCTION
   
!-------------------------------------------------------------------------------

END MODULE date_pack
