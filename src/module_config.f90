module config
implicit none

integer, parameter                             :: max_radar_file=100
character(len=1024)                            :: wrf_file, info_file
character(len=1024), dimension(max_radar_file) :: radar_file
integer                                        :: n_radar

logical :: if_use_wrf, if_radar_qc, if_output_qc, if_output_raw, if_output_obs, &
           if_terminal_velocity, if_output_model, if_debug=.false., &
           if_remove_lowref
real    :: min_valid_nyquist=20
real, parameter :: missing_r=-888888.

namelist /file/ wrf_file, info_file, radar_file

namelist /param/ if_use_wrf, if_radar_qc, if_output_qc, if_output_raw, if_output_obs, &
                 if_terminal_velocity, if_output_model, if_debug, min_valid_nyquist, if_remove_lowref

contains

   subroutine read_namelist(filename)
   implicit none

   character(len=*), intent(in) :: filename

   integer, parameter :: iunit=11
   integer :: i

   radar_file=""
   wrf_file  =""
   info_file ="radar_info.txt"

   if_use_wrf           =.false. 
   if_radar_qc          =.false. 
   if_output_qc         =.false.
   if_output_raw        =.false. 
   if_output_obs        =.false. 
   if_terminal_velocity =.false. 
   if_output_model      =.false. 
   if_remove_lowref     =.false. 
   if_debug             =.false.

   open(iunit, file=filename, status="old")
   write(*,"(2A)") "Reading ", trim(filename)
   read(iunit, nml=file )
!  write(*,*) "radar_file",radar_file
   n_radar=0
   do i=1, max_radar_file
      if(len_trim(radar_file(i))>0)then
         n_radar=i
      endif
   enddo

   if(len_trim(wrf_file)==0)then
      if_use_wrf=.false.
   endif
   if(.NOT.if_use_wrf)then
      if_output_model=.false.
   endif
   write(*,*) "&file"
   write(*,*) "wrf_file   = ", '"'//trim(wrf_file)//'"'
   write(*,*) "radar_file = ", ('"'//trim(radar_file(i))//'" , ',i=1, n_radar)
   write(*,*) "info_file  = ", '"'//trim(info_file)//'"'
   write(*,*) "/"

   read(iunit, nml=param)
   write(*,*) "&param"
   write(*,*) "if_use_wrf           = ", if_use_wrf
   write(*,*) "if_radar_qc          = ", if_radar_qc  
   write(*,*) "if_output_qc         = ", if_output_qc 
   write(*,*) "if_output_raw        = ", if_output_raw
   write(*,*) "if_output_model      = ", if_output_model
   write(*,*) "if_output_obs        = ", if_output_obs
   write(*,*) "if_terminal_velocity = ", if_terminal_velocity
   write(*,*) "if_remove_lowref     = ", if_remove_lowref
   write(*,*) "if_debug             = ", if_debug
   write(*,*) "min_valid_nyquist    = ", min_valid_nyquist
   write(*,*) "/"
   close(iunit)

   if(n_radar<1)then
      write(*,*) "no radar_file to read!!!"
      stop
   endif
   end subroutine
end module
