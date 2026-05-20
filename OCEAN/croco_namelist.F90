#include "cppdefs.h"

MODULE croco_namelist
   implicit none
   save
   public
   ! &croco_title
   character(80) :: title = "CROCO simulation"

   ! &croco_time_stepping
   real    :: dt = 0.d0
   integer :: ntimes = 0
   integer :: ndtfast = 20
   integer :: ninfo = 1

#ifdef USE_CALENDAR
   ! &croco_start_date
   character(len=19) :: start_date = '2000-01-01 00:00:00'

   ! &croco_end_date
   character(len=19) :: end_date = '2000-02-01 00:00:00'

   ! &output_time_steps:
   real :: dt_his = 1
   real :: dt_avg = 6
   real :: dt_rst = 12
#endif

   ! namelist filename (set via read_nml_fname from command-line arg 2)
   character(len=200) :: fname_nml = 'croco.nml'

contains
   subroutine read_nml_fname()
      implicit none
      integer :: nargs
      nargs = command_argument_count()

      call get_command_argument(2, fname_nml)
      if (len_trim(fname_nml) == 0) then
         fname_nml = 'croco.nml'
      end if
#ifdef AGRIF
      if (.Not. Agrif_Root()) then
# ifdef AGRIF_ADAPTIVE
         fname_nml = trim(fname_nml)//'.1'
# else
         fname_nml = trim(fname_nml)//'.'//Agrif_Cfixed()
# endif
      end if
#endif
   end subroutine read_nml_fname

   subroutine read_nml(ierr)
      use param, ONLY: stdout
      implicit none
      integer, intent(out) :: ierr

      integer :: nmlunit, ios

      namelist /croco_title/              title
      namelist /croco_time_stepping/      dt, ntimes, ndtfast, ninfo
#ifdef USE_CALENDAR
      namelist /croco_start_date/         start_date
      namelist /croco_end_date/           end_date
      namelist /croco_output_time_steps/  dt_his, dt_avg, dt_rst
#endif
      ierr = 0
      nmlunit = 10

      open (unit=nmlunit, file=trim(fname_nml), status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'WARNING: namelist file not found: ', trim(fname_nml)
         write (*, *) 'Using default values.'
         return
      end if

      read (nmlunit, nml=croco_title, iostat=ios)
      rewind (nmlunit)

      read (nmlunit, nml=croco_time_stepping, iostat=ios)
      rewind (nmlunit)
#ifdef USE_CALENDAR
      read (nmlunit, nml=croco_start_date, iostat=ios)
      rewind (nmlunit)
      WRITE (stdout, nml=croco_start_date) 

      read (nmlunit, nml=croco_end_date, iostat=ios)
      rewind (nmlunit)
      WRITE (stdout, nml=croco_end_date)
      
      read (nmlunit, nml=croco_output_time_steps, iostat=ios)
      rewind (nmlunit)
      WRITE (stdout, nml=croco_output_time_steps)
#endif

      call check_nml_croco_time_stepping(ierr)
      call init_time_stepping_param
#ifdef USE_CALENDAR
      call init_calendar_param
#endif

      close (nmlunit)

   end subroutine read_nml

   subroutine init_time_stepping_param
      use scalars, ONLY: dtfast
      implicit none
      ! TODO : place this in initialisation phase ??
      dtfast = dt/float(ndtfast)     ! set barotropic time step.
   end subroutine init_time_stepping_param

   subroutine check_nml_croco_time_stepping(ierr)
      use param, ONLY: stdout, NWEIGHT
      implicit none
      integer, intent(inout) :: ierr
      if (NWEIGHT < (2*ndtfast - 1)) then
         write (stdout, '(a,i0)') 'Error - Number of 2D timesteps (2*ndtfast-1): ', 2*ndtfast - 1
         write (stdout, '(a,i0)') 'exceeds barotopic weight dimension: ', NWEIGHT
         ierr = ierr + 1
      end if
      if (ntimes == 0) then
         write (stdout, '(a,i0)') 'Error - Null number timestep ntimes: ', ntimes
         ierr = ierr + 1
      end if
      if (dt == 0.d0) then
         write (stdout, '(a,f10.1)') 'Error - Null time step dt: ', dt
         ierr = ierr + 1
      end if
   end subroutine check_nml_croco_time_stepping

#ifdef USE_CALENDAR
   subroutine init_calendar_param
      use scalars, ONLY: nrrec,start_time
      use ncscrum, ONLY: origin_year,origin_month, origin_day, origin_hour, &
                         origin_minute,origin_second, origin_date,origin_date_in_sec
      implicit none
      ! TODO : place this in initialisation phase 
      real*8 :: tool_datosec 
# ifdef ANA_INITIAL
        if (nrrec.eq.0) then
          origin_date=start_date
          origin_date_in_sec=tool_datosec(origin_date)
          READ(origin_date(1:4),fmt='(i4)') origin_year
          READ(origin_date(6:7),fmt='(i2)') origin_month
          READ(origin_date(9:10),fmt='(i2)') origin_day
          READ(origin_date(12:13),fmt='(i2)') origin_hour
          READ(origin_date(15:16),fmt='(i2)') origin_minute
          READ(origin_date(18:19),fmt='(i2)') origin_second
          start_time = origin_date_in_sec
        endif
# endif
   end subroutine init_calendar_param
#endif

end module croco_namelist

