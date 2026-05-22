#include "cppdefs.h"

MODULE croco_namelist
   implicit none
   save
   public

   ! &croco_title
   character(len=80) :: title = "CROCO simulation"

   ! &croco_logfile
   character(len=180) :: logname = "croco.log"

   ! &croco_time_stepping
   real    :: dt = 0.d0
   integer :: ntimes = 0
   integer :: ndtfast = 20
   integer :: ninfo = 1

#ifdef NBQ
   ! &croco_time_stepping_nbq
   integer :: ndtnbq = 1
   real :: csound_nbq = 1000
   real :: visc2_nbq = 0.01d0
#endif

#ifdef SOLVE3D
   ! &croco_s_coord
   real :: theta_s = 7.0d0
   real :: theta_b = 2.0d0
   real :: Tcline = 200.0d0
#endif

#ifdef USE_CALENDAR
   ! &croco_use_calendar
   character(len=19) :: start_date = '2000-01-01 00:00:00'
   character(len=19) :: end_date = '2000-02-01 00:00:00'
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
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      implicit none
      integer, intent(out) :: ierr

      integer :: nmlunit, ios

      namelist /croco_title/ title
      namelist /croco_logfile/ logname
      namelist /croco_time_stepping/ dt, ntimes, ndtfast, ninfo
#ifdef NBQ
      namelist /croco_time_stepping_nbq/ ndtnbq, csound_nbq, visc2_nbq
#endif
#ifdef SOLVE3D
      namelist /croco_s_coord/ theta_s, theta_b, Tcline
#endif
#ifdef USE_CALENDAR
      namelist /croco_use_calendar/ start_date, end_date, dt_his, dt_avg, dt_rst
#endif

      ierr = 0
      nmlunit = 10

      MPI_master_only write (stdout, *) '*** READING NAMELIST ***'
      open (unit=nmlunit, file=trim(fname_nml), status='old', &
            action='read', iostat=ios)
      if (ios /= 0) then
         MPI_master_only write (stdout, *) 'WARNING: namelist file not found: ', trim(fname_nml)
         MPI_master_only write (stdout, *) 'Using default values.'
         return
      end if

      read (nmlunit, nml=croco_title, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_title")

      read (nmlunit, nml=croco_time_stepping, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_time_stepping")
      call check_nml_croco_time_stepping(ierr)
      call init_time_stepping_param
#ifdef NBQ
      read (nmlunit, nml=croco_time_stepping_nbq, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_time_stepping_nbq")
      call check_nml_croco_time_stepping_nbq(ierr)
      call init_time_stepping_nbq_param
# endif
#ifdef SOLVE3D
      read (nmlunit, nml=croco_s_coord, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_s_coord")
      call check_nml_croco_s_coord(ierr)
      ! TODO call init_s_coord_param
#endif

#ifdef USE_CALENDAR
      read (nmlunit, nml=croco_use_calendar, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_use_calendar")
      call init_calendar_param
#endif

      ! put LOGFILE at the end so that all the warning about missing nml
      ! are write in stdout before change on logfile
#ifdef LOGFILE
      read (nmlunit, nml=croco_logfile, iostat=ios); rewind (nmlunit)
      call warn_if_nml_missing(ios, "croco_logfile")
      call init_logfile_param(ierr)
#endif

      close (nmlunit)

      ! write all namelist in stdout
      MPI_master_only WRITE (stdout, nml=croco_title)
#ifdef LOGFILE
      MPI_master_only WRITE (stdout, nml=croco_logfile)
#endif
      MPI_master_only WRITE (stdout, nml=croco_time_stepping)
#ifdef NBQ
      MPI_master_only WRITE (stdout, nml=croco_time_stepping_nbq)
#endif
#ifdef SOLVE3D
      MPI_master_only WRITE (stdout, nml=croco_s_coord)
#endif
#ifdef USE_CALENDAR
      MPI_master_only WRITE (stdout, nml=croco_use_calendar)
#endif

   end subroutine read_nml

#ifdef LOGFILE
   subroutine init_logfile_param(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      use scalars, ONLY: mynode2, NNODES2
      implicit none
#ifdef ENSEMBLE
      ! needed for ENSEMBLE cmember
#include "mpi_cpl.h"
#endif
      integer, intent(inout) :: ierr
      integer :: ios
      character(len=len(logname) + 20) :: logfile_path

#if defined MPI && defined PARALLEL_FILES
      call insert_node(logname, len_trim(logname), mynode2, NNODES2, ierr)
#endif

      if (mynode == 0) then
#ifndef AGRIF
         stdout = 80
         logfile_path = trim(logname)
#else
         stdout = Agrif_Get_Unit()
         if (.Not. Agrif_Root()) then
            logfile_path = trim(logname)//'.'//Agrif_Cfixed()
         else
            logfile_path = trim(logname)
         end if
#endif

#ifdef ENSEMBLE
         logfile_path = trim(cmember)//trim(logname)
#endif
         open (unit=stdout, file=trim(logfile_path), &
               status='REPLACE', form='formatted', iostat=ios)
         if (ios /= 0) then
            write (6, '(/1x,3A/)') 'ERROR: Cannot open log file: ', trim(logfile_path)
            ierr = ierr + 1
            return
         end if

      else
         stdout = 6
      end if

   end subroutine init_logfile_param
#endif /* LOGFILE */

   subroutine init_time_stepping_param
      use scalars, ONLY: dtfast
      implicit none
      ! TODO : place this in initialisation phase ??
      dtfast = dt/float(ndtfast)     ! set barotropic time step.
   end subroutine init_time_stepping_param

   subroutine check_nml_croco_time_stepping(ierr)
      use param, ONLY: stdout, NWEIGHT
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      implicit none
      integer, intent(inout) :: ierr
      if (NWEIGHT < (2*ndtfast - 1)) then
         MPI_master_only write (stdout, '(a,i0)') 'Error - Number of 2D timesteps (2*ndtfast-1): ', 2*ndtfast - 1
         MPI_master_only write (stdout, '(a,i0)') 'exceeds barotopic weight dimension: ', NWEIGHT
         ierr = ierr + 1
      end if
      if (ntimes == 0) then
         MPI_master_only write (stdout, '(a,i0)') 'Error - Null number timestep ntimes: ', ntimes
         ierr = ierr + 1
      end if
      if (dt == 0.d0) then
         MPI_master_only write (stdout, '(a,f10.1)') 'Error - Null time step dt: ', dt
         ierr = ierr + 1
      end if
   end subroutine check_nml_croco_time_stepping

#if defined NBQ
   subroutine init_time_stepping_nbq_param
      use scalars, ONLY: dtfast
      implicit none
      real dtnbq
      common/time_nbq2/dtnbq
      ! TODO : place this in initialisation phase ??
      dtnbq = dtfast
      ndtnbq = 1
   end subroutine init_time_stepping_nbq_param

   subroutine check_nml_croco_time_stepping_nbq(ierr)
      use param, ONLY: stdout
      use scalars, ONLY: g
#  if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#  endif
      implicit none
      integer, intent(inout) :: ierr
      if (ndtnbq <= 0) then
         MPI_master_only write (stdout, '(a,i0)') 'Error - NBQ acoustic substep ratio ndtnbq must be strictly positive: ', ndtnbq
         ierr = ierr + 1
      end if
      ! TODO : à l'initialisation: (pas encore possible a la lecture de la namelist)
      !if (csound_nbq <= 5.d0 * sqrt(g* hmax)) then
      !   write (stdout, '(a,f12.4,a,f12.4)') 'Error - pseudo-acoustic speed csound_nbq = ', csound_nbq, &
      !                               ' must exceed 5*sqrt(g*hmax) = ', 5.d0 * sqrt(g * hmax)
      !   ierr = ierr + 1
      !end if
      if (csound_nbq > 1500.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') 'Error - NBQ pseudo-acoustic speed csound_nbq = ', csound_nbq, &
            ' must not exceed real acoustic speed (1500 m/s).'
         ierr = ierr + 1
      end if
      if (visc2_nbq < 0.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') 'Error - NBQ bulk viscosity visc2_nbq = ', visc2_nbq, &
            ' must be positive.'
         ierr = ierr + 1
      end if
   end subroutine check_nml_croco_time_stepping_nbq
#endif

#ifdef SOLVE3D
   subroutine check_nml_croco_s_coord(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      implicit none
      integer, intent(inout) :: ierr
      if (theta_s < 0.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') 'Error - S-coord surface control parameter theta_s = ', theta_s, &
            ' must be positive.'
         ierr = ierr + 1
      end if
      if (theta_s > 20.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') 'Error - S-coord surface control parameter theta_s = ', theta_s, &
            ' is too large. It should not exceed 20.'
         ierr = ierr + 1
      end if
      if (theta_b < 0.d0) then
         MPI_master_only write (stdout, '(a,f12.4,a)') 'Error - S-coord bottom control parameter theta_b = ', theta_b, &
            ' must be positive.'
         ierr = ierr + 1
      end if
      if (Tcline < 0.d0) then
         MPI_master_only write (stdout, *) 'Error - S-coordinate surface/bottom layer width used in', &
            'vertical coordinate stretching Tcline = ', Tcline, &
            ' must be positive.'
         ierr = ierr + 1
      end if
   end subroutine check_nml_croco_s_coord
#endif

#ifdef USE_CALENDAR
   subroutine init_calendar_param
      use scalars, ONLY: nrrec, start_time
      use ncscrum, ONLY: origin_year, origin_month, origin_day, origin_hour, &
                         origin_minute, origin_second, origin_date, origin_date_in_sec
      implicit none
      ! TODO : place this in initialisation phase
      real*8 :: tool_datosec
# ifdef ANA_INITIAL
      if (nrrec .eq. 0) then
         origin_date = start_date
         origin_date_in_sec = tool_datosec(origin_date)
         READ (origin_date(1:4), fmt='(i4)') origin_year
         READ (origin_date(6:7), fmt='(i2)') origin_month
         READ (origin_date(9:10), fmt='(i2)') origin_day
         READ (origin_date(12:13), fmt='(i2)') origin_hour
         READ (origin_date(15:16), fmt='(i2)') origin_minute
         READ (origin_date(18:19), fmt='(i2)') origin_second
         start_time = origin_date_in_sec
      end if
# endif
   end subroutine init_calendar_param
#endif

   subroutine warn_if_nml_missing(ios, nml_name)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode   ! needed for MPI_master_only
#endif
      implicit none
      integer, intent(in) :: ios
      character(len=*), intent(in) :: nml_name
      if (ios /= 0) then
         MPI_master_only write (stdout, *) "WARNING : "//trim(nml_name)// &
            " namelist not found, default values are used"
      end if
   end subroutine warn_if_nml_missing

end module croco_namelist

