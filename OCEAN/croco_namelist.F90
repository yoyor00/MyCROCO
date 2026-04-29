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
      implicit none
      integer, intent(out) :: ierr

      integer :: nmlunit, ios

      namelist /croco_title/ title
      namelist /croco_time_stepping/ dt, ntimes, ndtfast, ninfo

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
      call check_nml_croco_time_stepping(ierr)
      call init_time_stepping_param

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

end module croco_namelist
