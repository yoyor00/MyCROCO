!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
module tools_string
   ! Utilities for string manipulation in Fortran.
   ! All functions operate character by character via iachar/achar
   ! so they do not depend on any locale or compiler extension.

   implicit none
   private
   public :: to_uppercase, to_lowercase

contains

   !---------------------------------------------------------------------
   !  to_uppercase
   !  Return an upper-case copy of str.
   !---------------------------------------------------------------------
   function to_uppercase(str) result(upper_str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str))      :: upper_str
      integer :: i, ich

      do i = 1, len(str)
         ich = iachar(str(i:i))
         if (ich >= iachar('a') .and. ich <= iachar('z')) then
            upper_str(i:i) = achar(ich - 32)
         else
            upper_str(i:i) = str(i:i)
         end if
      end do
   end function to_uppercase

   !---------------------------------------------------------------------
   !  to_lowercase
   !  Return a lower-case copy of str.
   !---------------------------------------------------------------------
   function to_lowercase(str) result(lower_str)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str))      :: lower_str
      integer :: i, ich

      do i = 1, len(str)
         ich = iachar(str(i:i))
         if (ich >= iachar('A') .and. ich <= iachar('Z')) then
            lower_str(i:i) = achar(ich + 32)
         else
            lower_str(i:i) = str(i:i)
         end if
      end do
   end function to_lowercase

end module tools_string
