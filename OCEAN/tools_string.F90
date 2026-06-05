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
   ! set of utils to manage string in Fortran

   implicit none

   ! default
   private
   public to_uppercase

contains

   function to_uppercase(str) result(upper_str)
      ! Return uppercase of a string
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper_str
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

end module tools_string
