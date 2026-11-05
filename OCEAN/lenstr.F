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
      integer function lenstr (string)
!
! Returns the character position of the last nonblank character in
! input string after removing the leading blank characters, if any.
! Algorithm: Find length of input character string using intrinsic
! function "len", then find positions of  the first "is" and the
! last "ie" nonblank characters. Define the length of the string
! based on these bounds. In the case of there are any leading blank
! characters, shift the whole string to the beginning.
!
      implicit none
      integer is, ie
      character*(*) string
      ie=len(string)
   1   if (string(ie:ie).eq.' ') then
         ie=ie-1
         if (ie.gt.0) goto 1
       endif
      is=0
   2   is=is+1
       if (string(is:is).eq.' ' .and. is.lt.ie) goto 2
      lenstr=ie-is+1
      if (is.gt.1) string=string(is:ie)
      return
      end
