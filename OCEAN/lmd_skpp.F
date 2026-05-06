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
#include "cppdefs.h"

#ifdef LMD_SKPP
# if defined LMD_SKPP2005
#  include "lmd_skpp2005.F"
# else
#  include "lmd_skpp1994.F"
# endif
#else
!========================
      subroutine lmd_skpp_empty
      return
      end
!========================
#endif
