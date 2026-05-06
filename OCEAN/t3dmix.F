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
#if defined SOLVE3D && (defined TS_DIF2 || defined TS_DIF4) &&\
 (defined TS_MIX_ISO || defined TS_MIX_GEO || defined TS_MIX_S)
# if defined TS_MIX_ISO || defined TS_MIX_GEO
#  include "t3dmix_ISO.F"
# elif defined TS_MIX_S
#  include "t3dmix_S.F"
# endif
#else
      subroutine t3dmix_empty
      end
#endif

