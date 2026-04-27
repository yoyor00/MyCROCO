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
#if defined SOLVE3D && defined UV_VIS2 && defined UV_MIX_S
# include "uv3dmix_S.F"
#elif defined SOLVE3D && defined UV_VIS2 && defined UV_MIX_GEO
# include "uv3dmix_GP.F"
#elif defined SOLVE3D && defined UV_VIS4 && defined UV_MIX_S
# include "uv3dmix4_S.F"
#elif defined SOLVE3D && defined UV_VIS4 && defined UV_MIX_GEO
# include "uv3dmix4_GP.F"
#else
      subroutine uv3dmix_empty
      end
#endif

