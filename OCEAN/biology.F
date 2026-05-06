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
#if defined SOLVE3D && defined BIOLOGY && defined PISCES
# include "biology_pisces.F"
#elif defined SOLVE3D && defined BIOLOGY && defined BIO_NChlPZD
# include "bio_NChlPZD.F"
#elif defined SOLVE3D && defined BIOLOGY && defined BIO_N2ChlPZD2
# include "bio_N2ChlPZD2.F"
#elif defined SOLVE3D && defined BIOLOGY && defined BIO_BioEBUS
# include "bio_BioEBUS.F"
#else
      subroutine biology_empty
      end
#endif

