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
#ifdef AGRIF
AGRIFISDEFINED
#endif
#ifdef OPENMP
OPENMPISDEFINED
#endif
#ifdef MPI
MPIISDEFINED
#endif
#if defined OA_COUPLING || defined OW_COUPLING
OACPLISDEFINED
#endif
#ifdef XIOS
XIOSISDEFINED
#endif
#ifdef NBQ
NBQISDEFINED
#endif
#ifdef OPENACC
OPENACCISDEFINED
#endif
#ifdef NBQ2
NBQ2ISDEFINED
#endif
#ifdef key_pisces_npzd
P2ZISDEFINED
#endif
#ifdef key_pisces_quota
P5ZISDEFINED
#endif
