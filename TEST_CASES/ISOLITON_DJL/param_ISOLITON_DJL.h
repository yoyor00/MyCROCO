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
!----------------------------------------------------------------------
! Dimensions of Physical Grid and array dimensions
!----------------------------------------------------------------------
      integer, parameter :: LLm0 = 256
      integer, parameter :: MMm0 = 1
      integer, parameter :: N = 128

#ifdef SOLVE3D
!----------------------------------------------------------------------
! Number of tracers
!----------------------------------------------------------------------
# ifdef PASSIVE_TRACER
      integer, parameter :: ntrc_pas = 1
# else
      integer, parameter :: ntrc_pas = 0
# endif

# if defined SUBSTANCE
      integer, parameter :: ntrc_subs = 0
      integer, parameter :: ntfix = 0
      integer, parameter :: ntrc_substot = ntrc_subs+ntfix
# else
      integer, parameter :: ntrc_subs = 0
      integer, parameter :: ntrc_substot = 0
# endif

# ifdef SEDIMENT
      integer, parameter :: NSAND = 2
      integer, parameter :: NMUD = 0
      integer, parameter :: NGRAV = 0
      integer, parameter :: ntrc_sed = NSAND+NMUD+NGRAV
      integer, parameter :: NST = ntrc_sed
# else
      integer, parameter :: ntrc_sed = 0
# endif

#endif /* SOLVE3D */

!----------------------------------------------------------------------
! Number of layers in Sediment (SL)
!----------------------------------------------------------------------
      integer, parameter :: N_sl = 0

!----------------------------------------------------------------------
! MPI related variables
!----------------------------------------------------------------------
#ifdef MPI
      integer, parameter :: NP_XI = 1
      integer, parameter :: NP_ETA = 4
      integer, parameter :: NNODES = NP_XI*NP_ETA
      integer, parameter :: NPP = 1
      integer, parameter :: NSUB_X = 1
      integer, parameter :: NSUB_E = 1
#elif defined OPENMP
      integer, parameter :: NPP = 4
      integer, parameter :: NSUB_X = 1
      integer, parameter :: NSUB_E = NPP
#else
      integer, parameter :: NPP = 1
      integer, parameter :: NSUB_X = 1
      integer, parameter :: NSUB_E = NPP
#endif

#include "param_dev.h"
