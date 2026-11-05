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
/* This is include file "coupling.h":
  ----------------------------------------------------
  Variables responsible for communication between two-
  and three-dimensional parts of the model.
*/
#ifdef SOLVE3D
# ifdef VAR_RHO_2D
      real rhoA(GLOBAL_2D_ARRAY)
      real rhoS(GLOBAL_2D_ARRAY)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
# endif
      real rufrc(GLOBAL_2D_ARRAY)
      real rvfrc(GLOBAL_2D_ARRAY)
      real rufrc_bak(GLOBAL_2D_ARRAY,2)
      real rvfrc_bak(GLOBAL_2D_ARRAY,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak

      real Zt_avg1(GLOBAL_2D_ARRAY)
      real DU_avg1(GLOBAL_2D_ARRAY,5)
      real DV_avg1(GLOBAL_2D_ARRAY,5)
      real DU_avg2(GLOBAL_2D_ARRAY)
      real DV_avg2(GLOBAL_2D_ARRAY)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
#endif
