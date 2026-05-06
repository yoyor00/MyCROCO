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
!
! This is "work.h": declaration of utility work array.
!
#ifdef SOLVE3D
      real work(GLOBAL_2D_ARRAY,0:N)
      common /work3d/ work
      real workr(GLOBAL_2D_ARRAY,1:N)
      common /work3d_r/ workr
#endif

#ifdef ABL1D
      real work3dabl(GLOBAL_2D_ARRAY,N_abl)
      common /work3d_abl/ work3dabl
#endif

      real work2d(GLOBAL_2D_ARRAY)
      common /work2d/ work2d

      real work2d2(GLOBAL_2D_ARRAY)
      common /work2d2/ work2d2



