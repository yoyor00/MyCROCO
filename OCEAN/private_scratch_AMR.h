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
      real A2d(N2d,NiSA,0:NPP-1), A3d(N3d,4,0:NPP-1)
     &  ,A3dHz(N3dHz,4,0:NPP-1)	
# ifdef SEDIMENT
      integer B2d(N2d,0:NPP-1)
# endif

      common /private_scratch/ A2d,A3d,A3dHz
#  ifdef SEDIMENT
      common /private_scratch_bis/ B2d
#  endif


