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
#ifdef STATIONS
      subroutine init_arrays_sta

      implicit none
# include "param.h"
# include "sta.h"
# include "scalars.h"

      integer iflt, j, i, k
!
!
!  Initialize model control points
!
      do iflt=1,Msta
c       stagrd(iflt)=-1 ! -1 is the dead float value.
        do j=istagrd,istazgrd
           stainfo(j,iflt)=staspval
        enddo
        do j=1,NSTAVARS
          stadata(j,iflt)=0.
        enddo
# ifdef ALL_SIGMA
        do j=istadpt,istav
         do k=1,N
           staSigm(j,iflt,k)=staspval
         enddo
        enddo
# endif
      enddo

      return
      end

#else
      subroutine init_arrays_sta_empty
      return
      end

#endif /* STATIONS */

