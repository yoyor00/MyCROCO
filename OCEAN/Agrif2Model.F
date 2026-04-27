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
!$AGRIF_DO_NOT_TREAT
!     **************************************************************************
!!!   Subroutine   Agrif_Set_numberofcells
!     **************************************************************************
!
      Subroutine Agrif_Set_numberofcells(Agrif_Gr)
      USE Agrif_Grids
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
      if ( associated(Agrif_Curgrid) )then
# include "SetNumberofcells.h"
      endif
      End Subroutine Agrif_Set_numberofcells
!
!     **************************************************************************
!!!   Subroutine   Agrif_Get_numberofcells
!     **************************************************************************
      Subroutine Agrif_Get_numberofcells(Agrif_Gr)
      USE Agrif_Grids
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
# include "GetNumberofcells.h"
      End Subroutine Agrif_Get_numberofcells
!
!     **************************************************************************
!!!   Subroutine Agrif_Allocationcalls
!     **************************************************************************
      Subroutine Agrif_Allocationcalls(Agrif_Gr)
      USE Agrif_Grids
# include "include_use_Alloc_agrif.h"
      Implicit none
      Type(Agrif_Grid), Pointer :: Agrif_Gr
# include "allocations_calls_agrif.h"
      End Subroutine Agrif_Allocationcalls
!
!     **************************************************************************
!!!   Subroutine Agrif_probdim_modtype_def
!     **************************************************************************
      Subroutine Agrif_probdim_modtype_def()
      Use Agrif_Grids
      Implicit none
# include "modtype_agrif.h"
# include "probdim_agrif.h"
# include "keys_agrif.h"
      Return
      End Subroutine Agrif_probdim_modtype_def
!
!     **************************************************************************
!!!   Subroutine Agrif_clustering_def
!     **************************************************************************
      Subroutine Agrif_clustering_def()
      Use Agrif_Grids
      Implicit none
!!#include "clustering_agrif.h"
      Return
      End Subroutine Agrif_clustering_def
!$AGRIF_END_DO_NOT_TREAT
#else
      subroutine Agrif2Model
         write(*,*) 'Impossible to bet here'
      end subroutine Agrif2model
#endif
