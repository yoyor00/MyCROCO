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
#if defined OA_COUPLING || defined OW_COUPLING

       SUBROUTINE cpl_prism_init

          USE mod_prism             ! OASIS3 prism module
          IMPLICIT NONE

# include "cppdefs.h"
# include "param.h"
# include "mpi_cpl.h"

      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_prism_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!                exchange between WRF and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication
      !!--------------------------------------------------------------------

      !------------------------------------------------------------------
      ! Initialize the PRISM system for the application
      !------------------------------------------------------------------

      CALL prism_init_comp_proto ( comp_id, comp_name, comp_ierror)
      IF (comp_ierror .NE. 0) THEN
        CALL prism_abort_proto(comp_id,
     &                         'cpl_prism_ini', 'init_comp_proto')
      END IF

      !------------------------------------------------------------------
      ! Get an MPI communicator for CROCO/AGRIF local communication
      !------------------------------------------------------------------

      CALL prism_get_localcomm_proto ( ocean_grid_comm, comp_ierror )
      IF (comp_ierror .NE. 0) THEN
        CALL prism_abort_proto(comp_id, 'cpl_prism_ini',
     &                         'get_localcomm_proto')
      END IF
      return
      end
#else
!---------------------------------------------------------------------------------------
! Dummy subroutine                                                           NO Coupling
!---------------------------------------------------------------------------------------
      subroutine cpl_prism_init_empty
      end
#endif

