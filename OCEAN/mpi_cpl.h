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
!$AGRIF_DO_NOT_TREAT
      INTEGER :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
#if defined ENSEMBLE
! Parameter for ensemble simulation
      INTEGER :: kmember ! index of ensemble member computed by this processor
      CHARACTER(len = 3) :: cmember
      common /ens_comm/ kmember, cmember
#endif /* ENSEMBLE */
#if defined OA_COUPLING || defined OW_COUPLING
      INTEGER :: comp_id                       ! component identification
      CHARACTER(len=6)   :: comp_name = 'crocox'

      INTEGER :: comp_ierror
      INTEGER :: oasis_part_id
      INTEGER :: oasis_var_nodims(2)
      INTEGER :: oasis_var_shape(4)
      INTEGER :: oasis_var_type
      INTEGER , dimension(5) :: oasis_ig_paral ! Box partiton

      INTEGER, PARAMETER ::   nmaxfld = 60 ! Maximum number of coupling fields
      INTEGER, PARAMETER ::   nmaxatm =  5 ! Maximum number of atmospheric models

      ! Coupling fields
      CHARACTER(len = 64), DIMENSION(nmaxfld) :: srcv_clname 
      CHARACTER(len = 64), DIMENSION(nmaxfld) :: ssnd_clname
      INTEGER, DIMENSION(0:nmaxatm, nmaxfld) :: srcv_nid, ssnd_nid
      common /exchange_fields_oasis3/ srcv_clname, ssnd_clname
      common /exchange_fields_oasis3_id/ srcv_nid, ssnd_nid

      INTEGER :: oasis_time, oasis_runtime
      common /exchange_times_oasis3/ oasis_time, oasis_runtime

      REAL cplmsk(GLOBAL_2D_ARRAY,0:nmaxatm)
      common /coupling_mask/cplmsk
#endif /* OA_COUPLING */
