!$AGRIF_DO_NOT_TREAT
      INTEGER :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT 

#ifdef OA_COUPLING  
      INTEGER :: comp_id                       ! component identification	
      CHARACTER(len=6)   :: comp_name = 'roms3d'

      INTEGER :: comp_ierror
      INTEGER :: oasis_part_id      
      INTEGER :: oasis_var_nodims(2) 
      INTEGER :: oasis_var_shape(4) 
      INTEGER :: oasis_var_type
      INTEGER , dimension(5) :: oasis_ig_paral ! Box partiton
 
      integer nmaxfld
      parameter (nmaxfld = 50)

      INTEGER :: krcv, ksnd   ! number of coupling fields per agrif grids
      CHARACTER(len = 8), DIMENSION(nmaxfld) :: srcv_clname, ssnd_clname   ! Coupling fields
      INTEGER, DIMENSION(nmaxfld) :: srcv_nid, ssnd_nid
      common /exchange_fields_oasis3/ krcv, ksnd,
     &srcv_clname,ssnd_clname,srcv_nid,ssnd_nid

      INTEGER :: oasis_time,oasis_runtime
      common /exchange_times_oasis3/ oasis_time,oasis_runtime
#endif /* OA_COUPLING */

