#include "cppdefs.h"
! Modified by SINTEF 2024
! ---------------------
#ifdef NBQ || CFABM
      module module_param
# include "param.h"
#ifdef CFABM
      ! == Added by SINTEF 2024
      USE fabm, only: type_fabm_model, type_fabm_interior_variable_id

      TYPE type_fabm_model_array
        CLASS (type_fabm_model), pointer :: model
      END TYPE type_fabm_model_array
      TYPE T_FABM
        TYPE(type_fabm_model_array), dimension(:), allocatable :: f(:)
        TYPE(type_fabm_interior_variable_id), allocatable :: id_temp(:)
        TYPE(type_fabm_interior_variable_id), allocatable :: id_salt(:)
        logical, allocatable :: needsT(:)
        logical, allocatable :: needsS(:)
      END TYPE T_FABM
      TYPE (T_FABM), allocatable :: FMODEL(:)

      ! Not sure that this translates well from ROMS to CROCO, the multi nesting is not trivial...
      IF (.not.allocated(FMODEL)) THEN
            allocate ( FMODEL(gridlist%nitems) )
            DO ng=1, gridlist%nitems
                  Ntiles=NP_XI*NP_ETA-1
                  allocate ( FMODEL(ng) % f(-1:Ntiles) )
                  allocate ( FMODEL(ng) % id_temp(-1:Ntiles) )
                  allocate ( FMODEL(ng) % id_salt(-1:Ntiles) )
                  allocate ( FMODEL(ng) % needsT(-1:Ntiles) )
                  allocate ( FMODEL(ng) % needsS(-1:Ntiles) )
                  DO ivar=-1,Ntiles
                        FMODEL(ng) % needsT(ivar) = .false. ! Will be .true. for FABM models that require temperature
                        FMODEL(ng) % needsS(ivar) = .false. ! Will be .true. for FABM models that require salinity
                  END DO
            END DO
      END IF
      ! ==
#endif /* CFABM */
! ---------------------
      end module module_param
#else
      module module_param_empty
      end module module_param_empty
#endif/* NBQ || CFABM */