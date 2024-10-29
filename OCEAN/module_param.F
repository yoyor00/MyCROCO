#include "cppdefs.h"
#ifdef NBQ
      module module_param
# include "param.h"
! Modified by SINTEF 2024
! ---------------------
#ifdef CFABM
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

      IF (.not.allocated(FMODEL)) THEN
            allocate ( FMODEL(Ngrids) )
            DO ng=1,Ngrids
              Ntiles=NtileI(ng)*NtileJ(ng)-1
              allocate ( FMODEL(ng) % f(-1:Ntiles) )
              allocate ( FMODEL(ng) % id_temp(-1:Ntiles) )
              allocate ( FMODEL(ng) % id_salt(-1:Ntiles) )
              allocate ( FMODEL(ng) % needsT(-1:Ntiles) )
              allocate ( FMODEL(ng) % needsS(-1:Ntiles) )
              DO ivar=-1,Ntiles
                FMODEL(ng) % needsT(ivar) = .false.
                FMODEL(ng) % needsS(ivar) = .false.
              END DO
            END DO
          END IF
#endif
! ---------------------
      end module module_param
#else
      module module_param_empty
      end module module_param_empty
#endif