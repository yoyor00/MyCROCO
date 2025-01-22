      use Agrif_Util



#include "cppdefs.h"
! Modified by SINTEF 2024
! ---------------------
#ifdef NBQ || CFABM
      module module_param
# include "param.h"
#ifdef CFABM







#ifdef AGRIF
      Type(Agrif_pgrid),pointer  :: parcours
#endif

#ifdef AGRIF
!
!  Close the netcdf files also for the child grids
!
      parcours=>Agrif_Curgrid%child_list % first
      do while (associated(parcours))
        Call Agrif_Instance(parcours % gr)
        call closecdf
# ifdef XIOS
        call iom_context_finalize( "crocox") ! needed for XIOS+AGRIF
# endif
        parcours => parcours % next
      enddo
#endif






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
      TYPE (T_FABM), allocatable :: FMODEL

      ! Not sure that this translates well from ROMS to CROCO, the multi nesting is not trivial...
      IF (.not.allocated(FMODEL)) THEN
            allocate ( FMODEL )
            Ntiles=NP_XI*NP_ETA-1         !TODO: or is it: Ntiles=NSUB_X*NSUB_E-1
            allocate ( FMODEL % f(-1:Ntiles) )
            allocate ( FMODEL % id_temp(-1:Ntiles) )
            allocate ( FMODEL % id_salt(-1:Ntiles) )
            allocate ( FMODEL % needsT(-1:Ntiles) )
            allocate ( FMODEL % needsS(-1:Ntiles) )
            DO ivar=-1,Ntiles
                  FMODEL % needsT(ivar) = .false. ! Will be .true. for FABM models that require temperature
                  FMODEL % needsS(ivar) = .false. ! Will be .true. for FABM models that require salinity
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