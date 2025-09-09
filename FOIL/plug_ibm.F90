#include "cppdefs.h"

#if defined DEB_IBM 

 MODULE plug_ibm
   ! interface between croco and ibm module

   USE module_ibm
   USE ibm,        ONLY : ibm_init, ibm_3d

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ibm_init_main
   PUBLIC :: ibm_update_main

  ! ====================================================================
  CONTAINS


  SUBROUTINE ibm_init_main(tile)

   INTEGER :: tile
# include "compute_tile_bounds.h"

   CALL ibm_init(zeta,t(:,:,:,nstp,isalt),t(:,:,:,nstp,itemp),Istr,Iend,Jstr,Jend)

  END SUBROUTINE


  !======================================================================
  SUBROUTINE ibm_update_main(tile)
     
   INTEGER :: tile
# include "compute_tile_bounds.h"
   
   CALL ibm_3d(zeta,u,v,t(:,:,:,nstp,isalt),t(:,:,:,nstp,itemp),Istr,Iend,Jstr,Jend)

  END SUBROUTINE



  !=========================================================================
 END MODULE plug_ibm

#else

 MODULE plug_ibm_empty
 END MODULE plug_ibm_empty

#endif /* DEB_IBM */