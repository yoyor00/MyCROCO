!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
! This is include file "bottom_drag_2d.h" 
! 
! RDRG2_2D   bottom drag coefficient value (without dimension)
!
# ifdef BOTTOM_DRAG_2D
      real RDRG2_2D(GLOBAL_2D_ARRAY)
      common /bottom_drag_RDRG2_2D/RDRG2_2D 
# endif



