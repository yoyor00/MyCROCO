!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"

      subroutine sto_initial (tile)
      implicit none
#include "param.h"
       integer tile, trd
#include "compute_tile_bounds.h"
      call sto_initial_tile (Istr,Iend,Jstr,Jend)
      return
      end
!
      subroutine sto_initial_tile (Istr,Iend,Jstr,Jend)
!
!----------------------------------------------------------------------
!  Stochastically perturbe initial conditions 
!----------------------------------------------------------------------
!
#ifdef AGRIF
      use Agrif_Util
#endif
#ifdef SUBSTANCE
      USE comsubstance, ONLY: cini_wat
#endif
      USE stoarray        ! Stochastic arrays and options
      USE stoalloc        ! Stochastic working arrays
      USE stoics
      implicit none
#include "param.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "scalars.h"
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
!
#include "compute_auxiliary_bounds.h"
!
#if defined EW_PERIODIC && !defined MPI
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
#endif
#if defined NS_PERIODIC && !defined MPI
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#endif


!$acc kernels if(compute_on_device) default(present)
!
!----------------------------------------------------------------------
!  Initial conditions for free surface and 2D momentum components.
!----------------------------------------------------------------------
!

#ifdef SOLVE3D
!
!----------------------------------------------------------------------
!  Initial conditions for momentum components [m/s].
!----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
!  Set initial conditions for potential temperature [degC] and
!  salinity [PSU] and passive tracer
!----------------------------------------------------------------------
!
#   if defined TEMPERATURE
        call sto_ics( t(:,:,:,1,itemp), sto_xi3d(:,:,:) )
#   endif
!
#endif /*SOLVE3D */

#ifdef SOLVE3D
!
!--------------------------------------------------------------------
!  Initial conditions for passive tracer
!--------------------------------------------------------------------
!
!
!
!--------------------------------------------------------------------
!  Initial conditions for tracer-Substance
!--------------------------------------------------------------------
!
!----------------------------------------------------------------------
!  Initial conditions for sediment tracer variables.
!----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
!  Initial conditions for skin temperature
!----------------------------------------------------------------------
!
!
!--------------------------------------------------------------------
!  Initial conditions for biology tracer variables.
!--------------------------------------------------------------------
!
#endif /* SOLVE3D */
!
!--------------------------------------------------------------------
!  Exchange boundary information
!--------------------------------------------------------------------
!

#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
# ifdef SOLVE3D
#  ifdef TEMPERATURE
#   ifdef THREE_GHOST_POINTS_TS
      call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             t(START_2D_ARRAY,1,1,itemp))
      call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             t(START_2D_ARRAY,1,2,itemp))
#   else
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        t(START_2D_ARRAY,1,1,itemp))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        t(START_2D_ARRAY,1,2,itemp))
#   endif
#  endif /* TRACERS */
# endif
#endif
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE

      return
      end
