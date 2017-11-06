! $Id: hnbq 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#ifndef CHILD
!
# include "cppdefs.h"
# if defined NBQ && defined NBQ_IJK
      subroutine hnbq_bc_tile (Istr,Iend,Jstr,Jend)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      if (AGRIF_Root()) then
        call hnbq_bc_parent_tile (Istr,Iend,Jstr,Jend)
      else
        call hnbq_bc_child_tile (Istr,Iend,Jstr,Jend)
      endif
      return
      end
!
! PARENT
!
      subroutine hnbq_bc_parent_tile (Istr,Iend,Jstr,Jend)
#  endif
!
! Set lateral boundary conditions for U-component momentum
! qdm_nbq_a(U_index,2) for the parent grid.
!
# endif /* NBQ */
#else
# if defined NBQ && defined NBQ_IJK
!
! CHILD
!
      subroutine hnbq_bc_child_tile (Istr,Iend,Jstr,Jend)
!
! Set lateral boundary conditions for U-component momentum
! qdm_nbq_a(U_index,2) for the child grid.
!
# endif /* NBQ */
#endif /* CHILD */
#if defined NBQ && defined NBQ_IJK
!
! Common Code
!
# include "set_obc_definitions.h"
!
      USE module_nh
      USE module_nbq
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"  
# ifdef MRL_WCI
#  include "forces.h"
# endif
#include "nbq.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get unbqbry_east or unbqclm
!
# ifdef CHILD
!      call unbqijk_bc_interp_tile(Istr,Iend,Jstr,Jend)
# endif
!
# if defined NBQ_FRC_BRY || defined NBQNUDGING
      tau_in=dt*tauM_in
      tau_out=dt*tauM_out
# endif
!
!# define nnew 2
!# define nstp 1
!# define OBC_COM_NBQORLANSKI
!
# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
#  ifdef OBC_COM_WEST
        do k=1,N
          do j=Jstr,Jend
!            Hz (Istr-1,j,k)=Hz (Istr,j,k)
!#  ifdef MASKING
!     &                     *rmask(Istr-1,j)
!#  endif
            Hzr(Istr-1,j,k)=Hzr(Istr,j,k)
#  ifdef MASKING
     &                     *rmask(Istr-1,j)
#  endif
          enddo
        enddo
#  endif              /* OBC_COM_WEST */
      endif         !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
#  ifdef OBC_COM_EAST
        do k=1,N
          do j=Jstr,Jend
!            Hz (Iend+1,j,k)=Hz (Iend,j,k)
!#  ifdef MASKING
!     &                     *rmask(Iend+1,j)
!#  endif
            Hzr(Iend+1,j,k)=Hzr(Iend,j,k)
#  ifdef MASKING
     &                     *rmask(Iend+1,j)
#  endif
          enddo
        enddo
#  endif            /* OBC_COM_EAST */
      endif         !<-- EASTERN_EDGE
 
# endif             /* !EW_COM_PERIODIC */

# ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
#  ifdef OBC_COM_SOUTH
        do k=1,N
          do i=Istr,Iend
            Hz (i,Jstr-1,k)=Hz (i,Jstr,k)
!#  ifdef MASKING
!     &                     *rmask(i,Jstr-1)
!#  endif
            Hzr(i,Jstr-1,k)=Hzr(i,Jstr,k)
#  ifdef MASKING
     &                     *rmask(i,Jstr-1)
#  endif
          enddo
        enddo
#  endif
      endif         !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
#  ifdef OBC_COM_NORTH
        do k=1,N
          do i=Istr,Iend
!            Hz (i,Jend+1,k)=Hz (i,Jend,k)
!#  ifdef MASKING
!     &                     *rmask(i,Jend+1)
!#  endif
            Hzr(i,Jend+1,k)=Hzr(i,Jend,k)
#  ifdef MASKING
     &                     *rmask(i,Jend+1)
#  endif
          enddo
        enddo
#  endif
      endif         !<-- NORTHERN_EDGE

# endif        /* !NS_COM_PERIODIC */
 

                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
!          Hz(Istr,Jstr-1,k)=0.5*(Hz(Istr,Jstr,k)
!     &                                +Hz(Istr+1,Jstr-1,k))
!#  ifdef MASKING
!     &                                  *rmask(Istr,Jstr-1)
!#  endif
          Hzr(Istr,Jstr-1,k)=0.5*(Hzr(Istr,Jstr,k)
     &                                +Hzr(Istr+1,Jstr-1,k))
#  ifdef MASKING
     &                                  *rmask(Istr,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
!          Hz(Iend+1,Jstr-1,k)=0.5*(Hz(Iend+1,Jstr,k)
!     &                                  +Hz(Iend,Jstr-1,k))
!#  ifdef MASKING
!     &                                *rmask(Iend+1,Jstr-1)
!#  endif
          Hzr(Iend+1,Jstr-1,k)=0.5*(Hzr(Iend+1,Jstr,k)
     &                                  +Hzr(Iend,Jstr-1,k))
#  ifdef MASKING
     &                                *rmask(Iend+1,Jstr-1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
!          Hz(Istr,Jend+1,k)=0.5*(Hz(Istr+1,Jend+1,k)
!     &                                  +Hz(Istr,Jend,k))
!#  ifdef MASKING
!     &                                  *rmask(Istr,Jend+1)
!#  endif
          Hzr(Istr,Jend+1,k)=0.5*(Hzr(Istr+1,Jend+1,k)
     &                                  +Hzr(Istr,Jend,k))
#  ifdef MASKING
     &                                  *rmask(Istr,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
!          Hz(Iend+1,Jend+1,k)=0.5*(Hz(Iend,Jend+1,k)
!     &                                  +Hz(Iend+1,Jend,k)) 
!#  ifdef MASKING
!     &                                *rmask(Iend+1,Jend+1)
!#  endif
          Hzr(Iend+1,Jend+1,k)=0.5*(Hzr(Iend,Jend+1,k)
     &                                  +Hzr(Iend+1,Jend,k)) 
#  ifdef MASKING
     &                                *rmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif


      return
      end
#else
# ifndef CHILD
      subroutine hnbq_bc_parent_empty
      end
# else
      subroutine hnbq_bc_child_empty
      end
# endif
#endif /* NBQ */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "hnbq_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
