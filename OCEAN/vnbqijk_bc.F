! $Id: vnbqijk_bc.F 1458 2014-02-03 15:01:25Z gcambon $
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
      subroutine vnbqijk_bc_tile (Istr,Iend,Jstr,Jend,grad)
#  ifdef AGRIF      
      use AGRIF_Util
      integer Istr,Iend,Jstr,Jend
      real grad(PRIVATE_2D_SCRATCH_ARRAY)
      if (AGRIF_Root()) then
        call vnbqijk_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
      else
        call vnbqijk_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
c        call vnbqijk_bc_interp_tile(Istr,Iend,Jstr,Jend)
      endif
      return
      end
!
! PARENT
!
      subroutine vnbqijk_bc_parent_tile (Istr,Iend,Jstr,Jend,grad)
#  endif
!
! Set lateral boundary conditions for V-component momentum
! qdm_nbq_a(V_index,2) for the parent grid.
!
# endif /* NBQ */
#else
# if defined NBQ && defined NBQ_IJK
!
! CHILD
!
      subroutine vnbqijk_bc_child_tile (Istr,Iend,Jstr,Jend,grad)
!
! Set lateral boundary conditions for V-component momentum
! qdm_nbq_a(V_index,2) for the child grid.
!
# endif /* NBQ */
#endif /* CHILD */
# if defined NBQ && defined NBQ_IJK
!
! Common Code
!
# include "set_obc_definitions.h"
!
!
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "nbq.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"
      integer Istr,Iend,Jstr,Jend, i,j,k,
     &        l_b_j,l_b_jp1,l_bm1_j,l_bm1_jp1,l_bm2_j,
     &        l_b_i,l_b_im1,l_bm1_i,l_bm1_im1,l_bm2_i,
     &        l_c,l_c_im1,l_c_jm1
      real    grad(PRIVATE_2D_SCRATCH_ARRAY), cff,eps,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
      double precision :: mask1_v,mask2_v
      real sum_o(PRIVATE_2D_SCRATCH_ARRAY)
!
# include "compute_auxiliary_bounds.h"
!
! Interpolations of the parent values to get vnbqbry_east or vnbqclm
!
# ifdef CHILD
      call vnbqijk_bc_interp_tile(Istr,Iend,Jstr,Jend)
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
# ifndef NS_COM_PERIODIC
!
!====================================================================
!                            SOUTHERN BC
!====================================================================
      if (SOUTHERN_EDGE) then
#  ifdef OBC_COM_SOUTH
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========

            l_b_i     = ijk2lmom_nh(i  ,Jstr  ,k,2)
            l_b_im1   = ijk2lmom_nh(i-1,Jstr  ,k,2)
            if (l_b_im1.eq.0) then
               l_b_im1 = l_b_i
               mask1_v = 0.
            else
               mask1_v = 1.
            endif
            l_bm1_i   = ijk2lmom_nh(i  ,Jstr+1,k,2)
            l_bm1_im1 = ijk2lmom_nh(i-1,Jstr+1,k,2)
            if (l_bm1_im1.eq.0) then
               l_bm1_im1 = l_bm1_i
               mask2_v = 0.
            else
               mask2_v = 1.
            endif

            grad(i,Jstr  )=(qdm_nbq_a(l_b_i  )
     &                     -qdm_nbq_a(l_b_im1)*mask1_v)
#    ifdef MASKING
     &                              *pmask(i  ,Jstr)
#    endif
            grad(i,Jstr+1)=(qdm_nbq_a(l_bm1_i  )
     &                     -qdm_nbq_a(l_bm1_im1)*mask2_v)
#    ifdef MASKING
     &                                *pmask(i,Jstr+1)
#    endif
          enddo
          do i=Istr,Iend

            l_b_i     = ijk2lmom_nh(i  ,Jstr  ,k,2)
            l_bm1_i   = ijk2lmom_nh(i  ,Jstr+1,k,2)
            l_bm2_i   = ijk2lmom_nh(i  ,Jstr+2,k,2)

           dft=qdm_nbq_a(l_bm1_i)-qdm_nbq_a(l_bm1_i)
           dfx=qdm_nbq_a(l_bm1_i)-qdm_nbq_a(l_bm2_i)

            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i,Jstr+1)+grad(i+1,Jstr+1)) .gt. 0.) then
              dfy=grad(i,Jstr+1)
            else
              dfy=grad(i+1,Jstr+1)
            endif
 
#    ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_COM_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
           qdm_nbq_a(l_b_i)=( cff*qdm_nbq_a(l_b_i  )
     &                              +cx*qdm_nbq_a(l_bm1_i)
     &                               -max(cy,0.)*grad(i  ,Jstr)
     &                               -min(cy,0.)*grad(i+1,Jstr)
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQNUDGING 
          qdm_nbq_a(l_b_i)=(1.-tau)*qdm_nbq_a(l_b_i)
#     ifdef NBQ_FRC_BRY
     &                                   +tau*vnbqbry_south(i,k)
#     else     
!    &                                    +tau*vnbqclm(i,Jstr,k)
     &                                    +tau*(
     &                      vnbqclm(i,Jstr,k)    
     &   +vbclm(i,Jstr)*hz_half(i,Jstr,k)
     &                                         )
#     endif
#    endif
#    ifdef MASKING
            qdm_nbq_a(l_b_i)=qdm_nbq_a(l_b_i)
     &           *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
            qdm_nbq_a(l_b_i)=qdm_nbq_a(l_b_i)
     &                                             *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
#   elif defined OBC_COM_NBQSPECIFIED
!                                           Southern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
            l_b_i     = ijk2lmom_nh(i  ,Jstr  ,k,2)
#    ifdef NBQ_FRC_BRY
            qdm_nbq_a(l_b_i)=vnbqbry_south(i,k)        ! specified
#    else
            qdm_nbq_a(l_b_i)=vnbqclm(i,Jstr,k)
#    endif
#    ifdef MASKING
     &                               *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(i,Jstr)
#    endif
          enddo
        enddo
#   else
!                                           Southern edge gradient BC
!                                           ======== ==== ======== ==

        sum_o=0.
        do k=1,N
          do i=Istr,Iend
!            qdmv_nbq(i,Jstr,k)=qdmv_nbq(i,Jstr+1,k)
              qdmv_nbq(i,Jstr,k)=v(i,Jstr,k,nrhs)
     &                      *0.5*(Hz(i,Jstr,k)+ Hz(i,Jstr-1,k))
#    ifdef MASKING
     &                                     *vmask(i,Jstr)
#    endif
#    ifdef WET_DRY
     &                                 *vmask_wet(i,Jstr)
#    endif
              sum_o(i,Jstr)=sum_o(i,Jstr)+qdmv_nbq(i,Jstr,k)
          enddo
        enddo
      
        do k=1,N
          do i=Istr,Iend
               qdmv_nbq(i,Jstr,k)=qdmv_nbq(i,Jstr,k)-sum_o(i,Jstr)
     &                  /float(N)
     &                +DV_nbq(i,Jstr)/float(N)
!     &                 +vbar(i,Jstr,kstp)   
!     &                      *0.5*(Hz(i,Jstr,k)+ Hz(i,Jstr-1,k))                                             
          enddo
        enddo

#   endif
!#  else
!        do k=1,N                                     ! Southern edge closed
!          do i=Istr,Iend                             ! ======== ==== ======
!            l_b_i     = ijk2lmom_nh(i  ,Jstr  ,k,2)  !  (no-flux: default)
!            qdm_nbq_a(l_b_i)=0.
!          enddo
!        enddo
#  endif              /* OBC_COM_SOUTH */
      endif         !<-- SOUTHERN_EDGE
!
!====================================================================
!                            NORTHERN BC
!====================================================================
      if (NORTHERN_EDGE) then
#  ifdef OBC_COM_NORTH
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                            ! Northern edge radiation
          do i=Istr,Iend+1                  ! ======== ==== =========

            l_b_i     = ijk2lmom_nh(i  ,Jend+1,k,2)
            l_b_im1   = ijk2lmom_nh(i-1,Jend+1,k,2)
            if (l_b_im1.eq.0) then
               l_b_im1 = l_b_i
               mask1_v = 0.
            else
               mask1_v = 1.
            endif
            l_bm1_i   = ijk2lmom_nh(i  ,Jend  ,k,2)
            l_bm1_im1 = ijk2lmom_nh(i-1,Jend  ,k,2)
            if (l_bm1_im1.eq.0) then
               l_bm1_im1 = l_bm1_i
               mask1_v = 0.
            else
               mask1_v = 1.
            endif

            grad(i,Jend  )=(qdm_nbq_a(l_bm1_j  )
     &                     -qdm_nbq_a(l_bm1_im1)*mask1_v)
#    ifdef MASKING
     &                                  *pmask(i,Jend)
#    endif
            grad(i,Jend+1)=(qdm_nbq_a(l_b_j  )
     &                     -qdm_nbq_a(l_b_im1)*mask2_v)
#    ifdef MASKING
     &                                 *pmask(i,Jend+1)
#    endif
          enddo
          do i=Istr,Iend

            l_b_i     = ijk2lmom_nh(i  ,Jend+1,k,2)
            l_bm1_i   = ijk2lmom_nh(i  ,Jend  ,k,2)
            l_bm2_i   = ijk2lmom_nh(i  ,Jend-1,k,2)

           dft=qdm_nbq_a(l_bm1_i)-qdm_nbq_a(l_bm1_i)
           dfx=qdm_nbq_a(l_bm1_i)-qdm_nbq_a(l_bm2_i)

            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i,Jend)+grad(i+1,Jend)) .gt. 0.) then
              dfy=grad(i,Jend)
            else
              dfy=grad(i+1,Jend)
            endif
 
#    ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_COM_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
           qdm_nbq_a(l_b_i)=( cff*qdm_nbq_a(l_b_i  )
     &                              +cx*qdm_nbq_a(l_bm1_i)
     &                             -max(cy,0.)*grad(i  ,Jend+1)
     &                             -min(cy,0.)*grad(i+1,Jend+1)
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQNUDGING 
          qdm_nbq_a(l_b_i)=(1.-tau)*qdm_nbq_a(l_b_i)
#     ifdef NBQ_FRC_BRY
     &                                   +tau*vnbqbry_north(i,k)
#     else     
!    &                                  +tau*vnbqclm(i,Jend+1,k)
     &                                  +tau* (
     &              vnbqclm(i,Jend+1,k)  
     &            + vbclm(i,Jend+1)*hz_half(i,Jend+1,k)
     &                                        )       
#     endif
#    endif
#    ifdef MASKING
            qdm_nbq_a(l_b_i)=qdm_nbq_a(l_b_i)
     &        *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
            qdm_nbq_a(l_b_i)=qdm_nbq_a(l_b_i)
     &                                             *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                           Northern edge Specified BC
!                                           ======== ==== ========= ==
        do k=1,N
          do i=Istr,Iend
            l_b_i     = ijk2lmom_nh(i  ,Jend+1,k,2)
#    ifdef NBQ_FRC_BRY
            qdm_nbq_a(l_b_i)=vnbqbry_north(i,k)      ! specified
#    else
            qdm_nbq_a(l_b_i)=vnbqclm(i,Jend+1,k)
#    endif
#    ifdef MASKING
     &                               *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(i,Jend+1)
#    endif
          enddo
        enddo
#   else
        sum_o=0.
        do k=1,N
          do i=Istr,Iend
!                                           Northern edge gradient BC
!                                           ======== ==== ======== ==
!            qdmv_nbq(i,Jend+1,k)=qdmv_nbq(i,Jend,k)
             qdmv_nbq(i,Jend+1,k)=v(i,Jend+1,k,nrhs)
     &                      *0.5*(Hz(i,Jstr,k)+ Hz(i,Jstr-1,k))
#    ifdef MASKING
     &                                   *vmask(i,Jend+1)
#    endif
#    ifdef WET_DRY
     &                               *vmask_wet(i,Jend+1)
#    endif
              sum_o(i,Jend+1)=sum_o(i,Jend+1)+qdmv_nbq(i,Jend+1,k)
          enddo
        enddo

        do k=1,N
          do i=Istr,Iend
               qdmv_nbq(i,Jend+1,k)=qdmv_nbq(i,Jend+1,k)-sum_o(i,Jend+1)
     &                  /float(N)
     &                  +DV_nbq(i,Jend+1)/float(N)
!     &                 +vbar(i,Jend+1,kstp)   
!     &                      *0.5*(Hz(i,Jend+1,k)+ Hz(i,Jend,k))
!     &                                             
          enddo
        enddo
      
#   endif
!#  else
!        do k=1,N                                    ! Northern edge closed
!          do i=Istr,Iend                            ! ======== ==== ======
!            l_b_i     = ijk2lmom_nh(i  ,Jend+1,k,2) !   (no-flux: default)
!            qdm_nbq_a(l_b_i)=0.
!          enddo
!        enddo
#  endif
      endif     !<--  NORTHERN_EDGE
# endif          /* !NS_COM_PERIODIC */
 
# ifndef EW_COM_PERIODIC
!
!====================================================================
!                            WESTERN BC
!====================================================================
      if (WESTERN_EDGE) then
#  ifdef OBC_COM_WEST
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                             ! Western edge radiation
          do j=JstrV-1,Jend                  ! ======= ==== =========

            l_b_j     = ijk2lmom_nh(Istr-1,j  ,k,2)
            l_b_jp1   = ijk2lmom_nh(Istr-1,j+1,k,2)
            l_bm1_j   = ijk2lmom_nh(Istr  ,j  ,k,2)
            l_bm1_jp1 = ijk2lmom_nh(Istr  ,j+1,k,2)

            grad(Istr-1,j)=qdm_nbq_a(l_b_jp1  )
     &                    -qdm_nbq_a(l_b_j    )
            grad(Istr  ,j)=qdm_nbq_a(l_bm1_jp1)
     &                    -qdm_nbq_a(l_bm1_j  )
          enddo
          do j=JstrV,Jend

            l_b_j     = ijk2lmom_nh(Istr-1,j  ,k,2)
            l_bm1_j   = ijk2lmom_nh(Istr  ,j  ,k,2)
            l_bm2_j   = ijk2lmom_nh(Istr+1,j  ,k,2)

           dft=qdm_nbq_a(l_bm1_j)  -qdm_nbq_a(l_bm1_j)
           dfx=qdm_nbq_a(l_bm1_j)  -qdm_nbq_a(l_bm2_j)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(Istr,j-1)+grad(Istr,j)) .gt. 0.) then
              dfy=grad(Istr,j-1)
            else
              dfy=grad(Istr,j  )
            endif
 
#    ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_COM_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
          qdm_nbq_a(l_b_j)=( cff*qdm_nbq_a(l_b_j  )
     &                              +cx*qdm_nbq_a(l_bm1_j)
     &                             -max(cy,0.)*grad(Istr-1,j-1)
     &                             -min(cy,0.)*grad(Istr-1,j  )
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQNUDGING 
          qdm_nbq_a(l_b_j)=(1.-tau)*qdm_nbq_a(l_b_j)
#     ifdef NBQ_FRC_BRY
     &                                    +tau*vnbqbry_west(j,k)
#     else
!    &                                  +tau*vnbqclm(Istr-1,j,k)
     &         +tau*v(Istr-1,j,k,nnew)
#     endif
#    endif
#    ifdef MASKING
            qdm_nbq_a(l_b_j)=qdm_nbq_a(l_b_j)
     &              *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
            qdm_nbq_a(l_b_j)=qdm_nbq_a(l_b_j)
     &                                             *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Western edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=JstrV,Jend
            l_b_j     = ijk2lmom_nh(Istr-1,j  ,k,2)
#    ifdef NBQ_FRC_BRY
            qdm_nbq_a(l_b_j)=vnbqbry_west(j,k)       ! specified
#    else
            qdm_nbq_a(l_b_j)=vnbqclm(Istr-1,j,k)
#    endif
#    ifdef MASKING
     &                               *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(Istr-1,j)
#    endif
          enddo
        enddo
#   else
!                                            Western edge gradient BC
!                                            ======= ==== ======== ==
!       sum_o=0.
        do k=1,N
          do j=JstrV,Jend
!           qdmv_nbq(Istr-1,j,k)=v(Istr-1,j,k,nrhs)
!    &                      *0.5*(Hz(Istr,j,k)+ Hz(Istr,j-1,k))
            qdmv_nbq(Istr-1,j,k)=qdmv_nbq(Istr,j,k)
#    ifdef MASKING
     &                                   *vmask(Istr-1,j)
#    endif
#    ifdef WET_DRY
     &                               *vmask_wet(Istr-1,j)
#    endif

!            sum_o(Istr-1,j)=sum_o(Istr-1,j)+qdmv_nbq(Istr-1,j,k)
!          
          enddo
        enddo    
!        do k=1,N
!          do j=JstrV,Jend
!               qdmv_nbq(Istr-1,j,k)=qdmv_nbq(Istr-1,j,k)-sum_o(Istr-1,j)
!     &                  /float(N)
!     &                 +vbar(Istr-1,j,kstp)   
!     &                 *0.5*(Hz(Istr,j,k)+ Hz(Istr,j-1,k))
!     &                                             
!          enddo
!        enddo  
#   endif
!#  else
!#   ifdef NS_COM_PERIODIC
!#    define J_RANGE JstrV,Jend
!#   else
!#    define J_RANGE Jstr,JendR
!#   endif
!        do k=1,N                        ! Wall: free-slip (gamma2=+1)
!          do j=J_RANGE                  ! =====   no-slip (gamma2=-1)
!            l_b_j     = ijk2lmom_nh(Istr-1,j  ,k,2)
!            l_bm1_j   = ijk2lmom_nh(Istr  ,j  ,k,2)
!            qdm_nbq_a(l_b_j)=gamma2*qdm_nbq_a(l_bm1_j)
!#   ifdef MASKING
!     &                                          *vmask(Istr-1,j)
!#   endif
!          enddo
!        enddo
!#   undef J_RANGE
#  endif
      endif          !<-- WESTERN_EDGE
!
!====================================================================
!                            EASTERN BC
!====================================================================
      if (EASTERN_EDGE) then
#  ifdef OBC_COM_EAST
#   ifdef OBC_COM_NBQORLANSKI
        do k=1,N                             ! Eastern edge radiation
          do j=JstrV-1,Jend                  ! ======= ==== =========

            l_b_j     = ijk2lmom_nh(Iend+1,j  ,k,2)
            l_b_jp1   = ijk2lmom_nh(Iend+1,j+1,k,2)
            l_bm1_j   = ijk2lmom_nh(Iend  ,j  ,k,2)
            l_bm1_jp1 = ijk2lmom_nh(Iend  ,j+1,k,2)

            grad(Iend+1,j)=qdm_nbq_a(l_b_jp1  )
     &                    -qdm_nbq_a(l_b_j    )
            grad(Iend  ,j)=qdm_nbq_a(l_bm1_jp1)
     &                    -qdm_nbq_a(l_bm1_j  )
          enddo
          do j=JstrV,Jend

            l_b_j     = ijk2lmom_nh(Iend+1,j  ,k,2)
            l_bm1_j   = ijk2lmom_nh(Iend  ,j  ,k,2)
            l_bm2_j   = ijk2lmom_nh(Iend-1,j  ,k,2)

            dft=qdm_nbq_a(l_bm1_j)  -qdm_nbq_a(l_bm1_j)
            dfx=qdm_nbq_a(l_bm1_j)  -qdm_nbq_a(l_bm2_j)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined NBQ_FRC_BRY || defined NBQNUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(Iend,j-1)+grad(Iend,j)) .gt. 0.) then
              dfy=grad(Iend,j-1)
            else
              dfy=grad(Iend,j  )
            endif
 
#    ifdef OBC_COM_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_COM_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
           qdm_nbq_a(l_b_j)=( cff*qdm_nbq_a(l_b_j)
     &                              +cx*qdm_nbq_a(l_bm1_j)
     &                             -max(cy,0.)*grad(Iend+1,j-1)
     &                             -min(cy,0.)*grad(Iend+1,j  )
     &                                               )/(cff+cx)
#    if defined NBQ_FRC_BRY  || defined NBQNUDGING 
          qdm_nbq_a(l_b_j)=(1.-tau)*qdm_nbq_a(l_b_j)
#     ifdef NBQ_FRC_BRY
     &                                    +tau*vnbqbry_east(j,k)
#     else     
!    &                                  +tau*vnbqclm(Iend+1,j,k)
     &         +tau*v(Iend+1,j,k,nnew)
#     endif
#    endif
#    ifdef MASKING
            qdm_nbq_a(l_b_j)=qdm_nbq_a(l_b_j)
     &                     *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
           qdm_nbq_a(l_b_j)=qdm_nbq_a(l_b_j)
     &                                             *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
!
#   elif defined OBC_COM_NBQSPECIFIED
!                                            Eastern edge Specified BC
!                                            ======= ==== ========= ==
        do k=1,N
          do j=Jstr,Jend
            l_b_j     = ijk2lmom_nh(Iend+1,j  ,k,2)
#    ifdef NBQ_FRC_BRY
            qdm_nbq_a(l_b_j)=vnbqbry_east(j,k)       ! specified
#    else
            qdm_nbq_a(l_b_j)=vnbqclm(Iend+1,j,k)
#    endif
#    ifdef MASKING
     &                               *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                           *vmask_wet(Iend+1,j)
#    endif
          enddo
        enddo
#   else
!                                            Eastern edge gradient BC
!                                            ======= ==== ======== ==
        sum_o=0.
        do k=1,N
          do j=Jstr,Jend
!           qdmv_nbq(Iend+1,j,k)=v(Iend+1,j,k,nrhs)
!    &                      *0.5*(Hz(Iend+1,j,k)+ Hz(Iend+1,j-1,k))
            qdmv_nbq(Iend+1,j,k)=qdmv_nbq(Iend,j,k)
#    ifdef MASKING
     &                                  *vmask(Iend+1,j)
#    endif
#    ifdef WET_DRY
     &                              *vmask_wet(Iend+1,j)
#    endif
!            sum_o(Iend+1,j)=sum_o(Iend+1,j)+qdmv_nbq(Iend+1,j,k)
!          
          enddo
        enddo      
!        do k=1,N
!          do j=JstrV,Jend
!               qdmv_nbq(Iend+1,j,k)=qdmv_nbq(Iend+1,j,k)-sum_o(Iend+1,j)
!     &                  /float(N)
!     &                 +vbar(Iend+1,j,kstp)   
!     &                      *0.5*(Hz(Iend+1,j,k)+ Hz(Iend+1,j-1,k))
!     &                                             
!          enddo
!        enddo  
#   endif
!#  else
!#   ifdef NS_COM_PERIODIC
!#    define J_RANGE JstrV,Jend
!#   else
!#    define J_RANGE Jstr,JendR
!#   endif
!        do k=1,N                        ! Wall: free-slip (gamma2=+1)
!          do j=J_RANGE                  ! ====    no-slip (gamma2=-1)
!            l_b_j     = ijk2lmom_nh(Iend+1,j  ,k,2)
!            l_bm1_j   = ijk2lmom_nh(Iend  ,j  ,k,2)
!            qdm_nbq_a(l_b_j)=gamma2*qdm_nbq_a(l_bm1_j)
!#   ifdef MASKING
!     &                                          *vmask(Iend+1,j)
!#   endif
!          enddo
!        enddo
!#   undef J_RANGE
#  endif
      endif     !<-- EASTERN_EDGE
# endif          /* !EW_COM_PERIODIC */
 
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
      
# if defined OBC_COM_SOUTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Istr-1,Jstr,k)=0.5*(qdmv_nbq(Istr-1,Jstr+1,k)
     &                                +qdmv_nbq(Istr,Jstr,k))
#  ifdef MASKING
     &                                  *vmask(Istr-1,Jstr)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_SOUTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Iend+1,Jstr,k)=0.5*(qdmv_nbq(Iend,Jstr+1,k)
     &                                +qdmv_nbq(Iend,Jstr,k))
#  ifdef MASKING
     &                                  *vmask(Iend+1,Jstr)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Istr-1,Jend+1,k)=0.5*(qdmv_nbq(Istr-1,Jend,k)
     &                                  +qdmv_nbq(Istr  ,Jend+1,k))     
#  ifdef MASKING
     &                                *vmask(Istr-1,Jend+1)
#  endif
        enddo
      endif
# endif
# if defined OBC_COM_NORTH && defined OBC_COM_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          qdmv_nbq(Iend+1,Jend+1,k)=0.5*(qdmv_nbq(Iend+1,Jend,k)
     &                                  +qdmv_nbq(Iend  ,Jend+1,k))          
#  ifdef MASKING
     &                                *vmask(Iend+1,Jend+1)
#  endif
        enddo
      endif
# endif

!# undef nnew
!# undef nstp

      return
      end
#else
# ifndef CHILD
      subroutine vnbqijk_bc_parent_empty
      end
# else
      subroutine vnbqijk_bc_child_empty
      end
# endif
#endif /* NBQ */
#ifndef CHILD
# define CHILD
# ifdef AGRIF
#  include "vnbqijk_bc.F"
# endif
# undef CHILD
#endif  /* !CHILD */
 
