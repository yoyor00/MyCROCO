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
#ifdef BACKSCATTER
      subroutine uv3dmix_BS (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
      call uv3dmix_BS_tile (Istr,Iend,Jstr,Jend,
     &                      A2d(1,1,trd), A2d(1,2,trd),
     &                      A2d(1,3,trd), A2d(1,4,trd))
      return
      end
!
!---------------------------------------------------------------------
!
      subroutine uv3dmix_BS_tile (Istr,Iend,Jstr,Jend, 
     &                            UFx,UFe,VFx,VFe)
!
!--------------------------------------------------------------------
! Computes antidiffusion of momentum, permitting energy backscatter
! from unresolved to resolved scales. 
!
! This is done with harmonic mixing along constant S-surfaces as 
! horizontal divergence of the stress tensor of components:
!
!                         du     dv
!         s_xx = -s_yy = ---- - -----
!                         dx     dy
!
!                         du     dv
!         s_xy =  s_yx = ---- + ----
!                         dy     dx
!
! References:
!
! Wajsowicz, R.C, 1993: A consistent formulation of the anisotropic
! stress tensor for use in models of the large-scale ocean
! circulation, JCP, 105, 333-338.
!
! Jansen, M.F., Held, I.M., 2014. Parameterizing subgrid-scale eddy 
! effects using energetically consistent backscatter. Ocean Model, 
! 80, 36–48
!
! Juricke S., Danilov S., Kutsenkob A., Oliver M., 2019: Ocean kinetic 
! energy backscatter parametrizations on unstructured grids. Ocean Model,
! 138, 51-67
!---------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, indx
# include "param.h"
      real UFe(PRIVATE_2D_SCRATCH_ARRAY), cff,
     &     UFx(PRIVATE_2D_SCRATCH_ARRAY), cff1,
     &     VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &     VFx(PRIVATE_2D_SCRATCH_ARRAY)
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
# include "penalization.h"
# include "diagnostics.h"
# include "diags_ek.h"
!
# include "compute_auxiliary_bounds.h"
!
      indx=3-nstp    !--> time index for target arrays;
!
!$acc kernels if(compute_on_device) default(present)
!
!---------------------------------------------------------------------
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!---------------------------------------------------------------------
!
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            cff=Hz(i,j,k)*visc3d_BS_r(i,j,k)*
     &          ( pmon_r(i,j)*( pn_u(i+1,j)*u(i+1,j,k,nstp)
     &                         -pn_u(i  ,j)*u(i  ,j,k,nstp) )
     &           -pnom_r(i,j)*( pm_v(i,j+1)*v(i,j+1,k,nstp)
     &                         -pm_v(i,j  )*v(i,j  ,k,nstp) ) )
            UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=Istr,Iend+1
            cff=0.25*visc3d_BS_p(i,j,k)*
     &                        (Hz(i-1,j  ,k)+Hz(i,j  ,k)+
     &                         Hz(i-1,j-1,k)+Hz(i,j-1,k))*
     &      (  pmon_p(i,j)
     &                     *( pn_v(i  ,j)*v(i  ,j,k,nstp)
     &                       -pn_v(i-1,j)*v(i-1,j,k,nstp) )
     &        +pnom_p(i,j)
     &                     *( pm_u(i,j  )*u(i,j  ,k,nstp)
     &                       -pm_u(i,j-1)*u(i,j-1,k,nstp) )  )
     &                                       SWITCH pmask(i,j)
            UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          enddo
        enddo
!
!---------------------------------------------------------------------
!  Apply antiviscous terms.
!
!  At this stage arrays u,v(...,indx) contain Hz*U and Hz*V with units 
!  of [m2/s]. Also compute vertical integral of antiviscous terms and 
!  add it into coupling terms for the barotropic mode. 

!  Finally, store antiviscous terms ~U*Hz/(pm*pn) for use in the 
!  unresolved KE equation (backscatter.F).
!---------------------------------------------------------------------
!
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=pn_u(i,j)*(UFx(i,j  )-UFx(i-1,j))
     &         +pm_u(i,j)*(UFe(i,j+1)-UFe(i  ,j))
# ifdef PERMEABILITY
            cff=cff*(1.-pena_u(i,j,k))
# endif
            cff1=pm_u(i,j)*pn_u(i,j)*cff

            rufrc(i,j)=rufrc(i,j) + cff
            u(i,j,k,indx)=u(i,j,k,indx) + dt*cff1

            KE_Hback_U(i,j,k)=cff*u(i,j,k,nstp) SWITCH umask(i,j)

# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL \
                            || defined DIAGNOSTICS_PV
            MHmix(i,j,k,1,indx) = cff SWITCH umask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                            + MHmix(i,j,k,1,indx)
#  endif
# elif defined DIAGNOSTICS_EK
            if (k.eq.1) then
              ekwrkHmix(i,j,1,indx) =
     &                         cff*u(i,j,k,nstp) SWITCH umask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                                     + ekwrkHmix(i,j,1,indx)
#  endif
            else
              ekwrkHmix(i,j,1,indx) = ekwrkHmix(i,j,1,indx)
     &                       + cff*u(i,j,k,nstp) SWITCH umask(i,j)
            endif
#  if defined DIAGNOSTICS_EK_MLD
            if (k.eq.kbl(i,j)) then
              ekwrkHmix_mld(i,j,1,indx) =
     &                          cff*u(i,j,k,nstp) SWITCH umask(i,j)
#   if defined UV_VIS4 || defined SPONGE_VIS2
     &                                  + ekwrkHmix_mld(i,j,1,indx)
#   endif
            elseif (k.gt.kbl(i,j)) then
              ekwrkHmix_mld(i,j,1,indx) = ekwrkHmix_mld(i,j,1,indx)
     &                        + cff*u(i,j,k,nstp) SWITCH umask(i,j)
            endif
#  endif
# endif
# if defined DIAGNOSTICS_VRT && !defined DIAGNOSTICS_UV
            if (k.eq.1) then
              wrkHmix(i,j,1,indx) = cff SWITCH umask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                              + wrkHmix(i,j,1,indx)
#  endif
            else
              wrkHmix(i,j,1,indx) = wrkHmix(i,j,1,indx)
     &                              + cff SWITCH umask(i,j)
            endif
# endif
          enddo
        enddo

        do j=JstrV,Jend
          do i=Istr,Iend
            cff=pn_v(i,j)*(VFx(i+1,j)-VFx(i,j  ))
     &         -pm_v(i,j)*(VFe(i  ,j)-VFe(i,j-1))
# ifdef PERMEABILITY
            cff=cff*(1.-pena_v(i,j,k))
# endif
            cff1=pm_v(i,j)*pn_v(i,j)*cff

            rvfrc(i,j)=rvfrc(i,j) + cff
            v(i,j,k,indx)=v(i,j,k,indx) + dt*cff1

            KE_Hback_V(i,j,k)=cff*v(i,j,k,nstp) SWITCH vmask(i,j)

# if defined DIAGNOSTICS_UV || defined DIAGNOSTICS_EK_FULL \
                            || defined DIAGNOSTICS_PV
            MHmix(i,j,k,2,indx) = cff SWITCH vmask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                                      + MHmix(i,j,k,2,indx)
#  endif
# elif defined DIAGNOSTICS_EK
            if (k.eq.1) then
              ekwrkHmix(i,j,2,indx) =
     &                        cff*v(i,j,k,nstp) SWITCH vmask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                                    + ekwrkHmix(i,j,2,indx)
#  endif
            else
                ekwrkHmix(i,j,2,indx) = ekwrkHmix(i,j,2,indx)
     &                      + cff*v(i,j,k,nstp) SWITCH vmask(i,j)
            endif
#  if defined DIAGNOSTICS_EK_MLD
              if (k.eq.kbl(i,j)) then
                ekwrkHmix_mld(i,j,2,indx) =
     &                        cff*v(i,j,k,nstp) SWITCH vmask(i,j)
#   if defined UV_VIS4 || defined SPONGE_VIS2
     &                                + ekwrkHmix_mld(i,j,2,indx)
#   endif
              elseif (k.gt.kbl(i,j)) then
                ekwrkHmix_mld(i,j,2,indx) = ekwrkHmix_mld(i,j,2,indx)
     &                       +cff*v(i,j,k,nstp) SWITCH vmask(i,j)
              endif
#  endif
# endif
# if defined DIAGNOSTICS_VRT && !defined DIAGNOSTICS_UV
              if (k.eq.1) then
                wrkHmix(i,j,2,indx) = cff WITCH vmask(i,j)
#  if defined UV_VIS4 || defined SPONGE_VIS2
     &                               + wrkHmix(i,j,2,indx)
#  endif
              else
                wrkHmix(i,j,2,indx) = wrkHmix(i,j,2,indx)
     &                             + cff SWITCH vmask(i,j)
              endif
# endif
          enddo
        enddo
      enddo ! <-- k

!$acc end kernels

!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
!
# if defined BD_EXCHANGE
      call exchange_u3d_tile (istr,iend,jstr,jend,KE_Hback_U)
      call exchange_v3d_tile (istr,iend,jstr,jend,KE_Hback_V)
# endif
!======================================================================
      return
      end
!======================================================================
#else
      subroutine uv3dmix_BS_empty (tile)
      end
#endif  /* BACKSCATTER */
!
