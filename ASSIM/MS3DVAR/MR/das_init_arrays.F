#include "cppdefs.h"

      subroutine das_init_arrays (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_init_arrays_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_init_arrays_tile (Istr,Iend,Jstr,Jend)
!
! Initialize (first touch) globally accessable arrays. Most of them
! are assigned to zeros, vertical mixing coefficients are assinged
! to their background values. These will remain unchenged if no
! vertical mixing scheme is applied. Because of the "first touch"
! default data distribution policy, this operation actually performs
! distribution of the shared arrays accross the cluster, unless
! another distribution policy is specified to override the default.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc
#include "param.h"
#include "das_param.h"
# include "das_param_lr.h"
#include "scalars.h"
#include "grid.h"
#include "das_ocean.h"
#include "das_ocean9.h"
#include "das_ocean_lr.h"
#include "das_ocean_smooth.h"
#include "das_covar.h"
#include "das_innov.h"

      real init    !!!!  0xFFFA5A5A ==> NaN

      parameter (init=0.)

!
# include "compute_extended_bounds.h"
!
!
!  Initialize 2-D primitive variables.
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_das(i,j)=init
          zeta_das9(i,j)=init
          zeta_adj(i,j)=init
          zeta_s(i,j)=init
          zeta_lr_s(i,j)=init
          zeta_sm(i,j)=init
          nsr_das(i,j)=init
          nsu_das(i,j)=init
          nsv_das(i,j)=init
          nzr_das(i,j)=init
          nzu_das(i,j)=init
          nzv_das(i,j)=init
        enddo
      enddo
#ifdef SOLVE3D
!
!  Initialize 3-D primitive variables.
!

      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            u_das(i,j,k)=init
            v_das(i,j,k)=init
            u_s(i,j,k)=init
            v_s(i,j,k)=init
            u_lr_s(i,j,k)=init
            v_lr_s(i,j,k)=init
            u_s_adj(i,j,k)=init
            v_s_adj(i,j,k)=init
            psi_das(i,j,k)=init
            chi_das(i,j,k)=init
            psi_s(i,j,k)=init
            chi_s(i,j,k)=init
            psi_adj(i,j,k)=init
            chi_adj(i,j,k)=init
            psi_s_adj(i,j,k)=init
            chi_s_adj(i,j,k)=init
            rho_s(i,j,k)=init
            rho_das9(i,j,k)=init
            pp_s(i,j,k)=init
            p_sm(i,j,k)=init
            p_sm_adj(i,j,k)=init
# ifdef MASKING
            rmask_das(i,j,k)=init
            pmask_das(i,j,k)=init
            umask_das(i,j,k)=init
            vmask_das(i,j,k)=init
#  endif
          enddo
        enddo
      enddo
      do itrc=1,NT
        do k=1,NDAS
          do j=JstrR,JendR
            do i=IstrR,IendR
              t_das(i,j,k,itrc)=init
              t_das9(i,j,k,itrc)=init
              t_s(i,j,k,itrc)=init
              t_lr_s(i,j,k,itrc)=init
              t_sm(i,j,k,itrc)=init
              t_adj(i,j,k,itrc)=init
              t_s_adj(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
#endif /* SOLVE3D */
!
! initialize those not innitialized in ROMS
!
#ifdef MASKING
      do j=JstrR,JendR
        do i=IstrR,IendR
          rmask(i,j)=init
          pmask(i,j)=init
          umask(i,j)=init
          vmask(i,j)=init
        enddo
      enddo
#endif
!
! initialize error variances and correlations
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          bz_das(i,j)=init
        enddo
      enddo

      do k=1,NDAS
        do j=JstrR,JendR
          do i=IstrR,IendR
            bpsi_das(i,j,k)=init
            bchi_das(i,j,k)=init
            bt_das(i,j,k,itemp)=init
            bt_das(i,j,k,isalt)=init
          enddo
        enddo
      enddo
!
      do j=JstrR,JendR
        do i=JstrR,JendR
          cpsiE_das(i,j)=init
          cchiE_das(i,j)=init
          ctE_das(i,j,itemp)=init
          ctE_das(i,j,isalt)=init
        enddo
      enddo
      do j=IstrR,IendR
        do i=IstrR,IendR
          cpsiX_das(i,j)=init
          cchiX_das(i,j)=init
          ctX_das(i,j,itemp)=init
          ctX_das(i,j,isalt)=init
        enddo
      enddo
      do j=1,NDAS
        do i=1,NDAS
          cpsiZ_das(i,j)=init
          cchiZ_das(i,j)=init
          ctZ_das(i,j,itemp)=init
          ctZ_das(i,j,isalt)=init
        enddo
      enddo
!
! initialize innovection vectors
!
#ifdef DAS_TMISST
      do j=JstrR,JendR
        do i=IstrR,IendR
          sst_tmi(i,j)=init
          tmi_mask(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_MCSST
      do j=JstrR,JendR
        do i=IstrR,IendR
# if !defined DAS_FDN_MCSST
          sst_mc(i,j)=init
          mc_mask(i,j)=init
# else
          sst_mc (i,j,1:nfndlev)=init
          mc_mask(i,j,1:nfndlev)=init
# endif
        enddo
      enddo
#endif
#ifdef DAS_GOES_SST
      do j=JstrR,JendR
        do i=IstrR,IendR
          sst_goes(i,j)=init
          goes_mask(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_SWOTSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh_swot(i,j)=init
          swot_mask(i,j)=init
          swot_oin(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_JASONSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh_js1(i,j)=init
          js1_mask(i,j)=init
          js1_oin(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_TPSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh_tp(i,j)=init
          tp_mask(i,j)=init
          tp_oin(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_SATSSSS
      do j=JstrR,JendR
        do i=IstrR,IendR
          sss_sats(i,j)=init
          sats_mask(i,j)=init
          sats_oin(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_NPSFLIGHT
      do j=JstrR,JendR
        do i=IstrR,IendR
          sst_fl(i,j)=init
          fl_mask(i,j)=init
        enddo
      enddo
#endif
#ifdef DAS_HFRADAR
      do j=JstrR,JendR
        do i=IstrR,IendR
          u_hf(i,j)=init
          v_hf(i,j)=init
          hf_umask(i,j)=init
          hf_vmask(i,j)=init
          u_hf6(i,j)=init
          v_hf6(i,j)=init
          hf_uierr6(i,j)=init
          hf_vierr6(i,j)=init
          hf_umask6(i,j)=init
          hf_vmask6(i,j)=init
        enddo
      enddo
#endif
!
! TS profles
!
      do k=1,ndas-1
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_obs(i,j,k)=init
            s_obs(i,j,k)=init
            obs_t_mask(i,j,k)=init
            obs_s_mask(i,j,k)=init
            err_tobs(i,j,k)=init
            err_sobs(i,j,k)=init
            dep_tobs(i,j,k)=init
            dep_sobs(i,j,k)=init
          enddo
        enddo
      enddo
!
! SIO
!
      do k=1,max_prf
        do j=JstrR,JendR                  
          do i=IstrR,IendR               
            t_sio(i,j,k)=init           
            s_sio(i,j,k)=init          
            sio_t_mask(i,j,k)=init    
            sio_s_mask(i,j,k)=init   
          enddo                     
        enddo                      
      enddo
!
! WHOI
!
      do k=1,max_prf
        do j=JstrR,JendR                     
          do i=IstrR,IendR                  
            t_whoi(i,j,k)=init              
            s_whoi(i,j,k)=init             
            whoi_t_mask(i,j,k)=init       
            whoi_s_mask(i,j,k)=init      
          enddo                        
        enddo                         
      enddo
!
! PTSUR
!
      do k=1,max_prf
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_ptsur(i,j,k)=init
            s_ptsur(i,j,k)=init
            ptsur_t_mask(i,j,k)=init
            ptsur_s_mask(i,j,k)=init
          enddo
        enddo
      enddo
!           
! MAPPED
!           
      do k=1,max_prf
        do j=JstrR,JendR     
          do i=IstrR,IendR     
            t_map(i,j,k)=init  
            s_map(i,j,k)=init    
            map_t_mask(i,j,k)=init
            map_s_mask(i,j,k)=init
          enddo                 
        enddo              
      enddo
!
! CalPoly AUV
!
#ifdef DAS_CALPOLYAUV
      do k=1,max_prf_cal
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_cal(i,j,k)=init
            s_cal(i,j,k)=init
            cal_mask(i,j,k)=init
          enddo
        enddo
      enddo
#endif
!
! Dorado AUV
!
#ifdef DAS_DORADOAUV
      do k=1,max_prf_dor
        do j=JstrR,JendR
          do i=IstrR,IendR
            t_dor(i,j,k)=init
            s_dor(i,j,k)=init
            dor_mask(i,j,k)=init
          enddo
        enddo
      enddo
#endif
#if defined DAS_JASONSSH || defined DAS_TPSSH\
  || defined DAS_TPSSH || defined DAS_ERS2SSH\
  || defined DAS_GFOSSH
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh_ref(i,j)=init
        enddo
      enddo
#endif
!
      return
      end
