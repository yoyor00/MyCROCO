! $Id: copy_to_devices.F 1615 2014-12-17 13:27:07Z rblod $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!

#include "cppdefs.h"
#if defined OPENACC
       subroutine copy_to_devices(tile)
       integer tile
#include "param.h"
#include "compute_tile_bounds.h"
       call copy_to_devices_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_to_devices_tile(Istr,Iend,Jstr,Jend)
       implicit none
#include "param.h"
       integer Istr,Iend,Jstr,Jend
#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "grid.h"
#include "coupling.h"
#ifdef MPI
#include "mpi_buffers.h"
#endif

!$acc enter data if(compute_on_device) copyin(z_w,z_r,Hz_bak,Hz,
!$acc& cs_r,sc_r,cs_w,sc_w,h,hinv,zeta,Zt_avg1
!$acc& ,pm_u,pn_u,u,v,pm,pn,pn_v,pm_v
!$acc& ,om_v,on_u,on_r,om_r,fomn,pnom_p,pmon_p,om_p,on_p
!$acc& ,on_v,om_u
!$acc& ,pnom_r,pmon_r, pnom_v, pmon_u
!$acc& ,Zob
!$acc& ,f
#ifdef MASKING
!$acc& ,rmask, pmask, umask, vmask, pmask, pmask2
#endif
# if defined CURVGRID && defined UV_ADV
!$acc& ,dmde,dndx
#endif
#if defined CURVGRID
!$acc& ,angler
#endif
#ifdef SPHERICAL
!$acc&, latr, lonr, latu, lonu, latv, lonv
#else
!$acc& ,xp,xr,yp,yr
#endif
#ifdef MPI
!$acc&, buf_snd4,buf_snd2,buf_rev4,buf_rev2
!$acc&, buf_snd1,buf_snd3,buf_rev1,buf_rev3
!$acc&, ibuf_sndN, ibuf_revN,
!$acc&     ibuf_sndS, ibuf_revS,
!$acc&     jbuf_sndW, jbuf_sndE,
!$acc&     jbuf_revW, jbuf_revE
#ifdef SOLVE3D
!$acc&, buf_snd4_3D, buf_snd2_3D,
!$acc&                  buf_rev4_3D, buf_rev2_3D,
!$acc&                  buf_snd1_3D, buf_snd3_3D,
!$acc&                  buf_rev1_3D, buf_rev3_3D
!$acc&, ibuf_sndN_3D, ibuf_revN_3D,
!$acc&     ibuf_sndS_3D, ibuf_revS_3D,
!$acc&     jbuf_sndW_3D, jbuf_sndE_3D,
!$acc&     jbuf_revW_3D, jbuf_revE_3D
#endif
#endif
!$acc& )

       return
       end

       subroutine copy_from_devices(tile)
       integer tile
#include "param.h"
#include "compute_tile_bounds.h"
       call copy_from_devices_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_from_devices_tile(Istr,Iend,Jstr,Jend)
       implicit none
#include "param.h"
       integer Istr,Iend,Jstr,Jend
#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "grid.h"
#include "coupling.h"

#include "compute_auxiliary_bounds.h"

!$acc exit data if(compute_on_device)
!$acc&  copyout(z_w(IstrR:IendR,JstrR:JendR,:)
!$acc& ,z_r(IstrR:IendR,JstrR:JendR,:)
!$acc& ,Hz_bak(IstrR:IendR,JstrR:JendR,:)
!$acc& ,Hz(IstrR:IendR,JstrR:JendR,:)
!$acc& ,hinv(IstrR:IendR,JstrR:JendR),Zt_avg1(IstrR:IendR,JstrR:JendR)
!$acc& )

       return
       end

       subroutine copy_to_devices_2(tile)
       integer tile
#include "param.h"
#include "compute_tile_bounds.h"
       call copy_to_devices_2_tile(Istr,Iend,Jstr,Jend)
       return
       end

       subroutine copy_to_devices_2_tile(Istr,Iend,Jstr,Jend)
       implicit none
#include "param.h"
       integer Istr,Iend,Jstr,Jend
#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "grid.h"
#include "coupling.h"
#include "private_scratch.h"
#include "mixing.h"
#include "forces.h"
#include "climat.h"
#include "work.h"
#include "ncscrum.h"
#include "averages.h"
#include "lmd_kpp.h"


!$acc enter data if(compute_on_device) copyin(A2d,A3d
!$acc& ,u,v,t
!$acc& ,huon,hvom,We
!$acc& ,weight
!$acc& ,Hz_bak,Zt_avg1,zeta
!$acc& ,rho,Hz,rho1,z_w,z_r
!$acc& ,rufrc,rvfrc,sustr,svstr,bustr,bvstr
#if defined UV_VIS2
!$acc& ,visc2_r,visc2_p
#endif
!$acc& ,cs_r,sc_r,cs_w,sc_w,h,hinv
!$acc& , DU_avg2, DV_avg2, DU_avg1, DV_avg1, ubar, vbar
!$acc& , Akv, Akt
#if defined TS_DIF2
!$acc& , diff2
#endif
!$acc& , btflx, stflx
!$acc& , rufrc_bak, rvfrc_bak
!$acc& , srflx
#if defined TS_DIF4
!$acc& ,diff4_sponge, diff2_sponge
#endif
#ifdef VAR_RHO_2D
!$acc& ,rhoA,rhoS
#endif
#if defined ANA_VMIX || defined BVF_MIXING || defined GLS_MIXING \
     || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP
!$acc& ,bvf
#endif
#if defined NONLIN_EOS && defined SPLIT_EOS
!$acc& ,qp1
#endif
#if defined OBC_EAST || defined OBC_WEST || defined OBC_NORTH || defined \
     OBC_SOUTH
#if defined T_FRC_BRY
#else
#if defined SOLVE3D && (defined M3CLIMATOLOGY || \
                        (defined AGRIF && !defined M3_FRC_BRY))
!$acc&, uclm, vclm,
#endif                        
#if defined M2CLIMATOLOGY || (defined AGRIF && !defined M2_FRC_BRY)
!$acc&, vbclm, ubclm                        
#endif  
#if defined TRACERS && (defined TCLIMATOLOGY || (defined AGRIF && !defined T_FRC_BRY))
!$acc&, tclm,
#endif
#if defined ZCLIMATOLOGY || defined AGRIF                        
!$acc&, ssh,                  
#endif
#endif
#endif
#if defined ZNUDGING
!$acc&, Znudgcof
#endif
#if defined M2NUDGING
!$acc&, M2nudgcof
#endif
#if defined M3NUDGING
!$acc&, M3nudgcof
#endif
#if defined TNUDGING
!$acc&, Tnudgcof
#endif
#if defined TS_MIX_ISO || defined TS_MIX_GEO
!$acc&, dRdx,dRde, idRz
#endif
#if defined LMD_SKPP || defined LMD_BKPP
!$acc&, kbl, kbbl, hbbl
#ifdef LMD_SKPP2005
!$acc&, hbls
#endif
#ifdef LMD_NONLOCAL
!$acc&, ghats
#endif
#endif
#ifdef DIF_COEF_3D
!$acc&, diff3d_u, diff3d_v
#endif
#ifdef TS_DIF4
!$acc&, diff4
#endif
#if defined TS_DIF2 || defined SPONGE_DIF2
!$acc& ,diff2, diff2_sponge
#endif
#ifdef LMD_MIXING
!$acc& ,ustar
#endif
# if defined SEDIMENT || defined LMD_MIXING
!$acc& ,B2d
#endif
#if defined UV_VIS2 || defined SPONGE_VIS2
!$acc& , visc2_sponge_p, visc2_sponge_r
#endif
#ifdef AVERAGES
!$acc& ,wrtavg
!$acc& ,srflx_avg, stflx_avg, sustr_avg, svstr_avg, wstr_avg
!$acc& ,bustr_avg, bvstr_avg, bostr_avg
!$acc& ,omega_avg, rho_avg, t_avg, u_avg, v_avg, w_avg
!$acc& ,zeta_avg, ubar_avg, vbar_avg
!$acc& ,workr
#ifdef AVERAGES_K
!$acc& ,Akv_avg, Akt_avg
#endif
#if defined LMD_SKPP || defined GLS_MIXING
!$acc& ,hbl_avg
#endif
# ifdef LMD_BKPP
!$acc& ,hbbl_avg
#endif
#ifdef DIF_COEF_3D
!$acc& ,diff3d_avg
#endif
#if defined LMD_SKPP || defined LMD_BKPP || defined GLS_MIXING
!$acc& ,Jwtype,lmd_mu1,lmd_mu2,lmd_r1
#endif
# if defined ANA_VMIX || defined BVF_MIXING \
  || defined LMD_MIXING || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
!$acc& , bvf_avg
#endif
#endif
#if defined TCLIMATOLOGY && !defined ANA_TCLIMA
!$acc& ,tclima
#endif
#if defined M2CLIMATOLOGY && !defined ANA_M2CLIMA ||\
   (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
!$acc& ,ubclima,vbclima,uclima,vclima
#endif
#if defined ZCLIMATOLOGY && !defined ANA_SSH
!$acc& ,sshg
#endif
#ifndef ANA_SMFLUX
!$acc& ,sustrg,svstrg
#endif
#if !defined ANA_STFLUX || !defined ANA_SSFLUX
!$acc& ,stflxg
#endif
#ifndef ANA_SRFLUX
!$acc& ,srflxg
#endif
# ifdef ANA_DIURNAL_SW
!$acc& ,tan_phi, sin_phi, cos_phi
#endif
# if defined QCORRECTION && (defined TEMPERATURE || defined SALINITY)
!$acc& ,dqdt, sst
#ifndef ANA_SST
!$acc&, sstg, dqdtg
#endif
#endif
#if defined SALINITY && defined SFLX_CORR
!$acc& , sss
#ifndef ANA_SSS
!$acc& ,sssg
#endif
#endif
!$acc& )

       return
       end
#else
       subroutine copy_to_devices_empty
       return
       end
#endif
