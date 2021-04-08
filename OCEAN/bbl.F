! $Id: bbl.F 1588 2014-08-04 16:26:01Z marchesiello $
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
#ifdef BBL
      subroutine bblm (tile)
!====================================================================!
!  Compute bottom stresses for combined waves & currents             !
!  using the parametric approximation by Soulsby 1997:               !
!  t_cw = t_c[1+1.2(t_w/(t_c+t_w))^3.2]                              !
!  in which                                                          !
!     t_cw = the combined wave-averaged stress (in current dir)      !
!     t_c  = stress due to currents if waves would be absent         ! 
!     t_w  = amplitude of stress due to waves without currents       !
!  and                                                               !
!     t_cw_max = SQRT([t_cw+t_w cos(phi_cw)]^2 + [t_w sin(phi_cw)]^2)!
!     in which                                                       !
!       t_cw_max ~ the maximum combined wave-averaged stress         !
!       phi_cw   = the angle between current and waves               !
!--------------------------------------------------------------------!
! References:
!
!  - Dyer 1986, Coastal & Estuarine Sediment Dynamics, Wiley, 342 pp.
!  - Harris & Wiberg 2001, Comp. & Geosci. 27, 675-690
!  - Li & Amos 2001, Comp. & Geosci. 27, 619-645
!  - Soulsby 1997, Dynamics of Marine Sands, Telford  Publ., 249 pp.
!  - Soulsby 1995, Bed shear-stresses due to combined waves and currents,
!    in: Stive et al: Advances in Coastal Morphodynamics, Wiley, 4.20-4.23
!  - Wiberg & Harris 1994, J. Geophys. Res. 99(C4), 775-789
!
! First Implementation: Meinte Blaas 2002
! Changes: P. Marchesiello & R. Benshila 2013-2014:
!  --> WKB wave model; bottom stress stability limit; 2D application
!          P. Marchesiello 2020: updates with BEDLOAD_VANDERA
!                                and new SANDBAR test case
!
!====================================================================!
      implicit none
      integer  tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
!
# include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
      call bblm_tile (Istr,Iend,Jstr,Jend,
     &                  A2d(1, 1,trd),A2d(1, 2,trd),A2d(1, 3,trd),
     &                  A2d(1, 4,trd),A2d(1, 5,trd),A2d(1, 6,trd),
     &                  A2d(1, 7,trd),A2d(1, 8,trd),A2d(1, 9,trd))
      return
      end

!********************************************************************!
      SUBROUTINE bblm_tile (Istr,Iend,Jstr,Jend,Ub,Zr,Ur,Vr,Umag,
     &                      tauc,tauw,taucw,taucwmax)
!********************************************************************!
!  output: bu/vstr = tauc = t_cw 
!     effective wave-averaged bottom stress
!     (after applying apparent roughness)  
!     bu/vstrw = tauw = t_cw_max
!     maximum skin-frictional bottom stress
!     rheight = ripple height (m)
!     rlength = ripple length (m)
!                                
!  input: Awave, Pwave, Dwave, h, Zob, d50, u,v
!         NB! Dwave is assumed to be in Cartesian convention (CROCO): 
!             the direction to where the vector points, measured 
!             counterclockwise from the positive x-axis of the
!             model grid (must account for grid angle)
!********************************************************************!
!
      implicit none
# include "param.h"
# include "bbl.h"
# include "forces.h"
# include "grid.h"
# ifdef SOLVE3D
#  include "ocean3d.h"
# else
#  include "ocean2d.h"
# endif
# include "scalars.h"
# include "sediment.h"
# if defined M3FAST && defined BSTRESS_FAST
#  include "nbq.h"
# endif
      logical SAND, SILT
      integer Iend, Istr, Jend, Jstr, i, ised, j
      real    Ab, anglec, anglew, cff, cff1, cff2,
     &        Fwave, Kbh, Kbh2, Kdh, phic, phicw, 
     &        d50, rhosed, rhow, smgd, osmgd,
     &        tau_cb, tau_ex, tau_bf, tau_up, tau_en,
     &        tau_w, tau_ws, tau_c, tau_cs, tau_cw, tau_cws,
     &        thetw, twopi, Ucur, Vcur, Ubc, visk, 
     &        znots, znot, znota, znot_bl, znot_rip,
     &        eps, zomin, zomax
      parameter (Ubc   = 1.e-2,  ! minimum orbital velocity for waves
     &           eps   = 1.e-10, ! small number to avoid zero division
     &           zomin = 1.e-5,  ! minimum roughness
     &           zomax = 3.e-2)  ! maximum roughness
# ifdef Z0_RIP
      real rhgt, rlen, rhmax, rhmin, rlmin, psi,
     &     rhbio, rlbio, rhbiomax, rhfac
      parameter (
     &       rhmin = 1.e-4,      ! minimum ripple height [m]  
     &       rlmin = rhmin/0.12, ! minimum ripple length [m]  
     &       rhbiomax = 6.e-3,   ! maximum biogenic ripple height [m]  
     &       rhfac = 0.0164)     ! biogenic ripple factor
# endif
      real        tauc(PRIVATE_2D_SCRATCH_ARRAY),
     &            tauw(PRIVATE_2D_SCRATCH_ARRAY),
     &           taucw(PRIVATE_2D_SCRATCH_ARRAY),
     &        taucwmax(PRIVATE_2D_SCRATCH_ARRAY),
     &              Ub(PRIVATE_2D_SCRATCH_ARRAY),
     &            Umag(PRIVATE_2D_SCRATCH_ARRAY),
     &              Ur(PRIVATE_2D_SCRATCH_ARRAY),
     &              Vr(PRIVATE_2D_SCRATCH_ARRAY),
     &              Zr(PRIVATE_2D_SCRATCH_ARRAY)
      real K1, K2, K3, K4, K5, K6
      parameter (K1=0.6666666666, K2=0.3555555555,
     &           K3=0.1608465608, K4=0.0632098765,
     &           K5=0.0217540484, K6=0.0065407983)
      real scf1, scf2, scf3, scf4, scf5
      parameter (scf1 = 0.5 * 1.39, scf2 = 0.52, 
     &           scf3 = 2.0 - scf2, scf4 = 1.2,
     &           scf5 = 3.2)
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
# ifdef WET_DRY
#  define SWITCH_WET *
# else
#  define SWITCH_WET !
# endif
!
      twopi=2.*pi
!
!-----------------------------------------------------------------------
!  Initalize stresses due to currents and waves.
!-----------------------------------------------------------------------
!
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          tauc(i,j)=0.       ! [m^2/s^2]
          tauw(i,j)=0.
          taucw(i,j)=0.
          taucwmax(i,j)=0.
        enddo
      enddo
!
!-----------------------------------------------------------------------
!  Set currents above bed.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend+1
        DO i=IstrU-1,Iend+1
# ifdef SOLVE3D
          Zr(i,j)=max(z_r(i,j,1)-z_w(i,j,0),Zob(i,j)+1.E-4)
#  if defined M3FAST && defined BSTRESS_FAST
          Ur(i,j)=2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
          Vr(i,j)=2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
#  else
          Ur(i,j)=u(i,j,1,nrhs)
          Vr(i,j)=v(i,j,1,nrhs)
#  endif
# else
          Zr(i,j)=0.5*(h(i,j)+zeta(i,j,kstp))
          Ur(i,j)=ubar(i,j,kstp)
          Vr(i,j)=vbar(i,j,kstp)
# endif
        enddo
      enddo
!
!=======================================================================
!  Compute bottom stresses.
!=======================================================================
!
! start main j,i loop
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
# ifdef Z0_RIP
          rhbio = 0.
          rlbio = 0.
          rlen = Lripple(i,j)
          rhgt = Hripple(i,j)
# endif
# ifdef SOLVE3D
          rhow = rho(i,j,1)+rho0
          visk = 1.3e-3/rhow         ! kinem. viscosity
# else
          visk = 1.3e-3/rho0
# endif
!
!---------------------------------------------------------------------
! Compute bed wave orbital velocity (m/s), excursion amplitude (m)
! from wind-induced waves, and angle between current and waves.
!---------------------------------------------------------------------
!
! Use Dean & Dalrymple 1991 6th-degree polynomial to approximate 
! wavenumber on shoaling water.
!
          Fwave=twopi/Pwave(i,j)
          Kdh=(h(i,j)+z_w(i,j,N))*Fwave*Fwave/g
          Kbh2=Kdh*Kdh+Kdh/(1.+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+
     &                         Kdh*(K5+K6*Kdh))))))
          Kbh = SQRT(Kbh2)
!
!  Compute bed wave orbital velocity and excursion amplitude.
!
          Ab=Awave(i,j)/SINH(Kbh) + eps
          Ub(i,j)=Fwave*Ab SWITCH rmask(i,j)
     &                     SWITCH_WET rmask_wet(i,j)
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur=0.5*(Ur(i,j)+Ur(i+1,j))
          Vcur=0.5*(Vr(i,j)+Vr(i,j+1))
          Umag(i,j)=SQRT(Ucur*Ucur+Vcur*Vcur)
     &                             SWITCH rmask(i,j)
     &                             SWITCH_WET rmask_wet(i,j)
     &                             + eps
!
!  Compute angle between currents and waves (radians)
!
          if (Ucur .ne. 0.) then
            phic=ATAN2(Vcur,Ucur)
          else
            phic=0.5*pi*SIGN(1.,Vcur)
          endif
          phicw=Dwave(i,j)-phic  ! Dwave in CROCO grid convention
!
!---------------------------------------------------------------------
! Sediment parameters:
! Establish local median grain size for all calculations in this 
! subroutine. Since most parameterizations have been derived ignoring 
! multiple grain sizes, we apply this single d50 also in the case of 
! mixed beds. 
!---------------------------------------------------------------------
!
# if defined ANA_BSEDIM || defined SEDIMENT
          d50 = Ssize(i,j)             ! [m]
          tau_cb = taucb(i,j)          ! [m^2/s^2]
          rhosed = Sdens(i,j)/rhow     ! [ ] relative density
# else
          d50 = 0.16e-3       ! default values: fine sand
          tau_cb = 0.17/rhow
          rhosed = 2650./rhow
# endif /* ANA_BSEDIM || SEDIMENT */
!
          smgd=(rhosed-1.)*g*d50
          osmgd=1./smgd
!
!---------------------------------------------------------------------
! Determine grain roughness from sediment size
!---------------------------------------------------------------------
!
          znots = d50/12.0
!
!---------------------------------------------------------------------
! Determine critical stresses for ripple formation
!---------------------------------------------------------------------
!
# ifdef Z0_RIP
!
! Transition to sheet flow (Li & Amos, 2001)
!
          tau_up = 0.03*smgd*d50**(-0.376)                ! [m^2/s^2]
!
! Break off (Grant & Madsen,1982)
!
          tau_bf = 0.78*(d50/visk)**0.6*smgd**0.3*tau_cb  ! [m^2/s^2]
# endif
!
!---------------------------------------------------------------------
! If significant waves (Ub > Ubc= 0.01 m/s):
! Wave-current interaction case according to Soulsby 1995.  
! Otherwise: tauw = tauc for sediment purposes
!---------------------------------------------------------------------
!
          if(Ub(i,j).gt.Ubc) then   !<======
!
!---------------------------------------------------------------------
! Determine skin stresses for pure waves and combined flow (m^2/s^2)
! using Soulsby approximation of the wave friction factor fw
!---------------------------------------------------------------------

! Pure wave skin stress
!  tau_ws=0.5*fw*Ub**2 with fw=2*scf1*(znots/Ab)**scf2
            tau_ws = scf1*((znots*Fwave)**scf2)*Ub(i,j)**scf3

! current-only skin stress
            cff1 = vonKar/LOG(Zr(i,j)/znots)
            tau_cs = cff1*cff1*Umag(i,j)*Umag(i,j)     ! [m^2/s^2] 

!  combined wave-current skin stress 
            tau_cws = tau_cs*(1.+scf4*((tau_ws/(tau_ws+tau_cs))**scf5))

!  Maximum of combined wave-current skin stress
#  ifdef SED_TAUCWMAX
            tau_cws = SQRT((tau_cws + tau_ws*COS(phicw))**2
     &                              +(tau_ws*SIN(phicw))**2)
#  endif

!  Store for sediment model 
            taucwmax(i,j) = tau_cws      ! [m^2/s^2]
            tauw(i,j)     = tau_ws

# ifdef Z0_BL
!---------------------------------------------------------------------
! Compute apparent bedload roughness for ripple predictor (and sediment 
! resuspension if needed). At high transport stages, friction depends on 
! thickness of bedload layer (Li & Amos 2001).
!---------------------------------------------------------------------
!
            tau_ex=max(tau_cws-tau_cb,0.)
            znot_bl=17.4*d50*(tau_ex*osmgd)**0.75
            znot_bl= znots + znot_bl
!
!---------------------------------------------------------------------
!  Compute bedload stress for sediment purpose [m^2/s^2], 
!  using grain + bedload roughness.
!---------------------------------------------------------------------
!
            cff1=vonKar/LOG(Zr(i,j)/znot_bl)
            cff2 = MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
            tau_c = cff2*Umag(i,j)*Umag(i,j)    
            tau_w = scf1*((znot_bl*Fwave)**scf2)*(Ub(i,j)**scf3)
            tau_cw = tau_c*(1.+scf4*((tau_w/(tau_w+tau_c))**scf5))
#  ifdef SED_TAUCWMAX
            tau_cw = SQRT((tau_cw+tau_w*COS(phicw))**2 +
     &                           (tau_w*SIN(phicw))**2)
#  endif
!
            taucwmax(i,j) = tau_cw
            tauw(i,j)     = tau_w
            tauc(i,j)     = tau_c

#  ifdef Z0_RIP        
!---------------------------------------------------------------------
! Determine bedform roughness ripple height (m) and ripple length (m) 
! for sandy beds. 
!---------------------------------------------------------------------
!
#   ifdef Z0_RIP_VANDERA
!  O'Donoghue et al. (2006) - Van der A et al. (2013)
!
            psi=osmgd*(1.27*Ub(i,j))**2          ! Max mobility number
            CALL rip_dim(psi, d50, rhgt, rlen)   ! Ripple dimensions
            rhgt=max(rhmin,rhgt*Ab)              ! height (m)
            rlen=max(rlmin,rlen*Ab)              ! length (m)
#   else
!  Li & Amos (2001) - Blaas et al. (2007)
!
            if (d50 .ge. 0.063e-3) then
              rhmax=0.25*rlen
              rhgt=max(min(rhmax,rhgt),rhmin)
              tau_en=tau_cws*max(1.,(rlen/(rlen-pi*rhgt))**2) 
              if (tau_cws.lt.tau_cb .and. tau_en.ge.tau_cb) then 
                rhgt = (19.6*SQRT(tau_cws/tau_cb)+20.9)*d50
                rlen = rhgt/0.12       ! local transport
              elseif (tau_cws.ge.tau_cb .and. tau_cw.lt.tau_bf) then
                if (tau_ws.ge.1.56*tau_cs) then
                  rhgt = (27.14*SQRT(tau_cw/tau_cb)+16.36)*d50
                  rlen = rhgt/0.15     ! wave bedload regime
                else
                  rhgt = (22.15*SQRT(tau_cw/tau_cb)+6.38)*d50
                  rlen = rhgt/0.12     ! curr bedload regime
                endif
              elseif (tau_cw.ge.tau_bf .and. tau_cw.lt.tau_up) then
                rlen = 535.*d50        ! break-off regime
                rhgt = 0.15*rlen*(SQRT(tau_up)-SQRT(tau_cw))/
     &                           (SQRT(tau_up)-SQRT(tau_bf )) 
              elseif (tau_cw.ge.tau_up) then 
                rlen = 0.              ! sheet flow (plane bed)
                rhgt = 0.
              else 
                rhgt=Hripple(i,j)      ! tau_en < tau_cb: no transport, 
                rlen=Lripple(i,j)      ! pre-existing conditions
              endif !tau_cws
            endif !d50
#   endif
#  endif /* Z0_RIP */
   
#  ifdef Z0_BIO
!---------------------------------------------------------------------
! Determine (biogenic) bedform roughness, ripple height (m)
! and ripple length (m) for silty beds, using Harris & Wiberg 2001.
!---------------------------------------------------------------------
!
! Use 10 cm default biogenic ripple length rlbio (Wheatcroft 1994) 
!
            if (d50.lt.0.063e-3) then
              rlbio = 0.1
              thetw = tau_cws*osmgd
              rhbio = (thetw**(-1.67))*rlbio*rhfac
              rhgt = min(rhbio,rhbiomax)
              rlen = rlbio
            endif
#  endif /* Z0_BIO */

#  if defined Z0_RIP || defined Z0_BIO
!---------------------------------------------------------------------
! Total roughness znota consists of grain-scale roughness +
! bedload thickness + sandy ripple-enhanced roughness + silty 
! biogenic ripple roughness.
!---------------------------------------------------------------------
!
            znot_rip = 0.267*rhgt*rhgt/(max(rlen,rlmin)) ! Nielsen (1992)
            znota = znot_bl + znot_rip

#  endif /* Z0_RIP || Z0_BIO */
!
!---------------------------------------------------------------------
! Limit total roughness between zomin or Zob and zomax
!---------------------------------------------------------------------
!
            znota = max(zomin,min(zomax,max(znota,Zob(i,j))))
!
# endif /* Z0_BL */
!
!---------------------------------------------------------------------
! Compute bottom stress components based on total roughness [m/s]^2
! and store for use in computing flow drag and eddy diffusivity/viscosity
!---------------------------------------------------------------------
!
            cff1  = vonKar/LOG(Zr(i,j)/znota)
            cff2 = MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
            tau_c = cff2*Umag(i,j)*Umag(i,j)  
            tau_w = scf1*((znota*Fwave)**scf2)*(Ub(i,j)**scf3)
            tau_cw= tau_c*(1.+scf4*((tau_w/(tau_w+tau_c))**scf5))
!
            taucw(i,j) = tau_cw
!          
          else ! (Ub < Ubc)  <================
!
!---------------------------------------------------------------------
! If current-only: tauw = tauc(skin) for use in sediment.F (ifdef BBL)
! tauc for current still depending on roughness due to current ripples 
! (ifdef Z0_RIP)
!---------------------------------------------------------------------
!
# ifdef Z0_RIP    
            if(tau_cs.gt.tau_up) then
              rhgt=0.
              rlen=0.
            else if(tau_cs.lt.tau_cb) then
              rhgt=Hripple(i,j) 
              rlen=Lripple(i,j) 
            else
              rlen=1000.*d50            ! Yalin (1964)
              rhgt=0.117*rlen**1.19     ! Allen (1970)
            endif
            znota = znots + 0.267*rhgt*rhgt/(max(rlen,rlmin))
# else
            znota = znots
# endif
            znota = max(zomin,min(zomax,max(znota,Zob(i,j))))
            cff1=vonKar/LOG(Zr(i,j)/znota)
            cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
            tau_c=cff2*Umag(i,j)*Umag(i,j)
!
            tauc(i,j)     = tau_c
            taucw(i,j)    = tau_c
            taucwmax(i,j) = tau_c
            tauw(i,j)     = 0.

          endif ! (Ub < Ubc)  <================
!
!---------------------------------------------------------------------
!  Load variables for output purposes.
!---------------------------------------------------------------------
!
          Abed(i,j)    = Ab            SWITCH rmask(i,j)
         Zbnot(i,j)    = znots         SWITCH rmask(i,j)
         Zbapp(i,j)    = znota         SWITCH rmask(i,j)
# ifdef Z0_RIP
          Hripple(i,j) = rhgt          SWITCH rmask(i,j)
          Lripple(i,j) = rlen          SWITCH rmask(i,j)
# endif
       enddo
      enddo  ! end main loop
!
!---------------------------------------------------------------------
!  Compute kinematic bottom stress components (m2/s2) for flow
!  due to combined current and wind-induced waves.
!---------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1=0.5*(taucw(i-1,j)+taucw(i,j))   ! interp at u point
          anglec=Ur(i,j)/(0.5*(Umag(i-1,j)+Umag(i,j)))
          bustr(i,j)=cff1*anglec SWITCH umask(i,j)
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
          cff1=0.5*(taucw(i,j-1)+taucw(i,j))    ! interp at v point
          anglec=Vr(i,j)/(0.5*(Umag(i,j-1)+Umag(i,j)))
          bvstr(i,j)=cff1*anglec SWITCH vmask(i,j)
        enddo
      enddo
!
!  Compute skin stress components at rho points
!
      do j=Jstr,Jend
        do i=Istr,Iend
          Ucur=0.5*(Ur(i,j)+Ur(i+1,j))
          Vcur=0.5*(Vr(i,j)+Vr(i,j+1))
!
          anglec=Ucur/Umag(i,j)
          anglew=cos(Dwave(i,j))
          bustrc(i,j)=tauc(i,j)*anglec
          bustrw(i,j)=tauw(i,j)*anglew
          bustrcwmax(i,j)=taucwmax(i,j)*anglew
!          Ubot(i,j)=Ub(i,j)*anglew
!          Ur(i,j)=Ucur
!
          anglec=Vcur/Umag(i,j)
          anglew=sin(Dwave(i,j))
          bvstrc(i,j)=tauc(i,j)*anglec
          bvstrw(i,j)=tauw(i,j)*anglew
          bvstrcwmax(i,j)=taucwmax(i,j)*anglew
!          Vbot(i,j)=Ub(i,j)*anglew
!          Vr(i,j)=Vcur
!
        enddo
      enddo

# ifdef LIMIT_BSTRESS
!---------------------------------------------------------------------
!  From J. Warner's code:
!  Set limiting factor for bottom stress. The bottom stress is adjusted
!  to not change the direction of momentum.  It only should slow down
!  to zero.  The value of 0.75 is arbitrary limitation assigment.
!---------------------------------------------------------------------
!
      cff=0.75/dt
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1=cff*0.5*(Zr(i-1,j)+Zr(i,j))
          bustr(i,j)=SIGN(1.0, bustr(i,j))* 
     &               MIN(ABS(bustr(i,j)),
     &                   ABS(Ur(i,j))*cff1)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1=cff*0.5*(Zr(i,j-1)+Zr(i,j))
          bvstr(i,j)=SIGN(1.0, bvstr(i,j))*
     &               MIN(ABS(bvstr(i,j)),
     &                   ABS(Vr(i,j))*cff1)
        enddo
      enddo
# endif
!
!---------------------------------------------------------------------
! Set boundary conditions
!---------------------------------------------------------------------
!
#  if defined M3FAST && defined BSTRESS_FAST
      if (LAST_FAST_STEP) then
#  endif

# ifndef EW_PERIODIC
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          bustr(Iend+1,j)=bustr(Iend,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Iend+1,j)=bvstr(Iend,j)
        END DO
        !
        DO j=Jstr,Jend
          bustrc(Iend+1,j)=bustrc(Iend,j)
          bustrw(Iend+1,j)=bustrw(Iend,j)
          bustrcwmax(Iend+1,j)=bustrcwmax(Iend,j)
          bvstrc(Iend+1,j)=bvstrc(Iend,j)
          bvstrw(Iend+1,j)=bvstrw(Iend,j)
          bvstrcwmax(Iend+1,j)=bvstrcwmax(Iend,j)
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          bustr(IstrU-1,j)=bustr(IstrU,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Istr-1,j)=bvstr(Istr,j)
        END DO
        DO j=Jstr,Jend
          bustrc(Istr-1,j)=bustrc(Istr,j)
          bustrw(Istr-1,j)=bustrw(Istr,j)
          bustrcwmax(Istr-1,j)=bustrcwmax(Istr,j)
          bvstrc(Istr-1,j)=bvstrc(Istr,j)
          bvstrw(Istr-1,j)=bvstrw(Istr,j)
          bvstrcwmax(Istr-1,j)=bvstrcwmax(Istr,j)
        END DO        
      END IF
# endif
# ifndef NS_PERIODIC
      IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jend+1) =bustr(i,Jend)
        END DO
        DO i=Istr,Iend
          bvstr(i,Jend+1) =bvstr(i,Jend)
        END DO
        DO i=Istr,Iend
          bustrc(i,Jend+1) =bustrc(i,Jend)
          bustrw(i,Jend+1) =bustrw(i,Jend)
          bustrcwmax(i,Jend+1) =bustrcwmax(i,Jend)
          bvstrc(i,Jend+1) =bvstrc(i,Jend)
          bvstrw(i,Jend+1) =bvstrw(i,Jend)
          bvstrcwmax(i,Jend+1) =bvstrcwmax(i,Jend)
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          bustr(i,Jstr-1)=bustr(i,Jstr)
        END DO
        DO i=Istr,Iend
          bvstr(i,JstrV-1)=bvstr(i,JstrV)
        END DO
        DO i=Istr,Iend
          bustrc(i,Jstr-1)=bustrc(i,Jstr)
          bustrw(i,Jstr-1)=bustrw(i,Jstr)
          bustrcwmax(i,Jstr-1)=bustrcwmax(i,Jstr)
          bvstrc(i,Jstr-1)=bvstrc(i,Jstr)
          bvstrw(i,Jstr-1)=bvstrw(i,Jstr)
          bvstrcwmax(i,Jstr-1)=bvstrcwmax(i,Jstr)
        END DO
      END IF
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,Jstr-1)=0.5*(bustr(Istr+1,Jstr-1)+bustr(Istr,Jstr))
        bvstr(Istr-1,Jstr)=0.5*(bvstr(Istr,Jstr)+bvstr(Istr-1,Jstr+1))
        bustrc(Istr-1,Jstr-1)=0.5*(bustrc(Istr,Jstr-1)+bustrc(Istr-1,Jstr))
        bvstrc(Istr-1,Jstr-1)=0.5*(bvstrc(Istr,Jstr-1)+bvstrc(Istr-1,Jstr))
        bustrw(Istr-1,Jstr-1)=0.5*(bustrw(Istr,Jstr-1)+bustrw(Istr-1,Jstr))
        bvstrw(Istr-1,Jstr-1)=0.5*(bvstrw(Istr,Jstr-1)+bvstrw(Istr-1,Jstr))
        bustrcwmax(Istr-1,Jstr-1)=0.5*(bustrcwmax(Istr,Jstr-1)
     &   +bustrcwmax(Istr-1,Jstr))
        bvstrcwmax(Istr-1,Jstr-1)=0.5*(bvstrcwmax(Istr,Jstr-1)
     &   +bvstrcwmax(Istr-1,Jstr))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(Iend+1,Jstr-1)=0.5*(bustr(Iend+1,Jstr)+bustr(Iend,Jstr-1))
        bvstr(Iend+1,Jstr)=0.5*(bvstr(Iend+1,Jstr+1)+bvstr(Iend,Jstr))
        bustrc(Iend+1,Jstr-1)=0.5*(bustrc(Iend+1,Jstr)+bustrc(Iend,Jstr-1))
        bvstrc(Iend+1,Jstr-1)=0.5*(bvstrc(Iend+1,Jstr)+bvstrc(Iend,Jstr-1))
        bustrw(Iend+1,Jstr-1)=0.5*(bustrw(Iend+1,Jstr)+bustrw(Iend,Jstr-1))
        bvstrw(Iend+1,Jstr-1)=0.5*(bvstrw(Iend+1,Jstr)+bvstrw(Iend,Jstr-1))
        bustrcwmax(Iend+1,Jstr-1)=0.5*(bustrcwmax(Iend+1,Jstr)
     &   +bustrcwmax(Iend,Jstr-1))
        bvstrcwmax(Iend+1,Jstr-1)=0.5*(bvstrcwmax(Iend+1,Jstr)
     &   +bvstrcwmax(Iend,Jstr-1))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        bustr(Istr,Jend+1)=0.5*(bustr(Istr,Jend)+bustr(Istr+1,Jend+1))
        bvstr(Istr-1,Jend+1)=0.5*(bvstr(Istr-1,Jend)+bvstr(Istr,Jend+1))
        bustrc(Istr-1,Jend+1)=0.5*(bustrc(Istr-1,Jend)+bustrc(Istr,Jend+1))
        bvstrc(Istr-1,Jend+1)=0.5*(bvstrc(Istr-1,Jend)+bvstrc(Istr,Jend+1))
        bustrw(Istr-1,Jend+1)=0.5*(bustrw(Istr-1,Jend)+bustrw(Istr,Jend+1))
        bvstrw(Istr-1,Jend+1)=0.5*(bvstrw(Istr-1,Jend)+bvstrw(Istr,Jend+1))
        bustrcwmax(Istr-1,Jend+1)=0.5*(bustrcwmax(Istr-1,Jend)
     &   +bustrcwmax(Istr,Jend+1))
        bvstrcwmax(Istr-1,Jend+1)=0.5*(bvstrcwmax(Istr-1,Jend)
     &   +bvstrcwmax(Istr,Jend+1))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        bustr(Iend+1,Jend+1)=0.5*(bustr(Iend+1,Jend)+bustr(Iend,Jend+1))
        bvstr(Iend+1,Jend+1)=0.5*(bvstr(Iend+1,Jend)+bvstr(Iend,Jend+1))
        bustrc(Iend+1,Jend+1)=0.5*(bustrc(Iend+1,Jend)+bustrc(Iend,Jend+1))
        bvstrc(Iend+1,Jend+1)=0.5*(bvstrc(Iend+1,Jend)+bvstrc(Iend,Jend+1))
        bustrw(Iend+1,Jend+1)=0.5*(bustrw(Iend+1,Jend)+bustrw(Iend,Jend+1))
        bvstrw(Iend+1,Jend+1)=0.5*(bvstrw(Iend+1,Jend)+bvstrw(Iend,Jend+1))
        bustrcwmax(Iend+1,Jend+1)=0.5*(bustrcwmax(Iend+1,Jend)
     &   +bustrcwmax(Iend,Jend+1))
        bvstrcwmax(Iend+1,Jend+1)=0.5*(bvstrcwmax(Iend+1,Jend)
     &   +bvstrcwmax(Iend,Jend+1))
      END IF
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,Abed(START_2D_ARRAY))
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,Hripple(START_2D_ARRAY))
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,Lripple(START_2D_ARRAY))
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,Zbnot(START_2D_ARRAY))
!      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,Zbapp(START_2D_ARRAY))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,bustr(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,bvstr(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bustrc(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bvstrc(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bustrw(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bvstrw(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bustrcwmax(START_2D_ARRAY))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,bvstrcwmax(START_2D_ARRAY))
# endif

#  if defined M3FAST && defined BSTRESS_FAST
      endif
#  endif

      return
      end
!
!======================================================================
!
      SUBROUTINE rip_dim(psi, d50, rhgt, rlen)
! 
! This subroutine returns rhgt, rlen: ripple height and length
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
!        
      implicit none 
#  include "param.h"
#  include "scalars.h"
      real psi, d50, rhgt, rlen
      real m_eta, m_lam, n_eta, n_lam
      real d50_mm, eps
      parameter(eps=1.e-14)
!     
      d50_mm=0.001*d50
      IF(d50_mm.lt.0.22) THEN
        m_eta=0.55
        m_lam=0.73
      ELSEIF(d50_mm.ge.0.22.and.d50_mm.lt.0.30) THEN
        m_eta=0.55+(0.45*(d50_mm-0.22)/(0.30-0.22))
        m_lam=0.73+(0.27*(d50_mm-0.22)/(0.30-0.22))
      ELSE
        m_eta=1.
        m_lam=1.
      ENDIF
! 
! Smooth transition between ripple regime and bed sheet flow regime 
!
      IF(psi.le.190.) THEN
        n_eta=1.
      ELSEIF(psi.gt.190..and.psi.lt.240.) THEN
        n_eta=0.5*(1.+cos(pi*(psi-190.)/50.))
      ELSEIF(psi.ge.240.) THEN
        n_eta=0.
      ENDIF
      n_lam=n_eta
!
      rhgt=MAX(0., m_eta*n_eta*(0.275-0.022*psi**0.42))
      rlen=MAX(eps,m_lam*n_lam*(1.970-0.440*psi**0.21))
!
      RETURN
      END
!
!======================================================================
!
#else
      subroutine bblm_empty
      end
#endif

