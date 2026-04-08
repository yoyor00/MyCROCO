!------------------------------------------------------------------------------
 MODULE sed_MUSTANG_CROCO
!------------------------------------------------------------------------------

#include "cppdefs.h"
#if defined MUSTANG 

!&E=========================================================================
!&E                   ***  MODULE  sed_MUSTANG_CROCO  ***
!&E
!&E ** Purpose : concerns subroutines related to sediment dynamics link to 
!&E              hydrodynamic model to be used in CROCO system
!&E 
!&E ** Description :
!&E     subroutine sed_MUSTANG_settlveloc ! settling velocity in the water 
!&E                                         column
!&E     subroutine sed_gradvit            ! calcul gradient de vitesse, u*
!&E     subroutine sed_skinstress         ! computes the skin stress
!&E     subroutine sed_bottom_slope
!&E     subroutine sed_exchange_w2s ! MPI treatment of slip deposit fluxes
!&E     subroutine sed_exchange_flxbedload ! MPI treatment of bedload fluxes
!&E     subroutine sed_exchange_maskbedload ! MPI exchange of mask for 
!&E                                           bedload
!&E     subroutine sed_exchange_corflu ! MPI treatment of corflu fluxes
!&E     subroutine sed_obc_corflu ! corflu fluxes at boundaries
!&E     subroutine sed_meshedges_corflu ! corflu fluxes interpolation at 
!&E                                       mesh edges
!&E
!&E==========================================================================

    !! * Modules used
    USE comMUSTANG
    USE comsubstance
    USE module_substance
# if defined key_MUSTANG_flocmod
    USE flocmod, ONLY : f_ws
#endif
    IMPLICIT NONE

    !! * Accessibility 
    PUBLIC sed_skinstress
    PUBLIC sed_gradvit
    PUBLIC sed_MUSTANG_settlveloc
#ifdef key_MUSTANG_bedload
    PUBLIC sed_bottom_slope
#if defined MPI 
    PUBLIC sed_exchange_flxbedload
    PUBLIC sed_exchange_maskbedload
#endif
#endif
#if defined MPI  && defined key_MUSTANG_slipdeposit
    PUBLIC sed_exchange_w2s
#endif

    PUBLIC sed_obc_corflu
    PUBLIC sed_meshedges_corflu
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
    PUBLIC sed_exchange_corflu
#endif


PRIVATE

CONTAINS
 
!!=============================================================================
SUBROUTINE sed_MUSTANG_settlveloc(ifirst, ilast, jfirst, jlast,   &
                                    WATER_CONCENTRATION) 
!&E----------------------------------------------------------------------------
!&E                 ***  ROUTINE sed_MUSTANG_settlveloc  ***
!&E
!&E ** Purpose : settling velocity computation
!&E
!&E ** Description : use arguments and common variable 
!&E  arguments IN : 
!&E         WATER_CONCENTRATION = t : WATER_CONCENTRATION 
!&E  arguments OUT:
!&E         ws_part : settling velocities for CROCO
!&E         ws3_bottom_MUSTANG: settling velocities in  bottom cell
!&E
!&E  need to be know 
!&E         g : gravity
!&E         kmax=N : number of layer in water
!&E          
!&E  need to be know by code treated substance 
!&E  (if not ==> coupler_MUSTANG.F90)
!&E         imud1, nvpc, nvp, nv_adv, isand1, isand2
!&E         f_ws(iv) (if key_MUSTANG_flocmod)
!&E         ws_free_opt, ws_free_para, ws_free_min, ws_free_max,
!&E         ws_hind_opt, ws_hind_para   
!&E     
!&E  use module MUSTANG variables  :
!&E         ros(iv)
!&E         ws_sand(iv)
!&E        
!&E ** Called by :  MUSTANG_update
!&E
!&E----------------------------------------------------------------------------

!! * Arguments
INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,N,3,NT), INTENT(IN) :: WATER_CONCENTRATION  
   
!! * Local declarations
INTEGER                    :: iv, k, ivpc, i, j
REAL(KIND=rsh)             :: cmes, phi, phiv, De, WSfree, Hind
REAL(KIND=rsh), PARAMETER  :: nuw = 0.00000102_rsh

!! * Executable part
DO j = jfirst, jlast
DO i = ifirst, ilast
           
    IF (htot(i, j) > h0fond) THEN
        DO k = 1, N
            cmes = 0.0_rsh
            DO ivpc = imud1, nvpc
                cmes = cmes + WATER_CONCENTRATION(i, j, k, nstp, itemp + ntrc_salt + ivpc)
            ENDDO
            cmes = MAX(0.0_rsh, cmes)

            ! first calculating sand settling velocity
            DO iv = isand1, isand2
                ws_part(i, j, k, itemp + ntrc_salt + iv) = ws_sand(iv)
            ENDDO
            
            ! next mud settling velocity 
#ifdef key_MUSTANG_flocmod   
            ws_part(i, j, k, itemp + ntrc_salt + imud1 : itemp + ntrc_salt + imud2 ) = f_ws(1:nv_mud)  
#else
            DO iv = imud1, nvp
                ! Free settling velocity - flocculation
                IF(ws_free_opt(iv) == 0) THEN ! constant settling velocity
                    WSfree = ws_free_min(iv)
                ELSEIF (ws_free_opt(iv) == 1) THEN ! Van Leussen 1994
                    WSfree = ws_free_para(1, iv) * cmes**ws_free_para(2, iv) &
                        * (1._rsh + ws_free_para(3, iv) * gradvit(k, i, j)) / &
                        (1._rsh + ws_free_para(4, iv) * gradvit(k, i, j)**2._rsh)
                ELSEIF (ws_free_opt(iv) == 2) THEN ! Winterwerp 1999
                    De = ws_free_para(1, iv) &
                        + ws_free_para(2, iv) &
                        * cmes / (ws_free_para(3, iv) * sqrt(gradvit(k, i, j)))
                    IF (De .GT. sqrt(nuw / gradvit(k, i, j))) THEN 
                        De = sqrt(nuw / gradvit(k, i, j)) 
                        ! in case of large C/low G limit floc size to kolmogorov microscale
                    ENDIF
                    WSfree = (ros(iv) - rho0) * g / (18._rsh * rho0 * nuw)  &
                        * ws_free_para(1, iv)**(3._rsh - ws_free_para(4, iv))  &
                        * De**(ws_free_para(4, iv) - 1._rsh)
                ELSEIF (ws_free_opt(iv) == 3) THEN ! Wolanski et al., 1989
                    WSfree = ws_free_para(1, iv) * cmes**ws_free_para(2, iv)
                ENDIF

                ! Hindered settling
                ! if ws_hind_opt.EQ.0 : no hindered settling... Hind = 1
                Hind = 1._rsh
                IF (ws_hind_opt(iv) == 1) THEN ! Scott, 1984
                    phi = MIN(1.0_rsh, cmes / ws_hind_para(1, iv))
                    Hind = (1._rsh - phi)**ws_hind_para(2, iv)
                ELSEIF (ws_hind_opt(iv) == 2) THEN ! Winterwerp, 2002 
                    ! WARNING : ros(iv) must be the same for all MUDS variables
                    phi = cmes / ros(iv)
                    IF (ws_free_opt(iv) == 2) THEN
                        phiv = phi * (De / ws_free_para(1, iv))**(3._rsh - ws_free_para(4, iv))
                    ELSE
                        phiv = cmes / ws_hind_para(1, iv)
                    ENDIF
                    Hind = (1._rsh - phiv)**ws_hind_para(2, iv) * &
                        (1._rsh - phi) / (1._rsh + 2.5_rsh * phiv)
                ELSEIF (ws_hind_opt(iv) == 3) THEN ! wolanski et al., 1989
                    IF (ws_free_opt(iv) == 3) THEN
                        Hind = 1._rsh / (cmes**2._rsh + ws_hind_para(1, iv)**2._rsh)**ws_hind_para(2, iv)
                    ENDIF  ! ws_hind_opt(iv) == 3 only if ws_free_opt(iv) == 3
                ENDIF

                ! limiting ws with min/max values...
                ! necessary if ndt_part not updated during the simulation, 
                ! ndt_part calculated in t3dmix_tridiagonal_settling from max ws
                ws_part(i, j, k, itemp + ntrc_salt + iv) = max(ws_free_min(iv), &
                    min(ws_free_max(iv), WSfree * Hind))
            ENDDO

#endif  /* key_MUSTANG_flocmod */

            DO iv = nvpc+1, nvp
                IF(irkm_var_assoc(iv) < imud1 .AND. irkm_var_assoc(iv) > 0) THEN    
                    ! sorbed substances on sands
                    ws_part(i, j, k, itemp + ntrc_salt + iv) = &
                        ws_part(i, j, k, itemp + ntrc_salt + irkm_var_assoc(iv))
                ENDIF
            ENDDO
            
            DO iv = nvp+1, nv_adv
                ws_part(i, j, k, itemp + ntrc_salt + iv) = 0.0_rsh
            ENDDO

        ENDDO ! do loop on k
        ws3_bottom_MUSTANG(1:nvp, i, j) = ws_part(i, j, 1, itemp+ntrc_salt+1:itemp+ntrc_salt+nvp)

    ELSE  ! htot(i,j) <= h0fond
        ws_part(i, j, :, :) = 0.0_rsh
        ws3_bottom_MUSTANG(:, i, j) = 0.0_rsh
    ENDIF

ENDDO
ENDDO

END SUBROUTINE sed_MUSTANG_settlveloc     

!!=============================================================================
SUBROUTINE sed_gradvit(ifirst, ilast, jfirst, jlast)
!&E--------------------------------------------------------------------------
!&E                 ***  ROUTINE sed_gradvit  ***
!&E
!&E ** Purpose : calculation of the turbulence energy G  
!&E
!&E ** Description : G= sqrt(turbulence dissipation/viscosity)
!&E                 to be programmed using hydrodynamic knowledge
!&E           using htot, rho0, sig, epn, nz ..
!&E
!&E     output : gradvit (in comMUSTANG)
!&E
!&E ** Called by :  MUSTANG_update
!&E
!&E--------------------------------------------------------------------------
!! * Modules used
#if defined key_MUSTANG_flocmod && defined SED_TOY_FLOC_0D
    USE flocmod, ONLY : flocmod_comp_g
#endif
#  include "mixing.h"
#  include "ocean3d.h"

!! * Arguments 
INTEGER, INTENT(IN)       :: ifirst, ilast, jfirst, jlast

!! * Local declarations
INTEGER        :: i, j, k
REAL(KIND=rsh) :: dist_surf_on_bottom, nuw, diss

nuw = 1.0e-6

DO j = jfirst, jlast
DO i = ifirst, ilast
    IF(htot(i, j) .GT. h0fond)  THEN
    DO k = 1, N
#if defined key_MUSTANG_flocmod && defined SED_TOY_FLOC_0D
        call flocmod_comp_g(gradvit(k, i, j), time-time_start)
#else
#if defined GLS_MIXING
        !
        ! Dissipation from turbulence clossure
        !
        if (k.eq.1) then
            diss = Eps_gls(i,j,k)
        elseif (k.eq.N) then
            diss = Eps_gls(i,j,k-1)
        else
            diss = 0.5*(Eps_gls(i,j,k-1)+Eps_gls(i,j,k))
        endif
        gradvit(k, i, j) = sqrt(diss/nuw)
#else
    ! gradvit : G=sqrt( turbulence dissipation rate/ vertical viscosity coefficient)
    ! if  turbulence dissipation rate has not been already evaluated: 
    ! use empirical formula from   Nezu and Nakawaga (1993)
    ! turbulence dissipation_rate = ustarbot**3 /Karman/Htot * (distance from surface/distance from bottom)
        dist_surf_on_bottom = ((z_w(i, j, N) - z_r(i, j, k)) / (z_r(i, j, k) - z_w(i, j, 0)))
        gradvit(k, i, j) = sqrt(ustarbot(i, j)**3._rsh / 0.4_rsh / htot(i, j) / &
                        (nuw + epsilon_MUSTANG) * dist_surf_on_bottom) 
#endif
#endif
    END DO
    ENDIF
ENDDO
ENDDO

END SUBROUTINE sed_gradvit


SUBROUTINE sed_skinstress(ifirst, ilast, jfirst, jlast)
! Compute bottom skin-friction stress (tauskin, ustarbot)
! from combined wave and current interaction (Soulsby 1995).
! Called by : MUSTANG_update
!
! l_tauskin_center : compute tauskin_c at rho point directly
! l_tauskin_ubar   : use depth-averaged velocity instead of u(k=1)
! l_tauskin_upwind : use upwind interpolation for tauskin_x/y
# ifdef WAVE_OFFLINE
    USE module_substance, ONLY : Uwave, Dwave, Pwave
# endif
#  ifdef BBL
#  include "bbl.h"
#  endif

INTEGER, INTENT(IN) :: ifirst, ilast, jfirst, jlast

REAL(KIND=rsh), DIMENSION(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1) :: Zr, Zref
REAL(KIND=rsh), DIMENSION(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1) :: ustar2_u, ustar2_v
REAL(KIND=rsh), DIMENSION(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1) :: ux, vx
REAL(KIND=rsh) :: z0_factor

!--- Bottom cell height, clipped to z0sed + margin -----------------------
Zr = compute_zr(z_r(ifirst-1:ilast+1, jfirst-1:jlast+1, 1),  &
                z_w(ifirst-1:ilast+1, jfirst-1:jlast+1, 0),  &
                z0sed(ifirst-1:ilast+1, jfirst-1:jlast+1))

#ifdef BBL
!--- BBL branch : tauskin from wave-current bottom stress ----------------
!    d50 is constant (160 microns) in bustrw/bvstrw — see bbl.F
CALL compute_tauskin_bbl(ifirst, ilast, jfirst, jlast, &
    u(:, :, 1, nnew), v(:, :, 1, nnew), &
    bustrw, bvstrw, rho0, &
#  if defined WET_DRY && defined MASKING
    rmask_wet, &
#  endif
    tauskin, tauskin_c, &
    tauskin_x, tauskin_y)
#else
!--- Select velocity field and log-argument numerator --------------------
!    ubar / full depth if l_tauskin_ubar
!    u(k=1) / bottom cell height otherwise
IF (l_tauskin_ubar) THEN
  ux        = ubar(ifirst-1:ilast+1, jfirst-1:jlast+1, nnew)
  vx        = vbar(ifirst-1:ilast+1, jfirst-1:jlast+1, nnew)
  Zref      = z_w(ifirst-1:ilast+1, jfirst-1:jlast+1, N) &
            - z_w(ifirst-1:ilast+1, jfirst-1:jlast+1, 0) ! full water column
  z0_factor = 2.718_rsh ! Euler number — log-layer correction
ELSE
  ux        = u(ifirst-1:ilast+1, jfirst-1:jlast+1, 1, nnew)
  vx        = v(ifirst-1:ilast+1, jfirst-1:jlast+1, 1, nnew)
  Zref      = Zr ! bottom cell height
  z0_factor = 1.0_rsh
ENDIF

!--- Log-law skin friction ------------------------------------------------
IF (l_tauskin_center) THEN
  CALL compute_tauskin_c_center(ifirst, ilast, jfirst, jlast,             &
       ux(ifirst-1:ilast+1, jfirst-1:jlast+1),                            &
       vx(ifirst-1:ilast+1, jfirst-1:jlast+1), Zref,                      &
       z0sed  (ifirst-1:ilast+1, jfirst-1:jlast+1),                       &
       z0_factor,                                                         &
       rho    (ifirst-1:ilast+1, jfirst-1:jlast+1, :),                    &
       rho0,                                                              &
       tauskin_c(ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       tauskin_x(ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       tauskin_y(ifirst-1:ilast+1, jfirst-1:jlast+1))
ELSE
  CALL compute_ustar2_uv(ifirst, ilast, jfirst, jlast,                    &
       ux(ifirst-1:ilast+1, jfirst-1:jlast+1),                            &
       vx(ifirst-1:ilast+1, jfirst-1:jlast+1), Zref,                      &
       z0sed    (ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       z0_factor,                                                         &
       ustar2_u(ifirst-1:ilast+1, jfirst-1:jlast+1),                      &
       ustar2_v(ifirst-1:ilast+1, jfirst-1:jlast+1),                      &
       raphbx   (ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       raphby   (ifirst-1:ilast+1, jfirst-1:jlast+1))

  CALL interpolate_tauskin_c(ifirst, ilast, jfirst, jlast,                &
       l_tauskin_upwind,                                                  &
       ustar2_u(ifirst-1:ilast+1, jfirst-1:jlast+1),                      &
       ustar2_v(ifirst-1:ilast+1, jfirst-1:jlast+1),                      &
       raphbx   (ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       raphby   (ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       rho      (ifirst-1:ilast+1, jfirst-1:jlast+1, :),                  &
       rho0,                                                              &
       tauskin_x(ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       tauskin_y(ifirst-1:ilast+1, jfirst-1:jlast+1),                     &
       tauskin_c(ifirst-1:ilast+1, jfirst-1:jlast+1))
ENDIF

!--- Wave + current combination (Soulsby 1995) ----------------------------
#  ifdef WAVE_OFFLINE
CALL combine_wave_current(ifirst, ilast, jfirst, jlast,               &
     tauskin_c(ifirst:ilast, jfirst:jlast),                           &
     Uwave(ifirst:ilast, jfirst:jlast),                               &
     Dwave(ifirst:ilast, jfirst:jlast),                               &
     Pwave(ifirst:ilast, jfirst:jlast),                               &
     ubar(ifirst:ilast, jfirst:jlast, nnew),                          &
     vbar(ifirst:ilast, jfirst:jlast, nnew),                          &
     z0sed(ifirst:ilast, jfirst:jlast),                               &
     rho(ifirst:ilast, jfirst:jlast, 1),                              &
     rho0, l_fricwave, fws2, tauskin(ifirst:ilast, jfirst:jlast),     &
     tauskin_w(ifirst:ilast, jfirst:jlast))
#  else
tauskin(ifirst:ilast, jfirst:jlast) = tauskin_c(ifirst:ilast, jfirst:jlast)
#  endif

#  if defined WET_DRY && defined MASKING
tauskin(ifirst:ilast, jfirst:jlast) = tauskin(ifirst:ilast, jfirst:jlast)* &
                                      rmask_wet(ifirst:ilast, jfirst:jlast)
#  endif

#endif  /* BBL */

!--- Bottom friction velocity ---------------------------------------------
WHERE (htot(ifirst:ilast, jfirst:jlast) > h0fond)
ustarbot(ifirst:ilast, jfirst:jlast) = &
  SQRT(MAX(0.0_rsh, tauskin(ifirst:ilast, jfirst:jlast))/rho0)
ELSEWHERE
ustarbot(ifirst:ilast, jfirst:jlast) = 0.0_rsh
END WHERE

END SUBROUTINE sed_skinstress


ELEMENTAL PURE FUNCTION compute_zr(z_r1, z_w0, z0sed) RESULT(Zr)
! Bottom cell height clipped to z0sed + z0_margin to prevent log(0).
REAL(KIND=rsh), INTENT(IN) :: z_r1, z_w0, z0sed
REAL(KIND=rsh)             :: Zr
REAL(KIND=rsh), PARAMETER  :: z0_margin = 1.0e-4_rsh
Zr = MAX(z_r1 - z_w0 - z0sed, z0_margin) + z0sed
END FUNCTION compute_zr

ELEMENTAL PURE FUNCTION log_law_drag(numerator, log_arg) RESULT(ustar2)
! Log-law drag : ustar2 = 0.16 * log(Z/z0)^-2 * numerator
! Caller provides the numerator (speed*u or speed^2) and Z/z0.
! The 0.16 coefficient is the squared von Karman drag coefficient.
REAL(KIND=rsh), INTENT(IN) :: numerator   ! speed*u_component or speed^2
REAL(KIND=rsh), INTENT(IN) :: log_arg     ! Z/z0, must be > 0
REAL(KIND=rsh)             :: ustar2
ustar2 = 0.16_rsh*LOG(log_arg)**(-2)*numerator
END FUNCTION log_law_drag

ELEMENTAL PURE FUNCTION upwind_interp(f_i, f_ip1) RESULT(res)
! Upwind interpolation using SIGN — branchless equivalent of:
!   f_i > 0 and f_ip1 > 0 : take f_i
!   f_i < 0 and f_ip1 < 0 : take f_ip1
!   otherwise              : take average
REAL(KIND=rsh), INTENT(IN) :: f_i, f_ip1
REAL(KIND=rsh)             :: res, cff1, cff2, cff3, cff4
cff1 = 0.5_rsh*(1.0_rsh + SIGN(1.0_rsh, f_ip1))
cff2 = 0.5_rsh*(1.0_rsh - SIGN(1.0_rsh, f_ip1))
cff3 = 0.5_rsh*(1.0_rsh + SIGN(1.0_rsh, f_i))
cff4 = 0.5_rsh*(1.0_rsh - SIGN(1.0_rsh, f_i))
res = cff3*(cff1*f_i + cff2*0.5_rsh*(f_i + f_ip1)) &
+ cff4*(cff2*f_ip1 + cff1*0.5_rsh*(f_i + f_ip1))
END FUNCTION upwind_interp

ELEMENTAL PURE FUNCTION weighted_avg(f_i, f_ip1, w_i, w_ip1) RESULT(res)
! Weighted average of f at two adjacent points using weights w.
! w = |u| / (|u| + eps) — larger weight on the upwind side.
! eps guard prevents division by zero when both velocities vanish.
REAL(KIND=rsh), INTENT(IN) :: f_i, f_ip1, w_i, w_ip1
REAL(KIND=rsh)             :: res
res = (f_i*w_i + f_ip1*w_ip1)/(w_i + w_ip1 + epsilon_MUSTANG)
END FUNCTION weighted_avg

#ifdef WAVE_OFFLINE
ELEMENTAL PURE FUNCTION compute_fws2(Uwave, Pwave, z0sed, &
             l_fricwave, fws2_default) &
RESULT(fws2ij)
! Wave friction factor (Soulsby 1997).
! Adjusted for significant waves when l_fricwave is true.
! Falls back to fws2_default otherwise.
REAL(KIND=rsh), INTENT(IN) :: Uwave, Pwave, z0sed, fws2_default
LOGICAL, INTENT(IN) :: l_fricwave
REAL(KIND=rsh)             :: fws2ij
REAL(KIND=rsh), PARAMETER  :: uwave_min = 0.001_rsh
REAL(KIND=rsh), PARAMETER  :: pwave_min = 0.001_rsh
fws2ij = fws2_default
IF (l_fricwave .AND. &
Uwave*Pwave > 0.0_rsh .AND. &
Pwave > pwave_min .AND. &
Uwave > uwave_min) THEN
fws2ij = 0.5_rsh*1.39_rsh &
*(Uwave*Pwave/REAL(2.0_rlg*pi*z0sed, rsh))**(-0.52_rsh)
END IF
END FUNCTION compute_fws2

ELEMENTAL PURE FUNCTION compute_tauskin_wc(tauskin_c, tauskin_w, &
                   ubar, vbar, Dwave) &
RESULT(tauskin)
! Soulsby (1995) scalar wave+current stress at one rho point.
! Current stress amplified by wave presence, then combined
! with wave stress accounting for relative angle.
! Falls back to simple addition if either stress vanishes.
REAL(KIND=rsh), INTENT(IN) :: tauskin_c, tauskin_w
REAL(KIND=rsh), INTENT(IN) :: ubar, vbar, Dwave
REAL(KIND=rsh)             :: tauskin
REAL(KIND=rsh)             :: speedbar, tauskin_cw, alpha, beta
speedbar = SQRT(ubar**2 + vbar**2)
IF (tauskin_c > 0.0_rsh .AND. &
tauskin_w > 0.0_rsh .AND. &
speedbar > 0.0_rsh) THEN
tauskin_cw = tauskin_c*(1.0_rsh + 1.2_rsh &
      *(tauskin_w/(tauskin_w + tauskin_c))**3.2_rsh)
alpha = ACOS(vbar/speedbar)   ! current direction from north
beta = Dwave                    ! wave direction from north
tauskin = SQRT((tauskin_cw + tauskin_w*ABS(COS(alpha - beta)))**2 &
+ (tauskin_w*ABS(SIN(alpha - beta)))**2)
ELSE
tauskin = tauskin_w + tauskin_c
END IF
END FUNCTION compute_tauskin_wc
#endif

PURE SUBROUTINE compute_ustar2_uv(ifirst, ilast, jfirst, jlast, &
          ux, vx, Zref, z0sed, z0_factor, &
          ustar2_u, ustar2_v, &
          raphbx, raphby)
! Log-law drag at u and v staggered points via log_law_drag kernel.
!
! ustar2_u(i,j) = 0.16 * log(Zref_u / (z0_u * z0_factor))^-2 * speed_u * ux(i,j)
! ustar2_v(i,j) = 0.16 * log(Zref_v / (z0_v * z0_factor))^-2 * speed_v * vx(i,j)
!
! ux/vx      : velocity field (ubar or u(k=1)) selected by caller
! Zref       : log-argument numerator (full depth or bottom cell height)
! z0_factor  : 2.718 if ubar (log-layer correction), 1.0 if u(k=1)
! raphbx/by  : upwind direction weights |u|/(|u|+eps) at u,v points

INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), INTENT(IN)  :: ux(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: vx(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: Zref(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: z0sed(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: z0_factor
REAL(KIND=rsh), INTENT(OUT) :: ustar2_u(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: ustar2_v(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: raphbx(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: raphby(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)

INTEGER        :: i, j
REAL(KIND=rsh) :: speed_u, speed_v, z0_u, z0_v, Zref_u, Zref_v

!--- u points : (ifirst:ilast+1, jfirst:jlast) ---------------------------
DO j = jfirst, jlast
DO i = ifirst, ilast + 1
raphbx(i, j) = ABS(ux(i, j))/(ABS(ux(i, j)) + epsilon_MUSTANG)
z0_u = 0.5_rsh*(z0sed(i - 1, j) + z0sed(i, j))
Zref_u = 0.5_rsh*(Zref(i - 1, j) + Zref(i, j))
speed_u = SQRT(0.0625_rsh*(vx(i, j) + vx(i, j + 1) + &
         vx(i - 1, j) + vx(i - 1, j + 1))**2 &
+ ux(i, j)**2)*ux(i, j)
ustar2_u(i, j) = log_law_drag(speed_u, &
            Zref_u/MAX(z0_u*z0_factor, &
                       epsilon_MUSTANG))
END DO
END DO

!--- v points : (ifirst:ilast, jfirst:jlast+1) ---------------------------
DO j = jfirst, jlast + 1
DO i = ifirst, ilast
raphby(i, j) = ABS(vx(i, j))/(ABS(vx(i, j)) + epsilon_MUSTANG)
z0_v = 0.5_rsh*(z0sed(i, j - 1) + z0sed(i, j))
Zref_v = 0.5_rsh*(Zref(i, j - 1) + Zref(i, j))
speed_v = SQRT(0.0625_rsh*(ux(i, j) + ux(i + 1, j) + &
         ux(i, j - 1) + ux(i + 1, j - 1))**2 &
+ vx(i, j)**2)*vx(i, j)
ustar2_v(i, j) = log_law_drag(speed_v, &
            Zref_v/MAX(z0_v*z0_factor, &
                       epsilon_MUSTANG))
END DO
END DO

END SUBROUTINE compute_ustar2_uv

PURE SUBROUTINE compute_tauskin_c_center(ifirst, ilast, jfirst, jlast, &
                 ux, vx, Zref, z0sed, z0_factor, &
                 rho, rho0_val, &
                 tauskin_c, tauskin_x, tauskin_y)
! Log-law skin friction at rho point via log_law_drag kernel.
! tauskin_c   = rho * 0.16 * log(Zref/(z0*z0_factor))^-2 * speed^2
! tauskin_x/y = tauskin_c * velocity direction unit vector
!
! ux/vx      : velocity field (ubar or u(k=1)) selected by caller
! Zref       : log-argument numerator (full depth or bottom cell height)
! z0_factor  : 2.718 if ubar (log-layer correction), 1.0 if u(k=1)

INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), INTENT(IN)  :: ux(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: vx(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: Zref(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: z0sed(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: z0_factor
REAL(KIND=rsh), INTENT(IN)  :: rho(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1, N)
REAL(KIND=rsh), INTENT(IN)  :: rho0_val
REAL(KIND=rsh), INTENT(OUT) :: tauskin_c(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_x(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_y(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)

INTEGER        :: i, j
REAL(KIND=rsh) :: urho, vrho, speed, rhotot

DO j = jfirst, jlast
DO i = ifirst, ilast
urho = 0.5_rsh*(ux(i, j) + ux(i + 1, j))
vrho = 0.5_rsh*(vx(i, j) + vx(i, j + 1))
speed = SQRT(urho**2 + vrho**2)
rhotot = rho(i, j, 1) + rho0_val
tauskin_c(i, j) = log_law_drag(speed**2, &
             Zref(i, j)/MAX(z0sed(i, j)*z0_factor, &
                            epsilon_MUSTANG)) &
*rhotot
tauskin_x(i, j) = urho/(speed + epsilon_MUSTANG)*tauskin_c(i, j)
tauskin_y(i, j) = vrho/(speed + epsilon_MUSTANG)*tauskin_c(i, j)
END DO
END DO

END SUBROUTINE compute_tauskin_c_center

PURE SUBROUTINE interpolate_tauskin_c(ifirst, ilast, jfirst, jlast, &
              l_upwind, &
              ustar2_u, ustar2_v, &
              raphbx, raphby, &
              rho, rho0_val, &
              tauskin_x, tauskin_y, tauskin_c)
! Interpolate ustar2 from u,v staggered points to rho point.
! l_upwind=T : upwind scheme (SIGN-based, branchless)
! l_upwind=F : weighted centred average using raphbx/raphby

INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
LOGICAL, INTENT(IN)  :: l_upwind
REAL(KIND=rsh), INTENT(IN)  :: ustar2_u(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: ustar2_v(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: raphbx(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: raphby(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: rho(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1, N)
REAL(KIND=rsh), INTENT(IN)  :: rho0_val
REAL(KIND=rsh), INTENT(OUT) :: tauskin_x(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_y(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_c(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)

INTEGER        :: i, j
REAL(KIND=rsh) :: rhotot, tx, ty

IF (l_upwind) THEN
DO j = jfirst, jlast
DO i = ifirst, ilast
tx = upwind_interp(ustar2_u(i, j), ustar2_u(i + 1, j))
ty = upwind_interp(ustar2_v(i, j), ustar2_v(i, j + 1))
rhotot = rho(i, j, 1) + rho0_val
tauskin_c(i, j) = SQRT(tx**2 + ty**2)*rhotot
tauskin_x(i, j) = tx*rhotot
tauskin_y(i, j) = ty*rhotot
END DO
END DO
ELSE
DO j = jfirst, jlast
DO i = ifirst, ilast
tx = weighted_avg(ustar2_u(i, j), ustar2_u(i + 1, j), &
raphbx(i, j), raphbx(i + 1, j))
ty = weighted_avg(ustar2_v(i, j), ustar2_v(i, j + 1), &
raphby(i, j), raphby(i, j + 1))
rhotot = rho(i, j, 1) + rho0_val
tauskin_c(i, j) = SQRT(tx**2 + ty**2)*rhotot
tauskin_x(i, j) = tx*rhotot
tauskin_y(i, j) = ty*rhotot
END DO
END DO
END IF

END SUBROUTINE interpolate_tauskin_c

PURE SUBROUTINE compute_tauskin_bbl(ifirst, ilast, jfirst, jlast, &
            ux, vx, &
            bustrw, bvstrw, rho0_val, &
#if defined WET_DRY && defined MASKING
            rmask_wet, &
#endif
            tauskin, tauskin_c, &
            tauskin_x, tauskin_y)
! BBL variant : tauskin from precomputed wave-current stresses.
! tauskin   = rho0 * sqrt(bustrw^2 + bvstrw^2)
! tauskin_c = tauskin (BBL: no wave/current separation)
! tauskin_x/y = tauskin_c * velocity direction unit vector (key_MUSTANG_bedload only)

INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), INTENT(IN)  :: ux(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: vx(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: bustrw(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: bvstrw(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(IN)  :: rho0_val
#if defined WET_DRY && defined MASKING
REAL(KIND=rsh), INTENT(IN)  :: rmask_wet(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
#endif
REAL(KIND=rsh), INTENT(OUT) :: tauskin(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_c(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_x(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_y(ifirst - 1:ilast + 1, jfirst - 1:jlast + 1)

#ifdef key_MUSTANG_bedload
INTEGER        :: i, j
REAL(KIND=rsh) :: urho, vrho, speed
#endif

tauskin = SQRT(bustrw**2 + bvstrw**2)*rho0_val
#if defined WET_DRY && defined MASKING
tauskin = tauskin*rmask_wet
#endif
tauskin_c = tauskin

#ifdef key_MUSTANG_bedload
DO j = jfirst, jlast
DO i = ifirst, ilast
urho = 0.5_rsh*(ux(i, j) + ux(i + 1, j))
vrho = 0.5_rsh*(vx(i, j) + vx(i, j + 1))
speed = SQRT(urho**2 + vrho**2)
tauskin_x(i, j) = urho/(speed + epsilon_MUSTANG)*tauskin_c(i, j)
tauskin_y(i, j) = vrho/(speed + epsilon_MUSTANG)*tauskin_c(i, j)
END DO
END DO
#else
tauskin_x = 0.0_rsh
tauskin_y = 0.0_rsh
#endif

END SUBROUTINE compute_tauskin_bbl

#ifdef WAVE_OFFLINE
PURE SUBROUTINE combine_wave_current(ifirst, ilast, jfirst, jlast, &
             tauskin_c, Uwave, Dwave, Pwave, &
             ubar, vbar, z0sed, &
             rho, rho0_val, l_fricwave, &
             fws2_default, &
             tauskin, tauskin_w)
! Soulsby (1995) wave + current combined bed stress.
! tauskin direction follows current direction.
! If either current or wave stress vanishes, stresses are simply added.

INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
REAL(KIND=rsh), INTENT(IN)  :: tauskin_c(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: Uwave(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: Dwave(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: Pwave(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: ubar(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: vbar(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: z0sed(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: rho(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(IN)  :: rho0_val, fws2_default
LOGICAL, INTENT(IN)  :: l_fricwave
REAL(KIND=rsh), INTENT(OUT) :: tauskin(ifirst:ilast, jfirst:jlast)
REAL(KIND=rsh), INTENT(OUT) :: tauskin_w(ifirst:ilast, jfirst:jlast)

! tauskin_w : vectorised over full domain
tauskin_w = (rho + rho0_val)                                              &
            * compute_fws2(Uwave, Pwave, z0sed, l_fricwave, fws2_default) &
            * Uwave**2

! tauskin : wave + current combination at each rho point
tauskin = compute_tauskin_wc(tauskin_c, tauskin_w, ubar, vbar, Dwave)
END SUBROUTINE combine_wave_current
#endif


!!==============================================================================
#ifdef key_MUSTANG_bedload
  SUBROUTINE sed_bottom_slope(ifirst, ilast, jfirst, jlast, bathy)
   !&E--------------------------------------------------------------------------                         
   !&E                 ***  ROUTINE sed_bottom_slope  ***
   !&E
   !&E ** Purpose : evaluation of bottom slope
   !&E
   !&E ** Description : depend on host model grid, compute slope_dhdx and 
   !&E    slope_dhdy from bathy of neigbour cells if they are not masked
   !&E
   !&E ** Called by :  MUSTANG_update, MUSTANG_morpho 
   !&E      om_r = cell dx 
   !&E      on_r = cell dy
   !&E
   !&E--------------------------------------------------------------------------
   !! * Arguments
   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT(IN)  :: bathy  ! bathymetry (m)

   !! * Local declarations
   INTEGER :: i, j

      DO j = jfirst, jlast
        DO i = ifirst, ilast
          IF (bathy(i+1, j).LE. -valmanq .OR. bathy(i-1, j).LE. -valmanq) then
             slope_dhdx(i, j) = 0.0_rsh
          ELSE
             slope_dhdx(i, j) = -1.0_rsh*(-bathy(i+1, j)+bathy(i-1, j)) / (2.0_rsh * om_r(i, j))
          ENDIF
          IF (bathy(i, j+1).LE. -valmanq .OR. bathy(i, j-1).LE. -valmanq) then
             slope_dhdy(i, j) = 0.0_rsh
          ELSE
             slope_dhdy(i, j) = -1.0_rsh*(-bathy(i, j+1)+bathy(i, j-1)) / (2.0_rsh * on_r(i, j))
          ENDIF
        ENDDO
      ENDDO

  END SUBROUTINE sed_bottom_slope
!!=============================================================================
#if defined MPI
  SUBROUTINE sed_exchange_flxbedload(ifirst, ilast, jfirst, jlast)
    !&E--------------------------------------------------------------------------                         
    !&E                 ***  ROUTINE sed_exchange_flxbedload  ***
    !&E
    !&E ** Purpose : treatment of bedload fluxes
    !&E
    !&E ** Description : MPI exchange between processors
    !&E
    !&E ** Called by :  MUSTANG_update
    !&E
    !&E--------------------------------------------------------------------------
    ! * Arguments
    INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast

    !! * Local declarations
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv

    do iv=ibedload1,ibedload2
        workexch(:,:) = flx_bx(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
             &          workexch(START_2D_ARRAY))
        flx_bx(iv,:,:) = workexch(:,:)
   
        workexch(:,:) = flx_by(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
             &          workexch(START_2D_ARRAY))
        flx_by(iv,:,:) = workexch(:,:)
    enddo

  END SUBROUTINE sed_exchange_flxbedload
#endif /* MPI */
#endif /* key_MUSTANG_bedload */

!!=============================================================================

#if defined MPI && defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
    SUBROUTINE sed_exchange_maskbedload(ifirst, ilast, jfirst, jlast)
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_exchange_maskbedload ***
    !&E
    !&E ** Purpose : exchange MPI mask for bedload
    !&E
    !&E ** Description : exchange MPI mask for bedload
    !&E
    !&E ** Called by : MUSTANG_update
    !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    workexch(:,:) = sedimask_h0plusxe(:,:)
    call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
              &           workexch(START_2D_ARRAY))
    sedimask_h0plusxe(:,:) = workexch(:,:)

    END SUBROUTINE sed_exchange_maskbedload
#endif /* defined MPI && defined key_MUSTANG_V2 && defined key_MUSTANG_bedload */

!!=============================================================================

#if defined MPI && defined key_MUSTANG_slipdeposit
    SUBROUTINE sed_exchange_w2s(ifirst, ilast, jfirst, jlast)
    !&E-------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_exchange_w2s ***
    !&E
    !&E ** Purpose : MPI exchange of slip deposit flux between processors
    !&E
    !&E ** Description : MPI exchange between processors
    !&E      used only if slopefac .NE. 0 (slip deposit if steep slope)
    !&E
    !&E ** Called by : MUSTANG_update
    !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv
      
    do iv = isand2+1, nvp
        workexch(:,:) = flx_w2s_corim1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corim1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corip1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corip1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corjm1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corjm1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corjp1(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corjp1(iv,:,:) = workexch(:,:)

        workexch(:,:) = flx_w2s_corin(iv,:,:)
        call exchange_r2d_tile (ifirst,ilast,jfirst,jlast,  &
                &          workexch(START_2D_ARRAY))
        flx_w2s_corin(iv,:,:) = workexch(:,:)
    enddo
  
    END SUBROUTINE sed_exchange_w2s
#endif /* defined MPI && defined key_MUSTANG_slipdeposit */
!!=============================================================================

#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
  SUBROUTINE sed_exchange_corflu(ifirst, ilast, jfirst, jlast, type)
   !&E-------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_corflu ***
   !&E
   !&E ** Purpose : treatment of horizontal flow corrections for the transport 
   !&E              of sand in suspension
   !&E
   !&E ** Description : periodic borders and MPI exchange between processors
   !&E
   !&E ** Called by : MUSTANG_update
   !&E-------------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast, type

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch
    INTEGER :: iv
    
    if (type .eq. 0) then  ! corflux and corfluy are still at rho point
        do iv = isand1, isand2
            workexch(:, :) = corflux(iv, :, :)
            call exchange_r2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corflux(iv, :, :) = workexch(:, :)
    
            workexch(:, :) = corfluy(iv, :, :)
            call exchange_r2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corfluy(iv, :, :) = workexch(:, :)
        enddo
    else  ! corflux and corfluy are at u,v point
        do iv = isand1, isand2
            workexch(:, :) = corflux(iv, :, :)
            call exchange_u2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corflux(iv, :, :) = workexch(:, :)
    
            workexch(:, :) = corfluy(iv, :, :)
            call exchange_v2d_tile (ifirst, ilast, jfirst, jlast,  &
                &          workexch(START_2D_ARRAY))
            corfluy(iv, :, :) = workexch(:, :)
        enddo
    endif

   END SUBROUTINE sed_exchange_corflu
#endif /* defined EW_PERIODIC || defined NS_PERIODIC || defined MPI */
!!=============================================================================
   
   SUBROUTINE sed_obc_corflu(istr, iend, jstr, jend)
 
    !&E------------------------------------------------------------------------
    !&E                 ***  ROUTINE sed_obc_corflu ***
    !&E
    !&E ** Purpose : treatment of horizontal flow corrections for the transport 
    !&E              of sand in suspension
    !&E
    !&E ** Description : extrapolation at borders 
    !&E   Use SOUTHERN_EDGE, WESTERN_EDGE, NORTHERN_EDGE, EASTERN_EDGE defined 
    !&E   in OCEAN/set_global_definitions.h using 
    !&E   variables istr, iend, jstr, jend
    !&E
    !&E ** Called by : MUSTANG_update
    !&E--------------------------------------------------------------------------
 
    !! * Arguments
    INTEGER,INTENT(IN) :: istr, iend, jstr, jend
 
    !! * Local declarations
    INTEGER :: i, j, ivp

    do ivp = isand1, isand2
      if (WESTERN_EDGE) then
          corflux(ivp, istr, :) = corflux(ivp, istr+1, :)
          corfluy(ivp, istr, :) = corfluy(ivp, istr+1, :)
          corflux(ivp, istr-1, :) = corflux(ivp, istr+1, :)
          corfluy(ivp, istr-1, :) = corfluy(ivp, istr+1, :)
      endif
      if (EASTERN_EDGE) then
          corflux(ivp, iend, :) = corflux(ivp, iend-1, :)
          corfluy(ivp, iend, :) = corfluy(ivp, iend-1, :)
          corflux(ivp, iend+1, :) = corflux(ivp, iend-1, :)
          corfluy(ivp, iend+1, :) = corfluy(ivp, iend-1, :)
      endif
      if (SOUTHERN_EDGE) then
          corflux(ivp, :, jstr) = corflux(ivp, :, jstr+1)
          corfluy(ivp, :, jstr) = corfluy(ivp, :, jstr+1)
          corflux(ivp, :, jstr-1) = corflux(ivp, :, jstr+1)
          corfluy(ivp, :, jstr-1) = corfluy(ivp, :, jstr+1)
      endif
      if (NORTHERN_EDGE) then
          corflux(ivp, :, jend) = corflux(ivp, :, jend-1)
          corfluy(ivp, :, jend) = corfluy(ivp, :, jend-1)
          corflux(ivp, :, jend+1) = corflux(ivp, :, jend-1)
          corfluy(ivp, :, jend+1) = corfluy(ivp, :, jend-1)
      endif
! ! corners
      if ((SOUTHERN_EDGE) .and. (WESTERN_EDGE))then
          corflux(ivp, istr, jstr) = corflux(ivp, istr+1, jstr+1)
          corfluy(ivp, istr, jstr) = corfluy(ivp, istr+1, jstr+1)
      endif
      if ((NORTHERN_EDGE) .and. (WESTERN_EDGE))then
          corflux(ivp, istr, jend) = corflux(ivp, istr+1, jend-1)
          corfluy(ivp, istr, jend) = corfluy(ivp, istr+1, jend-1)
        endif
      if ((NORTHERN_EDGE) .and. (EASTERN_EDGE))then
          corflux(ivp, iend, jend) = corflux(ivp, iend-1, jend-1)
          corfluy(ivp, iend, jend) = corfluy(ivp, iend-1, jend-1)
        endif
      if ((SOUTHERN_EDGE) .and. (EASTERN_EDGE))then
          corflux(ivp, iend, jstr) = corflux(ivp, iend-1, jstr+1)
          corfluy(ivp, iend, jstr) = corfluy(ivp, iend-1, jstr+1)
      endif
    enddo ! ivp

    END SUBROUTINE sed_obc_corflu
!!=============================================================================
   
    SUBROUTINE sed_meshedges_corflu(ifirst, ilast, jfirst, jlast)
 
        !&E------------------------------------------------------------------------
        !&E       ***  ROUTINE sed_meshedges_corflu ***
        !&E
        !&E ** Purpose :  interpolate corflux and corfluy on mesh edges (in u & v) 
        !&E
        !&E ** Description :  make interpolation from rho point to u,v points
        !&E
        !&E ** Called by : MUSTANG_update
        !&E--------------------------------------------------------------------------
        !! * Arguments
        INTEGER,INTENT(IN) :: ifirst, ilast, jfirst, jlast

        !! * Local declarations
        INTEGER :: i,j,iv, ivp
        REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: tmpx
        REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: tmpy

        !! * Executable part
        DO iv = isand1, isand2           
            tmpy(:, :) = corfluy(iv, :, :)
            tmpx(:, :) = corflux(iv, :, :)
            DO j = jfirst, jlast+1
                DO i = ifirst, ilast+1
                    corflux(iv, i, j) = 0.5_rsh * (tmpx(i-1, j) + tmpx(i, j))
                    corfluy(iv, i, j) = 0.5_rsh * (tmpy(i, j-1) + tmpy(i, j))
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE sed_meshedges_corflu
!!=============================================================================

 
! **TODO** code for CROCO if needed
! SUBROUTINE bathy_actu_fromfile(h0)
! subroutine sed_exchange_hxe_MARS(iwhat,xh0,xssh) 
! SUBROUTINE sed_exchange_maskbedload_MARS
! SUBROUTINE sed_outsaverestart(h0)

#endif

END MODULE sed_MUSTANG_CROCO
