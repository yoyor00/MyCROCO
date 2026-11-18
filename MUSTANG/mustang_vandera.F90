! Copyright (C) 2022-2026 IFREMER
! License: CeCILL-C
! See LICENSES/LICENSE_MUSTANG.txt

#include "cppdefs.h"

MODULE mustang_vandera

#if defined MUSTANG

! No use, no commons, this module is independant.
!
   IMPLICIT NONE
   PRIVATE

! subroutines
   PUBLIC :: vandera_qbw

! parameters
   INTEGER, PARAMETER :: rsh = 8
   REAL(kind=rsh), PARAMETER   :: vandera_epsilon = 1.D-8
   REAL(kind=rsh), PARAMETER   :: pi = atan(1.D0)*4.D0 ! pi value in the module
   REAL(kind=rsh), PARAMETER   :: twopi = 2.D0*pi
   REAL(kind=rsh), PARAMETER   :: otwopi = 1.D0/twopi
   REAL(kind=rsh), PARAMETER   :: deg2rad = pi/180.D0 ! convert degree to radian
   REAL(kind=rsh), PARAMETER   :: rad2deg = 180.D0/pi ! convert radian to degree
   REAL(kind=rsh), PARAMETER   :: sqrt2 = SQRT(2.D0)
   REAL(kind=rsh), PARAMETER   :: g = 9.81D0          ! g
   REAL(kind=rsh), PARAMETER   :: rho0 = 1000.D0      ! density of water (kg/m3)

CONTAINS

!===============================================================================
   SUBROUTINE vandera_qbw(depth, Hs, Td, Dwave, rhos, d50, &
                          qbw_x, qbw_y)
!! Computes the bed-load transport flux due to wave skewness and asymmetry,
!! following the Van der A et al. (2013) / SANTOSS half-cycle formulation
!! (as implemented in COAWST by J. Warner), for the zero-current case
!! (BEDLOAD_VANDERA_ZEROCURR): the current contribution to bedload is
!! computed separately by the caller and simply added to qbw_x/qbw_y.

      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) ::  depth !! water depth (m)
      REAL(KIND=rsh), INTENT(IN) ::  Hs    !! wave height (m)
      REAL(KIND=rsh), INTENT(IN) ::  Td    !! wave period (s)
      REAL(KIND=rsh), INTENT(IN) ::  Dwave !! wave propagation direction (rad)
      REAL(KIND=rsh), INTENT(IN) ::  rhos  !! sediment grain density (kg/m3)
      REAL(KIND=rsh), INTENT(IN) ::  d50   !! median grain diameter (m)
      REAL(KIND=rsh), INTENT(INOUT) :: qbw_x !! bedload flux, grid x-direction (kg/m/s)
      REAL(KIND=rsh), INTENT(INOUT) :: qbw_y !! bedload flux, grid y-direction (kg/m/s)

      ! Local variables
      REAL(KIND=rsh) k_wn          ! wave number (rad/m)
      REAL(KIND=rsh) c_w           ! wave celerity (m/s)
      REAL(KIND=rsh) d90           ! 90th-percentile grain diameter (m), here = d50
      REAL(KIND=rsh) d50_mix       ! median grain diameter of local sediment mix (m), here = d50
      REAL(KIND=rsh) smgd          ! (rhos/rho0-1)*g*d50, Shields-parameter denominator
      REAL(KIND=rsh) smgdr         ! sqrt(smgd)*d50*rhos, dimensional bedload scaling factor
      REAL(KIND=rsh) Sd            ! grain diameter (m), here = d50, used in smgdr
      REAL(KIND=rsh) phi_curwave   ! angle between current and wave direction (rad); unused here (zero current)
      REAL(KIND=rsh) urms          ! root-mean-square bed wave orbital velocity (m/s)
      REAL(KIND=rsh) udelta_wbl    ! current velocity at top of wave boundary layer (m/s); 0 here
      REAL(KIND=rsh) ksd_wbl       ! current-related bed roughness at wave boundary layer (m); 0 here
      REAL(KIND=rsh) thck_wbl      ! wave boundary layer thickness (m); 0 here
      REAL(KIND=rsh) ustrc_wbl     ! current friction velocity at wave boundary layer (m/s); 0 here
      REAL(KIND=rsh) skew_asym_r   ! Ruessink/Malarkey-Davies skewness-asymmetry parameter r (-)
      REAL(KIND=rsh) skew_asym_phi ! Ruessink/Malarkey-Davies skewness-asymmetry phase parameter phi (rad)
      REAL(KIND=rsh) ursell_no     ! Ursell number, wave nonlinearity parameter (-)
      REAL(KIND=rsh) T_crest       ! duration of the wave crest half-cycle (s)
      REAL(KIND=rsh) T_trough      ! duration of the wave trough half-cycle (s)
      REAL(KIND=rsh) T_cu          ! time from zero up-crossing to velocity peak, crest half-cycle (s)
      REAL(KIND=rsh) T_tu          ! time from zero down-crossing to velocity trough, trough half-cycle (s)
      REAL(KIND=rsh) umax          ! maximum (crest) orbital velocity (m/s)
      REAL(KIND=rsh) umin          ! minimum (trough) orbital velocity (m/s)
      REAL(KIND=rsh) RR_asymwave   ! Abreu et al. velocity skewness parameter (-)
      REAL(KIND=rsh) beta_asymwave ! Abreu et al. velocity asymmetry/phase parameter (rad)
      REAL(KIND=rsh) osmgd         ! 1/smgd
      REAL(KIND=rsh) uhat          ! representative full-cycle orbital velocity amplitude (m/s)
      REAL(KIND=rsh) ahat          ! representative orbital excursion amplitude (m)
      REAL(KIND=rsh) mag_bstrc     ! magnitude of current bed shear stress (m2/s2); 0 here
      REAL(KIND=rsh) fd_wbl        ! current friction factor at wave boundary layer; 0 here
      REAL(KIND=rsh) alpha         ! current/wave friction-factor weighting coefficient; 0 here
      REAL(KIND=rsh) eta           ! ripple height (m)
      REAL(KIND=rsh) ksw           ! wave-related bed roughness (m)
      REAL(KIND=rsh) tau_wRe       ! wave Reynolds stress contribution (m2/s2); 0 here
      REAL(KIND=rsh) alphac        ! current contribution scale factor for bedload (-)
      REAL(KIND=rsh) alphaw        ! wave contribution scale factor for bedload (-)
      REAL(KIND=rsh) ucrest_r      ! representative crest half-cycle orbital velocity (m/s)
      REAL(KIND=rsh) uhat_c        ! crest half-cycle peak orbital velocity (m/s), = umax
      REAL(KIND=rsh) dsf_c         ! sheet-flow layer thickness, crest half-cycle (m)
      REAL(KIND=rsh) theta_cx      ! x-component of Shields parameter, crest half-cycle (-)
      REAL(KIND=rsh) theta_cy      ! y-component of Shields parameter, crest half-cycle (-)
      REAL(KIND=rsh) mag_theta_c   ! magnitude of Shields parameter, crest half-cycle (-)
      REAL(KIND=rsh) wavecycle     ! half-cycle selector for sandload_vandera: +1 crest, -1 trough
      REAL(KIND=rsh) om_cc         ! sediment load entrained and transported within the crest half-cycle (-)
      REAL(KIND=rsh) om_ct         ! sediment load entrained in crest, transported in following trough half-cycle (-)
      REAL(KIND=rsh) utrough_r     ! representative trough half-cycle orbital velocity (m/s)
      REAL(KIND=rsh) uhat_t        ! trough half-cycle peak orbital velocity (m/s), = umin
      REAL(KIND=rsh) dsf_t         ! sheet-flow layer thickness, trough half-cycle (m)
      REAL(KIND=rsh) theta_tx      ! x-component of Shields parameter, trough half-cycle (-)
      REAL(KIND=rsh) theta_ty      ! y-component of Shields parameter, trough half-cycle (-)
      REAL(KIND=rsh) mag_theta_t   ! magnitude of Shields parameter, trough half-cycle (-)
      REAL(KIND=rsh) om_tt         ! sediment load entrained and transported within the trough half-cycle (-)
      REAL(KIND=rsh) om_tc         ! sediment load entrained in trough, transported in following crest half-cycle (-)
      REAL(KIND=rsh) cff1          ! scratch coefficient for the final bedload summation
      REAL(KIND=rsh) cff2          ! scratch coefficient for the final bedload summation
      REAL(KIND=rsh) cff3          ! scratch coefficient for the final bedload summation
      REAL(KIND=rsh) bedld_cx      ! crest half-cycle bedload, wave-aligned x-axis
      REAL(KIND=rsh) bedld_cy      ! crest half-cycle bedload, wave-aligned y-axis
      REAL(KIND=rsh) bedld_tx      ! trough half-cycle bedload, wave-aligned x-axis
      REAL(KIND=rsh) bedld_ty      ! trough half-cycle bedload, wave-aligned y-axis
      REAL(KIND=rsh) bedld_x       ! total (crest+trough) bedload, wave-aligned x-axis (kg/m/s)
      REAL(KIND=rsh) bedld_y       ! total (crest+trough) bedload, wave-aligned y-axis (kg/m/s)
!
! Wave and sediment parameters
!

      k_wn = kh_calc(Td, depth)/depth ! Wave number
      c_w = twopi/(k_wn*Td)          ! Wave speed
      d50_mix = d50                  ! (m)
      d90 = d50                      ! (m)
      smgd = (rhos/rho0 - 1.D0)*g*d50
      osmgd = 1.D0/smgd
      Sd = d50
      smgdr = SQRT(smgd)*Sd*rhos
      alphac = 1.0D0   ! bedload scale fac for curr contrib
      alphaw = 1.0D0   ! bedload scale fac for wave contrib
!
! Compute bed wave orbital velocity (m/s) & excursion amplitude
!
      urms = twopi/Td*(Hs/2.D0)/(SINH(k_wn*depth) + vandera_epsilon)
      uhat = urms*sqrt2
      ahat = uhat*Td*otwopi
!
! Magnitude of stress for computing current velocity
! at the wave boundary layer
      mag_bstrc = 0.D0
!
! Compute near bottom current at wave boundary layer (wbl) in the zero current case
!
      udelta_wbl = 0.D0
      ksd_wbl = 0.D0
      thck_wbl = 0.D0
      ustrc_wbl = 0.D0
! phi_curwave (angle between current and wave direction) is not used by
! full_wave_cycle_stress_factors / half_wave_cycle_stress_factors in this
! zero-current implementation (BEDLOAD_VANDERA_ZEROCURR), but is set here
! to avoid passing it on to them unset.
      phi_curwave = 0.D0
!
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "r" and "phi"
! common to both crest and trough.
!
      call skewness_params(Hs, Td, depth, skew_asym_r, skew_asym_phi, ursell_no)
!
! Abreu et al. use skewness params (r, phi) to get representative critical
! orbital velocity for crest and trough cycles.
!
      call abreu_points(skew_asym_r, skew_asym_phi, uhat, Td, &
                        T_crest, T_trough, T_cu, T_tu, umax, &
                        umin, RR_asymwave, beta_asymwave)
!
! Representative crest half cycle water particle velocity
! as well as full cycle orbital velocity and excursion.
! VA2013 Equation 10, 11.
!
      uhat_c = umax
      uhat_t = umin
      ucrest_r = 0.5D0*sqrt2*uhat_c
      utrough_r = 0.5D0*sqrt2*uhat_t
!
!
! Full wave cycle
!
      call full_wave_cycle_stress_factors(d50, d90, osmgd, Td, c_w, depth, &
                                          phi_curwave, RR_asymwave, uhat, &
                                          ahat, umax, umin, mag_bstrc, &
                                          T_crest, T_trough, T_cu, T_tu, &
                                          ksd_wbl, ustrc_wbl, thck_wbl, &
                                          udelta_wbl, fd_wbl, alpha, eta, ksw, &
                                          tau_wRe)
!
! Bed shear stress (Shields parameter) for Crest half cycle
! alpha VA2013 Eqn. 19
!
      call half_wave_cycle_stress_factors(T_cu, T_crest, ahat, ksw, fd_wbl, &
                                          alpha, alphac, alphaw, d50, osmgd, &
                                          ucrest_r, uhat_c, udelta_wbl, &
                                          phi_curwave, tau_wRe, dsf_c, &
                                          theta_cx, theta_cy, mag_theta_c)
!
! Compute sediment load entrained during each crest half cycle
!
      wavecycle = 1.D0
      call sandload_vandera(wavecycle, Hs, Td, depth, RR_asymwave, d50, &
                            d50_mix, rhos, c_w, eta, dsf_c, T_crest, T_cu, &
                            uhat_c, mag_theta_c, om_cc, om_ct)
!
!
! Bed shear stress (Shields parameter) for Trough half cycle
! alpha VA2013 Eqn. 19
!
      call half_wave_cycle_stress_factors(T_tu, T_trough, ahat, ksw, fd_wbl, &
                                          alpha, alphac, alphaw, d50, osmgd, &
                                          utrough_r, uhat_t, udelta_wbl, &
                                          phi_curwave, tau_wRe, dsf_t, &
                                          theta_tx, theta_ty, mag_theta_t)
!
! Compute sediment load entrained during each trough half cycle
!
      wavecycle = -1.D0
      call sandload_vandera(wavecycle, Hs, Td, depth, RR_asymwave, d50, &
                            d50_mix, rhos, c_w, eta, dsf_t, T_trough, T_tu, &
                            uhat_t, mag_theta_t, om_tt, om_tc)

      !
! Compute sediment load entrained during crest and trough half cycles
!
      cff1 = MAX(0.5D0*T_crest/T_cu, 0.D0)
      cff2 = sqrt(mag_theta_c)*(om_cc + cff1*om_tc)
      cff3 = (theta_cx/mag_theta_c)
      bedld_cx = cff2*cff3
      cff3 = (theta_cy/mag_theta_c)
      bedld_cy = cff2*cff3
!
      cff1 = MAX(0.5D0*T_trough/T_tu, 0.D0)
      cff2 = sqrt(mag_theta_t)*(om_tt + cff1*om_ct)
      cff3 = (theta_tx/mag_theta_t)
      bedld_tx = cff2*cff3
      cff3 = (theta_ty/mag_theta_t)
      bedld_ty = cff2*cff3
!
! Total bedload (kg m-1 sec-1). Use velocity-load equation 1 of VA2013
!
      bedld_x = smgdr*(bedld_cx*T_crest + bedld_tx*T_trough)/Td
      bedld_y = smgdr*(bedld_cy*T_crest + bedld_ty*T_trough)/Td
!
! Partition bedload into xi and eta directions, still at rho points
! (FX_r and FE_r have dimensions of kg).
! --> project on CROCO grid (bedld is wave-aligned in Van Der A)
!
      qbw_x = (bedld_x*COS(Dwave) - bedld_y*SIN(Dwave))
      qbw_y = (bedld_x*SIN(Dwave) + bedld_y*COS(Dwave))

!
   END SUBROUTINE vandera_qbw

!===============================================================================
   SUBROUTINE skewness_params(H_s, Td, depth, skew_asym_r, skew_asym_phi, Ur)
!
! Ruessink et al. (2012) provides equations for calculating skewness
! parameters. Uses Malarkey and Davies equations to get "bb" and "r",
! given input of H_s, T and depth
!
! skew_asym_r     - skewness/asymmetry parameter r=2b/(1+b^2)      [value]
! skew_asym_phi   - skewness/asymmetry parameter                   [value]
! Su    - umax/(umax-umin)                               [value]
! Au    - amax/(amax-amin)                               [value]
! alpha - tmax/pi                                        [value]
!
      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: H_s, Td, depth
      REAL(KIND=rsh), INTENT(INOUT) :: Ur
      REAL(KIND=rsh), INTENT(INOUT) :: skew_asym_r, skew_asym_phi
! Local variables
      REAL(KIND=rsh), PARAMETER :: p1 = 0.D0
      REAL(KIND=rsh), PARAMETER :: p2 = 0.857D0
      REAL(KIND=rsh), PARAMETER :: p3 = -0.471D0
      REAL(KIND=rsh), PARAMETER :: p4 = 0.297D0
      REAL(KIND=rsh), PARAMETER :: p5 = 0.815D0
      REAL(KIND=rsh), PARAMETER :: p6 = 0.672D0
      REAL(KIND=rsh) :: B, psi, bb
      REAL(KIND=rsh) :: k_wn, cff
      REAL(KIND=rsh) :: Su, Au
!
! Ruessink et al. (2012), Coastal Engineering 65:56-63.
!
! k is the local wave number computed with linear wave theory.
!
      k_wn = kh_calc(Td, depth)/depth
!
! Ursell number defined as 3/(32*pi^2) * Hs*L^2/h^3
!               (Fourier coeff for NL terms)
!
      Ur = max(0.375D0*H_s/(k_wn*k_wn*depth**3.D0), 1.D-12)
!
! Ruessink et al. (2012), Eq 9
!
      cff = EXP((p3 - log10(Ur))/p4)
      B = p1 + ((p2 - p1)/(1.D0 + cff))
!
      psi = -90.D0*deg2rad*(1.D0 - TANH(p5/Ur**p6))
!
!
! Markaley and Davies, equation provides bb which is "b" in paper
! Check from where CRS found these equations
!
      bb = sqrt(2.D0)*B/(sqrt(2.D0*B**2.D0 + 9.D0))
      skew_asym_r = 2.D0*bb/(bb**2.D0 + 1.D0)
!
! Ruessink et al., 2012 under Equation 12.
!
      skew_asym_phi = -psi - 0.5D0*pi
!
! Skewness Su and asymmetry Au
!
      Su = B*cos(psi)
      Au = B*sin(psi)
!
      RETURN
   END SUBROUTINE skewness_params
!
!===============================================================================
   SUBROUTINE abreu_points(r, phi, Uw, T, T_c, T_t, T_cu, T_tu, umax, umin, &
                           RR, beta)
!
!  Calculate umax, umin, and phases of asymmetrical wave orbital velocity
!
!  Use the asymmetry parameters from Ruessink et al, 2012
!  to get the umax, umin and phases of asymettrical wave
!  orbital velocity to be used by Van Der A.
!  T_c is duration of crest
!  T_cu Duration of accerating flow within crest half cycle
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: r, phi, Uw, T
      REAL(KIND=rsh), INTENT(INOUT) ::  T_c, T_t, T_cu, T_tu
      REAL(KIND=rsh), INTENT(INOUT) ::  umax, umin, RR, beta
!
! Local variables
      REAL(KIND=rsh) :: b, c, ratio, tmt, tmc, tzd, tzu
      REAL(KIND=rsh) :: omega, w, phi_new
      REAL(KIND=rsh) :: P, F0, betar_0
      REAL(KIND=rsh) :: cff1, cff2, cff
      REAL(KIND=rsh) :: Sk, Ak
!
      omega = 2.D0*pi/T
!
      phi_new = -phi

! Malarkey and Davies (Under equation 16b)
      P = SQRT(1.D0 - r*r)
!
! Malarkey and Davies (Under equation 16b)
!
      b = r/(1.D0 + P)
!
! Appendix E of Malarkey and Davies
!
      c = b*SIN(phi_new)
!
      cff1 = 4.D0*c*(b*b - c*c) + (1.D0 - b*b)*(1.D0 + b*b - 2.D0*c*c)
      cff2 = (1.D0 + b*b)**2.D0 - 4.D0*c*c
      ratio = cff1/cff2
!
! These "if" conditionals prevent ASIN to be between [-1,1] and prevent NaNs
! Not a problem in the MATLAB code
!
      IF (ratio .gt. 1.D0) THEN
         ratio = 1.D0
      END IF
      IF (ratio .lt. -1.D0) THEN
         ratio = -1.D0
      END IF
      tmc = ASIN(ratio)
!
      cff1 = 4.D0*c*(b*b - c*c) - (1.D0 - b*b)*(1.D0 + b*b - 2.D0*c*c)
      cff2 = (1.D0 + b*b)**2.D0 - 4.D0*c*c
      ratio = cff1/cff2
      IF (ratio .gt. 1.D0) THEN
         ratio = 1.D0
      END IF
      IF (ratio .lt. -1.D0) THEN
         ratio = -1.D0
      END IF
      tmt = ASIN(ratio)
!
      IF (tmc .lt. 0.D0) THEN
         tmc = tmc + 2.D0*pi
      END IF
      IF (tmt .lt. 0.D0) THEN
         tmt = tmt + 2.D0*pi
      END IF
!
! Non dimensional umax and umin, under E5 in Malarkey and Davies
!
      umax = 1.D0 + c
      umin = umax - 2.D0
!
!     Dimensionalize
!
      umax = umax*Uw
      umin = umin*Uw
!
! phase of zero upcrossing and downcrossing (radians)
!
      tzu = ASIN(b*SIN(phi_new))
      tzd = 2.D0*ACOS(c) + tzu
!
! MD, equation 17
!
      RR = 0.5D0*(1.D0 + b*SIN(phi_new))
!
! MD, under equation 18
!
      IF (r .le. 0.5D0) THEN
         F0 = 1.0D0 - 0.27D0*(2.D0*r)**(2.1D0)
      ELSE
         F0 = 0.59D0 + 0.14D0*(2.D0*r)**(-6.2D0)
      END IF
!
! MD, Equation 15a,b
!
      IF (r .ge. 0.D0 .and. r .lt. 0.5D0) THEN
         betar_0 = 0.5D0*(1.D0 + r)
      ELSEIF (r .gt. 0.5D0 .and. r .lt. 1.D0) THEN
         cff1 = 4.D0*r*(1.D0 + r)
         cff2 = cff1 + 1.D0
         betar_0 = cff1/cff2
      END IF
!
! MD, Equation 18
!
      cff = SIN((0.5D0*pi - ABS(phi_new))*F0)/SIN(0.5D0*pi*F0)
      beta = 0.5D0 + (betar_0 - 0.5D0)*cff
!
! MD, Table 1, get asymmetry parameterization
! using GSSO (10a,b)
!
      cff = SQRT(2.D0*(1.D0 + b*b)**3.D0)
      Sk = 3.D0*SIN(phi_new)/cff
      Ak = -3.D0*COS(phi_new)/cff
!
! These are the dimensional fractions of wave periods needed by Van der A eqn.
!
      w = 1.D0/omega
      T_c = (tzd - tzu)*w
      T_t = T - T_c
      T_cu = (tmc - tzu)*w
      T_tu = (tmt - tzd)*w
!
      RETURN
   END SUBROUTINE abreu_points

!===============================================================================
   SUBROUTINE full_wave_cycle_stress_factors(d50, d90, osmgd, Td, c_w, depth, &
                                             phi_curwave, RR, uhat, ahat, &
                                             umax, umin, mag_bstrc, &
                                             T_c, T_t, T_cu, T_tu, ksd, ustrc, &
                                             delw, udelta, fd, alpha, eta, &
                                             ksw, tau_wRe)
!
! Iterative solution to obtain current and wave related bed roughness,
! computed from converged time-averaged Shields parameter (VA2013
! Apendix A). Also, compute full-cycle friction factor and wave Reynolds
! stress
!
!  This subroutine returns the following:
!  eta                 : ripple height
!  udelta              : current velocity at the wave boundary layer
!  fd                  : current friction factor
!  tau_wRe             : Wave averaged Reynolds stress
!  T_c, T_t, T_cu, T_tu: Updated time periods in half cycles
!                        based on current velocity
!
      IMPLICIT NONE
!  Imported variable declarations.
!
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles
!
      REAL(KIND=rsh), INTENT(IN) :: d50, d90, osmgd
      REAL(KIND=rsh), INTENT(IN) ::  Td, c_w, depth
      REAL(KIND=rsh), INTENT(IN) ::  phi_curwave
      REAL(KIND=rsh), INTENT(IN) ::  RR, uhat, ahat
      REAL(KIND=rsh), INTENT(IN) ::  umax, umin
      REAL(KIND=rsh), INTENT(IN) ::  mag_bstrc
      REAL(KIND=rsh), INTENT(IN) ::  udelta, delw
      REAL(KIND=rsh), INTENT(INOUT) ::  T_c, T_t, T_cu, T_tu !not used here as we are in BEDLOAD_VANDERA_ZEROCURR
      REAL(KIND=rsh), INTENT(INOUT) ::  fd, alpha, eta, ksw, tau_wRe
!
!  Local variables
      INTEGER iter
      INTEGER, PARAMETER :: total_iters = 15
      REAL(KIND=rsh), PARAMETER :: tol = 0.001D0
      REAL(KIND=rsh), PARAMETER :: von_k = 0.41D0
      REAL(KIND=rsh), PARAMETER :: eps = 1.0D-14
      REAL(KIND=rsh), PARAMETER :: crs_fac = 1.D0
      REAL(KIND=rsh) theta_avg_old, theta_avg, theta_hat_i
      REAL(KIND=rsh) psi  ! maximum mobility number
      REAL(KIND=rsh) rlen ! ripple length
      REAL(KIND=rsh) ksd, ustrc
      REAL(KIND=rsh) fw
      REAL(KIND=rsh) alpha_w, fwd
      REAL(KIND=rsh) ustarw
!
! Maximum mobility number at crest and trough
! For irregular waves, use Rayleigh distributed maximum value
! - VA2013, text under equation Appendix B.4
!
      psi = (1.27D0*uhat)**2.D0*osmgd
!
! Use Appendix B eqn B.1 and B.2 to get ripple height and length
!
      call ripple_dim(psi, d50, eta, rlen)
!
      eta = eta*ahat
      rlen = rlen*ahat
!
! Initiliaze with theta_avg=0 and theta_hat_i=theta_avg
!
      theta_avg = 0.D0
      theta_avg_old = 0.D0
!
! Calculate friction factor fd
!
      fd = 0.D0 ! (fd =0 in the case VDA no current)
!
      do iter = 1, total_iters !---- ! Iterative solution
!
! Wave related bed roughness - VA2013 A.5
!
         ksw = ksw_calc(d50, mu_calc(d50), theta_avg, eta, rlen)
!
! Full-cycle wave friction factor - VA2013 Appendix Eq A.4
!
         fw = fw_calc(ahat, ksw)
!
!
! Current-related bed roughness - VA2013 Appendix A.1
!
         ksd = ksd_calc(d50, d90, mu_calc(d50), theta_avg, eta, rlen)
!
! Friction factor from ksd, udelta, delta
!
         fd = fd_calc_santoss(udelta, delw, ksd)
!
! Time-averaged absolute Shields stress - VA2013 Appendix Eq. A.3
!
         theta_avg = osmgd*0.25D0*fw*uhat**2.D0
!
         IF (ABS(theta_avg - theta_avg_old) .lt. tol) THEN
            EXIT
         END IF
         theta_avg_old = theta_avg
!
      END DO  ! end iterations ----
!
! Full-cycle current friction factor - VA2013 Eq 20
!
      fwd = fw
      alpha = 0.D0
!
! Wave Reynolds stress - VA2013 Eq 22
! (contribution associated with progressive surface waves)
!
      tau_wRe = 0.D0
!
      RETURN
   END SUBROUTINE full_wave_cycle_stress_factors
!

!===============================================================================
   SUBROUTINE half_wave_cycle_stress_factors(T_iu, T_i, ahat, ksw, fd, alpha, &
                                             alphac, alphaw, d50, osmgd, ui_r, &
                                             uhat_i, udelta, phi_curwave, &
                                             tau_wRe, dsf, theta_ix, theta_iy, &
                                             mag_theta_i)
!
!  This subroutine returns the following:
!  dsf                 : sheetflow thickness
!  theta_ix, theta_iy  : Shields parameter in x and y dir.
!  mag_theta_i         : Magnitude of Shields parameter for half cycle
!
      IMPLICIT NONE
!
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles
!
      REAL(KIND=rsh), INTENT(IN) :: T_iu, T_i, ahat, ksw
      REAL(KIND=rsh), INTENT(IN) :: fd, alpha
      REAL(KIND=rsh), INTENT(IN) :: alphac, alphaw
      REAL(KIND=rsh), INTENT(IN) :: d50, osmgd
      REAL(KIND=rsh), INTENT(IN) :: ui_r, uhat_i, udelta, phi_curwave
      REAL(KIND=rsh), INTENT(IN) :: tau_wRe
      REAL(KIND=rsh), INTENT(INOUT) :: dsf, theta_ix, theta_iy, mag_theta_i
!  Local variables
!
      REAL(KIND=rsh), PARAMETER :: eps = 1.0D-14
      REAL(KIND=rsh) :: fw_i, fwd_i
      REAL(KIND=rsh) :: alpha_w, fwd, k, c_w
      REAL(KIND=rsh) :: theta_hat_i
      REAL(KIND=rsh) :: ui_rx, ui_ry, mag_ui
!
! Wave friction factor for wave and crest half cycle - VA2013 Eq 21
!
      fw_i = fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave current friction factor (Madsen and Grant) - VA2013 Eq 18
! Different for crest and trough
!
      fwd_i = fw_i
!
! Magnitude of Shields parameter - VA2013 Eq 17
!
      theta_hat_i = 0.5D0*fwd_i*uhat_i**2.D0*osmgd
!
! Sheet flow thickness - VA2013 Appendix C C.1
!
      dsf = dsf_calc(d50, theta_hat_i) ! this dsf is in m
!
! Compute velocity magnitude based on representative velocities
!  - VA2013 Eq 12
!
! Get the representative trough half cycle water particle velocity
! as well as full cycle orbital velocity and excursion
!
      ui_rx = ui_r*alphaw
      ui_ry = 0.D0
!
! mag_ui is set to a min value to avoid non-zero division
!
      mag_ui = MAX(SQRT(ui_rx*ui_rx), eps)
!
! Magnitude of Shields parameter - VA2013 Eq 17
!
      mag_theta_i = MAX(0.5D0*fwd_i*osmgd*(mag_ui**2.D0), eps)
!
! Shields parameter in crest cycle
! rho0 is required for non-dimensionalizing
!
      theta_ix = ABS(mag_theta_i)*ui_rx/(mag_ui) + tau_wRe*osmgd/rho0
      theta_iy = 0.D0
!
! mag_theta_i is set to a min value to avoid non-zero division
!
      mag_theta_i = MAX(sqrt(theta_ix*theta_ix + theta_iy*theta_iy), eps)
!
      RETURN
   END SUBROUTINE half_wave_cycle_stress_factors
!
!===============================================================================
   SUBROUTINE sandload_vandera(wavecycle, Hs, Td, depth, RR, d50, d50_mix, &
                               rhos, c_w, eta, dsf, T_i, T_iu, &
                               uhat_i, mag_theta_i, &
                               om_ii, om_iy)

      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: wavecycle
      REAL(KIND=rsh), INTENT(IN) ::  Hs, Td, depth, RR
      REAL(KIND=rsh), INTENT(IN) ::  d50, d50_mix, rhos, c_w
      REAL(KIND=rsh), INTENT(IN) ::  eta, dsf
      REAL(KIND=rsh), INTENT(IN) ::  T_i, T_iu
      REAL(KIND=rsh), INTENT(IN) ::  uhat_i, mag_theta_i
      REAL(KIND=rsh), INTENT(INOUT) ::  om_ii, om_iy
      ! Local variables
      REAL(KIND=rsh), PARAMETER :: m_fac = 11.D0
      REAL(KIND=rsh), PARAMETER :: n_fac = 1.2D0
      REAL(KIND=rsh), PARAMETER :: alpha_fac = 8.2D0
      REAL(KIND=rsh), PARAMETER :: xi = 1.7D0 ! Based on Santoss_core.m
      REAL(KIND=rsh), PARAMETER :: eps = 1.0D-14
      REAL(KIND=rsh) :: eps_eff, P
      REAL(KIND=rsh) :: om_i, theta_diff, theta_ieff, theta_cr
      REAL(KIND=rsh) :: w_s, what, w_sc, delta
      REAL(KIND=rsh) :: cff, cff1, d50_8
!
! Layer thickness
!
      IF (eta .gt. 0.D0) THEN
         delta = eta   ! For ripple regime
      ELSE
         delta = dsf   ! For sheet flow regime
      END IF
!
! Find settling velocity based on Soulsby (1997).
! VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
!
      d50_8 = 0.8D0*d50
      w_s = w_s_calc(d50_8, rhos)
!
! Correction by vertical advection (Stokes 2nd-order)
!    wavecycle=+1 => VA2013 Eq 29, for crest cycle
!    wavecycle=-1 => VA2013 Eq 30, for trough cycle
!
      what = what_calc(Hs, Td, depth, RR, delta)
      w_sc = MAX(w_s + wavecycle*what, eps)
!
! VA2013 Eq 27-28, Phase lag parameter
!
      cff = (1.D0 - wavecycle*xi*uhat_i)/c_w
      cff1 = 1.D0/MAX(2.D0*(T_i - T_iu)*w_sc, eps)
      P = alpha_fac*delta*cff*cff1
!
! Hiding/exposure correction factor
!
!      eps_eff=(d50/d50_mix)**0.25
      eps_eff = 1.D0
!
! CRS for multiple sed types
!
      theta_ieff = eps_eff*mag_theta_i
!
! Find critical Shields parameters based on Soulsby (1997).
!
      theta_cr = theta_cr_calc(d50, rhos)
!
! VA2013 Eq 2: Sandload entrained during each half-cycle
!
      theta_diff = MAX((theta_ieff - theta_cr), 0.D0)
      om_i = m_fac*(theta_diff)**n_fac
!
! VA2013 Eq 23-26: Sandload transported within the
!        half-cycle (om_ii) or in the next (om_iy)
!
      IF (P .lt. 1.D0) THEN
         om_ii = om_i
         om_iy = 0.D0
      ELSE
         om_ii = om_i/P
         om_iy = om_i - om_ii
      END IF
!
      RETURN
   END SUBROUTINE sandload_vandera
!
!===============================================================================
   SUBROUTINE ripple_dim(psi, d50, eta, rlen)
! This subroutine returns the following:
! eta, rlen : Ripple dimensions: (height and length)
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: psi, d50
      REAL(KIND=rsh), INTENT(INOUT) :: eta, rlen
      ! Local variables
      REAL(KIND=rsh) :: d50_mm
      REAL(KIND=rsh) :: m_eta, m_lambda, n_eta, n_lambda
      REAL(KIND=rsh), PARAMETER :: eps = 1.D-14
!
      d50_mm = 1000.D0*d50
      IF (d50_mm .lt. 0.22D0) THEN
         m_eta = 0.55D0
         m_lambda = 0.73D0
      ELSEIF (d50_mm .ge. 0.22D0 .and. d50_mm .lt. 0.30D0) THEN
         m_eta = 0.55D0 + (0.45D0*(d50_mm - 0.22D0)/(0.30D0 - 0.22D0))
         m_lambda = 0.73D0 + (0.27D0*(d50_mm - 0.22D0)/(0.30D0 - 0.22D0))
      ELSE
         m_eta = 1.D0
         m_lambda = 1.D0
      END IF
!
! Smooth transition between ripple regime and bed sheet flow regime
!
      IF (psi .le. 190.D0) THEN
         n_eta = 1.D0
      ELSEIF (psi .gt. 190.D0 .and. psi .lt. 240.D0) THEN
         n_eta = 0.5D0*(1.D0 + cos(pi*(psi - 190.D0)/(50.D0)))
      ELSEIF (psi .ge. 240.D0) THEN
         n_eta = 0.D0
      END IF
      n_lambda = n_eta
!
      eta = MAX(0.D0, m_eta*n_eta*(0.275D0 - 0.022D0*psi**0.42D0))
      rlen = MAX(eps, m_lambda*n_lambda*(1.97D0 - 0.44D0*psi**0.21D0))
!
      RETURN
   END SUBROUTINE ripple_dim
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION kh_calc(Td, depth)
! Calculate wave number from Wave period and depth
!
      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: Td, depth
      !Local variables
      REAL(KIND=rsh) :: omega, Khd
      REAL(KIND=rsh) :: cff, x, y, t
!
      omega = 2.D0*pi/Td
!
! Soulsby (2006) "Simplified calculation of wave orbital velocities"
!
      IF (depth .lt. 0.D0) THEN  ! avoid negative wavenumber
         x = 0.D0
      ELSE
         x = omega**2.D0*depth/g
      END IF
!
      IF (x .lt. 1.D0) THEN
         y = SQRT(x)
      ELSE
         y = x
      END IF
!
      t = TANH(y)   ! Iteratively solving 3 times (Soulsby 1997)
      cff = (y*t - x)/(t + y*(1.D0 - t*t))
      y = y - cff
!
      t = TANH(y)
      cff = (y*t - x)/(t + y*(1.D0 - t*t))
      y = y - cff
!
      t = TANH(y)
      cff = (y*t - x)/(t + y*(1.D0 - t*t))
      y = y - cff
!
      kh_calc = y
!
      RETURN
   END FUNCTION kh_calc
!
!===============================================================================
   REAL(kind=rsh) FUNCTION w_s_calc(d50, rhos)
! Critical Shields parameter from Soulsby (1997).
! Dynamics of Marine Sands
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50, rhos
      ! Local variables
      REAL(KIND=rsh), PARAMETER :: nu = 1.36D-6
      REAL(KIND=rsh) :: s, dstar
      REAL(KIND=rsh) :: cff, cff1
!
      s = rhos/rho0
      dstar = (g*(s - 1.D0)/(nu*nu))**(1.D0/3.D0)*d50
      cff = nu/d50
      cff1 = 10.36D0
      w_s_calc = cff*(SQRT(cff1*cff1 + 1.049D0*dstar**3.D0) - cff1)
!
      RETURN
   END FUNCTION w_s_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION what_calc(Hs, Td, depth, RR, zws)
! Second order Stokes theory to get vertical velocity of water particle
! at a given elevation based on santoss_core.m
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) ::  Hs, Td, depth, RR, zws
      ! Local variables
      REAL(KIND=rsh), PARAMETER :: eps_inv = 1.D14
      REAL(KIND=rsh), PARAMETER :: eps = 1.D-14
      REAL(KIND=rsh) :: cff, worb1, worb2, worb
!
      worb1 = pi*Hs*zws/(Td*depth)
      worb2 = max(worb1*2.D0*(RR + RR - 1.D0), eps)
!
!  Using the SANTOSS model formulation
!
      cff = 0.125D0
      worb = cff*worb1* &
             SQRT(64.D0 - &
                  (-worb1 + &
                   SQRT(worb1**2.D0 + 32.D0*worb2**2.D0))**2.D0/ &
                  (worb2**2.D0)) + &
             worb2*SIN(2.*ACOS( &
                       cff*(-worb1 + &
                            SQRT(worb1**2.D0 + 32.D0*worb2**2.D0))/ &
                       worb2))
!
! Prevent worb from going to Infinity when worb2=0.0
!
      worb = MIN(worb, eps_inv)
      what_calc = worb
!
      RETURN
   END FUNCTION what_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION theta_cr_calc(d50, rhos)
! Critical Shields parameter from Soulsby (1997).
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50, rhos
      ! Local variables
      REAL(kind=rsh), PARAMETER :: nu = 1.36D-6
      REAL(KIND=rsh) :: s, dstar
      REAL(KIND=rsh) :: cff1, cff2
!
      s = rhos/rho0
      dstar = (g*(s - 1.D0)/(nu*nu))**(1.D0/3.D0)*d50
      cff1 = 0.3D0/(1.D0 + 1.2D0*dstar)
      cff2 = 0.055D0*(1.-EXP(-0.02D0*dstar))
      theta_cr_calc = cff1 + cff2
!
      RETURN
   END FUNCTION theta_cr_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION fwi_calc(T_iu, T_i, ahat, ksw)
! Wave friction factor for half cycles - VA2013 Eq 21
!
      IMPLICIT NONE

      REAL(KIND=rsh), INTENT(IN) :: T_iu, T_i, ahat, ksw
      ! Local variables
      REAL(kind=rsh), PARAMETER :: eps = 1.0D-14
      REAL(KIND=rsh) :: c1, ratio, fwi
      REAL(KIND=rsh) :: cff1, cff2, cff3
!
      fwi_calc = 0.3D0
!
      c1 = 2.6D0
      ratio = ahat/ksw
      IF (ratio .gt. 1.587D0) THEN
         cff1 = MAX((T_iu/T_i), 0.D0)
         cff2 = (2.D0*cff1)**c1
         cff3 = cff2*ratio
!
! These "if" condition prevents arithematic overflow error
! when 0.0**-0.19
!
         IF (cff3 .le. 0.D0) THEN
            fwi_calc = 0.D0
         ELSE
            fwi_calc = 0.00251D0*EXP(5.21D0*(cff3)**(-0.19D0))
         END IF
      END IF
!
      RETURN
   END FUNCTION fwi_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION dsf_calc(d50, theta_i)
! Sheet flow thickness - VA2013 Appendix C.1.
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50, theta_i
      ! Local variables
      REAL(KIND=rsh) :: d50_mm, cff
!
      d50_mm = d50*1000.D0
      IF (d50_mm .le. 0.15D0) THEN
         cff = 25.D0*theta_i
      ELSEIF (d50_mm .gt. 0.15 .and. d50_mm .lt. 0.2D0) THEN
         cff = 25.D0 - (12.D0*(d50_mm - 0.15D0)/0.05D0)
      ELSEIF (d50_mm .ge. 0.2D0) THEN
         cff = 13.D0*theta_i
      END IF
      dsf_calc = MAX(d50*cff, d50)
!
      RETURN
   END FUNCTION dsf_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION mu_calc(d50)
! Bed roughness factor based on grain size - VA2013 Appendix A
! required for current and wave related bed roughness.
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50
      ! Local variables
      REAL(KIND=rsh) :: d50_mm
!
      d50_mm = d50*1000.D0
!
      if (d50_mm .le. 0.15D0) THEN
         mu_calc = 6.D0
      ELSEIF (d50_mm .gt. 0.15D0 .and. d50_mm .lt. 0.2D0) THEN
         mu_calc = 6.D0 - 5.D0*((d50_mm - 0.15D0)/(0.2D0 - 0.15D0))
      ELSEIF (d50_mm .ge. 0.2D0) THEN
         mu_calc = 1.D0
      END IF
!
      RETURN
   END FUNCTION mu_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION ksw_calc(d50, mu, theta_avg, eta, rlen)
! Wave related bed roughness - VA2013 Eq A.5
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50, mu, theta_avg, eta, rlen
      ! Local variables
      REAL(KIND=rsh) :: ripple_fac, ksw, rlenmax
!
      rlenmax = MAX(rlen, d50)
      ripple_fac = 0.4D0*eta**2.D0/rlenmax
      ksw_calc = MAX(d50, d50*(mu + 6.D0*(theta_avg - 1.D0))) + ripple_fac
!
      RETURN
   END FUNCTION ksw_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION fw_calc(ahat, ksw)
! Full-cycle wave friction factor - VA2013 Eq A.4
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: ahat, ksw
      ! Local variables
      REAL(KIND=rsh) :: ratio, fw
!
      ratio = ahat/ksw
      IF (ratio .gt. 1.587D0) THEN
         fw_calc = 0.00251D0*EXP(5.21D0*(ratio)**(-0.19D0))
      ELSE
         fw_calc = 0.3D0
      END IF
!
      RETURN
   END FUNCTION fw_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION ksd_calc(d50, d90, mu, theta_avg, eta, rlen)
! Current-related bed roughness - VA2013 Appendix A.1
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: d50, d90, mu, theta_avg, eta, rlen
      ! Local variables
      REAL(KIND=rsh) :: ripple_fac, rlenmax
!
      rlenmax = MAX(rlen, d50)
      ripple_fac = 0.4D0*eta**2.D0/rlenmax
      ksd_calc = MAX(3.D0*d90, d50*(mu + 6.D0*(theta_avg - 1.D0))) + ripple_fac
!
      RETURN
   END FUNCTION ksd_calc
!
!===============================================================================
   REAL(KIND=rsh) FUNCTION fd_calc_santoss(udelta, delta, ksd)
! Current friction factor Assuming logarithmic velocity profile. - VA2013 Eq 20
!
      IMPLICIT NONE
      REAL(KIND=rsh), INTENT(IN) :: udelta, delta, ksd
      ! Local variables
      REAL(kind=rsh), PARAMETER :: min_udelta = 0.0001D0
      REAL(kind=rsh), PARAMETER :: von_k = 0.41D0

      IF (udelta .lt. min_udelta) THEN
         fd_calc_santoss = 0.D0
      ELSE
         fd_calc_santoss = 2.D0*(von_k/LOG(30.D0*delta/ksd))**2.D0
      END IF

      RETURN
   END FUNCTION fd_calc_santoss

#endif /* ifdef MUSTANG */

END MODULE mustang_vandera
