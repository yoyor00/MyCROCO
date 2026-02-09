      integer NDAS,NDASp
!      parameter (NDAS=32)
! Temporary test with CMEMS MED REA geopot z levels decimation 1:2, 1:3
!     parameter (NDAS=62)
      parameter (NDAS=42)

      parameter (NDASp=NDAS+1)
!
! minimization lbfgs
!
      integer NDIM
      PARAMETER(NDIM=
     &        +NDAS*(Lm+2)*(Mm+2)     !T
     &        +NDAS*(Lm+2)*(Mm+2)     !S
#if !defined DAS_GEOS_STRONG || defined DAS_HFRADAR
     &        +NDAS*(Lm+1)*(Mm+1)     !psi
     &        +NDAS*(Lm+2)*(Mm+2)     !chi
#endif
#if !defined DAS_HYDRO_STRONG
     &        +(Lm+2)*(Mm+2)
#endif
     &                                  )
      integer MSAVE,NXM,NWORK
      PARAMETER(MSAVE=5,NXM=NDIM*MSAVE,NWORK=NDIM+2*MSAVE)
!
! in-situ observation
!
      integer max_prf
      PARAMETER( max_prf=NDAS)   !!! fixed
!
      integer max_sio, max_whoi, max_flight,
     &        max_hfradar,max_hfradar6,max_js1,max_tp
      PARAMETER(max_sio=3000,max_whoi=100,max_flight=30000,
     &         max_hfradar=19000,max_hfradar6=10000,max_js1=30000,
     &         max_tp=35000)
      integer max_TS_obs
      PARAMETER(max_TS_obs=1200)
! ctd
!
      integer max_ptsur
      PARAMETER(max_ptsur=100)
      integer max_martn
      PARAMETER(max_martn=100)
      integer max_map
      PARAMETER(max_map=5000)
!
! inversion of observational error variance
! hcmin, h<hcmin, no ssh assimilated
!
      real sio_ot_raw,sio_os_raw,whoi_ot_raw,whoi_os_raw,
     &                   ptsur_ot_raw,ptsur_os_raw,
     &                   martn_ot,martn_os,map_ot_raw,map_os_raw,
     &                   js1_ossh,tp_ossh,swot_ossh,hcmin,sats_osss,
     &                   hf_ouv,hf_ouv6
      PARAMETER( sio_ot_raw=1.0/(0.6*0.6),sio_os_raw=1.0/(0.06*0.06),
     &        whoi_ot_raw=1.0/(0.6*0.6),whoi_os_raw=1.0/(0.06*0.06),
     &        map_ot_raw=1.0/(3.6*3.6),map_os_raw=1.0/(0.38*0.38),
     &        ptsur_ot_raw=1.0/(0.60*0.60),ptsur_os_raw=1.0/(0.05*0.05),
     &        martn_ot=1.0/(1.00*1.00),martn_os=1.0/(0.15*0.15),
     &        js1_ossh=1.0/(0.027*0.027), swot_ossh=1.0/(0.06*0.06),
     &        tp_ossh=1.0/(0.03*0.03), 
     &        hcmin=400.0, sats_osss=1.0/(0.13*0.13),
     &        hf_ouv=1.0/(0.01*0.01),hf_ouv6=1.0/(0.10*0.10)
     &         )
!
! auv
!
! CalPoly
      integer max_prf_cal,max_cal
      PARAMETER(max_prf_cal=NDAS,max_cal=100)
! Dorado
      integer max_prf_dor,max_dor
      PARAMETER(max_prf_dor=NDAS,max_dor=180)
!
! inversion of observational error variance
!
      real cal_ot,cal_os,dor_ot,dor_os
      PARAMETER( cal_ot=1.0/(1.5*1.5),
     &           cal_os=1.0/(0.25*0.25),
     &           dor_ot=1.0/(1.5*1.5),
     &           dor_os=1.0/(0.25*0.25)
     &         )
!
! geostrophic ratio
!
      real geo_ratio, cross_tsp
      PARAMETER( geo_ratio = 1.0, cross_tsp=0.)
      integer sm_rad   ! for smoothing radius
      PARAMETER( sm_rad=3)

      integer Local_len
      PARAMETER( Local_len=24)   ! correlation set zero beyond this length
!
      real sz_s_rad                       ! smoothing total SSH for smoothed SSH obs, both mapped and along track
      PARAMETER( sz_s_rad=5.0)            ! sz_s_rad_len should be 3 times larger than sz_s_rad
      integer sz_s_rad_len                ! used in das_zeta_das_smooth.F and das_zeta_s_smooth.F
      PARAMETER(sz_s_rad_len=14)           

      real sts_s_rad                       ! smoothing TS for gridded TS obs
      PARAMETER( sts_s_rad=3.0)            ! sts_s_rad_len should be 3 times larger than sts_s_rad
      integer sts_s_rad_len                ! used das_ts_das_smooth.F and das_ts_s_smooth.F 
      PARAMETER(sts_s_rad_len=8)          !  be integer

      real sp_rad                       ! length for smoothing pressure  
      PARAMETER( sp_rad=5.0)            ! before the calculation of geostrophic velocity
      integer sp_rad_len               
      PARAMETER( sp_rad_len=13)        

      real sz_rad                       ! after the minimization, protective smoothing
      PARAMETER( sz_rad=1.5)            ! used das_zeta_smooth.F
      integer sz_rad_len               
      PARAMETER( sz_rad_len=4)        

      real regp, reguv               !  Regular parameter                                                       
      PARAMETER( regp=0.4/(10.*10.) )
      PARAMETER( reguv=0.5/(0.6*0.6)  )   ! uv error 0.55 m/s

!
! increment of the surface pressure anomaly
!
! (1-alph) * zeta_h + beta * ( zeta_sm - zeta_h_sm ) 
!   + gama * zeta_h_sm
      real alph0, beta0, gama0
      PARAMETER( alph0 = 0.0, beta0=0.6, gama0=0.0 )

! Adimensional grid space WMED only from avg pn,pm 0.00049
#if defined DAS_ADIMPSICHI
      real adp
      PARAMETER( adp = 2025. )
#endif
! test with huge velocities hfradar in m/s 
#if defined DAS_HUGEHRADAR
      real fhuge
!     cr_amphuge, fvmaskhuge
      PARAMETER( fhuge = 1000.0 )
#endif
! Foundation temperature => number of impacted DAS levels + weights
! fndlev : array of weights oriented upward from NDAS-(nfndlev-1) to NDAS
! In das_innov.h sst_mc, mc_mask, mc_ion have a 3rd dimension
! which is also oriented upward from 1 to nfndlev
#if defined DAS_FDN_MCSST    
      integer nfndlev
      PARAMETER( nfndlev = 3 )
      real fndlev(nfndlev)
      PARAMETER( fndlev = (/0.45,0.40,0.15/) )
#endif
