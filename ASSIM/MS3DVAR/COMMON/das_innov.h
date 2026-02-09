/* This is include file "das_innov.h". 
  --------------------------------------------
*/
#ifdef DAS_TMISST
      real sst_tmi(GLOBAL_2D_ARRAY)
      real tmi_mask(GLOBAL_2D_ARRAY)
      real tmi_oin(GLOBAL_2D_ARRAY)
!      real dist_coast(GLOBAL_2D_ARRAY)
      common /innov_sst_tmi/sst_tmi
     &       /innov_mask_tmi/tmi_mask 
     &       /innov_oin_tmi/tmi_oin
!      common /innov_dist_coast/dist_coast
#endif
      real dist_coast(GLOBAL_2D_ARRAY)
      common /innov_dist_coast/dist_coast
#ifdef DAS_MCSST
# if !defined DAS_FDN_MCSST
      real sst_mc(GLOBAL_2D_ARRAY)
      real mc_mask(GLOBAL_2D_ARRAY)
      real mc_oin(GLOBAL_2D_ARRAY)
      common /innov_sst_mc/sst_mc
     &       /innov_mask_mc/mc_mask
     &       /innov_oin_mc/mc_oin
# else
      real sst_mc(GLOBAL_2D_ARRAY,nfndlev)
      real mc_mask(GLOBAL_2D_ARRAY,nfndlev)
      real mc_oin(GLOBAL_2D_ARRAY,nfndlev)
      common /innov_sst_mc/sst_mc
     &       /innov_mask_mc/mc_mask
     &       /innov_oin_mc/mc_oin
# endif
#endif
#ifdef DAS_GOES_SST
      real sst_goes(GLOBAL_2D_ARRAY)
      real goes_mask(GLOBAL_2D_ARRAY)
      real goes_oin(GLOBAL_2D_ARRAY)
      common /innov_sst_goes/sst_goes
     &       /innov_mask_goes/goes_mask
     &       /innov_oin_goes/goes_oin
#endif
#ifdef DAS_SWOTSSH
      real ssh_swot(GLOBAL_2D_ARRAY)
      real swot_mask(GLOBAL_2D_ARRAY)
      real swot_oin(GLOBAL_2D_ARRAY)
      common /innov_ssh_swot/ssh_swot
     &       /innov_mask_swot/swot_mask
     &       /innov_oin_swot/swot_oin
#endif
#ifdef DAS_SATSSSS
      real sss_sats(GLOBAL_2D_ARRAY)
      real sats_mask(GLOBAL_2D_ARRAY)
      real sats_oin(GLOBAL_2D_ARRAY)
      common /innov_sss_sats/sss_sats
     &       /innov_mask_sats/sats_mask
     &       /innov_oin_sats/sats_oin
#endif
#ifdef DAS_JASONSSH
      real ssh_js1(GLOBAL_2D_ARRAY)
      real js1_mask(GLOBAL_2D_ARRAY)
      real js1_oin(GLOBAL_2D_ARRAY)
      common /innov_ssh_js1/ssh_js1
     &       /innov_mask_js1/js1_mask
     &       /innov_oin_js1/js1_oin
!
! raw data
!
      real lon_js1(max_js1),lat_js1(max_js1),
     &      ssh_js1_raw(max_js1)
      common /innov_ssh_js1_raw/lon_js1,lat_js1,ssh_js1_raw
#endif
#ifdef DAS_TPSSH
      real ssh_tp(GLOBAL_2D_ARRAY)
      real tp_mask(GLOBAL_2D_ARRAY)
      real tp_oin(GLOBAL_2D_ARRAY)
      common /innov_ssh_tp/ssh_tp
     &       /innov_mask_tp/tp_mask
     &       /innov_oin_tp/tp_oin
!
! raw data
!
      real lon_tp(max_tp),lat_tp(max_tp),
     &      ssh_tp_raw(max_tp)
      common /innov_ssh_tp_raw/lon_tp,lat_tp,ssh_tp_raw
#endif
#if defined DAS_JASONSSH || defined DAS_TPSSH\
  || defined DAS_TPSSH || defined DAS_ERS2SSH\
  || defined DAS_GFOSSH || defined DAS_SWOTSSH
      real ssh_ref(GLOBAL_2D_ARRAY)
      real ref_mask(GLOBAL_2D_ARRAY)
      common /ssh_ref_yavg/ssh_ref
     &       /ssh_ref_mask/ref_mask
#endif
!
! IN SITU
!
!
! flight
!
      real sst_fl(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE sst_fl(BLOCK_PATTERN) BLOCK_CLAUSE
      real fl_mask(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE fl_mask(BLOCK_PATTERN) BLOCK_CLAUSE
      common /innov_sst_fl/sst_fl
     &       /innov_mask_fl/fl_mask
!
! raw data
!
      real lon_fl(max_flight),lat_fl(max_flight),
     &      sst_fl_raw(max_flight)
      common /innov_sst_fl_raw/lon_fl,lat_fl,sst_fl_raw
!
! HF radar
!
      real u_hf(GLOBAL_2D_ARRAY)
      real v_hf(GLOBAL_2D_ARRAY)
      real hf_umask(GLOBAL_2D_ARRAY)
      real hf_vmask(GLOBAL_2D_ARRAY)
      common /innov_u_hf/u_hf
     &       /innov_v_hf/v_hf
     &       /innov_umask_hf/hf_umask
     &       /innov_vmask_hf/hf_vmask
!
! raw radar data
!
      real lon_hf(max_hfradar),lat_hf(max_hfradar),
     &      amperr_hf(max_hfradar),fitdif_hf(max_hfradar),
     &      u_hf_raw(max_hfradar),v_hf_raw(max_hfradar)
      common /innov_uv_hf_raw/lon_hf,lat_hf,
     &       amperr_hf,fitdif_hf,u_hf_raw,v_hf_raw
!
! 6KM
! HF radar
!
      real u_hf6(GLOBAL_2D_ARRAY)
      real v_hf6(GLOBAL_2D_ARRAY)
      real hf_umask6(GLOBAL_2D_ARRAY)
      real hf_vmask6(GLOBAL_2D_ARRAY)
      real hf_uierr6(GLOBAL_2D_ARRAY)  !inverse of squared error
      real hf_vierr6(GLOBAL_2D_ARRAY)  !inverse of squared error
      common /innov_u_hf6/u_hf6
     &       /innov_v_hf6/v_hf6
     &       /innov_umask_hf6/hf_umask6
     &       /innov_vmask_hf6/hf_vmask6
     &       /innov_uierr_hf6/hf_uierr6
     &       /innov_vierr_hf6/hf_vierr6
!
! raw radar data
!
      real lon_hf6(max_hfradar6),lat_hf6(max_hfradar6),
     &      amperr_hf6(max_hfradar6),fitdif_hf6(max_hfradar6),
     &      u_hf_raw6(max_hfradar6),v_hf_raw6(max_hfradar6)
      common /innov_uv_hf_raw6/lon_hf6,lat_hf6,
     &       amperr_hf6,fitdif_hf6,u_hf_raw6,v_hf_raw6
!
! cal-poly AUV
! 0-300m in depth
!
      real t_cal(GLOBAL_2D_ARRAY,max_prf_cal),
     &     s_cal(GLOBAL_2D_ARRAY,max_prf_cal)
      real cal_mask(GLOBAL_2D_ARRAY,max_prf_cal)
      common /innov_temp_cal/t_cal
     &       /innov_salt_cal/s_cal
     &       /innov_mask_cal/cal_mask
!
! raw data
!
      real lon_cal(max_cal,max_prf_cal),
     &     lat_cal(max_cal,max_prf_cal),
     &     t_cal_raw(max_cal,max_prf_cal),
     &     s_cal_raw(max_cal,max_prf_cal)
      common /innov_cal_raw/
     &     lon_cal,lat_cal,t_cal_raw,s_cal_raw
!
! Dorado AUV
! 0-300m in depth
!
      real t_dor(GLOBAL_2D_ARRAY,max_prf_dor),
     &     s_dor(GLOBAL_2D_ARRAY,max_prf_dor)
      real dor_mask(GLOBAL_2D_ARRAY,max_prf_dor)
      common /innov_temp_dor/t_dor
     &       /innov_salt_dor/s_dor
     &       /innov_mask_dor/dor_mask
!
! raw data
!
      real lon_dor(max_dor,max_prf_dor),
     &     lat_dor(max_dor,max_prf_dor),
     &     t_dor_raw(max_dor,max_prf_dor),
     &     s_dor_raw(max_dor,max_prf_dor)
      common /innov_dor_raw/
     &     lon_dor,lat_dor,t_dor_raw,s_dor_raw
!
! NEW New new profiles
!
      real t_obs(GLOBAL_2D_ARRAY,ndas-1),
     &    s_obs(GLOBAL_2D_ARRAY,ndas-1),
     &    obs_t_mask(GLOBAL_2D_ARRAY,ndas-1),
     &    obs_s_mask(GLOBAL_2D_ARRAY,ndas-1),
     &    err_tobs(GLOBAL_2D_ARRAY,ndas-1),
     &    err_sobs(GLOBAL_2D_ARRAY,ndas-1),
     &    dep_tobs(GLOBAL_2D_ARRAY,ndas-1),
     &    dep_sobs(GLOBAL_2D_ARRAY,ndas-1)
      common /innov_t_obs/t_obs
     &       /innov_s_obs/s_obs
     &       /innov_mask_tobs/obs_t_mask
     &       /innov_mask_sobs/obs_s_mask
     &       /innov_err_t_obs/err_tobs
     &       /innov_err_s_obs/err_sobs
     &       /innov_dep_t_obs/dep_tobs
     &       /innov_dep_s_obs/dep_sobs
!
! TS obs raw  data
!
      real   t_obs_raw(max_TS_obs,ndas-1),
     &       s_obs_raw(max_TS_obs,ndas-1),
     &       err_tobs_raw(max_TS_obs,ndas-1),
     &       err_sobs_raw(max_TS_obs,ndas-1),
     &       dep_tobs_raw(max_TS_obs,ndas-1),
     &       dep_sobs_raw(max_TS_obs,ndas-1)
      real lon_obs(max_TS_obs), lat_obs(max_TS_obs),
     &     time_obs(max_TS_obs)
      common /obs_prof_st_raw/t_obs_raw, s_obs_raw
     &       /obs_err_st_raw/err_tobs_raw,err_sobs_raw
     &       /dep_dep_st_raw/dep_tobs_raw,dep_sobs_raw
     &       /obs_ts_lon_lat_time/lon_obs,lat_obs,time_obs
      integer I_obs(max_TS_obs), J_obs(max_TS_obs),
     &        K_obs(max_TS_obs), inst_obs(max_TS_obs)
      common /obs_index_st_raw/I_obs,J_obs,K_obs,inst_obs
!
! gliders: sio
!
      real t_sio(GLOBAL_2D_ARRAY,max_prf)
      real s_sio(GLOBAL_2D_ARRAY,max_prf)
      real sio_t_mask(GLOBAL_2D_ARRAY,max_prf)
      real sio_s_mask(GLOBAL_2D_ARRAY,max_prf)
      common /innov_t_sio/t_sio
     &       /innov_s_sio/s_sio
     &       /innov_mask_tsio/sio_t_mask
     &       /innov_mask_ssio/sio_s_mask
!
! raw radar data
!
      real   t_sio_raw(max_sio,max_prf), 
     &       s_sio_raw(max_sio,max_prf)
      real lon_sio(max_sio), lat_sio(max_sio)   
      real sio_ot(max_prf), sio_os(max_prf)
      common /sio_glides_st_raw/t_sio_raw, s_sio_raw, 
     &                   lon_sio, lat_sio, sio_ot, sio_os

!
! gliders: whoi
!
      real t_whoi(GLOBAL_2D_ARRAY,max_prf)
      real s_whoi(GLOBAL_2D_ARRAY,max_prf)
      real whoi_t_mask(GLOBAL_2D_ARRAY,max_prf)
      real whoi_s_mask(GLOBAL_2D_ARRAY,max_prf)
      common /innov_t_whoi/t_whoi
     &       /innov_s_whoi/s_whoi
     &       /innov_mask_twhoi/whoi_t_mask
     &       /innov_mask_swhoi/whoi_s_mask
!
! raw  data
!
      real   t_whoi_raw(max_whoi,max_prf),
     &       s_whoi_raw(max_whoi,max_prf)
      real lon_whoi(max_whoi), lat_whoi(max_whoi)
      real whoi_ot(max_prf), whoi_os(max_prf)
      common /whoi_glides_st_raw/t_whoi_raw, s_whoi_raw,
     &             lon_whoi, lat_whoi,whoi_ot, whoi_os
!      real s_whoi(max_whoi,max_prf), t_whoi(max_whoi,max_prf),
!     &     mask_whoi(max_whoi,max_prf)
!      real lon_whoi(max_whoi), lat_whoi(max_whoi)
!      common /whoi_glides_st/s_whoi,t_whoi,lon_whoi,lat_whoi,
!     &                      mask_whoi
!      integer Iwhoi(max_whoi),Jwhoi(max_whoi)
!      common /whoi_glides_loc/Iwhoi,Jwhoi
!
! pt sur ctd
      real t_ptsur(GLOBAL_2D_ARRAY,max_prf)
      real s_ptsur(GLOBAL_2D_ARRAY,max_prf)
      real ptsur_t_mask(GLOBAL_2D_ARRAY,max_prf)
      real ptsur_s_mask(GLOBAL_2D_ARRAY,max_prf)
      common /innov_t_ptsur/t_ptsur
     &       /innov_s_ptsur/s_ptsur
     &       /innov_mask_tptsur/ptsur_t_mask
     &       /innov_mask_sptsur/ptsur_s_mask
!
! raw  data
!
      real   t_ptsur_raw(max_ptsur,max_prf),
     &       s_ptsur_raw(max_ptsur,max_prf)
      real lon_ptsur(max_ptsur), lat_ptsur(max_ptsur)
      real ptsur_ot(max_prf), ptsur_os(max_prf)
      common /ptsur_glides_st_raw/t_ptsur_raw, s_ptsur_raw,
     &             lon_ptsur, lat_ptsur,ptsur_ot, ptsur_os
!      real s_ptsur(max_ptsur,max_prf), 
!     &     t_ptsur(max_ptsur,max_prf),
!     &     mask_ptsur(max_ptsur,max_prf)
!      real lon_ptsur(max_ptsur), lat_ptsur(max_ptsur)
!      common /ptsur_cdt_st/s_ptsur,t_ptsur,
!     &                      lon_ptsur,lat_ptsur,
!     &                      mask_ptsur
!      integer Iptsur(max_ptsur),Jptsur(max_ptsur)
!      common /ptsur_cdt_loc/Iptsur,Jptsur
! martin ctd
      real s_martn(max_martn,max_prf),
     &     t_martn(max_martn,max_prf),
     &     mask_martn(max_martn,max_prf)
      real lon_martn(max_martn), lat_martn(max_martn)
      common /martn_cdt_st/s_martn,t_martn,
     &                      lon_martn,lat_martn,
     &                      mask_martn
      integer Imartn(max_martn),Jmartn(max_martn)
      common /martn_cdt_loc/Imartn,Jmartn
!                                                                                                                               
! mapped TS data
!                                                                                                                               
      real t_map(GLOBAL_2D_ARRAY,max_prf)                                                                                      
      real s_map(GLOBAL_2D_ARRAY,max_prf)                                                                                      
      real map_t_mask(GLOBAL_2D_ARRAY,max_prf)                                                                                 
      real map_s_mask(GLOBAL_2D_ARRAY,max_prf)                                                                                 
      common /innov_t_map/t_map                                                                                               
     &       /innov_s_map/s_map                                                                                               
     &       /innov_mask_tmap/map_t_mask                                                                                      
     &       /innov_mask_smap/map_s_mask                                                                                      
!                                                                                                                               
! raw map data                                                                                                                 
!                                                                                                                               
      real   t_map_raw(max_map,max_prf),                                                                                      
     &       s_map_raw(max_map,max_prf)                                                                                       
      real lon_map(max_map), lat_map(max_map)                                                                                  
      real map_ot(max_prf), map_os(max_prf)
      common /map_ctd_raw/t_map_raw, s_map_raw,                                                                              
     &         lon_map, lat_map,map_ot,map_os
!
!  numbers
!
      integer prf_num_sio,prf_num_whoi,prf_num_map,
     &            prf_num_TS_obs,
     &            prf_num_ptsur,prf_num_martn,num_flight,
     &            num_hfradar,num_hfradar6,num_js1,num_tp
      common /insitu_num/prf_num_sio,prf_num_whoi,prf_num_map,
     &            prf_num_TS_obs,
     &            prf_num_ptsur,prf_num_martn,num_flight,
     &            num_hfradar,num_hfradar6,num_js1,num_tp
      integer num_cal(max_prf_cal),num_dor(max_prf_dor)
      common /insitu_num_cal/num_cal,num_dor
!
!
! NOTE: tmi_oin and mc_oin are the inverse of their
!       error variances
!
! BLXD ADDED flag_d2c
      logical flag_tmi, flag_mc, flag_goes,flag_js1,
     &        flag_swot,flag_tp,flag_sio,flag_whoi,flag_ptsur,
     &        flag_TS_obs,flag_martn,flag_flight,flag_cal,flag_dor,
     &        flag_map,flag_prof,flag_ship,flag_hfradar,
     &        flag_hfradar6, flag_ssh_ref,
     &        flag_sats,flag_d2c
      common /innov_flags/flag_tmi, flag_mc, flag_goes, 
     &       flag_js1,flag_swot,flag_tp,flag_sio,flag_whoi,flag_ptsur,
     &        flag_TS_obs,flag_martn,flag_flight,flag_cal,flag_dor,
     &        flag_map,flag_prof,flag_ship,flag_hfradar,
     &        flag_hfradar6, flag_ssh_ref,
     &        flag_sats,flag_d2c
!
! observation file: sst
!
      character*250 file_tmi,file_mc,file_goes
      common /innov_file_sst/file_tmi,file_mc,file_goes
!
! observation file: ssh
!
      character*250 file_js1, file_swot,file_tp, file_sats
      common /innov_file_ssh/file_js1,file_swot, file_tp, file_sats
!
! observation file: gliders
!
      character*250 file_sio_glider,
     &             file_whoi_glider
      common /innov_file_glider/file_sio_glider, file_whoi_glider
!
! observation file: gliders
!
      character*250 file_ptsur_ctd,file_martn_ctd,
     &           file_map_ctd,file_prof_ctd, file_prf_TS_obs
      common /innov_file_ctd/file_ptsur_ctd,
     &        file_martn_ctd,file_map_ctd,file_prof_ctd,
     &        file_prf_TS_obs
!
! flight sst
!
      character*250 file_flight_sst,file_cal_auv,
     &            file_dor_auv,file_ship_sst
      common /innov_file_flight/file_flight_sst,
     &              file_cal_auv,file_dor_auv,file_ship_sst
!
! HF FGAT
!
      character*250 file_hfradar_uv,file_hfradar_uv6
      common /innov_file_HF/file_hfradar_uv,file_hfradar_uv6

!
! distance to coast
!
      character*250 file_dist_coast
      common /dist_coast_file/file_dist_coast
!
! reference ssh
!
      character*250 file_ssh_ref
      common /ssh_ref_file/file_ssh_ref
      character*250 file_ssh_ref_js1
      common /ssh_ref_file_js1/file_ssh_ref_js1
!
! out put innovations
!
      character*250 file_innov_sio, file_innov_whoi
      common /innov_print_file/file_innov_sio,
     &                         file_innov_whoi
      character*250 file_innov_ana_sio, file_innov_ana_whoi
      common /innov_print_ana_file/file_innov_ana_sio,
     &                         file_innov_ana_whoi
!
! temp nc array
!
      real anc(Lm+2,Mm+2)
      common /temp_nc/anc
