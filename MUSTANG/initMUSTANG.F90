!------------------------------------------------------------------------------
MODULE initMUSTANG
!------------------------------------------------------------------------------

#include "cppdefs.h"

#ifdef MUSTANG

!&E============================================================================
!&E                   ***  MODULE  initMUSTANG  ***
!&E
!&E
!&E ** Purpose : concerns all subroutines related to sediment initialization
!&E 
!&E ** Description :
!&E     subroutine MUSTANG_init_sediment  ! initialize the 
!&E                                         parameters of sediment module
!&E                                       ! and the sediment layers
!&E     subroutine MUSTANG_param          ! defines parameterization of 
!&E                                         sediment
!&E     subroutine MUSTANG_compatibility  ! verification of compatibility 
!&E                                         between the various parameters
!&E
!&E     subroutine MUSTANG_sedinit        ! initialize the sediment if not 
!&E                                         from file
!&E
!&E     subroutine MUSTANG_init_hsed      ! initialize the sediment thickness 
!&E                                         function sediment parameters 
!&E                                         (MUSTANG_mixsed)
!&E     subroutine MUSTANG_init_param     ! initialize the settling velocity 
!&E                                         of sand variables and erosion 
!&E                                         parameters
!&E     subroutine MUSTANG_init_output    ! initialize dimensions of output 
!&E                                         tables
!&E     subroutine MUSTANG_morphoinit      ! initialize 
!&E     subroutine MUSTANG_morphoinit_mesh ! initialize 
!&E
!&E   If key_MUSTANG_flocmod (floculation module, Verney et all, 2011):
!&E     subroutine flocmod_init           ! initialize flocs characteristics
!&E     subroutine flocmod_kernels        ! computations of 
!&E                                         agregation/fragmentation kernels 
!&E                                         for FLOCMOD
!&E     subroutine flocmod_agregation_statistics ! computation of shear / 
!&E                                         differential settling statistics
!&E
!&E============================================================================
   !! * Modules used
#include "coupler_define_MUSTANG.h"

   USE comMUSTANG
   USE sed_MUSTANG_HOST,  ONLY : sedinit_fromfile
   USE sed_MUSTANG,  ONLY : MUSTANG_E0sand

   !USE sed_MUSTANG_HOST,  ONLY : MUSTANG_morphoinit_mesh,sed_exchange_hxe_HOST
   USE comsubstance
   USE module_substance

   IMPLICIT NONE

   !! * Accessibility
   PUBLIC MUSTANG_init_sediment

   PRIVATE

 CONTAINS
  
!!=============================================================================
   SUBROUTINE MUSTANG_init_sediment(ifirst, ilast, jfirst, jlast,  &
            WATER_ELEVATION,                                       &
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN  
            dhsed,                                                 &
#endif
# ifdef key_MUSTANG_flocmod
            TRANSPORT_TIME_STEP,                                   &
# endif
            h0fondin, z0hydro, WATER_CONCENTRATION )
 
    !&E------------------------------------------------------------------------
    !&E                 ***  ROUTINE MUSTANG_init_sediment  ***
    !&E
    !&E ** Purpose : initialize the sediment layers (only inside the domain,
    !&E              not at boundaries)
    !&E
    !&E ** Description : called at the beginning of the simulation
    !&E                  call MUSTANG routines for initialization
    !&E
    !&E ** Called by : mustang_init_sediment_main
    !&E
    !&E------------------------------------------------------------------------
   !! * Modules used
   USE coupler_MUSTANG,  ONLY : coupl_conv2MUSTANG
   USE sed_MUSTANG,  ONLY : sed_MUSTANG_comp_z0hydro
#ifdef key_MUSTANG_splitlayersurf
   USE sed_MUSTANG,  ONLY : sed_MUSTANG_split_surflayer
#endif
#ifdef key_MUSTANG_V2
   USE sed_MUSTANG,  ONLY : MUSTANGV2_comp_poro_mixsed
#ifdef key_MUSTANG_bedload
   USE sed_MUSTANG_HOST,  ONLY : sed_bottom_slope
#endif
#endif

   !! * Arguments
   INTEGER, INTENT(IN)                    :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),INTENT(IN)              :: h0fondin
   REAL(KIND=rsh),DIMENSION(ARRAY_Z0HYDRO),INTENT(INOUT)        :: z0hydro                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(INOUT):: WATER_ELEVATION                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(IN)       :: WATER_CONCENTRATION  
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN_MUSTANG_byHYDRO  
   REAL(KIND=rsh),DIMENSION(ARRAY_DHSED),INTENT(INOUT)          :: dhsed                       
#endif
#ifdef key_MUSTANG_flocmod
   REAL(KIND=rlg),INTENT(IN)          :: TRANSPORT_TIME_STEP    
#endif 
   !! * Local declarations
    INTEGER   :: i,j,k,iv,isplit
#ifdef key_MUSTANG_V2
   REAL(KIND=rsh)                  :: mass_tot,poro_kij,crel_mud_ini
   REAL(KIND=rsh),DIMENSION(1:nvpc):: frac_sed
   REAL(KIND=rsh),DIMENSION(1:nvp) :: mass_sed
#else
   REAL(KIND=rsh)                  :: somalp
#endif

   !!--------------------------------------------------------------------------
   !! * Executable part

    CALL MUSTANG_param('r')
    CALL MUSTANG_alloc
    CALL MUSTANG_init_param()
#ifdef key_MUSTANG_flocmod
    CALL flocmod_init(TRANSPORT_TIME_STEP)
#endif
#if ! defined key_noTSdiss_insed
    cp_s(:,:)=cp_suni
    emissivity_s(:,:)=emissivity_sed
#endif
    CALL MUSTANG_param('w')  

   h0fond = h0fondin
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! recovery of concentrations at the bottom layer
     ! which could be used to initiate concentrations in the sediment
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! here iappel = 0 only WATER_CONCENTRATION, SALINITY_MOD, TEMPERATURE_MOD are used     
      CALL coupl_conv2MUSTANG(ifirst,ilast,jfirst,jlast,0,BATHY_H0,WATER_ELEVATION,       &
                              WATER_CONCENTRATION )

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! reading new bathy issued from a previous run with morphocoupl 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(l_bathy_actu) THEN
        ! this routine depends on hydrodynamic hos model because reading save file from a revious run
        ! CALL bathy_actu_fromfile(BATHY_H0) **TODO** has not been program for CROCO

      ENDIF

#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluation of slope for bedload
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined MORPHODYN_MUSTANG_byHYDRO
      it_morphoYes=0
#endif
      CALL sed_bottom_slope(ifirst, ilast, jfirst, jlast, BATHY_H0)
#endif


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! definition of initiales conditions in sediment
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (l_repsed) THEN
        CALL sedinit_fromfile(BATHY_H0)
      ELSE
        CALL MUSTANG_sedinit(ifirst,ilast,jfirst,jlast,BATHY_H0)
      END IF
      
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! porosity estimation
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j=jfirst,jlast
        DO i=ifirst,ilast
          IF(ksma(i,j).NE.0)THEN

#ifdef key_MUSTANG_splitlayersurf
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !! Splitting surface layers if too thick
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             isplit=0
             DO k=ksma(i,j),ksma(i,j)-nlayer_surf_sed+1,-1
                 IF(k > ksmi(i,j)) THEN
                    IF(dzs(k,i,j) > dzsmax(i,j) + 5.0_rsh* dzsmin) isplit=1
                 ENDIF
             ENDDO
             IF(isplit==1 ) then            
                 CALL sed_MUSTANG_split_surflayer(i,j,ksma(i,j))
             ENDIF
#endif

            DO k=ksmi(i,j),ksma(i,j)
#ifdef key_MUSTANG_V2
             mass_tot = 0.0_rsh
             DO iv=igrav1,imud2
               mass_sed(iv)=cv_sed(iv,k,i,j)*dzs(k,i,j) 
               mass_tot=mass_tot+mass_sed(iv)
             END DO
             DO iv=igrav1,imud2
               frac_sed(iv)=mass_sed(iv)/mass_tot
             END DO
             crel_mud_ini=(1.0_rsh-poro_mud_ini)*ros(1)
             CALL MUSTANGV2_comp_poro_mixsed(frac_sed,poro_mud_ini,  &
                                crel_mud_ini,poro_kij)
             poro(k,i,j)=poro_kij
             crel_mud(k,i,j)=crel_mud_ini
             poro_mud(k,i,j)=poro_mud_ini
             c_sedtot(k,i,j)=(1.0_rsh-poro(k,i,j))*ros(1)
             DO iv=igrav1,imud2
               cv_sed(iv,k,i,j)=c_sedtot(k,i,j)*frac_sed(iv)
             END DO
#else
             somalp=0.0_rsh
             DO iv=1,nvpc
                somalp=somalp+cv_sed(iv,k,i,j)/ros(iv)
             ENDDO
             poro(k,i,j)=1.0_rsh-somalp
#endif

            ENDDO
          ELSE
            poro(:,i,j)=0.0_rsh
#ifdef key_MUSTANG_V2
            crel_mud(:,i,j)=0.0_rsh
            poro_mud(:,i,j)=0.0_rsh
#endif 
          ENDIF
        ENDDO
      ENDDO
      
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! estimation of sediment heights and bathy bedrock for morphodynamic
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      !IF(l_morphocoupl .AND. .NOT.l_initfromfile)THEN
      ! appelle meme si pas morpho pour si dredging et initialisation de hsed dans tous les cas, meme si initfromfile
      CALL MUSTANG_morphoinit(ifirst, ilast, jfirst, jlast, BATHY_H0, WATER_ELEVATION   &
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN_MUSTANG_byHYDRO  
                                                   ,dhsed               &
#endif
             )
      !ENDIF

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluation of Z0_hydro if l_z0_coupl_init
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      IF(l_z0hydro_coupl_init)CALL sed_MUSTANG_comp_z0hydro(ifirst,ilast,jfirst,jlast,z0hydro)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! initialization of output tables
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      CALL MUSTANG_init_output


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! verification of CPP keys for compatiblity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     CALL MUSTANG_compatibility()

    
#ifdef key_MUSTANG_debug
#if defined SPHERICAL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! initialization of position of point where we want informations for debug
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     IF(lat_debug .NE. 0._rlg .OR. lon_debug .NE. 0._rlg) THEN
       ! locating i,j from latitude, longitude 
        i_MUSTANG_debug=0
        j_MUSTANG_debug=0            
        DO j=jfirst,jlast
        DO i=ifirst,ilast
           IF(LATITUDE(i,j) == lat_debug .AND. LONGITUDE(i,j) == lon_debug) THEN
             i_MUSTANG_debug=i
             j_MUSTANG_debug=j             
           ENDIF
        ENDDO
        ENDDO
        IF(i_MUSTANG_debug==0) THEN
          MPI_master_only write(ierrorlog,*) 'CAUTION : the point for debug is not well known  '
          MPI_master_only write(ierrorlog,*) '*****************************************************'
          MPI_master_only write(ierrorlog,*) 'lat_debug and lon_debug do not correspond to one point in the grid'
          MPI_master_only write(ierrorlog,*) 'choose i_MUSTANG_debug et j_MUSTANG_debug directly in paraMUSTANGV2.txt'
          MPI_master_only write(ierrorlog,*) 'or change lat_debug and lon_debug'
          MPI_master_only write(ierrorlog,*) 'simulation stopped'
          stop    
        ENDIF
     ENDIF
#endif
#endif
    
    PRINT_DBG*, 'END MUSTANG_init_sediment'

  END SUBROUTINE MUSTANG_init_sediment
  
   !!===========================================================================
 
   SUBROUTINE MUSTANG_param(rw)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_param  ***
   !&E
   !&E ** Purpose : reads parameterization of sediment dynamics
   !&E
   !&E ** Description : read namelist file paraMUSTANGV1.txt or paraMUSTANGV2.txt)
   !&E
   !&E ** Note : GRAVITY must be known as a parameters transmitted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E
   !&E ** Called by :  MUSTANG_initialization
   !&E
   !&E ** External calls : tool_datosec(), tool_sectodat()
   !&E
   !&E ** Reference :
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   CHARACTER(LEN=1), INTENT( in ) :: rw

   !! * Local declarations
   LOGICAL               :: l_varassoc
   CHARACTER(len=lchain) :: filepc
   CHARACTER(len=lchain) :: filepc_user
   CHARACTER(len=19)     :: tool_sectodat
   INTEGER               :: iv,ivpc,isubs,IERR_MPI
#ifdef key_CROCO
   INTEGER               :: lstr,lenstr
#endif
   REAL(KIND=rlg)        :: tool_datosec,dtsedc,dtsedd,dtsedb

   NAMELIST/namsedim_init/ l_repsed,filrepsed,l_unised,fileinised,          &
                       date_start_morpho,date_start_dyninsed,               &
                       l_z0seduni,z0seduni,z0sedmud,z0sedbedrock,hseduni,   &
                       cseduni,ksmiuni,ksmauni,sini_sed,tini_sed,           &
                       l_init_hsed,csed_mud_ini,l_initsed_vardiss,poro_mud_ini
   NAMELIST/namsedim_layer/ l_dzsminuni,dzsminuni,l_dzsmaxuni,dzsmaxuni,dzsmax_bottom,  &
                       dzsmin,nlayer_surf_sed,k1HW97,k2HW97,fusion_para_activlayer
   NAMELIST/namsedim_erosion/ activlayer,frmudcr2,coef_frmudcr1,         &
                       x1toce_mud,x2toce_mud,E0_sand_option,E0_sand_para,&
                       n_eros_sand,E0_mud,n_eros_mud,ero_option,xexp_ero,&
                       htncrit_eros,E0_sand_Cst,     &
                       tau_cri_option,tau_cri_mud_option_eroindep,       &
                       l_peph_suspension,l_xexp_ero_cst,l_eroindep_mud,  &
                       l_eroindep_noncoh,E0_mud_para_indep
#ifdef key_MUSTANG_V2
   NAMELIST/namsedim_poro/poro_option,poro_min,Awooster,Bwooster,Bmax_wu 
   NAMELIST/namsedim_bedload/ l_peph_bedload,l_slope_effect_bedload,     &
                              alphabs,alphabn,hmin_bedload,l_fsusp
#endif

   NAMELIST/namsedim_bottomstress/ l_fricwave,fricwav,                   &
                                   l_z0hydro_coupl_init,l_z0hydro_coupl, &
                                   coef_z0_coupl,z0_hydro_mud,z0_hydro_bed
   NAMELIST/namsedim_deposition/ cfreshmud,csedmin,cmudcr,aref_sand,     &
                       cvolmaxsort,cvolmaxmel,slopefac
   NAMELIST/namsedim_lateral_erosion/coef_erolat,coef_tauskin_lat,        &
                       l_erolat_wet_cell
   NAMELIST/namsedim_consolidation/ l_consolid,xperm1,xperm2,xsigma1,    &
                       xsigma2,csegreg,csandseg,dt_consolid,subdt_consol
   NAMELIST/namsedim_diffusion/ l_diffused, choice_flxdiss_diffsed,      &
                       xdifs1,xdifs2,xdifsi1,xdifsi2,epdifi,fexcs,dt_diffused
   NAMELIST/namsedim_bioturb/l_bioturb,l_biodiffs,xbioturbmax_part,xbioturbk_part, &
                       dbiotu0_part,dbiotum_part,xbioturbmax_diss,xbioturbk_diss, &
                       dbiotu0_diss,dbiotum_diss,frmud_db_min,frmud_db_max,        &
                       dt_bioturb,subdt_bioturb
   NAMELIST/namsedim_morpho/l_morphocoupl,l_morphomesh,MF,l_bathy_smoothing, &
                            l_dredging,l_MF_dhsed,l_bathy_actu,dt_morpho, &
                            l_transfer2hydro_dhsed                                  

   NAMELIST/namsedoutput/ name_out_hsed,name_out_nblaysed,name_out_dzs,name_out_tauskin,      &
                          name_out_tauskin_c,name_out_tauskin_w,choice_nivsed_out,            &
                          riog_valid_min_hsed,riog_valid_max_hsed,                        &
                          riog_valid_min_nblaysed,riog_valid_max_nblaysed,riog_valid_min_dzs,   &
                          riog_valid_max_dzs,riog_valid_min_tauskin,riog_valid_max_tauskin, &
                          nk_nivsed_out,ep_nivsed_out,epmax_nivsed_out,                   &
                          l_outsed_flx_Bload_all,l_outsed_flx_WS_all,                     &
                          l_outsed_poro,l_outsed_activlayer,l_outsed_surf,                &
                          l_outsed_saltemp,l_outsed_toce,l_outsed_frmudsup,               &
                          l_outsed_bil_Bload_all,l_outsed_fsusp,l_outsed_dzs_ksmax,      &
                          l_outsed_peph,l_outsed_eroiter,l_outsed_z0sed,                  &
                          l_outsed_flx_WS_int,l_outsed_flx_Bload_int,l_outsed_bil_Bload_int
#ifdef key_MUSTANG_V2
#ifdef key_MUSTANG_debug
   NAMELIST/namsedim_debug/ l_debug_effdep,lon_debug,lat_debug, l_debug_erosion,          &
                              i_MUSTANG_debug,j_MUSTANG_debug ,date_start_debug
#endif
#endif

#ifdef key_MUSTANG_flocmod
  NAMELIST/namflocmod/ l_ADS,l_ASH,l_COLLFRAG,f_dp0,f_nf,f_nb_frag,f_alpha,f_beta,f_ater,   &
                             f_ero_frac,f_ero_nbfrag,f_ero_iv,f_mneg_param,f_collfragparam, &
                             f_dmin_frag,f_cfcst,f_fp,f_fy,l_out_MUDtot,l_out_f_dtmin,      &
                             l_out_G,l_out_f_d90,l_out_f_d10
#endif
#if ! defined key_noTSdiss_insed
  NAMELIST/namtempsed/mu_tempsed1,mu_tempsed2,mu_tempsed3,epsedmin_tempsed,  &
                      epsedmax_tempsed,alb,cp_suni,emissivity_sed
#endif

   !!--------------------------------------------------------------------------
   !! * Executable part

IF(rw == 'r')THEN

! namelists reading
!==================

#ifdef key_CROCO
    lstr=lenstr(sedname_must)
    filepc_user=sedname_must(1:lstr)
    filepc=REPFICNAMELIST//'/paraMUSTANG_default.txt'
#else
# ifdef key_MUSTANG_V2
    filepc=REPFICNAMELIST//'/paraMUSTANGV2.txt'
# else
    filepc=REPFICNAMELIST//'/paraMUSTANGV1.txt'
# endif
#endif

    OPEN(50,file=filepc,status='old',form='formatted',access='sequential')
    OPEN(51,file=filepc_user,status='old',form='formatted',access='sequential')

    MPI_master_only write(*,*) '*****************************************************'
    MPI_master_only write(*,*) 'READING DEFAULT AND USER-DEFINED  MUSTANG input file'
    MPI_master_only write(*,*) TRIM(filepc),' AND ',TRIM(filepc_user)
    MPI_master_only write(*,*) '*****************************************************'

    READ(50,namsedim_init)
    READ(50,namsedim_layer)
    READ(50,namsedim_bottomstress)
    READ(50,namsedim_deposition)
    READ(50,namsedim_erosion)
#ifdef key_MUSTANG_V2
    READ(50,namsedim_poro)
    READ(50,namsedim_bedload)
#endif    
    READ(50,namsedim_lateral_erosion)
    READ(50,namsedim_consolidation)
    READ(50,namsedim_diffusion)
    READ(50,namsedim_bioturb)
    READ(50,namsedim_morpho)
#if ! defined key_noTSdiss_insed
    READ(50,namtempsed)
#endif


    READ(51,namsedim_init)
    READ(51,namsedim_layer)
    READ(51,namsedim_bottomstress)
    READ(51,namsedim_deposition)
    READ(51,namsedim_erosion)
#ifdef key_MUSTANG_V2
    READ(51,namsedim_poro)
    READ(51,namsedim_bedload)
#endif
    READ(51,namsedim_lateral_erosion)
    READ(51,namsedim_consolidation)
    READ(51,namsedim_diffusion)
    READ(51,namsedim_bioturb)
    READ(51,namsedim_morpho)
#if ! defined key_noTSdiss_insed
    READ(51,namtempsed)
#endif


    xsigma1sg=xsigma1/GRAVITY
    IF(csed_mud_ini==0.0_rsh)csed_mud_ini=cfreshmud
    tstart_dyninsed=tool_datosec(date_start_dyninsed)
    t_dyninsed=tstart_dyninsed
    tstart_morpho=tool_datosec(date_start_morpho)
    dtsedc=365._rlg*86400._rlg
    MPI_master_only write(*,*) 'init dtsedc=',dtsedc
    dtsedd=dtsedc
    dtsedb=dtsedc
    IF(l_consolid) THEN
       dtsedc=dt_consolid
       MPI_master_only write(*,*) 'init dt_consolid=',dt_consolid
    ELSE
       dt_consolid=dtsedc
       MPI_master_only write(*,*) 'init dt_consolid=',dt_consolid
    ENDIF
    IF(l_diffused) THEN
       dtsedd=dt_diffused
    ELSE
       dt_diffused=dtsedd
    ENDIF
    IF(l_bioturb .OR. l_biodiffs) THEN
       dtsedb=dt_bioturb
    ELSE
       dt_bioturb=dtsedb
    ENDIF
    dt_dyninsed=MIN(dtsedc,dtsedd,dtsedb)
    dtsedc=365._rlg*86400._rlg
    dtsedb=dtsedc
    IF(l_consolid) THEN
      dtsedc=subdt_consol
    ELSE
      subdt_consol=dtsedc
    ENDIF
   IF(l_bioturb) THEN
      dtsedb=subdt_bioturb
    ELSE
      subdt_bioturb=dtsedb
    ENDIF    
    subdt_consol= MIN(dtsedb,dtsedc,dt_dyninsed)
#ifdef key_MUSTANG_V2    
    IF(l_dzsminuni) THEN 
       coeff_dzsmin=0._rsh        
    ELSE           
       coeff_dzsmin=1._rsh  
    ENDIF   
#endif  
    IF(l_morphocoupl) THEN
    t_morpho=tstart_morpho
       IF(l_MF_dhsed) THEN
          MF_dhsed=MF
          MF=1._rsh
       ELSE
          MF_dhsed=0._rsh  ! unused
         ! MF used, read from namelist
       ENDIF
    ELSE
      MF=1.0_rsh
    ENDIF
    cexcs=1.0_rsh-fexcs
    fws2=fricwav*0.5_rsh
    l_dyn_insed=.FALSE.
    IF(l_consolid .OR. l_bioturb .OR. l_diffused .OR. l_biodiffs)l_dyn_insed=.TRUE.

#if defined key_oasis && defined key_oasis_mars_ww3  
    IF(l_ww3morpho)l_transfer2hydro_dhsed = .TRUE.    
#endif

    READ(50,namsedoutput)
#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
    READ(50,namsedim_debug)
    t_start_debug=tool_datosec(date_start_debug)
#endif

    READ(51,namsedoutput)
#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
    READ(51,namsedim_debug)
    t_start_debug=tool_datosec(date_start_debug)
#endif




    ! conversion of thickness reading in mm to m for computation
    epmax_nivsed_out=epmax_nivsed_out/1000.0_rsh
    ep_nivsed_out(:)= ep_nivsed_out(:)/1000.0_rsh
    IF(choice_nivsed_out==1) nk_nivsed_out=ksdmax

#ifdef key_MUSTANG_flocmod
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!! module FLOCULATION  !!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !filepc='./paraflocmod.txt'
    !OPEN(50,file=filepc,status='old',form='formatted',access='sequential')
    READ(50,namflocmod)
    READ(51,namflocmod)

#endif
    CLOSE(50)
    CLOSE(51)                                  
       MPI_master_only write(*,*) 'initEND0 dt_consolid=',dt_consolid

ELSE
       MPI_master_only write(*,*) 'initEND1 dt_consolid=',dt_consolid


! namelists writing for information
!==================================
     IF_MPI (MASTER) THEN
    
        MPI_master_only WRITE(iscreenlog,namsedim_init)
        MPI_master_only WRITE(iscreenlog,namsedim_layer)
        MPI_master_only WRITE(iscreenlog,namsedim_bottomstress)
        MPI_master_only WRITE(iscreenlog,namsedim_deposition)
#ifdef key_MUSTANG_V2
        MPI_master_only WRITE(iscreenlog,namsedim_poro)
        MPI_master_only WRITE(iscreenlog,namsedim_bedload)
#else
        MPI_master_only WRITE(iscreenlog,namsedim_erosion)
#endif
        MPI_master_only WRITE(iscreenlog,namsedim_lateral_erosion)
        MPI_master_only WRITE(iscreenlog,namsedim_consolidation)
        MPI_master_only WRITE(iscreenlog,namsedim_diffusion)
        MPI_master_only WRITE(iscreenlog,namsedim_bioturb)
        MPI_master_only WRITE(iscreenlog,namsedim_morpho)
#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
        MPI_master_only WRITE(iscreenlog,namsedim_debug)
#endif
#if ! defined key_noTSdiss_insed
        MPI_master_only WRITE(iscreenlog,namtempsed)
#endif
   
  !!--------------------------------------------------
  
#ifdef key_MUSTANG_flocmod
   !! module floculation
       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) '    FLOCMOD'
       MPI_master_only WRITE(iscreenlog,*) '***********************'
       MPI_master_only WRITE(iscreenlog,*) 'class  diameter  volume  density  mass Ws'
       DO iv=1,nv_mud
         MPI_master_only WRITE(iscreenlog,*) iv,f_diam(iv),f_vol(iv),f_rho(iv),f_mass(iv),f_ws(iv)
       ENDDO
       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) ' *** PARAMETERS ***'
       MPI_master_only WRITE(iscreenlog,*) 'Primary particle size (f_dp0)                                : ',f_dp0
       MPI_master_only WRITE(iscreenlog,*) 'Fractal dimension (f_nf)                                     : ',f_nf
       MPI_master_only WRITE(iscreenlog,*) 'Flocculation efficiency (f_alpha)                            : ',f_alpha
       MPI_master_only WRITE(iscreenlog,*) 'Floc break up parameter (f_beta)                             : ',f_beta
       MPI_master_only WRITE(iscreenlog,*) 'Nb of fragments (f_nb_frag)                                  : ',f_nb_frag
       MPI_master_only WRITE(iscreenlog,*) 'Ternary fragmentation (f_ater)                               : ',f_ater
       MPI_master_only WRITE(iscreenlog,*) 'Floc erosion (% of mass) (f_ero_frac)                        : ',f_ero_frac
       MPI_master_only WRITE(iscreenlog,*) 'Nb of fragments by erosion (f_ero_nbfrag)                    : ',f_ero_nbfrag
       MPI_master_only WRITE(iscreenlog,*) 'fragment class (f_ero_iv)                                    : ',f_ero_iv
       MPI_master_only WRITE(iscreenlog,*) 'negative mass tolerated before redistribution (f_mneg_param) : ',f_mneg_param
       MPI_master_only WRITE(iscreenlog,*) 'Boolean for differential settling aggregation (L_ADS)        : ',l_ADS
       MPI_master_only WRITE(iscreenlog,*) 'Boolean for shear aggregation (L_ASH)                        : ',l_ASH
       MPI_master_only WRITE(iscreenlog,*) 'Boolean for collision fragmenation (L_COLLFRAG)              : ',l_COLLFRAG
       MPI_master_only WRITE(iscreenlog,*) 'Collision fragmentation parameter (f_collfragparam)          : ',f_collfragparam
       MPI_master_only WRITE(iscreenlog,*) ' '
       MPI_master_only WRITE(iscreenlog,*) '*** END FLOCMOD INIT *** '    
    
       IF (.not.l_ADS .and. .not.l_ASH) THEN
        MPI_master_only write(ierrorlog,*) 'CAUTION : incompatible flocculation kernel options : '
        MPI_master_only write(ierrorlog,*) '*****************************************************'
        MPI_master_only write(ierrorlog,*) 'l_ADS=',l_ADS
        MPI_master_only write(ierrorlog,*) 'l_ASH=',l_ASH
        MPI_master_only write(ierrorlog,*) 'simulation stopped'
        stop
      ENDIF
#endif

    ENDIF_MPI
    
ENDIF
     PRINT_DBG*, 'END MUSTANG_param'   
   
  END SUBROUTINE MUSTANG_param

  !!===========================================================================
 
   SUBROUTINE MUSTANG_compatibility()
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_compatibility  ***
   !&E
   !&E ** Purpose : verifies compatibility between the various parameters (sedim)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : MUSTANG_init_sedimentcas1DV apres cas_init_sedim
   !&E                (in MARS also by cas_init_sed for testcases)
   !&E
   !&E ** External calls : tool_datosec(), tool_sectodat()
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   
   !! * Local declarations
   INTEGER               :: iv,ivpc,isubs,IERR_MPI,iv2
   REAL(rsh)             :: rohomog
 
 
    ! test dissolved variables in sediment must be used without  key_nofluxwat_IWS
#if defined key_nofluxwat_IWS && ! defined key_noTSdiss_insed
       MPI_master_only WRITE(iwarnlog,*)
       MPI_master_only WRITE(iwarnlog,*)'****************************'
       MPI_master_only WRITE(iwarnlog,*)'        WARNING           ' 
       MPI_master_only WRITE(iwarnlog,*)'****************************'
       MPI_master_only WRITE(iwarnlog,*)' You are using CPP key : key_nofluxwat_IWS ' 
       MPI_master_only WRITE(iwarnlog,*)' So you are not taking into account water fluxes threw water/sediment interface ' 
       MPI_master_only WRITE(iwarnlog,*)' which are generated by erosion/ settling/ consolidation..'
       MPI_master_only WRITE(iwarnlog,*)' AND you are NOT using CPP key : key_noTSdiss_insed ' 
       MPI_master_only WRITE(iwarnlog,*)' so you want to simulate dissolved variables in sediment'
       MPI_master_only WRITE(iwarnlog,*)' You will have probably conservativity problems '
       MPI_master_only WRITE(iwarnlog,*)' If you simulate dissolved variables in water AND in sediment, don t use key_nofluxwat_IWS '
#endif


    ! test compatibility tocd /consolidation
    IF_MPI (MASTER) THEN
    IF (l_consolid) THEN
       MPI_master_only WRITE(iwarnlog,*)
       MPI_master_only WRITE(iwarnlog,*)'****************************'
       MPI_master_only WRITE(iwarnlog,*)' You are taking into account CONSOLIDATION PROCESS'
       MPI_master_only WRITE(iwarnlog,*)' and fresh deposit concentration (cfreshmud) = ',cfreshmud
       MPI_master_only WRITE(iwarnlog,*)' and critical stress of deposition (tocd) of each sand and mud variable are equal to :'
       DO isubs=imud1,imud2
         MPI_master_only WRITE(iwarnlog,*)' variable',NAME_SUBS,'  tocd=',tocd(isubs)
       ENDDO
       MPI_master_only WRITE(iwarnlog,*)' If cfreshmud is big (> 100), you have to choose a tocd between 1 and 10'
       MPI_master_only WRITE(iwarnlog,*)' otherwise you have to choose a tocd > 10 or 20'
       
    ELSE
       DO isubs=isand2+1,nvp
         IF (tocd(isubs) > 1.5 ) THEN
           MPI_master_only WRITE(iwarnlog,*)'****************************'
           MPI_master_only WRITE(iwarnlog,*)'  WARNING    '
           MPI_master_only WRITE(iwarnlog,*)' You are not taking into account CONSOLIDATION PROCESS'
           MPI_master_only WRITE(iwarnlog,*)' fresh deposit concentration (cfreshmud) = ',cfreshmud
           MPI_master_only WRITE(iwarnlog,*)' and critical stress of deposition (tocd) of isusb=',isubs,' exceed 1 Pa :',tocd(isubs)
           MPI_master_only WRITE(iwarnlog,*)' it should be smaller (0.5 to 1 ?) '
         !  CALL_MPI MPI_FINALIZE(IERR_MPI)
         !  STOP
         ENDIF
       ENDDO
    ENDIF
    ENDIF_MPI
    
    ! test ros homogen for sand in all cases

    IF_MPI (MASTER) THEN
    IF (l_consolid) THEN
       rohomog=ros(irk_fil(1))
       DO iv=2,nvpc
        IF(ros(irk_fil(iv)) .NE. rohomog) THEN
            MPI_master_only WRITE(ierrorlog,*)
            MPI_master_only WRITE(ierrorlog,*)'****************************'
            MPI_master_only WRITE(ierrorlog,*)' You are taking into account CONSOLIDATION PROCESS'
            MPI_master_only WRITE(ierrorlog,*)' and you must have same density for all constitutive particles'
            MPI_master_only WRITE(ierrorlog,*)' It is not the case : (see in variable.dat) '
            MPI_master_only WRITE(ierrorlog,*)' ros(1)= ',ros(irk_fil(1))
            MPI_master_only WRITE(ierrorlog,*)' ros(',irk_fil(iv),')= ',ros(irk_fil(iv))

            CALL_MPI MPI_FINALIZE(IERR_MPI)
            STOP
         ENDIF
       ENDDO
     ELSE
       IF(isand2 > isand1) THEN
        rohomog=ros(irk_fil(isand1))
        DO iv=isand1+1,isand2
         IF(ros(irk_fil(iv)) .NE. rohomog) THEN
            MPI_master_only WRITE(ierrorlog,*)
            MPI_master_only WRITE(ierrorlog,*)'****************************'
            MPI_master_only WRITE(ierrorlog,*)' you must have same density for all sand particles'
            MPI_master_only WRITE(ierrorlog,*)' It is not the case : (see in variable.dat) '
            MPI_master_only WRITE(ierrorlog,*)' ros(1)= ',ros(irk_fil(isand1))
            MPI_master_only WRITE(ierrorlog,*)' ros(',irk_fil(iv),')= ',ros(irk_fil(iv))

            CALL_MPI MPI_FINALIZE(IERR_MPI)
            STOP
          ENDIF
          DO iv2=isand1,iv
           IF(diam_sed(irk_fil(iv)) > diam_sed(irk_fil(iv2))) THEN
             MPI_master_only WRITE(ierrorlog,*)
             MPI_master_only WRITE(ierrorlog,*)'****************************'
             MPI_master_only WRITE(ierrorlog,*)' The "SAND" variables should be stored in substance file '
             MPI_master_only WRITE(ierrorlog,*)'in order of decreasing diameters (the coarsest sand to the finest sand)'
             MPI_master_only WRITE(ierrorlog,*)' It is not the case : (see in substance file) '

             CALL_MPI MPI_FINALIZE(IERR_MPI)
             STOP
           ENDIF
          ENDDO
        ENDDO
       ENDIF
     ENDIF

    ENDIF_MPI


    ! test compatibility keys sedimento+biolo
#if defined key_biolo && defined key_benthos && ! defined key_Pconstitonly_insed
    IF_MPI (MASTER) THEN
       MPI_master_only WRITE(ierrorlog,*)
       MPI_master_only WRITE(ierrorlog,*)'****************************'
       MPI_master_only WRITE(ierrorlog,*)' You are using key_MUSTANG + key_biolo + key_benthos'
       MPI_master_only WRITE(ierrorlog,*)' It is not compatible'
       MPI_master_only WRITE(ierrorlog,*)' You must add the key key_Pconstitonly_insed and biological variables '
       MPI_master_only WRITE(ierrorlog,*)' will not be included in sedim module but in benthos module -in biolo)'
       MPI_master_only WRITE(ierrorlog,*)
       MPI_master_only WRITE(ierrorlog,*)' OR you want simulate biological variables in sediment module'
       MPI_master_only WRITE(ierrorlog,*)' and you have to use : key_MUSTANG + key_biolo without key_benthos'
        CALL_MPI MPI_FINALIZE(IERR_MPI)
        STOP
    ENDIF_MPI
#endif

   ! test compatibilite key_BLOOM_insed
#if defined key_BLOOM_insed && ! defined key_oxygen
    IF_MPI (MASTER) THEN
       MPI_master_only WRITE(iscreenlog,*)
       MPI_master_only WRITE(iscreenlog,*)'****************************'
       MPI_master_only WRITE(iscreenlog,*)' It is not possible to use key_BLOOM_insed  without the key_oxygen'
       MPI_master_only WRITE(iscreenlog,*)' ATTENTION : use CPP key : key_oxygen in Makefile'
        CALL_MPI MPI_FINALIZE(IERR_MPI)
    ENDIF_MPI
#endif

   ! test compatibility bedload with MUSTANG Version 2 
#if defined key_MUSTANG_bedload && ! defined key_MUSTANG_V2
    IF_MPI (MASTER) THEN
       MPI_master_only WRITE(ierrorlog,*)
       MPI_master_only WRITE(ierrorlog,*)'****************************'
       MPI_master_only WRITE(ierrorlog,*)' You are using key_MUSTANG_bedload but '
       MPI_master_only WRITE(ierrorlog,*)' You must add the key key_MUSTANG_V2 '
       MPI_master_only WRITE(ierrorlog,*)' in order to use the new version V2 of MUSTANG which include bedload behavior'
       MPI_master_only WRITE(ierrorlog,*)
        CALL_MPI MPI_FINALIZE(IERR_MPI)
        STOP
    ENDIF_MPI
#endif

#ifdef key_MUSTANG_V2
  ! test version V2 consolidation/poro_option
    IF_MPI (MASTER) THEN
        IF (l_consolid .AND. poro_option .NE. 2) THEN
           MPI_master_only WRITE(iwarnlog,*)'****************************'
           MPI_master_only WRITE(iwarnlog,*)'  WARNING    '
           MPI_master_only WRITE(iwarnlog,*)' You are  taking into account CONSOLIDATION PROCESS'
           MPI_master_only WRITE(iwarnlog,*)' and have choose poro_option=1'
           MPI_master_only WRITE(iwarnlog,*)' poro_option must be =2. it has been changed automatically'
           MPI_master_only WRITE(iwarnlog,*)
           poro_option=2
        ENDIF
    ENDIF_MPI
#endif

    PRINT_DBG*, 'END MUSTANG_compatibility'   
   
  END SUBROUTINE MUSTANG_compatibility
 
   !!===========================================================================
 
  SUBROUTINE MUSTANG_init_param()
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_init_ero_param  ***
   !&E
   !&E ** Purpose : initialize sand settlings velocities 
   !&E
   !&E ** Description : 
   !&E    - compute for each sand and gravel (for 1 to isand2)  
   !&E       * diamstar
   !&E       * ws_sand
   !&E       * tetacri0
   !&E       * stresscri0
   !&E       * xnielsen
   !&E       * psi_sed (in V2)
   !&E    - compute for each sand (for isand1 to isand2) 
   !&E       * E0_sand (in V2)
   !&E    - compute for each gravel, dans and mud (for 1 to imud2) 
   !&E       * rosmrowsros
   !&E    - ros_sand_homogen = ros(isand1) if there is sand (isand1 > 0)
   !&E
   !&E ** Note : GRAVITY and RHOREF must be known 
   !&E
   !&E ** Called by :  MUSTANG_initialization
   !&E 
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------

   !! * Local declarations
   INTEGER        :: iv 
   REAL(KIND=rsh),PARAMETER :: shield_cri_wu = 0.03_rsh          

   !!--------------------------------------------------------------------------
   !! * Executable part

! valable aussi pour les graviers ???
    DO iv = 1, isand2

      ! here, for simplification, water density = RHOREF 
      ! here, for simplification, viscosity = 1.e-6 constant
      ! 10000 = (1/viscosite**2)**(1/3)
      rosmrowsros(iv) = (ros(iv) - RHOREF) / ros(iv)
      diamstar(iv) = diam_sed(iv) * 10000.0_rsh * (GRAVITY * (ros(iv) / RHOREF - 1.0_rsh))**0.33_rsh
      ! according to Soulsby, 1997, and if viscosity = 10-6 m/s :
      ws_sand(iv)=.000001_rsh*((107.33_rsh+1.049_rsh*diamstar(iv)**3)**0.5_rsh-10.36_rsh)/diam_sed(iv)
      xnielsen(iv)=12.0_rsh*1000.0_rsh/GRAVITY/(1.0_rsh-RHOREF/ros(iv))

       ! Critical shear stress in erosion law (N/m2)
      IF (tau_cri_option == 0) THEN
        tetacri0(iv)=0.3_rsh/(1.0_rsh+1.2_rsh*diamstar(iv))+0.055_rsh*(1.0_rsh-EXP(-0.02_rsh*diamstar(iv)))
        stresscri0(iv)=tetacri0(iv)*GRAVITY*(ros(iv)-RHOREF)*diam_sed(iv)
      ELSE IF (tau_cri_option == 1) THEN ! Wu and Lin (2014)
        stresscri0(iv)=GRAVITY*(ros(iv)-RHOREF)*diam_sed(iv)*shield_cri_wu
      END IF

#ifdef key_MUSTANG_V2
      psi_sed(iv)=LOG(diam_sed(iv)*1000.0_rsh)/LOG(2.0_rsh)
      IF (iv .ge. isand1) THEN
            E0_sand(iv) = MUSTANG_E0sand(diam_sed(iv), stresscri0(iv), ros(iv), ws_sand(iv))
      ENDIF
#endif
    ENDDO ! iv = 1, isand2

    DO iv=imud1,imud2
      rosmrowsros(iv)=(ros(iv)-RHOREF)/ros(iv)
    ENDDO

    ! homogeneous sand density
    IF (isand1 > 0) THEN
       iv=isand1
       ros_sand_homogen=ros(iv)
    ENDIF

  END SUBROUTINE MUSTANG_init_param

   !!===========================================================================
!
  SUBROUTINE MUSTANG_init_hsed(cv_sedini)
!
!     initialize sediment thickness, in order to be coherent with sediment parameters
!     function of cseduni, cvolmax, csed_ini of each sediment,
!      valid only if cv_sed  and dzs not variable
! **********************************************************************
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   REAL(KIND=rsh), INTENT(IN),DIMENSION(nvpc) :: cv_sedini

   !! * Local declarations
   INTEGER        :: iv                
   REAL(KIND=rsh) :: dzsgrv,dzssan,dzsmud,voldepgrv,voldepsan,masdepmud,     &
                     dvolsan,dmasmud,dmasmudgrav,dmasmudsand

   !!---------------------------------------------------------------------------
   !! * Executable part

    
      !!! Volume a deposer
      voldepgrv=0.0_rsh
      voldepsan=0.0_rsh        
      DO iv=igrav1,igrav2
            voldepgrv=voldepgrv+cv_sedini(iv)/ros(iv)*hseduni  ! kg/m3 * m /(kg/m3)  (vol_g)
      ENDDO
      DO iv=isand1,isand2
            voldepsan=voldepsan+cv_sedini(iv)/ros(iv)*hseduni  ! kg/m3 * m /(kg/m3)  (vol_s)
      ENDDO

       !!! Masse a deposer
       masdepmud=0.0_rsh
       DO iv=imud1,imud2
            masdepmud=masdepmud+cv_sedini(iv)*hseduni  ! kg/m2  (m_v)
       ENDDO

       !!! nouvelle Couche de gravier
       IF(voldepgrv.GE.0.0_rsh)THEN
            dzsgrv=voldepgrv/cvolmaxsort
       ENDIF
    
       !!! Couche de sable
       IF(voldepsan.GE.0.0_rsh)THEN
         dvolsan=MIN(voldepsan,dzsgrv*(cvolmaxmel-cvolmaxsort)) ! volume de sable qui s insere dans le gravier (cvolmaxmel-cvolmaxsort)
         dzssan=MAX(0.,(voldepsan-dvolsan)/cvolmaxsort)
       ENDIF
    
       !!! Couche de vase
       IF(masdepmud.GE.0.0_rsh) THEN
            IF(dzsgrv >0.0_rsh .and. dzssan > 0.0_rsh) THEN
               ! melange sable+gravier  
               dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxmel)  ! masse de vase qui s insere dans le gravier (1-cvolmaxmel)
               IF (isand2>isand1 ) THEN
                  dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxmel) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
               ELSE
                  dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxsort) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
               ENDIF
            ELSE IF (dzssan > 0.0_rsh .AND. isand2>isand1 ) THEN
               ! melange sables  
               dmasmudgrav=0.
               dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxmel) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
            ELSE IF(dzsgrv >0.0_rsh .AND. igrav2>igrav1 ) THEN
               ! melange graviers  
               dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxmel)  ! masse de vase qui s insere dans le gravier (1-cvolmaxmel)
               dmasmudsand=0.
            ELSE
               ! un seul type de sediment 
               dmasmudgrav=dzsgrv*csed_mud_ini*(1.0_rsh-cvolmaxsort)  ! masse de vase qui s insere dans le gravier ou sable (1-cvolmaxsort)
               dmasmudsand=dzssan*csed_mud_ini*(1.0_rsh-cvolmaxsort) ! masse de vase qui s insere dans le sable (1-cvolmaxsort)
            ENDIF
         ! nouvelle couche de vase
         dzsmud=(masdepmud-dmasmudgrav-dmasmudsand)/csed_mud_ini
         IF(dzsmud <0.0_rsh) THEN
          dzsmud=0.0_rsh
          IF (masdepmud-dmasmudgrav < 0.0_rsh) THEN
             dmasmudgrav = masdepmud
             dmasmudsand = 0.0_rsh
          ELSE
             dmasmudsand = masdepmud-dmasmudgrav
          ENDIF
        ENDIF     
        dmasmud = dmasmudgrav +dmasmudsand
      ENDIF

      !!! Nouvelle Hauteur totale
      hsed_new = dzsgrv +dzssan +dzsmud

  END SUBROUTINE MUSTANG_init_hsed

!!===========================================================================================
   SUBROUTINE MUSTANG_sedinit(ifirst,ilast,jfirst,jlast,BATHY_H0)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_sedinit  ***
   !&E
   !&E ** Purpose : initialize the sediment geometry & content
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by :  MUSTANG_init_sediment
   !&E                 only if not.l_repsed
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER,INTENT(IN)                                        ::   ifirst,ilast,jfirst,jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)       :: BATHY_H0      

   !! * Local declarations
   CHARACTER*50   :: filesand,filemud
   INTEGER        :: i,j,k,iv,ivg,imud_nobio,IERR_MPI
   REAL(KIND=rsh) :: somalp,dzsf,dzsgrv,dzssan,dzsmud,voldepgrv,voldepsan,masdepmud,h0max_def,     &
                     dvolsan,dmasmud,dmasmudgrav,dmasmudsand,cini_modif,tot_csed_init,dzsmaxmin

   !!--------------------------------------------------------------------------
   !! * Executable part

    dzsmax(PROC_IN_ARRAY)=dzsmaxuni
    z0sed(PROC_IN_ARRAY)=z0seduni
    ksmi(PROC_IN_ARRAY)=ksmiuni
    ksma(PROC_IN_ARRAY)=0
    hsed(PROC_IN_ARRAY)=-valmanq
    dzs(ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq
    cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq
    c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq

    IF (l_unised)THEN
   
      ! uniform bed :
      WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) ksmi(PROC_IN_ARRAY)=ksmiuni
      WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) ksma(PROC_IN_ARRAY)=ksmauni
      DO k=ksdmin,ksdmax
        WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) dzs(k,PROC_IN_ARRAY)=0.0_rsh
        WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) c_sedtot(k,PROC_IN_ARRAY)=0.0_rsh
        DO iv=-1,nv_tot
          WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) cv_sed(iv,k,PROC_IN_ARRAY)=0.0_rsh
        ENDDO
      ENDDO
 
 ! if biological variable (detrital organic matter) are constitutive particulaite variable
 !     readjustment of cini_sed (fraction of csedtot : somme=1)
      cini_modif=0.0_rsh
      imud_nobio=1
      DO iv=1,isand2
        cini_modif=cini_modif+cini_sed(iv)
      ENDDO
      DO iv=imud1,imud2
        IF(unit_modif_mudbio_N2dw(irk_fil(iv)) .NE. 1.0_rsh) THEN
            ! variable bio (Ndet) mis en constitutive mais exprime dans la donnee initiale en mmole/kg de sediment total
            !cini_sed(iv)=cini_sed(iv)/cseduni
            cini_modif=cini_modif+cini_sed(iv)
        ELSE 
            ! vase non detritique
            IF(cini_sed(iv) > (1.0_rsh-cini_modif)) imud_nobio=iv
            cini_modif=cini_modif+cini_sed(iv)
        ENDIF    
      ENDDO
      IF(cini_modif > 1.0_rsh) cini_sed(imud_nobio)=cini_sed(imud_nobio)-cini_modif+1.0_rsh
   
      tot_csed_init=0.0_rsh
      DO iv=1,nvpc
        ! constitutive variables : cini_sed in fraction of cseduni
        tot_csed_init=tot_csed_init+cini_sed(iv)              
      ENDDO
      IF(tot_csed_init .NE. 1.0_rsh) THEN
        MPI_master_only WRITE(ierrorlog,*)' '
        MPI_master_only WRITE(ierrorlog,*)' inital concentrations of constitutive particulate matter'
        MPI_master_only WRITE(ierrorlog,*)' have not been well defined'
        MPI_master_only WRITE(ierrorlog,*)' See variable.dat to give fraction of total concentration '
        MPI_master_only WRITE(ierrorlog,*)' for each constitutive variable (GRAV, SAND, MUD) (but not MUDB) '
        MPI_master_only WRITE(ierrorlog,*)' : the sum of all must be =1 '
        CALL_MPI MPI_FINALIZE(IERR_MPI)
        STOP         
      ENDIF
      
      ! ajustement du hsed et cv_sed
      ! ==============================
      IF(l_init_hsed .AND. hseduni > 0.0_rsh) THEN
         ALLOCATE(cv_sedini(1:nvpc))
         !!! Concentration massique
         DO iv=1,nvpc
            ! constitutive variables : cini_sed in % of cseduni
            cv_sedini(iv)=cseduni*cini_sed(iv)
         ENDDO

         CALL MUSTANG_init_hsed(cv_sedini)
         IF_MPI (MASTER) THEN
           MPI_master_only WRITE(iscreenlog,*)'NEW TOTAL SEDIMENT THICKNESS AFTER ADJUSTEMENT = ',hsed_new
         ENDIF_MPI
         cseduni=cseduni*hseduni/hsed_new         
         hseduni=hsed_new
         WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) hsed(PROC_IN_ARRAY)=hseduni
         DEALLOCATE(cv_sedini)
      ELSE
         WHERE (BATHY_H0(PROC_IN_ARRAY) /= -valmanq) hsed(PROC_IN_ARRAY)=hseduni
      ENDIF

      DO j=jfirst,jlast
        DO i=ifirst,ilast
          IF(BATHY_H0(i,j) /= -valmanq) THEN
            DO k=ksmi(i,j),ksma(i,j)
              c_sedtot(k,i,j)=0.0_rsh
              DO iv=1,nvpc
                ! constitutive variables : cini_sed in % of cseduni
                cv_sed(iv,k,i,j)=cseduni*cini_sed(iv)        
                c_sedtot(k,i,j)=c_sedtot(k,i,j)+cv_sed(iv,k,i,j)
              ENDDO
              ! porosity estimation moved after dzs evaluation below
#ifdef key_Pconstitonly_insed
              cv_sed(nvpc+1:nv_state,k,i,j)=0.0_rsh 
#else
              DO iv=nvpc+1,nvp
                ! non constitutive variables : cini_sed in conc/kg of sediment
                IF (irkm_var_assoc(iv) >0) THEN
                  cv_sed(iv,k,i,j)=cini_sed(iv)*cv_sed(irkm_var_assoc(iv),k,i,j)
                ELSE
                  cv_sed(iv,k,i,j)=cini_sed(iv)*c_sedtot(k,i,j)
                END IF
              ENDDO
              IF(l_initsed_vardiss) THEN
              ! initalisation avec cv_wat (CW) du fond
                DO iv=nvp+1,nv_adv
                  ! dissolved variables : cini_sed in conc/m3 of pore water
                  cv_sed(iv,k,i,j)=cw_bottom_MUSTANG(iv,i,j)
                ENDDO
              ELSE
                ! initalisation avec cini (variable.dat)
                DO iv=nvp+1,nv_adv
                  ! dissolved variables : cini_sed in conc/m3 of pore water
                  cv_sed(iv,k,i,j)=cini_sed(iv)
                ENDDO
              ENDIF
              DO iv=nv_adv+1,nv_state
                ! + fixed variables : cini_sed in conc/m3 of sediment
                cv_sed(iv,k,i,j)=cini_sed(iv)
              ENDDO
#endif
              IF(l_initsed_vardiss) THEN
              ! initalisation avec valeurs du fond
                cv_sed(-1,k,i,j)=temp_bottom_MUSTANG(i,j)
                cv_sed( 0,k,i,j)=sal_bottom_MUSTANG(i,j)
              ELSE
                cv_sed(-1,k,i,j)=tini_sed
                cv_sed( 0,k,i,j)=sini_sed
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO

    ELSE

      ! non uniform bed :
      ! -----------------------
      PRINT*, 'NON UNIFORM BED COVERAGE'
      PRINT*, 'you have to read a netcdf file describing for each grid cell (file_inised)'
      PRINT*, 'the values for hsed(i,j),z0sed(i,j),ksmi(i,j),ksma(i,j)'
      PRINT*, 'c_sedtot(i,j),cini_sed(iv,iv=1:nvstate)'
      STOP
      
    ENDIF
!
!   END NON UNIFORM BED COVERAGE L_UNISED
!   -------------------------------------------------
!
!   regular vertical mesh
!   ---------------------

    DO j=jfirst,jlast
      DO i=ifirst,ilast
        IF(ksma(i,j).NE.0)THEN
          dzsf=hsed(i,j)/float(ksma(i,j)-ksmi(i,j)+1)
          DO k=ksmi(i,j),ksma(i,j)
            dzs(k,i,j)=dzsf
          ENDDO
        ELSE
          dzs(:,i,j)=0.0_rsh
        ENDIF
      ENDDO
    ENDDO

! DZSmax non uniform
    IF(.NOT. l_dzsmaxuni) THEN
       IF_MPI (MASTER) THEN
            IF(dzsmaxuni == 0.0_rsh) THEN
               MPI_master_only write(ierrorlog,*)' ERROR in paraMUSTANG : dzsmaxuni =0. and l_repsed=.FALSE.'
               MPI_master_only write(ierrorlog,*)' if dzsmaxuni =0. dzsmax must be read in an initial file (l_repsed=T)'
               MPI_master_only write(ierrorlog,*)' if l_repsed=.FALSE. dzsmax must be evaluated from dzsmaxuni not null'

               CALL_MPI MPI_FINALIZE(IERR_MPI)
               STOP
            ENDIF
       ENDIF_MPI

        ! definition of dzsmax
        ! estimated as a function of depth and dzsmaxuni considered as the max value of dzsmax
        ! here minimum value of dzsmax in depth = dzsmaxuni / 100
        ! max depth giving the lowest value of dzsma to inform, setting by default = 50m
       dzsmaxmin=dzsmaxuni/100._rsh
       h0max_def=50.0_rsh
       MPI_master_only write(iscreenlog,*)
       MPI_master_only write(iscreenlog,*)'DZSMAX NON UNIFORM'
       MPI_master_only write(iscreenlog,*)'by default linearly computed in MUSTANG_sedinit (sedim.F90)'
       MPI_master_only write(iscreenlog,*)'from dzsmaxuni (',dzsmaxuni,'m) to dzsmaxuni/100 (',dzsmaxmin,'m) depending on water depth'
       MPI_master_only write(iscreenlog,*)'dzsmax minimum for depth > ',h0max_def
       MPI_master_only write(iscreenlog,*)'PROGRAM your own initialization of dzsmax (or READ a file) in MUSTANG_sedinit'
       dzsmax(PROC_IN_ARRAY)=(h0max_def-BATHY_H0(PROC_IN_ARRAY))/h0max_def*dzsmaxuni+dzsmaxmin
       WHERE (dzsmax(PROC_IN_ARRAY) > dzsmaxuni ) dzsmax(PROC_IN_ARRAY)=dzsmaxuni
       WHERE (dzsmax(PROC_IN_ARRAY) < dzsmaxmin ) dzsmax(PROC_IN_ARRAY)=dzsmaxmin

       ! or to be read in a data file (to be programmed) 
        
    ENDIF 
    
    PRINT_DBG*, 'END sed_MUSTANG_sedINIT'

  END SUBROUTINE MUSTANG_sedinit
  
   !!===========================================================================

   SUBROUTINE MUSTANG_init_output
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE MUSTANG_init_output  ***
   !&E
   !&E ** Purpose : define output arrays in sediment
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by :  MUSTANG_init_sediment
   !&E                 (by main in MARS if testcase)
   !&E 
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#if defined key_BLOOM_insed
   USE bioloinit,  ONLY : ndiag_tot,ndiag_3d_sed,ndiag_2d_sed,ndiag_1d,ndiag_2d
#endif

   !! * Arguments

   !! * Local declarations
   INTEGER        :: k,nk_nivsed_outlu,nv_out
   REAL(KIND=rsh) :: dzs_estim
                     
   !!--------------------------------------------------------------------------
   !! * Executable part

     IF_MPI (MASTER) THEN
           MPI_master_only WRITE(iscreenlog,*)
           MPI_master_only WRITE(iscreenlog,*)'***************************************************************'
           MPI_master_only WRITE(iscreenlog,*)'*****************     SEDIMENT OUTPUT  ************************'
           MPI_master_only WRITE(iscreenlog,*)'***************************************************************'
           MPI_master_only WRITE(iscreenlog,*)
    ENDIF_MPI
    nk_nivsed_outlu=nk_nivsed_out
    IF(choice_nivsed_out == 1 ) THEN
       ALLOCATE(nivsed_out(0:ksdmax+1))
       nivsed_out(0)=0
       DO k=ksdmin,ksdmax
             nivsed_out(k)=ksdmax+1-k
       END DO
       nivsed_out(ksdmax+1)=0
       IF_MPI (MASTER)THEN
           MPI_master_only WRITE(iscreenlog,*)'results in sediment will be save on all the sediment layers (ksdmax)'
       ENDIF_MPI
    ELSE IF(choice_nivsed_out == 2 ) THEN
       IF(nk_nivsed_outlu > ksdmax) THEN
          nk_nivsed_out=ksdmax
       ENDIF
       ALLOCATE(nivsed_out(0:nk_nivsed_out+1))
       nivsed_out(0)=0
       DO k=1,nk_nivsed_out
         nivsed_out(k)=nk_nivsed_out+1-k
       END DO
       nivsed_out(nk_nivsed_out+1)=0       
       IF_MPI (MASTER) THEN
          MPI_master_only WRITE(iscreenlog,*)'results in sediment will be save only for the first ', nk_nivsed_out, &
                              ' layers from sediment surface, one average for the last'
       ENDIF_MPI
    ELSE IF(choice_nivsed_out == 3 ) THEN
       dzs_estim=MIN(dzsmaxuni,hseduni/ksmauni)
       nk_nivsed_out = MIN(ksdmax, INT(epmax_nivsed_out / dzs_estim) +3) 
       ALLOCATE(nivsed_out(0:nk_nivsed_out+1))
       nivsed_out(0)=0
       DO k=1,nk_nivsed_out
         nivsed_out(k)=nk_nivsed_out+1-k
       END DO
       nivsed_out(nk_nivsed_out+1)=0
       IF_MPI (MASTER) THEN
           MPI_master_only WRITE(iscreenlog,*)'results in sediment will be save from sediment surface till ',  &
                                epmax_nivsed_out*1000,'mmeters (integration below)'
           MPI_master_only WRITE(iscreenlog,*)'the last layer saved in the sediment output file will be that  ',  &
                              'the bottom of which exceeds the desired maximum thickness '          
!           MPI_master_only WRITE(iscreenlog,*)'number of sediment layers in output file :',nk_nivsed_out,'INT(',epmax_nivsed_out,'/',dzs_estim,'+3)'
           IF (.NOT. l_dzsmaxuni) THEN
             MPI_master_only WRITE(iscreenlog,*)
!             MPI_master_only WRITE(iscreenlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!             MPI_master_only WRITE(iscreenlog,*)'WARNING : the maximum number of saved sed. layers has been evaluated from dzsmaxuni=',dzsmaxuni
             MPI_master_only WRITE(iscreenlog,*)'            if dzsmax is not uniform, it could be smaller than dzsmaxuni and  &
                           then there may be points where '
             MPI_master_only WRITE(iscreenlog,*)'             this number is not sufficient to describe the maximum thickness'
             MPI_master_only WRITE(iscreenlog,*)
             MPI_master_only WRITE(iscreenlog,*)'This option is not recommended with l_dzsmaxuni=.False.'
!             MPI_master_only WRITE(iscreenlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           ENDIF
       ENDIF_MPI
     ELSE IF(choice_nivsed_out == 4 ) THEN
       IF(nk_nivsed_outlu > 5 ) THEN
            IF_MPI (MASTER) THEN
               MPI_master_only WRITE(iscreenlog,*)'ERROR in namelist namsedoutput (paraMUSTANG.txt) '
               MPI_master_only WRITE(iscreenlog,*)'nk_nivsed_out must be <=5'
               MPI_master_only WRITE(iscreenlog,*)'nk_nivsed_out is automatically set to 5 '
               MPI_master_only WRITE(iscreenlog,*)'and an latter layer (6th) is added to integrate till the bottom sediment  '
            ENDIF_MPI
            nk_nivsed_outlu = 5
       ENDIF        
       nk_nivsed_out =  nk_nivsed_outlu+1   
       ALLOCATE(ep_nivsed_outp1(nk_nivsed_out))
       ep_nivsed_outp1(1:nk_nivsed_outlu)=ep_nivsed_out(1:nk_nivsed_outlu)
       ep_nivsed_outp1(nk_nivsed_out)=10.0_rsh  
       IF_MPI (MASTER) THEN
             MPI_master_only WRITE(iscreenlog,*)'results in sediment will be save on ',nk_nivsed_out-1, &
                 'integrated layers whom thickness are constant and given by user - first is  sediment surface'
             MPI_master_only WRITE(iscreenlog,*)'the ',nk_nivsed_out,'the layer will be an integrated layer till the bottom' 
       ENDIF_MPI  
       ALLOCATE(nivsed_out(0:nk_nivsed_out+1))
       nivsed_out(0)=0
       DO k=1,nk_nivsed_out
           nivsed_out(k)=nk_nivsed_out+1-k
       END DO
       nivsed_out(nk_nivsed_out+1)=0
      ENDIF

       ALLOCATE(var3D_dzs(nk_nivsed_out,PROC_IN_ARRAY))     
       ALLOCATE(var3D_TEMP(nk_nivsed_out,PROC_IN_ARRAY))
       ALLOCATE(var3D_SAL(nk_nivsed_out,PROC_IN_ARRAY))
#ifdef key_Pconstitonly_insed
       nv_out=nvpc
#else
       nv_out=nv_adv
#endif   
       ALLOCATE(var3D_cvsed(nk_nivsed_out,PROC_IN_ARRAY,nv_out))
#ifdef key_BLOOM_insed
       ALLOCATE(var2D_diagsed(PROC_IN_ARRAY,ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d))
       ALLOCATE(var3D_diagsed(nk_nivsed_out,PROC_IN_ARRAY,ndiag_tot-ndiag_3d_sed+1:ndiag_tot))
#endif
#if defined key_MUSTANG_specif_outputs               
       ALLOCATE(var3D_specifout(nk_nivsed_out,PROC_IN_ARRAY,nv_out3Dk_specif))
#endif
   
    PRINT_DBG*, 'END sed_MUSTANG_init_output'

  END SUBROUTINE MUSTANG_init_output
   !!===========================================================================

 SUBROUTINE MUSTANG_morphoinit(ifirst, ilast, jfirst, jlast, BATHY_H0, WATER_ELEVATION  &
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN_MUSTANG_byHYDRO  
                  ,dhsed                                                &
#endif
                  )

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE morphoinit  ***
   !&E
   !&E ** Purpose : initialization of hsed, hsed0, h0_bedrock, hsed_previous,
   !&E                 morpho0 or morphox & morphoy (MARS)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :  MUSTANG_init_sediment 
   !&E                 (in MARS by casinit and some cas tests as casestuar, casobsthydroflume, castflat2dv
   !&E                     and AgrifUser)
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments 
   INTEGER, INTENT(IN)                    :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(INOUT)        :: BATHY_H0                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(INOUT) :: WATER_ELEVATION
#if (defined key_oasis && defined key_oasis_mars_ww3) || defined MORPHODYN_MUSTANG_byHYDRO  
   REAL(KIND=rsh),DIMENSION(ARRAY_DHSED),INTENT(INOUT)           :: dhsed                       
#endif
   !! * Local declarations
   INTEGER                  :: i, j, k

   !!--------------------------------------------------------------------------  
   !! * Executable part

      PRINT_DBG*, 'DANS SED_MORPHO_INIT'

      IF(l_morphocoupl) THEN
           !!  *1*  initialisation of morphox, morphoy  !!!
                  !!! ATTENTION si pas de hx,hy pas de morphox, morphoy
           ! morpho = 1 if morphodynamic effective
           ! morpho = 0 if depth not vary
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        morpho0(ARRAY_morpho)=1.0_rsh
       ! morpho=0 aux limites du domaine

     !  **TODO** To Program
     !   morpho0(boundaries)=0.0_rsh
       
        ! call a special routine if morpho must be intialize at 0 in some meshes
        !     routine in sed_MUSTANG_HOST
      
        IF(l_morphomesh) CALL MUSTANG_morphoinit_mesh
     
      !!  *2*  deallocate unused array   !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if !defined MORPHODYN_MUSTANG_byHYDRO
        IF(.NOT. l_MF_dhsed )  DEALLOCATE (hsed_previous)
#endif
#if defined key_oasis && defined key_oasis_mars_ww3          
        IF (.NOT. l_transfer2hydro_dhsed) DEALLOCATE (dhsed_save)
#endif   
      ENDIF
 
 
       !!  *3*  initialisation of hsed  inside the domain  !!!
      !!        eliminate the meshes at open boundaries
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO j=jfirst,jlast
       DO i=ifirst,ilast
      !  WARNING : not in first and last mesh on i and j axis (boundaries)
          hsed(i,j)=0.0_rsh
          DO k=ksmi(i,j),ksma(i,j)
            hsed(i,j)=hsed(i,j)+dzs(k,i,j)
          ENDDO
       ENDDO
      ENDDO       

      IF(l_morphocoupl .OR. l_dredging)THEN
#if !defined MORPHODYN_MUSTANG_byHYDRO
        IF (l_MF_dhsed) THEN
#endif
         DO j=jfirst,jlast
           DO i=ifirst,ilast
           !  WARNING: not in first and last mesh on i and j axis (boundaries)
             hsed_previous(i,j)=hsed(i,j)
           ENDDO
         ENDDO       
#if !defined MORPHODYN_MUSTANG_byHYDRO
        END IF
#endif
        
        IF(.NOT.l_repsed)THEN   
          hsed0(:,:)=0.0_rsh
          DO j=jfirst,jlast
            DO i=ifirst,ilast
           !  WARNING : not in first and last mesh on i and j axis (boundaries)
              hsed0(i,j)=hsed(i,j)
            ENDDO
          ENDDO
        ENDIF                       
        
        !!  *4*  definition de la cote du socle (sous le sediment) total depth invariable
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO j=jfirst,jlast
          DO i=ifirst,ilast
           IF(BATHY_H0(i,j) .NE. -valmanq)  THEN
              h0_bedrock(i,j)=BATHY_H0(i,j)+hsed(i,j)         ! here BATHY_H0 is the inital bathy if l_bathy_actu=F 
                                             ! or a new bathy readinf from a previous morphodynamic run if l_bathy_actu=T
           ELSE
              h0_bedrock(i,j)=BATHY_H0(i,j)
           ENDIF
          ENDDO
        ENDDO
      ENDIF   ! endif l_morphocoupl (OR l_dredging)

      IF(l_morphocoupl)THEN    

       IF(l_repsed)THEN          

           !!  *5*  update bathy and water elevation if initalisation from file !!!
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO j=jfirst,jlast
          DO i=ifirst,ilast
             !ssh(i,j)=MAX(ssh(i,j),-BATHY_H0(i,j))         
             SURF_ELEVATION_ij=MAX(SURF_ELEVATION_ij,-BATHY_H0(i,j))         
          ENDDO
         ENDDO
         
           !!  *6* echange MPI of BATHY_H0 and WATER_ELEVATION for neighboring cells
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
           !! To Program, has not been program for CROCO **TODO** CALL sed_exchange_hxe_MARS(1,xh0=BATHY_H0,xssh=WATER_ELEVATION)

       ENDIF   ! endif l_repsed


          !!  *9* transfer dhsed (variation from initial time) if need for coupling WW3 with oasis
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
#if defined key_oasis && defined key_oasis_mars_ww3       
         IF (l_transfer2hydro_dhsed) THEN  
           DO j=jfirst,jlast
             DO i=ifirst,ilast
               dhsed(i,j)=hsed0(i,j)-hsed(i,j)
               dhsed_save(i,j)=dhsed(i,j)
             ENDDO
           ENDDO
         END IF
#endif

#if defined MORPHODYN_MUSTANG_byHYDRO
           DO j=jfirst,jlast
             DO i=ifirst,ilast
               dhsed(i,j)=hsed0(i,j)-hsed(i,j)
              ! dhsed_savedd(i,j)=dhsed(i,j)
#ifdef key_MUSTANG_debug
               IF (i==i_MUSTANG_debug .AND. j==j_MUSTANG_debug ) THEN
             ! if (i.eq.20.and.j.eq.2) then
                MPI_master_only write(*,*) 'dhsed(',i,',',j,') initial:',dhsed(i,j),hsed0(i,j),hsed(i,j)
               ENDIF
#endif
             ENDDO
           ENDDO
#endif

      ENDIF        ! endif l_morphocoupl

      PRINT_DBG*, 'END MUSTANG_MORPHOINIT'

  END SUBROUTINE MUSTANG_morphoinit
   !!===========================================================================

  subroutine MUSTANG_morphoinit_mesh 
     !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE MUSTANG_morphoinit_mesh  ***
  !&E
  !&E ** Purpose : initalize morpho0 in the domain
  !&E 
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E --------------------------------------------------------------------------
  !! * Executable part

   !! morpho0(i,j)=1._rsh   !! mesh where morphodynamic coupling is effective
   !! morpho0(i,j)=0._rsh   !! mesh where bathy cannot change (i.e. dam ..)

  PRINT_DBG*, 'END MUSTANG_morphoinit_mesh'
  
  END SUBROUTINE MUSTANG_morphoinit_mesh
   !!===========================================================================  

SUBROUTINE MUSTANG_alloc(l_filesubs)
 
  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE MUSTANG_alloc  ***
  !&E
  !&E ** Purpose : allocation of arrays relative to sediment
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : MUSTANG_initialization 
  !&E                (and in MARS by cas_init_subs_disspart for test cases
  !&E                        and  by Agrif_User for Agrif runs)
  !&E
  !&E
  !&E--------------------------------------------------------------------------

  USE comsubstance

  !! * Arguments
  LOGICAL, INTENT(IN), OPTIONAL   :: l_filesubs
  !! * Local declarations
  INTEGER               :: iv
  CHARACTER(len=lchain) :: filepc
  LOGICAL               :: l_filesubsr

  !!--------------------------------------------------------------------------
  !! * Executable part
  IF(PRESENT(l_filesubs)) THEN
    l_filesubsr=l_filesubs
  ELSE
    l_filesubsr=.FALSE.
  ENDIF
  
  ALLOCATE(ws_sand(nvp))        
  ALLOCATE(diamstar(nvp))
  ALLOCATE(rosmrowsros(nvp))
  ALLOCATE(stresscri0(nvp))
  ALLOCATE(tetacri0(nvp))
  ALLOCATE(xnielsen(nvp))
  ALLOCATE (typart(-1:nv_adv))
  typart(-1:0)=0.0_rsh
  typart(1:nvpc)=1.0_rsh
  typart(nvpc+1:nv_adv)=0.0_rsh
#ifdef key_MUSTANG_V2
  ALLOCATE(E0_sand(nvp))
  E0_sand(1:nvp)=0.0_rsh
#endif   
#if  ! defined key_noTSdiss_insed
  ALLOCATE(ivdiss(-1:nv_adv-nvp))
#endif

     
#ifdef key_MUSTANG_flocmod
  ALLOCATE(f_ws(1:nv_mud))
  f_ws(1:nv_mud)=0.0_rsh
#endif   

#ifdef key_MUSTANG_V2
  ALLOCATE(psi_sed(nvp))
  psi_sed(1:nvp)=0.0_rsh
#endif   

!  allocation of spatial variables  
!  dimensions defined dans coupler_define_MUSTANG.h
!  dimensions in MARS :  PROC_IN_ARRAY       = limin:limax,ljmin:ljmax
!                        PROC_IN_ARRAY_m1p2  = liminm1:limaxp2,ljminm1:ljmaxp2
!                        PROC_IN_ARRAY_m1p1  = liminm1:limaxp1,ljminm1:ljmaxp1
!                        PROC_IN_ARRAY_0p1   = limin:limaxp1,ljmin:ljmaxp1

  ALLOCATE(ksmi(PROC_IN_ARRAY))
  ALLOCATE(ksma(PROC_IN_ARRAY))
  ALLOCATE(hsed(PROC_IN_ARRAY))  
  ALLOCATE(z0sed(PROC_IN_ARRAY))
  ALLOCATE(tauskin(PROC_IN_ARRAY))
  ALLOCATE(tauskin_c(PROC_IN_ARRAY))
  ALLOCATE(tauskin_w(PROC_IN_ARRAY))
  ALLOCATE(ustarbot(PROC_IN_ARRAY))
  ALLOCATE(dzsmax(PROC_IN_ARRAY))
#if ! defined key_noTSdiss_insed
  ALLOCATE(phitemp_s(PROC_IN_ARRAY))
  ALLOCATE(phitemp_sout(PROC_IN_ARRAY))
  ALLOCATE(cp_s(PROC_IN_ARRAY))
  !ALLOCATE(mu_tempsedsurf(PROC_IN_ARRAY))
  ALLOCATE(poro_sedsurf(PROC_IN_ARRAY))
#endif
  ALLOCATE(emissivity_s(PROC_IN_ARRAY))
  ALLOCATE(htot(PROC_IN_ARRAY_m2p2))
  ALLOCATE(alt_cw1(PROC_IN_ARRAY))
  ALLOCATE(epn_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))  
  ALLOCATE(sal_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))  
  ALLOCATE(temp_bottom_MUSTANG(PROC_IN_ARRAY_m1p2))
  ALLOCATE(cw_bottom_MUSTANG(nv_tot,PROC_IN_ARRAY_m1p2))
  ALLOCATE(ws3_bottom_MUSTANG(nvp,PROC_IN_ARRAY_m1p2))  
  ALLOCATE(roswat_bot(PROC_IN_ARRAY))  
#ifdef key_MUSTANG_V2
  ALLOCATE(sigmapsg(ksdmin:ksdmax))
  ALLOCATE(stateconsol(ksdmin:ksdmax))
  ALLOCATE(permeab(ksdmin:ksdmax))
#endif
#ifdef key_sand2D
  ALLOCATE(rouse2D(nv_adv,PROC_IN_ARRAY))
  ALLOCATE(sum_tmp(nv_adv,PROC_IN_ARRAY))
  rouse2D(1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  sum_tmp(1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
#endif

  ksmi(PROC_IN_ARRAY)=0
  ksma(PROC_IN_ARRAY)=0
  hsed(PROC_IN_ARRAY)=0.0_rsh
  z0sed(PROC_IN_ARRAY)=0.0_rsh
  tauskin(PROC_IN_ARRAY)=0.0_rsh
  tauskin_c(PROC_IN_ARRAY)=0.0_rsh
  tauskin_w(PROC_IN_ARRAY)=0.0_rsh
  ustarbot(PROC_IN_ARRAY)=0.0_rsh
  
  dzsmax(PROC_IN_ARRAY)=0.0_rsh
  htot(PROC_IN_ARRAY_m2p2)=0.0_rsh
  alt_cw1(PROC_IN_ARRAY)=0.0_rsh
  emissivity_s(PROC_IN_ARRAY)=0.0_rsh
#if ! defined key_noTSdiss_insed
  phitemp_s(PROC_IN_ARRAY)=0.0_rsh
  phitemp_sout(PROC_IN_ARRAY)=0.0_rsh
  cp_s(PROC_IN_ARRAY)=0.0_rsh
  !mu_tempsedsurf(PROC_IN_ARRAY)=0.0_rsh
  poro_sedsurf(PROC_IN_ARRAY)=0.0_rsh
#endif
  epn_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
  sal_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
  temp_bottom_MUSTANG(PROC_IN_ARRAY_m1p2)=0.0_rsh
  cw_bottom_MUSTANG(nv_tot,PROC_IN_ARRAY_m1p2)=0.0_rsh
  ws3_bottom_MUSTANG(nvp,PROC_IN_ARRAY_m1p2)=0.0_rsh
  roswat_bot(PROC_IN_ARRAY)=0.0_rsh

  ALLOCATE(cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(poro(ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(dzs(ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(flx_s2w(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(flx_w2s(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(flx_w2s_sum(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(corflux(nv_adv, PROC_IN_ARRAY_m1p1))
  ALLOCATE(corfluy(nv_adv, PROC_IN_ARRAY_m1p1))
  ALLOCATE(fludif(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(fluconsol(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(fluconsol_drycell(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(flu_dyninsed(-1:nv_adv,PROC_IN_ARRAY))
  ALLOCATE(gradvit(NB_LAYER_WAT,PROC_IN_ARRAY))

  cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  poro(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  dzs(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  flx_s2w(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  flx_w2s(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  flx_w2s_sum(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  corflux(1:nv_adv, PROC_IN_ARRAY_m1p1) = 1.0_rsh
  corfluy(1:nv_adv, PROC_IN_ARRAY_m1p1) = 1.0_rsh
  fludif(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  fluconsol(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  fluconsol_drycell(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  flu_dyninsed(-1:nv_adv,PROC_IN_ARRAY)=0.0_rsh
  gradvit(NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh

#ifdef key_MUSTANG_specif_outputs
 ! outputs
 ! variables 3D /k
  nv_out3Dk_specif=1
     ! 1 : poro_save  
#if defined key_MUSTANG_add_consol_outputs && defined key_MUSTANG_V2
  nv_out3Dk_specif=nv_out3Dk_specif+8
     ! 2 : loadograv_save
     ! 3 : permeab_save
     ! 4 : sigmapsg_save
     ! 5 : dtsdzs_save
     ! 6 : hinder_save
     ! 7 : sed_rate_save
     ! 8 : sigmadjge_save
     ! 9 : stateconsol_save
#endif
  ALLOCATE(varspecif3Dk_save(nv_out3Dk_specif,ksdmin:ksdmax,PROC_IN_ARRAY))
  varspecif3Dk_save(1:nv_out3Dk_specif,ksdmin:ksdmax,PROC_IN_ARRAY)= 0.0_rsh

   nv_out3Dnv_specif=3
     ! 1 : toce_save
     ! 2 : flx_s2w_save
     ! 3 : flx_w2s_save
#ifdef key_MUSTANG_V2
   nv_out3Dnv_specif=nv_out3Dnv_specif+1
     ! 4 : pephm_fcor_save  
#ifdef key_MUSTANG_bedload
   nv_out3Dnv_specif=nv_out3Dnv_specif+4
     ! 5 : flx_bx
     ! 6 : flx_by
     ! 7 : bil_bedload
     ! 8 : fsusp
#endif
#endif
  ALLOCATE(varspecif3Dnv_save(nv_out3Dnv_specif,nvpc,PROC_IN_ARRAY))
  ALLOCATE(varspecif3Dnv_out(nv_out3Dnv_specif,nvpc,PROC_IN_ARRAY))
  varspecif3Dnv_save(1:nv_out3Dnv_specif,1:nvpc,PROC_IN_ARRAY)= 0.0_rsh
  
   nv_out2D_specif=2
     ! 1 : frmudsup 
     ! 2 : dzs_ksmax 
#ifdef key_MUSTANG_V2
      nv_out2D_specif=nv_out2D_specif+13
     ! 3 : dzs_aclay_comp_save
     ! 4 : dzs_aclay_kept_save
     ! 5 : tero_noncoh (cumulated time (in hours) elapsed in non cohesive regime)
     ! 6 : tero_coh (cumulated time (in hours) elapsed in cohesive regime)
     ! 7 : pct_iter_noncoh
     ! 8 : pct_iter_coh
     ! 9 : niter_ero
     ! 10: z0sed
     ! 11 : flx_s2w_noncoh
     ! 12 : flx_w2s_noncoh
     ! 13 : flx_s2w_coh
     ! 14 : flx_w2s_coh
     ! 15: z0hydro (if l_z0hydro_coupl)
#ifdef key_MUSTANG_bedload
   nv_out2D_specif=nv_out2D_specif+3
     ! 16 : flx_bx_int
     ! 17 : flx_by_int
     ! 18 : bil_bedload_int
#endif
!   end version MUSTANG V2
#endif

   ALLOCATE(varspecif2D_save(nv_out2D_specif,PROC_IN_ARRAY))
   ALLOCATE(varspecif2D_out(nv_out2D_specif,PROC_IN_ARRAY))
   varspecif2D_save(1:nv_out2D_specif,PROC_IN_ARRAY)= 0.0_rsh
   
!   end specifs output
#endif

#ifdef key_MUSTANG_V2
  ALLOCATE(poro_mud(ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(crel_mud(ksdmin:ksdmax,PROC_IN_ARRAY))
  ALLOCATE(l_isitcohesive(PROC_IN_ARRAY))
  poro_mud(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  crel_mud(ksdmin:ksdmax,PROC_IN_ARRAY)=0.0_rsh
  l_isitcohesive(PROC_IN_ARRAY)=.FALSE.

#ifdef key_MUSTANG_bedload
  ALLOCATE( flx_bx(1:nvp,PROC_IN_ARRAY_m1p1)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1
  ALLOCATE( flx_by(1:nvp,PROC_IN_ARRAY_m1p1)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1
  ALLOCATE( slope_dhdx(PROC_IN_ARRAY))
  ALLOCATE( slope_dhdy(PROC_IN_ARRAY))
  ALLOCATE( sedimask_h0plusxe(PROC_IN_ARRAY_m1p1)) ! Warning /Baptiste : m1p1 sur les 2 indices au lieu d1

  flx_bx(1:nvp,PROC_IN_ARRAY_m1p1)=0.0_rsh 
  flx_by(1:nvp,PROC_IN_ARRAY_m1p1)=0.0_rsh
  slope_dhdx(PROC_IN_ARRAY)=0.0_rsh
  slope_dhdy(PROC_IN_ARRAY)=0.0_rsh
  sedimask_h0plusxe(PROC_IN_ARRAY_m1p1)=0.0_rsh
#endif
#endif

  ALLOCATE(phieau_s2w(PROC_IN_ARRAY))
  ALLOCATE(phieau_s2w_drycell(PROC_IN_ARRAY))
  ALLOCATE(phieau_s2w_consol(PROC_IN_ARRAY))
  phieau_s2w(PROC_IN_ARRAY)=0.0_rlg
  phieau_s2w_drycell(PROC_IN_ARRAY)=0.0_rlg
  phieau_s2w_consol(PROC_IN_ARRAY)=0.0_rlg

  ALLOCATE( raphbx(PROC_IN_ARRAY_m1p1), raphby(PROC_IN_ARRAY_m1p1) )
  ALLOCATE( tauskin_x(PROC_IN_ARRAY), tauskin_y(PROC_IN_ARRAY) )
  ALLOCATE( dry_cell(PROC_IN_ARRAY))
  raphbx(PROC_IN_ARRAY_m1p1)=0.0_rsh
  raphby(PROC_IN_ARRAY_m1p1)=0.0_rsh
  dry_cell(PROC_IN_ARRAY)=0

  ALLOCATE(flx_s2w_corim1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_s2w_corip1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_s2w_corjm1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_s2w_corjp1(-1:nv_adv,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_w2s_corin(nvp,PROC_IN_ARRAY))
  ALLOCATE(flx_w2s_corim1(nvp,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_w2s_corip1(nvp,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_w2s_corjm1(nvp,PROC_IN_ARRAY_m1p1))
  ALLOCATE(flx_w2s_corjp1(nvp,PROC_IN_ARRAY_m1p1))
#if ! defined key_nofluxwat_IWS
  ALLOCATE(phieau_s2w_corim1(PROC_IN_ARRAY_m1p1))
  ALLOCATE(phieau_s2w_corip1(PROC_IN_ARRAY_m1p1))
  ALLOCATE(phieau_s2w_corjm1(PROC_IN_ARRAY_m1p1))
  ALLOCATE(phieau_s2w_corjp1(PROC_IN_ARRAY_m1p1))
#endif
  
! Initialization
  ALLOCATE(cini_sed(nv_state))
  
     DO iv=1,nv_adv
       cini_sed(iv)=cini_sed_r(irk_fil(iv))*unit_modif_mudbio_N2dw(irk_fil(iv))
     ENDDO   

#ifdef key_Pconstitonly_insed 
     nv_use=nvpc
#else
     nv_use=nvp
#endif  

#if ! defined key_noTSdiss_insed
! counting of dissolved variables for diffusion in the sediment
  ivdiss(:)=0
  ivdiss(-1)=-1
  ivdiss(0)=0
#if ! defined key_Pconstitonly_insed
  DO iv=1,nv_adv-nvp
     ivdiss(iv)=iv+nvp
  ENDDO  
#endif
#endif

!  option morpho
!!!!!!!!!!!!!!!!!!!!
  IF(l_morphocoupl) THEN
!! Warning : no hx, hy in other model than MARS 
      ALLOCATE(morpho0(ARRAY_morpho))
      ALLOCATE(h0_bedrock(ARRAY_h0_bedrock))
      ALLOCATE(hsed0(PROC_IN_ARRAY))
      ALLOCATE(hsed_previous(PROC_IN_ARRAY))
#if defined MORPHODYN_MUSTANG_byHYDRO
      ALLOCATE(dhsed_save(PROC_IN_ARRAY))
#endif
      hsed0(PROC_IN_ARRAY)=0.0_rsh
      hsed_previous(PROC_IN_ARRAY)=0.0_rsh

  ENDIF

 !! declaration of the MUSTANG variables needed in the hydro model 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(EROS_FLUX_s2w(ARRAY_EROS_FLUX_s2w))
  ALLOCATE(SETTL_FLUX_w2s(ARRAY_SETTL_FLUX_w2s))
  ALLOCATE(SETTL_FLUXSUM_w2s(ARRAY_SETTL_FLUXSUM_w2s))
  ALLOCATE(WATER_FLUX_INPUTS(ARRAY_WATER_FLUX_INPUTS))
  EROS_FLUX_s2w(ARRAY_EROS_FLUX_s2w)=0.0_rsh
  SETTL_FLUX_w2s(ARRAY_SETTL_FLUX_w2s)=0.0_rsh
  SETTL_FLUXSUM_w2s(ARRAY_SETTL_FLUXSUM_w2s)=0.0_rsh
  WATER_FLUX_INPUTS(ARRAY_WATER_FLUX_INPUTS)=0.0_rsh
  
  ALLOCATE(fwet(PROC_IN_ARRAY))
  fwet(:,:)=1.0_rsh
  
#ifdef key_MUSTANG_flocmod
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! module FLOCULATION  !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 ! floc characteristics
 ALLOCATE(f_diam(1:nv_mud))     ! floc diameter
 ALLOCATE(f_vol(1:nv_mud))      ! floc volume
 ALLOCATE(f_rho(1:nv_mud))      ! floc density
 ALLOCATE(f_mass(0:nv_mud+1))     ! floc mass
 
 ! mass concentration 
 !ALLOCATE(f_cv(1:nv_mud))       ! extracted from cv_wat(:,k,j,j) mass concentration for every mud variables

 ! agregation kernels
 ALLOCATE(f_coll_prob_sh(1:nv_mud,1:nv_mud)) !  shear agregation collision probability
 ALLOCATE(f_coll_prob_ds(1:nv_mud,1:nv_mud)) ! differential settling collision probability
 
 ALLOCATE(f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud)) ! shear agregation gain term
 ALLOCATE(f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud)) ! differential settling agregation gain term
 ALLOCATE(f_l1_sh(1:nv_mud,1:nv_mud)) ! shear agregation loss term
 ALLOCATE(f_l1_ds(1:nv_mud,1:nv_mud)) ! differential settling agregation loss term  
 ALLOCATE(f_g3(1:nv_mud,1:nv_mud)) ! fragmentation gain term     
 ALLOCATE(f_l3(1:nv_mud)) ! fragmentation loss term
   
 ALLOCATE(f_davg(1:NB_LAYER_WAT,PROC_IN_ARRAY))
 ALLOCATE(f_d50(1:NB_LAYER_WAT,PROC_IN_ARRAY))
 ALLOCATE(f_d90(1:NB_LAYER_WAT,PROC_IN_ARRAY))
 ALLOCATE(f_d10(1:NB_LAYER_WAT,PROC_IN_ARRAY))
 ALLOCATE(f_dtmin(1:NB_LAYER_WAT,PROC_IN_ARRAY))
 
 f_diam(1:nv_mud)=0.0_rsh
 f_vol(1:nv_mud)=0.0_rsh
 f_rho(1:nv_mud)=0.0_rsh
 f_mass(0:nv_mud+1)=0.0_rsh
 
 !f_cv(1:nv_mud)=0.0_rsh
 
 f_coll_prob_sh(1:nv_mud,1:nv_mud)=0.0_rsh
 f_coll_prob_ds(1:nv_mud,1:nv_mud)=0.0_rsh
 
 f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
 f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
 f_l1_sh(1:nv_mud,1:nv_mud)=0.0_rsh
 f_l1_ds(1:nv_mud,1:nv_mud)=0.0_rsh
 f_g3(1:nv_mud,1:nv_mud)=0.0_rsh
 f_l3(1:nv_mud)=0.0_rsh

 f_davg(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
 f_d50(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
 f_d10(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
 f_d90(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh
 f_dtmin(1:NB_LAYER_WAT,PROC_IN_ARRAY)=0.0_rsh  

#endif
  
  PRINT_DBG*, 'END MUSTANG_ALLOC'
  
  END SUBROUTINE MUSTANG_alloc



#ifdef key_MUSTANG_flocmod
   !!===========================================================================

   SUBROUTINE flocmod_init(TRANSPORT_TIME_STEP)
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE flocmod_init  ***
   !&E
   !&E ** Purpose : initalization of flocs characteristics
   !&E
   !&E ** Description : 
   !&E
   !&E ** Note : GRAVITY must be known as a parameters transmtted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E ** Called by :  MUSTANG_initialization
   !&E 
   !&E ** External calls : 
   !&E
   !&E ** Note : NUMBER_PI must be known as a parameters transmtted by coupleur 
   !&E           in MARS : coupleur_dimhydro.h (USE ..)
   !&E           in CROCO : module_MUSTANG.F (include..)
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   !! * Arguments
   REAL(KIND=rlg),INTENT(IN)          :: TRANSPORT_TIME_STEP 

   !! * Local declarations
                     
   !!--------------------------------------------------------------------------
   !! * Executable part
   !!--------------------------------------------------
   !! floc characteristics
   
    f_diam(1:nv_mud)=diam_sed(imud1:nvpc) ! floc size from variable.dat
    f_vol(1:nv_mud)=NUMBER_PI/6._rsh*(f_diam(1:nv_mud))**3
    f_rho(1:nv_mud)=RHOREF+(ros(1:nv_mud)-RHOREF)*(f_dp0/f_diam(1:nv_mud))**(3._rsh-f_nf)
    f_mass(1:nv_mud)=f_vol(1:nv_mud)*(f_rho(1:nv_mud)-RHOREF)
    f_mass(nv_mud+1)=f_mass(nv_mud)*2_rsh+1._rsh  
    IF (f_diam(1).eq.f_dp0)  THEN
      f_mass(1)=f_vol(1)*ros(1)
    ENDIF
    f_ws(1:nv_mud)=GRAVITY*(f_rho(1:nv_mud)-RHOREF)*f_diam(1:nv_mud)**2._rsh/(18._rsh*0.001_rsh)
    f_dtmin(1:NB_LAYER_WAT,PROC_IN_ARRAY)=REAL(TRANSPORT_TIME_STEP,rsh)

    ! kernels computation
    CALL flocmod_kernels()    
 
    PRINT_DBG*, 'END flocmod_init'

  END SUBROUTINE flocmod_init
  
   !!===========================================================================
   
  SUBROUTINE flocmod_kernels()

  !&E--------------------------------------------------------------------------
  !&E                 ***  flocmod_kernels  ***
  !&E
  !&E ** Purpose : computations of agregation/fragmentation kernels for FLOCMOD
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_init
  !&E
  !&E ** External calls : 
  !&E
  !&E--------------------------------------------------------------------------

  !! * Local declarations

   LOGICAL        ::  f_test
   REAL(KIND=rsh) ::  f_weight,mult
   INTEGER        ::  iv1,iv2,iv3
   
  !!--------------------------------------------------------------------------
  !! * Executable part

  f_test=.true.

  ! compute collision probability
  
  CALL flocmod_agregation_statistics

  !********************************************************************************
  ! agregation : GAIN : f_g1_sh and f_g1_ds
  !********************************************************************************
  
   DO iv1=1,nv_mud-1
    DO iv2=1,nv_mud
      DO iv3=iv2,nv_mud
        IF((f_mass(iv2)+f_mass(iv3)) .GT. f_mass(iv1-1) &
           .AND. ((f_mass(iv2)+f_mass(iv3)) .LE. f_mass(iv1))) THEN

           f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

        ELSE IF ((f_mass(iv2)+f_mass(iv3)) .GT. f_mass(iv1) &
          .AND. ((f_mass(iv2)+f_mass(iv3)) .LT. f_mass(iv1+1))) THEN
          
           f_weight=1._rsh-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

       ELSE
           f_weight=0.0_rsh
       ENDIF

       f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)
       f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)

      ENDDO
    ENDDO
   ENDDO
   iv1=nv_mud
    DO iv2=1,nv_mud
      DO iv3=iv2,nv_mud
        IF((f_mass(iv2)+f_mass(iv3)) .GT. f_mass(iv1-1) &
           .AND. ((f_mass(iv2)+f_mass(iv3)) .LE. f_mass(iv1))) THEN

           f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

        ELSE IF ((f_mass(iv2)+f_mass(iv3)) .GT. f_mass(iv1) &
          .AND. ((f_mass(iv2)+f_mass(iv3)) .LT. f_mass(iv1+1))) THEN

            f_weight=1._rsh

       ELSE
           f_weight=0.0_rsh
       ENDIF

       f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)
       f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)

      ENDDO
    ENDDO
  

  
  !********************************************************************************
  ! Shear fragmentation : GAIN : f_g3
  !********************************************************************************
  
   DO iv1=1,nv_mud
    DO iv2=iv1,nv_mud
    
      IF (f_diam(iv2) > f_dmin_frag) THEN
      ! binary fragmentation
    
        IF (f_mass(iv2)/f_nb_frag .GT. f_mass(iv1-1) &
           .AND. f_mass(iv2)/f_nb_frag .LE. f_mass(iv1)) THEN

          IF (iv1 == 1) THEN 
            f_weight=1._rsh
          ELSE
            f_weight=(f_mass(iv2)/f_nb_frag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
          ENDIF

        ELSEIF (f_mass(iv2)/f_nb_frag .GT. f_mass(iv1) &
                .AND. f_mass(iv2)/f_nb_frag .LT. f_mass(iv1+1)) THEN

          f_weight=1._rsh-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

        ELSE
      
          f_weight=0._rsh
      
        ENDIF
    
      ELSE
        f_weight=0.0_rsh
      ENDIF
    
      f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1._rsh-f_ero_frac)*(1._rsh-f_ater)*f_weight*f_beta &
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3._rsh-f_nf)           &
                  *f_mass(iv2)/f_mass(iv1)
 
   ! ternary fragmentation
      IF (f_diam(iv2) .GT. f_dmin_frag) THEN
        IF (f_mass(iv2)/(2._rsh*f_nb_frag) .GT. f_mass(iv1-1) &
         .AND. f_mass(iv2)/(2._rsh*f_nb_frag) .LE. f_mass(iv1)) THEN

          IF (iv1 == 1) THEN 
            f_weight=1._rsh
          ELSE
            f_weight=(f_mass(iv2)/(2._rsh*f_nb_frag)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
          ENDIF

        ELSE IF (f_mass(iv2)/(2._rsh*f_nb_frag) .GT. f_mass(iv1) &
          .and. f_mass(iv2)/(2._rsh*f_nb_frag) .LT. f_mass(iv1+1)) THEN

          f_weight=1._rsh-(f_mass(iv2)/(2._rsh*f_nb_frag)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))

        ELSE
          f_weight=0._rsh
      
        ENDIF
    ! update for ternary fragments
          f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1._rsh-f_ero_frac)*(f_ater)*f_weight*f_beta &
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3._rsh-f_nf)           &
                  *f_mass(iv2)/f_mass(iv1)   
 
   ! Floc erosion
   
        IF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .GT. f_mass(f_ero_iv)) THEN

          IF (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .GT. f_mass(iv1-1)) &
          .AND. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .LE. f_mass(iv1)) THEN
       
            IF (iv1 == 1) THEN
              f_weight=1._rsh
            ELSE
              f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
            ENDIF
 
          ELSE IF ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .GT. f_mass(iv1) &
          .AND. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .LT. f_mass(iv1+1)) THEN
       
            f_weight=1._rsh-(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
 
          ELSE
            f_weight=0._rsh
          ENDIF
     
      ! update for eroded floc masses 
      
          f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta                    &
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3._rsh-f_nf)           &
!		  *f_mass(iv2)/f_mass(iv1)                                          &
                  *(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag)/f_mass(iv1)
    
          IF (iv1 == f_ero_iv) THEN
    
            f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta                           &
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)**(3._rsh-f_nf)           &
!		  *f_mass(iv2)/f_mass(iv1)                                          &
                  *f_ero_nbfrag*f_mass(f_ero_iv)/f_mass(iv1)
          ENDIF
        ENDIF
      ENDIF ! condition on f_dmin_frag
    ENDDO
  ENDDO  
 
!********************************************************************************
!  Shear agregation : LOSS : f_l1
!********************************************************************************
 
  DO iv1=1,nv_mud
    DO iv2=1,nv_mud
    
      IF(iv2 == iv1) THEN
        mult=2._rsh
      ELSE
        mult=1._rsh
      ENDIF
      
      f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1) 
      f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1) 
       
       
    ENDDO
  ENDDO
 
!********************************************************************************
!  Shear fragmentation : LOSS : f_l2
!********************************************************************************
 
 
  DO iv1=1,nv_mud
    IF (f_diam(iv1) > f_dmin_frag) THEN
 ! shear fragmentation
      f_l3(iv1)=f_l3(iv1)+(1._rsh-f_ero_frac)*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3._rsh-f_nf)

 ! shear erosion
      IF ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv)) THEN
        f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3._rsh-f_nf)
      ENDIF
    ENDIF    
  ENDDO

  PRINT_DBG*, 'END FLOCMOD_KERNELS'

  END SUBROUTINE flocmod_kernels

!!===========================================================================

  SUBROUTINE flocmod_agregation_statistics

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_agregation_statistics  ***
  !&E
  !&E ** Purpose : computation of shear / differential settling statistics
  !&E
  !&E ** Description :
  !&E
  !&E ** Note : NUMBER_PI must be known as a parameters transmtted by coupleur 
  !&E           in MARS : coupleur_dimhydro.h (USE ..)
  !&E           in CROCO : module_MUSTANG.F (include..)
  !&E
  !&E ** Called by : flocmod_kernels
  !&E
  !&E ** External calls : 
  !&E
  !&E--------------------------------------------------------------------------

  !! * Local declarations
  INTEGER      :: iv1,iv2
  !REAL(KIND=rsh), PARAMETER :: mu=0.001
  !!--------------------------------------------------------------------------
  !! * Executable part


  DO iv1=1,nv_mud
    DO iv2=1,nv_mud

      f_coll_prob_sh(iv1,iv2)=1._rsh/6._rsh*(f_diam(iv1)+f_diam(iv2))**3._rsh
      
      !f_coll_prob_ds(iv1,iv2)=0.25_rsh*NUMBER_PI*(f_diam(iv1)+f_diam(iv2))**2._rsh &
                       !*GRAVITY/mu*abs((f_rho(iv1)-RHOREF)*f_diam(iv1)**2._rsh &
                       !-(f_rho(iv2)-RHOREF)*f_diam(iv2)**2._rsh)
      f_coll_prob_ds(iv1,iv2)=0.25_rsh*NUMBER_PI*(f_diam(iv1)+f_diam(iv2))**2._rsh &
                       * ABS(f_ws(iv1)-f_ws(iv2))
      
    ENDDO
  ENDDO
  
  PRINT_DBG*, 'END flocmod_agregation_statistics'

  END SUBROUTINE flocmod_agregation_statistics
 
 !!===========================================================================
#endif /* key_MUSTANG_flocmod */

#endif /* MUSTANG */

END MODULE initMUSTANG
