#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE coupleur_BIOLink_helping_functions
!
!---------------------------------------------------------------------------
   

#if defined BIOLink

   !&E======================================================================
   !&E                   ***  MODULE  BIOLink_helping_functions  ***
   !&E
   !&E ** Purpose : Adding some functions in BIOLink to complete
   !&E              missing features in hydrodynamics or tracer models
   !&E            
   !&E
   !&E ** Description :
   !&E     subroutine BIOLink_read_vardiag     ! reading diagnostics variables - 
   !&E                                          called by BIOLink_initialization
   !&E
   !&E     subroutine BIOLink_sinking_rate     ! update and limiting sinking rate 
   !&E                                           for each variables - called by 
   !&E                                           BIOLink_update
   !&E
   !&E     subroutine BIOLink_eval_PAR         ! evaluation of solar radiation 
   !&E                                          extinction, attenuation and PAR - 
   !&E                                          called by BIOLink_update
   !&E
   !&E     subroutine BIOLink_substance        ! definition number et type of 
   !&E                                           variables (if key_nosubstmodule)
   !&E
   !&E     subroutine BIOLink_SPMsat_file       ! special reading SPM satellite file 
   !&E                                               for key_messat
   !&E
   !&E   History :
   !&E    !  2022-02 (G. Koenig) - Creation from coupleur_BIOLink.F90
   !&E
   !&E=========================================================================================

!*****************************************************************************************!
!                                                                                         !
!                                      Interface                                          !
!                                                                                         !
!*****************************************************************************************!

#include "coupleur_define_BIOLink.h" ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

#ifdef key_MARS
#include "coupleur_dimhydro_BIOLink.h"
   USE sflxatm,      ONLY : rad
#else
   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)
#endif /* key_MARS */

   USE comBIOLink     ! Common variables of the BIOLink coupler


   IMPLICIT NONE
  
!*****************************************************************************************!
!                                                                                         !
!                                  Accessibility                                          !
!                                                                                         !
!
!                                                                                         !
!*****************************************************************************************!

     !====================================================================
     ! List of functions and routines
     !====================================================================

    PUBLIC BIOLink_read_vardiag   ! Reading of diagnostic variables

    PUBLIC BIOLink_sinking_rate   ! Handles the sinking rate in the absence
                                  ! of MUSTANG

    PUBLIC BIOLink_eval_PAR       ! Evaluates the PAR in the case where no 
                                  ! PAR modules are available

    PUBLIC BIOLink_substance      ! Handles some water settling
                                  ! features of tracers if MUSTANG
                                  ! is not avaialble

    PUBLIC BIOLink_SPMsat_file    ! Determines the concentration of 
                                  ! suspended matter by satellite 

#if defined key_nosubstmodule
    PUBLIC BIOLink_substance      ! Declaration of the characteristics of the tracers
                                  ! If the module substance is not declared
#endif /* BIOLink_substance */


CONTAINS

   SUBROUTINE BIOLink_sinking_rate(ifirst,ilast,jfirst,jlast)

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_sinking_rate  ***
   !&E
   !&E ** Biologic dynamics:  - limiting sinking rate for each variables
   !&E
   !&E ** Purpose :  Updating sinking rates
   !&E              
   !&E ** Description :  It first computes the settling velocity from 
   !&E                   the provided files. If not, it computes it
   !&E                   from the ratio of the cell thickness and BIOLink 
   !&E                   time step. In the case we use PEPTIC, it computes
   !&E                   it from PEPTIC 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : None
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from peptic_dynzwat - verti_quota 
   !&E                  (A. Menesguen, M. Sourrisseau)
   !&E       !  2022-02 (G. Koenig) commenting and editing
   !&E
   !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined key_MARS

   USE comvars2d,    ONLY : ig,id,hm
   USE obccombine, ONLY : l_obc_south, l_obc_north, l_obc_west, l_obc_east

#endif /* key_MARS */

     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN) :: ifirst,ilast,jfirst,jlast ! Limits of the MPI
                                                    ! subdomain
   
     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER                :: i,j,k,iv ! Counter variables

#if defined PEPTIC

   INTEGER                :: i_plkt,i_quota_loc,ind

#endif /* PEPTIC */

     !====================================================================
     ! Execution of the function
     !====================================================================

      !**************** Computation of sinking velocity ********************!
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,iv)

   DO j=jfirst,jlast ! Meridional loop

#if defined key_MARS

      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)

#else

      DO i=ifirst,ilast ! Zonal loop

#endif /* key_MARS */

        IF (TOTAL_WATER_HEIGHT(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN

#if ! defined MUSTANG

           DO iv=1,nvp

             DO k = NB_LAYER_WAT,1,-1

              WS_BIOLink(k,iv,i,j)=(ws_free_min(iv)+ws_free_max(iv))/2.0_rsh

             END DO ! k

           END DO ! iv
#endif /* MUSTANG */
 
#if defined PEPTIC
           ! I do not understand what happens
           DO i_plkt = 1 , bd_fp%nb_plct

             DO i_quota_loc = 1 , plct(i_plkt)%nb_quota

               IF (plct(i_plkt)%is_quota_var(i_quota_loc)== 1) THEN

                 iv=plct(i_plkt)%num_mod_mars(i_quota_loc) 

                 DO k = NB_LAYER_WAT,1,-1

                   WS_BIOLink(k,iv,i,j) = plct(i_plkt)%sink_max  !m.d-1

                 END DO ! k

               END IF ! Mysterious mystery

             END DO ! i_quota_loc

           END DO ! i_plkt

#endif /* PEPTIC */

           DO iv=nvpc+1,nvp

                 DO k = NB_LAYER_WAT,1,-1

                   WS_BIOLink(k,iv,i,j)=sign(MIN(0.95_rlg*THICKLAYERWC(k,i,j)/TRANSPORT_TIME_STEP,REAL(ABS(WS_BIOLink(k,iv,i,j)),rlg)),WS_BIOLink(k,iv,i,j))
                 
                 END DO ! k

           END DO  !iv

        ENDIF ! TOTAL WATER HEIGHT GT RESIDUAL WATER THICKNESS

      END DO  !i

    END DO  !j
!$OMP END DO
         

  END SUBROUTINE BIOLink_sinking_rate

#if defined key_nosubstmodule

  SUBROUTINE BIOLink_substance(icall)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE coupl_BIOLink_dim  ***
   !&E
   !&E ** Purpose : Remplacement of substance module for the settlement
   !&E              and interaction with the sediment
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : None
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2016-11  (B. Thouvenin) 
   !&E       !  2022-02  (G. Koenig)
   !&E--------------------------------------------------------------------------

     !====================================================================
     ! External arguments
     !====================================================================
 
   INTEGER,INTENT(IN)  :: icall ! if = 0 we just count the number of 
                                ! variables. If not we initialize the array

     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER    :: iv ! Tracer counter

     !====================================================================
     ! Execution of the function
     !====================================================================

   IF(icall==0) THEN


     ! definition of the number of simulated variables (substances)  

     !==========================================================================
     ! definition here ??? 
     ! or read a data file or a namelist file 
     !==========================================================================

     ! number of particulate variables  type gravel 
      nv_grav=0
     ! number de particulate variables type sand 
      nv_sand=1
     ! number de particulate variables type mud 
      nv_mud=2
     ! number de particulate variables type sorbed on another partculate variable 
      nv_sorb=0
     ! number of dissolved  variables 
      nv_dis=1


      ! initialize the number of variables according to their type
      nvpc=nv_mud+nv_sand+nv_grav ! Constitutive particulate variables
      nvp=nvpc+nv_ncp+nv_sorb ! Total particulate variables
      nv_adv=nvp+nv_dis ! Advected variables
      nv_state=nv_adv ! Variables whose quantity changes
      nv_tot=nv_state ! Total number of variables
    
    ELSE
    ! if icall = 1

      ALLOCATE(name_var(1,nv_tot)) ! Name array
      ALLOCATE(typart(nv_state)) ! Particulate variables array
      ALLOCATE(typdiss(nv_adv)) ! Dissolved variables array
      ALLOCATE(itypv(nv_tot) ! Total variables array
      ALLOCATE(irk_fil(nv_tot)) ! Order in the module array

      IF (nvp > 0) THEN

        ALLOCATE(ws_free_opt(nvp)) ! Settling
        ws_free_opt(:)=0 ! velocity 

        ALLOCATE(ws_hind_opt(nvp)) ! Settling velocity
        ws_hind_opt(:)=0 ! of last time step

        ALLOCATE(ws_free_para(4,nvp)) ! Parameters of 
        ws_free_para(:,:)=0.0_rsh ! water settling

        ALLOCATE(ws_free_min(nvp)) ! Minimal water
                                   ! settling velocity
        ALLOCATE(ws_free_max(nvp)) ! Maximal water 
        ws_free_max(:)=0.0_rsh     ! settling velocity

        ALLOCATE(ws_hind_para(2,nvp)) ! Parameters of the
        ws_hind_para(:,:)=0.0_rsh ! the water settling velocity
                                  ! at the last time step
        ALLOCATE(tocd(nvp)) ! I do not know
        tocd(:)=0.0_rsh

      ENDIF

      ALLOCATE(l_subs2D(nv_adv))
      ALLOCATE(ws3(ARRAY_WAT_SETTL)) ! Water settling 3d array
      ALLOCATE(cv_wat(ARRAY_WATER_CONC)) ! Concentration of tracer in water

! and info on  variables
    ! variable names (name_var), type (itypv), settling velocities  (ws_free..), 
    ! itypv=1 for gravels ; =2 pour sand ; =3 pour mud, 
    ! itypv=4 for non constitutives particulate variables, 
    ! itypv=5 for sorbed variables  
    ! itypv=6 for dissoolved
    ! for particulate, give diameter, tocd, ros 
    ! for muds and non constit particulates, give parameters to define settling velocities 
    ! initial concentrations  in water and  sediments
    ! for sand : treatment in 2D or 3D (l_subs2D)
    ! for sorbed particulate, give the indices or numbers of the constitutives 
    !                         part. variables on which they are sorbed (irkm_var_assoc)
    ! identification of igrav1,igrav2,isand1,isand2,imud1,imud2
    ! ---------------------------------------------------------
      igrav1=1
      igrav2=nv_grav
      isand1=igrav2+1
      isand2=igrav2+nv_sand
      imud1=isand2+1
      imud2=isand2+nv_mud 

      typart(1:nv_state)=0.
      typart(1:nvpc)=1.
      typdiss(nvp+1:nv_adv)=1.
      typdiss(1:nvp)=0.

      l_subs2D(:)=.false.

     ! we arrange here to give the variables in the right order for MUSTANG 
     ! so that they are well ranked in the table of concentrations.
     ! so :
      DO iv=1,nv_tot
        irk_fil(iv)=iv
        unit_modif_mudbio_N2dw(iv)=1.0_rsh
      ENDDO
  

   END SUBROUTINE BIOLink_substance

#endif /* key_nosubstmodule */


#if defined BIOLink_PAR_eval

  SUBROUTINE BIOLink_eval_PAR(ifirst,ilast,jfirst,jlast,cdate)

  !&E---------------------------------------------------------------------
  !&E                 ***  ROUTINE BIOLink_eval_PAR ***
  !&E
  !&E ** Purpose : Evaluation of solar radiation extinction, attenuation and PAR 
  !&E
  !&E ** Description : 
  !&E
  !&E ** Called by : BIOLink_update
  !&E
  !&E ** External calls : None
  !&E
  !&E ** Reference :
  !&E
  !&E ** History : B. Thouvenin ( date unknown), creation from verti_quota
  !&E              G. Koenig ( February 2022), commenting
  !&E
  !&E---------------------------------------------------------------------

     !====================================================================
     ! External arguments
     !====================================================================
   
   INTEGER, INTENT(IN)       :: ifirst,ilast,jfirst,jlast ! Limits of MPI
                                                          ! subdomains
   CHARACTER(LEN=19),INTENT(IN)  :: cdate                 ! Date

     !====================================================================
     ! Local declarations of variables
     !====================================================================

    INTEGER                  ::  i,j,k,iv,kmaxmod        ! loop indexes

#  if defined PEPTIC

    INTEGER                  ::  i_plkt,i_pom,ind,num_cell,i_quota_loc ! loop
                                                                ! indexes
    REAL(KIND=rsh),DIMENSION(bd_fp%nb_plct) :: conc ! Concentration ?

#  endif /* PEPTIC */

    REAL(KIND=rsh),DIMENSION(NB_LAYER_WAT)     :: attenuation ! Attenuation
                                                              ! of light

   INTEGER                  :: iday,jhour,numday_lum ! Variables for
   INTEGER                  :: numhour_lum,klum      ! averaging light
   INTEGER                  :: numday_extinction     ! and extinction over
   INTEGER                  :: numhour_extinction    ! 24 hours

     !====================================================================
     ! Execution of the function
     !====================================================================


!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k,kmaxmod,attenuation)

      !******************* Extinction du to water ************************!

   DO j=jfirst,jlast

#  if defined key_MARS

     DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)

       IF(TOTAL_WATER_HEIGHT(i,j) < hm ) THEN

               kmaxmod=1

       ELSE

               kmaxmod=NB_LAYER_WAT

       ENDIF ! TOTAL_WATER_HEIGHT

#  else

     DO i=ifirst,ilast

        kmaxmod=NB_LAYER_WAT ! We only compute at boundary layers
                             ! where MUSTANG is not applied 

#  endif /* key_MARS */

       PAR_top_layer(:,i,j)=0.0_rsh ! Initialization of tables
       EXTINCTION_RAD(:,i,j)=0.0_rsh ! for the extinction and the 
       attenuation(:)=0.0_rsh ! PAR 

#  if defined METeOR

       Flimrad_layer(:,i,j)=0.0_rsh

#  endif /* METeOR */

       IF (TOTAL_WATER_HEIGHT(i,j) .GT. RESIDUAL_THICKNESS_WAT) THEN

         DO k = LOOPK_SURF_TO_BOTTOM_WAT   ! kmaxmod,1,-1

           EXTINCTION_RAD(k,i,j) = PARAM_WAT_EXTINCT

#  if defined PEPTIC

           EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j)         & 
                                   + bd_fp%extincspim * cmes_3dmgl(k,i,j)

#  elif defined BLOOM || (defined METeOR && ! defined PEPTIC)

           EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j)         &
                                   + p_extincspim                & 
                                   *(cmes_3dmgl(k,i,j)+epsilon_BIOLink)
#  endif /* BLOOM || (METeOR && PEPTIC ) */

      !***************** Extinction due to chlorophyll *********************!

#  if defined PEPTIC

           IF (BIOLink_chloro(k,i,j) > 1.0e-10_rsh) THEN

#  endif /* PEPTIC */

#  if defined PEPTIC || defined BLOOM

             EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j)       &
                                     + PARAM_CHLORO1_EXTINCT     &
                                     * ( BIOLink_chloro(k,i,j)   & 
                                     ** PARAM_CHLORO2_EXTINCT )

#    if defined PEPTIC

           ENDIF ! BIOLink_chloro

#    endif /* PEPTIC */

#  endif /* PEPTIC/BLOOM */
             
      !****************** Extinction due to macroalgae *********************!

#  if defined BLOOM && defined key_zostera

           IF(k==1 .and. FIXCONCPOS(iv_zost_LB-nv_adv,1,i,j).gt.0.0_rsh) THEN

             EXTINCTION_RAD(k,i,j) = EXTINCTION_RAD(k,i,j)        &
                                   + p_zost_leafabscoef           &
                                   * FIXCONCPOS(iv_zost_LB-nv_adv,1,i,j) &
                                   * p_zost_klai

           ENDIF ! k==1 and FIXCONCPOS

#  endif /* BLOOM && key_zostera */
       
       !************* Estimation of attenuation at each cell ***************! 

           attenuation(k) = EXP( -EXTINCTION_RAD(k,i,j) * THICKLAYERWC(k,i,j))

         END DO !k

       !*********** Estimation of PAR at the top of each layer *************! 
          
#  if defined PEPTIC

#    if defined key_growth_diurne 

         PAR_top_layer(kmaxmod,i,j)=SOLAR_RAD(i,j)*bd_fp%parradratio / RAD_SRFSCALE

#    endif /* key_growth_diurne */

#  elif defined BLOOM

         PAR_top_layer(kmaxmod,i,j)=SOLAR_RAD(i,j)*p_parradratio / RAD_SRFSCALE 

#  elif defined METeOR

         Flimrad_layer(kmaxmod,i,j)=attenuation(kmaxmod)

#  endif /* PEPTIC/BLOOM/METeOR */

         DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT 

            PAR_top_layer(k,i,j) = PAR_top_layer(ABOVE_K,i,j) * attenuation(ABOVE_K)       

#  if defined METeOR

            Flimrad_layer(k,i,j)=Flimrad_layer(ABOVE_K,i,j)* attenuation(k)

#  endif /* METeOR */
                 
         END DO !k

         k=0 ! at bottom

         PAR_top_layer(k,i,j) = PAR_top_layer(ABOVE_K,i,j) * attenuation(ABOVE_K)

       !************* Estimation of PAR averaged in each layer *************! 

#  if defined PEPTIC

         DO k=LOOPK_SUBSURF_TO_BOTTOM_WAT   ! kmaxmod-1,1,-1   
            
            klum = k - 1

            IF (k == 1) klum = 1
              IF (EXTINCTION_RAD(klum,i,j) /= 0.0_rsh) THEN

#    if defined key_growth_diurne

               ! Conversion in microEinstein of the PAR fraction

               PAR_avg_layer_phyto(1,k,i,j) = PAR_top_layer(k,i,j) &
                                              / EXTINCTION_RAD(klum,i,j) &  
                                              *( 1.0_rsh -               &
                                              attenuation(klum))         &
                                              / THICKLAYERWC(k,i,j)      &
                                              * 4.6_rsh *86400.0_rsh 
#    endif /* key_growth_diurne */

            ENDIF ! EXTINCTION_RAD */

         END DO !k

#  endif /* PEPTIC */         
     
       ELSE

          PAR_top_layer(:,i,j)=0.0_rsh 

#  if defined PEPTIC

          PAR_avg_layer_phyto(1,:,i,j) = 0.0_rsh

#  endif /* PEPTIC */

       ENDIF ! d>RESIDUAL_THICKNESS_WAT

     END DO !i

   END DO !j
!$OMP END DO

       !****** Estimation of PAR averaged on 24 hours on the top layer *****! 

#  if defined PEPTIC

    !integration light / 24h ! Warning : Must restart at 00.00  and one day lag between realistic forcing and application

   IF ((iheure_BIOLINK == 0) .and. (iminu_BIOLINK == 0) .and. (isec_BIOLINK <= BIO_TIME_STEP/2._rsh)) then 

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k)

      DO j = jfirst,jlast

#    if defined key_MARS

        DO i = MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)

#    else 

        DO i=ifirst,ilast

#    endif /* key_MARS */

          IF ((i==i_BIOLink_verif) .and. (j==j_BIOLink_verif)) THEN

           MPI_master_only WRITE(iscreenlog,*) 'new daily_aver_PAR_W_m2',cdate,BIO_TIME_STEP,PAR_top_layer_day(:,i,j)/86400.0_rsh ! error?*bd_fp%parradratio !light_ave_daily(i,j) 

          ENDIF 

          PAR_top_layer_day(:,i,j)=0.0_rsh 

        END DO

      END DO

!$OMP END DO

   ELSE

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j,k)

      DO j = jfirst,jlast

#    if defined key_MARS

        DO i = MAX0(limin,ig(j)+1),MIN0(limax,id(j)-1)

#    else 

        DO i=ifirst,ilast

#    endif /* key_MARS */

           PAR_top_layer_day(:,i,j)=PAR_top_layer_day(:,i,j) + PAR_top_layer(:,i,j)*BIO_TIME_STEP

        END DO ! i

      END DO ! j

!$OMP END DO

    ENDIF

#  endif /* PEPTIC */


END SUBROUTINE  BIOLink_eval_PAR

#endif /* BIOLink_PAR_eval */

#endif /* BIOLink */

END MODULE coupleur_BIOLink_helping_functions

