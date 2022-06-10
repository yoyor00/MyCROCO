#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE coupleur_BIOLink_helping
!
!---------------------------------------------------------------------------
   

#if defined BIOLink

   !&E======================================================================
   !&E                   ***  MODULE  BIOLink_helping  ***
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

#include "coupleur_define_BIOLink.h" 
                                     ! Equivalence of names between the hydro and the different
                                     ! Tracer and biological model. Also contains a function
                                     ! For the allocation of the main tables

   USE module_BIOLink ! script that groups all the files of BIOLink together and allows the 
                      ! access to their functions/subroutines/variables
   USE comsubstance   ! Access to the functions/variables of SUBSTANCE ( from the MUSTANG
                      ! sediment model)

   USE comBIOLink     ! Common variables of the BIOLink coupler

   USE comBIOLink_helping ! Common variables of the BIOLink coupler
                          ! for helping functions

   USE comBIOLink_physics ! Common variables of the BIOLink coupler
                          ! For physics functions 

#if defined MUSTANG

   USE comMUSTANG ,  ONLY : htot ! Height of the water column

#endif /*MUSTANG*/

#if defined ECO3M
   USE mod_eco3m, ONLY : THICKLAYERWC
# endif

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
#if defined key_messat
    PUBLIC BIOLink_SPMsat_file    ! Determines the concentration of 
                                  ! suspended matter by satellite 
#endif /* key_messat */

#if defined key_nosubstmodule
    PUBLIC BIOLink_substance      ! Handles some water settling
                                  ! features of tracers if MUSTANG
                                  ! is not avaialble
#endif /* BIOLink_substance */


CONTAINS

  SUBROUTINE BIOLink_read_vardiag

   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_read_vardiag  ***
   !&E
   !&E ** Purpose : Reading of the file describing the diagnostic variables
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : bloom_init_id,bloom_create_var_diagtracer 
   !&E                     from bloom_initdefine
   !&E                     nb_var_tracerN,nb_var_tracerP from parameters
   !&E                    
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2019-08 (B. Thouvenin) issued from sub_read_vardiag 
   !&E       !  2022-02 (G. Koenig) commenting 
   !&E
   !&E---------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================

#if defined BLOOM

   USE bloom_initdefine, ONLY : bloom_init_id

#  if defined key_N_tracer

   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer

#  endif /* key_N_tracer */

#  if defined key_P_tracer

   USE bloom_initdefine, ONLY : bloom_create_vardiagtracer

#  endif /* key_P_tracer */

#endif /* BLOOM */

     !====================================================================
     ! Local declarations of variables
     !====================================================================

   LOGICAL               :: ex,l_diag_wat,l_diag_sed
   INTEGER               :: eof,isubs,isubs_r,dimvar,it,ind_white,IERR_MPI
   CHARACTER(LEN=lchain) :: namvar_r,long_name_var_r,standard_name_var_r,unitvar_r
   CHARACTER(LEN=5)      :: comment

#if defined key_N_tracer || defined key_P_tracer

   INTEGER               :: is,id,ivtra

#endif /* key_N_tracer || key_P_tracer */

     !====================================================================
     ! Execution of the function
     !====================================================================

      !********************** Saving into simu.log *************************!

   IF_MPI (MASTER) THEN
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
     MPI_master_only WRITE(iscreenlog,*) '**************** VAR_READ_DIAG.F90 ***************'
     MPI_master_only WRITE(iscreenlog,*) '**************************************************'
     MPI_master_only WRITE(iscreenlog,*) ' '
     MPI_master_only WRITE(iscreenlog,*) 'file defining usefull diagnostic variables : ',TRIM(filevardiag)

   ENDIF_MPI

      ! Initialize number of diagnostic variables according to their dimensions
      ndiag_1d = 0
      ndiag_2d = 0
      ndiag_3d = 0
      ndiag_3d_wat = 0
      ndiag_3d_sed = 0
      ndiag_2d_sed = 0
      ndiag_tot = 0
      isubs = 0
      l_out_subs_diag=.false.  ! no saving of diagnoses variable by default

#if defined MUSTANG

      l_out_subs_diag_sed=.false.  ! no saving of diagnoses variable by default

#endif /* MUSTANG */

      ! read total number of diagnostic variables

      eof = 0

      INQUIRE(file=filevardiag,exist=ex)

      IF (ex) THEN
        
        OPEN(49,file = filevardiag,form='formatted')
     
        comment='debut'
       
        DO

          READ(49,'(a)',iostat=eof) comment

          IF (comment.EQ.'*****') EXIT

        END DO

        DO WHILE(eof==0)

          READ(49,*,iostat=eof) ! variable number
          READ(49,'(a)',iostat=eof) ! variable name
          READ(49,'(a)',iostat=eof) !variable long_name
          READ(49,'(a)',iostat=eof) !variable standard_name
          READ(49,'(a)',iostat=eof) !unit
          READ(49,*,iostat=eof)     !variable valid_min value
          READ(49,*,iostat=eof)     !variable valid_max value
          READ(49,*,iostat=eof)     !dimension (1=1D ; 2=2D ; 3=3D)
          READ(49,*,iostat=eof)     !variable in water    (k,i,j)
          READ(49,*,iostat=eof)     !variable in sediment (i,j,k)
          READ(49,*,iostat=eof)     !saving in file

          IF (eof==0) ndiag_tot=ndiag_tot+1
            READ(49,'(a)',iostat=eof)

        END DO

     CLOSE(49)

#if defined BLOOM

#  if defined key_N_tracer

       ndiag_tracerN=0

       DO is=1,nb_source_tracerN

            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1

#    if defined key_age_tracer

            ndiag_tot=ndiag_tot+nb_var_tracerN+1
            ndiag_tracerN=ndiag_tracerN+nb_var_tracerN+1

#    endif /* key_age_tracer */

       END DO ! is

#  endif /* key_N_tracer */

#  if defined key_P_tracer

       ndiag_tracerP=0

       DO is=1,nb_source_tracerP

            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1

#    if defined key_age_tracer

            ndiag_tot=ndiag_tot+nb_var_tracerP+1
            ndiag_tracerP=ndiag_tracerP+nb_var_tracerP+1

#    endif /* key_age_tracer */

       END DO ! is

#  endif /* key_P_tracer */

#endif /* BLOOM */

     IF_MPI (MASTER) THEN

        MPI_master_only WRITE(iscreenlog,*) 'Number of diagnostic variables = ',ndiag_tot

     ENDIF_MPI

   ! allocate arrays for diagnostic variables
     ALLOCATE( idimv_r(ndiag_tot) )
     ALLOCATE( l_diagBIOLink_out(ndiag_tot) )
     ALLOCATE( name_vardiag(ndiag_tot) )
     ALLOCATE( long_name_vardiag(ndiag_tot) )
     ALLOCATE( standard_name_vardiag(ndiag_tot) )
     ALLOCATE( unit_vardiag(ndiag_tot) )
     ALLOCATE( valid_min_vardiag(ndiag_tot) )
     ALLOCATE( valid_max_vardiag(ndiag_tot) )
     ALLOCATE( irk_diag(ndiag_tot) )


   ! read diagnostic variables and order them according to their matrix dimensions
     OPEN(49,file = filevardiag,form='formatted')

     comment='debut'

     DO WHILE (comment /= '*****')

       READ(49,'(a)',iostat=eof) comment

     END DO

     DO WHILE(eof==0)

       READ(49,*,iostat=eof) isubs_r

       isubs=isubs+1

       IF (eof==0) THEN

         READ(49,'(a)',iostat=eof) namvar_r

         ind_white=INDEX(namvar_r,' ')
         name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(namvar_r(1:ind_white))))

         READ(49,'(a)',iostat=eof) long_name_var_r
         ind_white=INDEX(long_name_var_r,' ') 
         long_name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(long_name_var_r(1:ind_white))))

         READ(49,'(a)',iostat=eof) standard_name_var_r
         ind_white=INDEX(standard_name_var_r,' ')
         standard_name_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(standard_name_var_r(1:ind_white))))

         READ(49,'(a)',iostat=eof) unitvar_r
         ind_white=INDEX(unitvar_r,' ')
         unit_vardiag(isubs)=TRIM(ADJUSTL(ADJUSTR(unitvar_r(1:ind_white))))

         READ(49,*,iostat=eof) valid_min_vardiag(isubs)
         READ(49,*,iostat=eof) valid_max_vardiag(isubs)
         READ(49,*,iostat=eof) dimvar

         IF (dimvar==1) THEN

           idimv_r(isubs)=1
           ndiag_1d=ndiag_1d+1

           vname(1,indxBLMdiag1D+ndiag_1d) = namvar_r
           vname(2,indxBLMdiag1D+ndiag_1d) = long_name_var_r
           vname(3,indxBLMdiag1d+ndiag_1d) = standard_name_var_r
           vname(4,indxBLMdiag1D+ndiag_1d) = unitvar_r

         ELSE IF (dimvar==2) THEN

           idimv_r(isubs)=2
           ndiag_2d=ndiag_2d+1

           vname(1,indxBLMdiag2D+ndiag_2d) = namvar_r
           vname(2,indxBLMdiag2D+ndiag_2d) = long_name_var_r
           vname(3,indxBLMdiag2d+ndiag_2d) = standard_name_var_r
           vname(4,indxBLMdiag2D+ndiag_2d) = unitvar_r

         ELSE IF (dimvar==3) THEN

           idimv_r(isubs)=3
           ndiag_3d=ndiag_3d+1

           vname(1,indxBLMdiag3D+ndiag_3d) = namvar_r
           vname(2,indxBLMdiag3D+ndiag_3d) = long_name_var_r
           vname(3,indxBLMdiag3d+ndiag_3d) = standard_name_var_r
           vname(4,indxBLMdiag3D+ndiag_3d) = unitvar_r


         END IF

         READ(49,*,iostat=eof) l_diag_wat

         IF (l_diag_wat .and. dimvar==3) ndiag_3d_wat=ndiag_3d_wat+1

         READ(49,*,iostat=eof) l_diag_sed

         IF (l_diag_sed .and. dimvar==2) THEN

           ndiag_2d_sed=ndiag_2d_sed+1
           idimv_r(isubs)=6

         END IF

         IF (l_diag_sed .and. dimvar==3) THEN

           ndiag_3d_sed=ndiag_3d_sed+1

           IF (l_diag_wat) THEN

             idimv_r(isubs)=4

           ELSE

             idimv_r(isubs)=5

           ENDIF

         END IF

         READ(49,*,iostat=eof) l_diagBIOLink_out(isubs)
         READ(49,*,iostat=eof)

         IF (l_diag_sed) THEN

#if ! defined MUSTANG

          IF_MPI (MASTER) THEN

           MPI_master_only WRITE(iscreenlog,*)' '
           MPI_master_only WRITE(iscreenlog,*)' WARNING : diagnostic variable in sediment '
           MPI_master_only WRITE(iscreenlog,*)'           without CPP key MUSTANG '
           MPI_master_only WRITE(iscreenlog,*)' DIAG. VAR. NAME : ',TRIM(name_vardiag(isubs))

          ENDIF_MPI
#else

          IF (l_diagBIOLink_out(isubs)) THEN

             l_out_subs_diag_sed=.true.

          ENDIF

#endif /* MUSTANG */

         ELSE

           IF (l_diagBIOLink_out(isubs))  l_out_subs_diag=.true.

         ENDIF
         
         IF_MPI (MASTER) THEN

           MPI_master_only WRITE(iscreenlog,*)' '
           MPI_master_only WRITE(iscreenlog,*)' DIAG. VAR. NAME : ',TRIM(name_vardiag(isubs))
           MPI_master_only WRITE(iscreenlog,*)' UNIT            : ',TRIM(unit_vardiag(isubs))
           MPI_master_only WRITE(iscreenlog,*)' DIMENSION       : ', dimvar,'D'
           MPI_master_only WRITE(iscreenlog,*)' DIAG SAVE ?     : ', l_diagBIOLink_out(isubs)

           IF (l_diag_wat)THEN

              MPI_master_only WRITE(iscreenlog,*)' DIAGNOSTIC INTO THE SEA'

           ENDIF

           IF (l_diag_sed) THEN

              MPI_master_only WRITE(iscreenlog,*)' DIAGNOSTIC INTO THE SEDIMENT'

           ENDIF

         ENDIF_MPI

       END IF

     END DO

     CLOSE(49)

   ELSE

     MPI_master_only WRITE(iscreenlog,*)' filevardiag', TRIM(filevardiag), ' not found  !! '

   END IF

#if defined BLOOM && (defined key_N_tracer || defined key_P_tracer)

       call bloom_create_vardiagtracer

#endif /* BLOOM && ( key_N_tracer || key_P_tracer ) */

   IF_MPI (MASTER) THEN

     MPI_master_only WRITE(iscreenlog,*)' '
     MPI_master_only WRITE(iscreenlog,*)' STOCK OF DIAGNOSTIC VARIABLES'
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables total :',ndiag_tot
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 1 dim :',ndiag_1d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 2 dim :',ndiag_2d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3 dim :',ndiag_3d
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3d wat:',ndiag_3d_wat
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 2d sed:',ndiag_2d_sed
     MPI_master_only WRITE(iscreenlog,*)' number of diagnostic variables 3d sed:',ndiag_3d_sed

   ENDIF_MPI

   IF (ndiag_tot /= ndiag_1d+ndiag_2d+ndiag_3d) THEN

      MPI_master_only PRINT*,'WARNING number of diagnostic variables is incoherent.'
      MPI_master_only PRINT*,'Check in file :',filevardiag

   END IF

#if ! defined MUSTANG

   IF (ndiag_3d_sed /= 0 .OR. ndiag_2d_sed /= 0 ) THEN

      MPI_master_only PRINT*,'WARNING no MUSTANG, you should not have any diagnostic variable for sediment'
      MPI_master_only PRINT*,'simulation stopped'
      CALL_MPI MPI_FINALIZE(IERR_MPI)
      STOP

   END IF
#endif /* MUSTANG */

   ! allocate diagnostic variables
   ALLOCATE( diag_1d(1:ndiag_1d) )
   diag_1d(:)=0.0_rsh

   ALLOCATE( diag_2d(ndiag_1d+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   diag_2d(:,:,:)=0.0_rsh

   ALLOCATE( diag_3d_wat(ndiag_2d+1:ndiag_2d+ndiag_3d_wat,NB_LAYER_WAT,PROC_IN_ARRAY) )
   diag_3d_wat(:,:,:,:)=0.0_rsh

   ALLOCATE( diag_3d_CROCO(ndiag_2d+1:ndiag_2d+ndiag_3d_wat,PROC_IN_ARRAY,NB_LAYER_WAT) )
   diag_3d_CROCO(:,:,:,:)=0.0_rsh

#if defined MUSTANG && defined key_BLOOM_insed

   ALLOCATE( diag_3D_sed(ndiag_tot-ndiag_3d_sed+1:ndiag_tot,ksdmin:ksdmax,PROC_IN_ARRAY) ) 
   diag_3d_sed(:,:,:,:)=0.0_rsh

   ALLOCATE( diag_2D_sed(ndiag_1d+ndiag_2d-ndiag_2d_sed+1:ndiag_1d+ndiag_2d,PROC_IN_ARRAY) )
   diag_2D_sed(:,:,:)=0.0_rsh

#endif /* MUSTANG && key_BLOOM_insed */

   ! Storage of diagnostic variables within reading order

#if defined BLOOM

   it = 0
   ! dimvar=1 : diag1D, =2 diag2D
   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only
   ! dimvar=6 : diag2D in sed only

#  if defined key_BLOOM_insed
   ! dimvar=1 : diag1D, =2 diag2D

   DO dimvar = 1,2

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
!         irk_diag(it)=isubs
         irk_diag(isubs)=it

#    if defined key_N_tracer

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_tracer */

       END IF

     END DO

   END DO

   ! dimvar=6 : diag2D in sed only
   dimvar=6

   DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else
         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_P_tracer && BLOOM */

       END IF

   END DO

   ! dimvar=3 : diag3D in wat only
   ! dimvar=4 : diag 3D in wat and in sed
   ! dimvar=5 : diag3D in sed only

   DO dimvar = 3,5

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1

         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif

       END IF

     END DO

   END DO

#  else

   DO dimvar = 1,6

     DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

#    if defined key_N_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerN .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    elif defined key_P_tracer && defined BLOOM

         IF(isubs <= ndiag_tot-ndiag_tracerP .or. isubs==ndiag_tot) CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    else

         CALL bloom_init_id(isubs,standard_name_vardiag(isubs))

#    endif /* key_N_P_tracer && BLOOM */

       END IF

     END DO

   END DO
#  endif /* BLOOM */

#endif /* key_BLOOM_insed */

#if defined PEPTIC

   it = 0

   DO dimvar = 1,4

    DO isubs = 1,ndiag_tot

       IF (idimv_r(isubs) == dimvar) THEN

         it=it+1
         irk_diag(isubs)=it

       END IF

    END DO

   END DO

#endif /* PEPTIC */

   IF_MPI (MASTER) THEN

     DO isubs = 1,ndiag_tot

        MPI_master_only WRITE(iscreenlog,*) isubs,TRIM(name_vardiag(isubs)),idimv_r(isubs),irk_diag(isubs)

     END DO

   ENDIF_MPI

  END SUBROUTINE BIOLink_read_vardiag


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

      DO i=ifirst,ilast ! Zonal loop

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
#if !defined ECO3M
                   WS_BIOLink(k,iv,i,j)=sign(MIN(0.95_rlg*THICKLAYERWC(k,i,j)/TRANSPORT_TIME_STEP,REAL(ABS(WS_BIOLink(k,iv,i,j)),rlg)),WS_BIOLink(k,iv,i,j))
#else
! THICKLAYERWC is a (i,j,k) array 
                   WS_BIOLink(k,iv,i,j)=sign(MIN(0.95_rlg*THICKLAYERWC(i,j,k)/TRANSPORT_TIME_STEP,REAL(ABS(WS_BIOLink(k,iv,i,j)),rlg)),WS_BIOLink(k,iv,i,j))
#endif
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

      !******************* Extinction due to water ************************!


   DO j=jfirst,jlast

     DO i=ifirst,ilast

        kmaxmod=NB_LAYER_WAT ! We only compute at boundary layers
                             ! where MUSTANG is not applied 

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

        DO i=ifirst,ilast

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

        DO i=ifirst,ilast

           PAR_top_layer_day(:,i,j)=PAR_top_layer_day(:,i,j) + PAR_top_layer(:,i,j)*BIO_TIME_STEP

        END DO ! i

      END DO ! j

!$OMP END DO

    ENDIF

#  endif /* PEPTIC */
END SUBROUTINE  BIOLink_eval_PAR

#endif /* BIOLink_PAR_eval */

#if defined key_messat
     SUBROUTINE BIOLink_SPMsat_file(ifirst,ilast,jfirst,jlast,icall,forcSPM)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE BIOLink_SPMsat_file  ***
   !&E
   !&E ** Purpose : special reading SPM satellite file for key_messat
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : BIOLink_update
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History : !Creation by B. Thouvenin 
   !&E              ! Commenting by G. Koenig ( February 2022)
   !&E
   !&E--------------------------------------------------------------------------

     !====================================================================
     ! Routines from external models
     !====================================================================


   USE ionc4,       ONLY : ionc4_openr,ionc4_read_dimt,ionc4_read_time, &
                           ionc4_read_xyt,ionc4_close

     !====================================================================
     ! External arguments
     !====================================================================

   INTEGER, INTENT(IN)  :: ifirst,ilast,jfirst,jlast ! Limits of the 
                                                     ! horizontal domain
   INTEGER, INTENT(IN)  :: icall
   REAL(KIND=rsh), DIMENSION(PROC_IN_ARRAY),INTENT(OUT),OPTIONAL  :: forcSPM

     !====================================================================
     ! Local declarations of variables
     !====================================================================

   INTEGER    :: k,i,j,iv,it ! Spatial, tracer and temporal counters
   INTEGER    :: ijour,imois,ian,iheure,iminu,isec
   INTEGER    :: jjulien_init,tool_julien ! Date and time tool
   CHARACTER(LEN=19) :: date_start,tool_sectodat
   REAL(kind=rlg)    :: tbid
   REAL(kind=riosh), DIMENSION(COMPLETE_ARRAY) :: tab_mes ! Table of 
                                                          ! suspended matter
   INTEGER, DIMENSION(51)        :: t_clim_messat_int ! Weekly climatology
                                                      ! of satellite 
                                                      ! suspended matter
   CHARACTER(LEN=lchain)         :: nomfic ! Name of the climatology file
   INTEGER                       :: IERR_MPI ! Error flag for MPI
   REAL(kind=rlg)                :: torigin,tool_datosec

     !====================================================================
     ! Execution of the function
     !====================================================================

      !******************* Opening of the climatology file *****************!

   IF (icall == 0) THEN
     !! first call - initialization - open and read files
        
    DO it=1,idimt_messat

       CALL ionc4_read_xyt(filemessat_clim,'MES',tab_mes,imin,imax,jmin,jmax,it)

       DO j=jfirst,jlast

         DO i=ifirst,ilast

          IF(h0(i,j).GT.-valmanq) THEN

            IF (tab_mes(i,j).LT.valmanq) THEN

              messat(i,j,it)=tab_mes(i,j)

              IF(messat(i,j,it) .LT. 0.0_rsh) THEN

                MPI_master_only WRITE(iscreenlog,*)'pb messat maille ',i,j,it,messat(i,j,it)

              END IF

            ELSE

                   MPI_master_only WRITE(iscreenlog,*)'pas de MES climato en ',i,j,it,tab_mes(i,j)

            END IF

          END IF

        END DO

      END DO

      MPI_master_only WRITE(iscreenlog,*) t_clim_messat(it)

    END DO

    CALL ionc4_close(filemessat_clim)

#  if defined key_agrif

     IF (.NOT. l_initfromfile) THEN

       date_start=tool_sectodat(tdebagrif_messat)

       CALL tool_decompdate(date_start,ijour,imois,ian,iheure,iminu,isec)

       jjulien_init=tool_julien(ijour,imois,ian)-tool_julien(1,1,ian)+1

       date_s_annee_messat=jjulien_init*24*60*60

     ELSE

       date_s_annee_messat=jjulien_BIOLINK*24*60*60

     END IF

#  else

     date_s_annee_messat=jjulien_BIOLINK*24*60*60

#  endif /* key_agrif */

     it=1
     idateinf_messat=1
     t_clim_messat_int(1:idimt_messat)=INT(t_clim_messat(1:idimt_messat))

     IF (date_s_annee_messat .LE. t_clim_messat_int(1)) THEN

        idatesup_messat=1
        idateinf_messat=idimt_messat

     ELSE

        DO WHILE ((it .LT. idimt_messat) .AND.      &
             (t_clim_messat_int(it+1) .LT. date_s_annee_messat))

          it=it+1
          idateinf_messat=it

        END DO

        IF (idateinf_messat .EQ. idimt_messat) THEN

          idatesup_messat=1

        ELSE

          idatesup_messat=idateinf_messat+1

        END IF

     END IF

     IF_MPI (MASTER) THEN

       MPI_master_only WRITE(iscreenlog,*) 'Fin de lecture fichier climato satellite MES'
       MPI_master_only WRITE(iscreenlog,*) it,' champs MES lus'
       MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
       MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',date_s_annee_messat
       MPI_master_only WRITE(iscreenlog,*)'----------------------------------------------------'
       MPI_master_only WRITE(iscreenlog,*)''

     ENDIF_MPI

   END IF 

   IF (l_messat_obs) THEN

     IF_MPI (MASTER) THEN

       MPI_master_only WRITE(iscreenlog,*)'------------------------------------------------'
       MPI_master_only WRITE(iscreenlog,*)'lecture fichier obs moyennes quizaines pour mes-sat  '
       MPI_master_only WRITE(iscreenlog,*)'     donnees en g.l-1', TRIM(filemessat_obs)

     ENDIF_MPI

     CALL ionc4_openr(filemessat_obs)
     idimt_messat_obs=ionc4_read_dimt(filemessat_obs)
     ALLOCATE(messat_obs(limin:limax,ljmin:ljmax,idimt_messat_obs))
     ALLOCATE(t_obs_messat(idimt_messat_obs))

    DO it=1,idimt_messat_obs

      CALL ionc4_read_time(filemessat_obs,it,t_obs_messat(it))

    END DO

    t_obs_messat(:)=t_obs_messat(:)-torigin

    MPI_master_only WRITE(*,*)'tobd_messat debut et fin=',t_obs_messat(1),t_obs_messat(idimt_messat_obs)

    DO it=1,idimt_messat_obs

      CALL ionc4_read_xyt(filemessat_obs,'MES_SAT',tab_mes,imin,imax,jmin,jmax,it)
      DO j=jfirst,jlast

       DO i=ifirst,ilast

         IF(h0(i,j).GT.-valmanq) THEN

           IF (tab_mes(i,j).LT.valmanq) THEN

             messat_obs(i,j,it)=tab_mes(i,j)

             IF(messat_obs(i,j,it) .LT. 0.0_rsh) THEN

               MPI_master_only WRITE(iscreenlog,*)'pb messat_obs maille ',i,j,it,messat_obs(i,j,it)
 
             END IF

           ELSE

             MPI_master_only WRITE(iscreenlog,*)'pas de MES climato en ',i,j,it,tab_mes(i,j)

           END IF

         END IF

       END DO

     END DO

       MPI_master_only WRITE(iscreenlog,*) t_obs_messat(it)

    END DO

    CALL ionc4_close(filemessat_obs)

    it=1
    idateinf_messat=1
    tbid=CURRENT_TIME

    IF (l_initfromfile .and. l_init_rtime) then

       CALL ionc4_openr(file_init)
       CALL ionc4_read_time(file_init,1,tbid)
       CALL ionc4_close(file_init)

    END IF

    date_start=tool_sectodat(tbid)
    MPI_master_only WRITE(iscreenlog,*)'date_start=',date_start

    IF (tbid .LE. t_obs_messat(1)) THEN

      IF_MPI (MASTER) THEN

        MPI_master_only WRITE(ierrorlog,*)'la date de depart est anterieure a la premiere date du fichier messat'
        CALL_MPI MPI_FINALIZE(IERR_MPI)
        STOP

       ENDIF_MPI

    ELSE

      DO WHILE ((it .LT. idimt_messat_obs) .AND.(t_obs_messat(it+1) .LT. tbid))

          it=it+1
          idateinf_messat=it

      END DO

      IF (idateinf_messat .EQ. idimt_messat_obs) THEN

        IF_MPI (MASTER) THEN

          MPI_master_only WRITE(ierrorlog,*)'la date de depart est posterieure a la derniere date du fichier messat'    
          CALL_MPI MPI_FINALIZE(IERR_MPI)
          STOP

        ENDIF_MPI

      ELSE

         idatesup_messat=idateinf_messat+1

      END IF

    END IF
   
    IF_MPI (MASTER) THEN

      MPI_master_only WRITE(iscreenlog,*) 'Fin de lecture fichier observation quinzaine satellite MES'
      MPI_master_only WRITE(iscreenlog,*) it,' champs MES lus'
      MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
      MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',tbid
      MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', t_obs_messat(idateinf_messat),' et ',t_obs_messat(idatesup_messat)
      MPI_master_only WRITE(iscreenlog,*)'----------------------------------------------------'
      MPI_master_only WRITE(iscreenlog,*)''
      MPI_master_only WRITE(iscreenlog,*) 'date du run  = ',tbid
      MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', idateinf_messat,' et ',idatesup_messat
      MPI_master_only WRITE(iscreenlog,*) ' champs MES interp entre', t_obs_messat(idateinf_messat),' et ',t_obs_messat(idatesup_messat)

    ENDIF_MPI

  END IF !endif sur l_messat_obs

      !******************* Evaluation of spm at time t *********************!

  ELSE IF (icall == 1) THEN

      !***************** Interpolation factor computation ******************!

    IF (l_messat_clim) THEN

      IF (date_s_annee_messat .NE. jjulien_BIOLINK*24*60*60) THEN ! Definit des dates inf et sup pour l interp
 
        date_s_annee_messat= jjulien_BIOLINK*24*60*60              ! de la climato

        IF (idateinf_messat .LT. idatesup_messat) THEN

          IF (date_s_annee_messat .GT. INT(t_clim_messat(idatesup_messat))) THEN

            idateinf_messat=idatesup_messat

            IF (idatesup_messat .NE. idimt_messat) THEN

              idatesup_messat=idatesup_messat+1

            ELSE

              idatesup_messat=1

            END IF

          END IF

        ELSE

          IF ((date_s_annee_messat .GT. INT(t_clim_messat(idatesup_messat))) .AND.   &
               (date_s_annee_messat .LT. INT(t_clim_messat(idateinf_messat)))) THEN

             idateinf_messat=1
             idatesup_messat=idatesup_messat+1
          END IF

        END IF

     END IF !endif on date_s_annee_messat

     IF (idatesup_messat .GT. idateinf_messat)  THEN

         interp= (date_s_annee_messat -t_clim_messat(idateinf_messat))/        &
                 (t_clim_messat(idatesup_messat)-t_clim_messat(idateinf_messat))

     ELSE

        IF (date_s_annee_messat .GE. INT(t_clim_messat(idateinf_messat))) THEN

           interp= (date_s_annee_messat-t_clim_messat(idateinf_messat))/&
                  ((366*24*60*60+t_clim_messat(idatesup_messat))-t_clim_messat(idateinf_messat))

        ELSE

           interp=(date_s_annee_messat+366*24*60*60-t_clim_messat(idateinf_messat))/ &
                  ((366*24*60*60+t_clim_messat(idatesup_messat))-t_clim_messat(idateinf_messat))

        END IF

     END IF

   END IF  !endif on l_messat_clim

   IF (l_messat_obs) THEN

     IF (t .GT. t_obs_messat(idatesup_messat)) THEN

       idateinf_messat=idatesup_messat

       IF (idatesup_messat .NE. idimt_messat_obs) THEN

         idatesup_messat=idatesup_messat+1

       ELSE

         MPI_master_only WRITE(*,*)'derniere date du fichier messat depassee'

         STOP

       END IF

     END IF

     IF (idatesup_messat .GT. idateinf_messat)  THEN

        interp= (t -t_obs_messat(idateinf_messat))/        &
                 (t_obs_messat(idatesup_messat)-t_obs_messat(idateinf_messat))

     END IF
       
   END IF

!$OMP DO SCHEDULE(RUNTIME)

     DO j=jfirst,jlast

       DO i=ifirst,ilast

         IF (TOTAL_WATER_HEIGHT .LE. RESIDUAL_THICKNESS_WAT) THEN

            forcSPM(i,j)=valmanq

         ELSE

           IF (l_messat_obs) THEN

               forcSPM(i,j)=messat_obs(i,j,idateinf_messat)+                                &
               (messat_obs(i,j,idatesup_messat)-messat_obs(i,j,idateinf_messat))*interp

           END IF !endifon l_messat_obs

           IF (l_messat_clim) THEN

             forcSPM(i,j)=messat(i,j,idateinf_messat)+                                &
                 (messat(i,j,idatesup_messat)-messat(i,j,idateinf_messat))*interp

           ENDIF  !endif on l_messat_clim

         ENDIF  !endif on htot

        END DO

      END DO
!$OMP END DO
  
  ENDIF ! icall

   END SUBROUTINE BIOLink_SPMsat_file

#endif /* key_messat */

#endif /* BIOLink */

END MODULE coupleur_BIOLink_helping

