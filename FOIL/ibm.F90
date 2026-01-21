MODULE ibm

!!======================================================================
!!                   ***  MODULE ibm   ***
!! 
!!         M. Huret 01/2011
!!         IBM (Individual Based model)
!!         To be adapted for each case for which a CPP key has to be defined
!&E
!&E ** History :
!&E       !  2011-10 (M. Huret) Introduction of MPI
!&E       !  2012-02 (M. Huret) Introduction of NAMELIST
!&E       !  2024    (M. Caillaud, M. Huret, D. Gourves) Coupled with CROCO   
!&E
!!======================================================================

#include "cppdefs.h"
#include "toolcpp.h"

#ifdef MPI
    use mpi
#endif

#if defined  DEB_IBM

    USE module_ibm
    USE comtraj,        ONLY : kmax,rsh,rlg,lchain,type_position,    &
                               rg_valmanq_io,dg_valmanq_io,valmanq,   &
                               type_particle,type_patch,patches
    USE comtraj,        ONLY : nb_species, duration


    IMPLICIT NONE
    PRIVATE

    !! * Accessibility
    PUBLIC :: ibm_3d
    PUBLIC :: ibm_init

    !! * Shared module variables

    !! * Private variables                  
    REAL(kind=rsh)                                  :: w_max,alpha_w                    ! From paraibm namibmbio namelist


    LOGICAL                                         :: repro                            ! From paraibm namibmpop, activate or not repro
    REAL(KIND=rlg)                                  :: dt_spawn, dt_save                ! From paraibm namibmpop, spawning interval in hours of patch
                                                                                        ! From paraibm namibmpop namemist, time between two savings of data        
    INTEGER                                         :: max_part                         ! From paraibm namibmpop, max nb of particles in a nc file

    LOGICAL, DIMENSION(nb_species)                  :: spawn, first_spawn, newseason    ! Logicals to manage spawn
    INTEGER                                         :: current_year, current_day        ! Integers to save current year and day for repro and age of fish
    INTEGER                                         :: yearclass                        ! Logical to update fish AgeClass
    INTEGER, PARAMETER                              :: out = 4                          ! Type of reals in IBM output

#ifdef IBM_SPECIES
    ! From paraibm, namibmpop namelist arguments
    LOGICAL                                         :: fish_mort, density_dependent
    REAL(KIND=rsh)                                  :: multiplier_tac

    ! Number of eggs and stade 1 for larvae mortality
    REAL(KIND=rsh), DIMENSION(nb_species)           :: number_eggs

    ! Variables pour la 2e methode de repro
    REAL(KIND=rsh),DIMENSION(nb_species)            :: nb_indv_ponte

    REAL(KIND=rlg)                                  :: slope     = 0.000001722786_rlg
    REAL(KIND=rsh)                                  :: alpha_sel = 0.876931491863804146725_rsh
    REAL(KIND=rsh)                                  :: beta_sel  = 11.725820182044399686561_rsh
#endif





 !!===================================================================================================================================
 !!===================================================================================================================================
 !!===================================================================================================================================
   
 CONTAINS


 SUBROUTINE ibm_init(xe,sal,temp,Istr,Iend,Jstr,Jend)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE ibm_init  ***
    !&E
    !&E ** Purpose : Initialize properties of particles with call to traj_init3d
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_init_main
    !&E ** External calls : LAGRANGIAN_init,ionc4_openr,ionc4_read_trajt,ionc4_read_traj
    !&E                     ionc4_close,ionc4_read_time,ionc4_read_dimt,deb_init,fish_move_init
    !&E                     ibm_parameter_init
    !&E ** Reference      : Huret et al. (2016)
    !&E
    !&E ** History :
    !&E       !  2010-11 (M. Huret) Original code
    !&E       !  2011-10 (M. Huret) Introduction of MPI
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves, M. Huret) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used

    ! Import subroutines
    USE trajinitsave,   ONLY : LAGRANGIAN_init
    USE ionc4,          ONLY : ionc4_openr, ionc4_read_trajt, ionc4_read_traj,   &
                               ionc4_close, ionc4_read_time, ionc4_read_dimt,    &
                               ionc4_gatt_char_read, ionc4_read_dimtraj
#ifdef IBM_SPECIES
    USE debmodel,       ONLY : deb_init
    USE ibmmove,        ONLY : fish_move_init
    USE ibmtools,       ONLY : ibm_parameter_init
#endif

    ! Import variables
    USE comtraj,        ONLY : iscreenlog
    USE comtraj,        ONLY : type_particle,type_patch,patches,file_trajec,                    &
                               file_pathout,itypetraj,ndtz,ibm_restart
#ifdef IBM_SPECIES
    USE comtraj,        ONLY : debuse,F_Fix,ffix,file_food,file_NBSS,frac_deb_death,            &
                               fileanchovy,filesardine,fileprobadistrib_anc,nbSizeClass_anc,    &
                               sizemin_anc,fileprobadistrib_sar,nbSizeClass_sar,sizemin_sar,    &
                               catch_anc_bob,catch_sar_bob,fishing_strategy
#endif

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY     ),        INTENT( in ) :: xe
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax),        INTENT( in ) :: sal,temp
    INTEGER,                                                INTENT( in ) :: Istr,Iend,Jstr,Jend    

    !! * Local declarations
    CHARACTER(LEN=lchain)                           :: file_inp                 ! Name of netcdf restart file with data
    INTEGER                                         :: lstr, lenstr             ! To read paraibm file
    INTEGER                                         :: idimt                    ! Read last time in restart file 
    INTEGER                                         :: num                      ! For restart loop to keep good num info

    LOGICAL                                         :: ibm_l_time               ! From paraibm in namibmrestart namelist, restart info

    INTEGER                                         :: nb_part_nc               ! Number of particles in netcdf for patch
    INTEGER                                         :: duration_ibm_anc,duration_ibm_sar ! Life time of anchovy and sardine

    INTEGER                                         :: i,j,n,m, il                 ! Integers for loops
    INTEGER                                         :: index_num                          ! Integers for indexing in restart
    ! To convert date to seconds or seconds to date
    CHARACTER(len=19)                               :: tool_sectodat 
    INTEGER                                         :: mm_clock,hh,minu,sec           ! jj and aaaa are saved as current_year/day for later

#ifdef IBM_SPECIES 
    ! From paraibm, spatial size for eggs merging in a superindividual

#endif     

    TYPE(type_patch),    POINTER                    :: patch                        ! Temporary shortcut for a patch in loops  
    TYPE(type_particle), POINTER                    :: particle                     ! Temporary shortcut for a particle in loops
    ! Temporary arrays to read data from netcdf if restart
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:)       :: flag_nc, temp_nc, super_nc   
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:)       :: dens_nc, size_nc, drate_nc
    REAL(KIND=rlg), ALLOCATABLE, DIMENSION(:)       :: dayb_nc
    INTEGER,        ALLOCATABLE, DIMENSION(:)       :: stage_nc, age_nc, ageClass_nc, num_nc 

    ! Definition of namelists in paraibm
    NAMELIST/namibmin/      file_trajec,file_pathout,itypetraj,ndtz
    NAMELIST/namibmrestart/ ibm_restart,ibm_l_time
    NAMELIST/namibmbio/     w_max, alpha_w
#ifdef IBM_SPECIES
    NAMELIST/namibmpop/     repro,dt_spawn,dt_save,max_part,duration_ibm_anc,duration_ibm_sar, &
                            fish_mort,fishing_strategy,multiplier_tac,density_dependent
    NAMELIST/namibmdeb/     debuse,F_Fix,ffix,file_NBSS,file_food,frac_deb_death
    NAMELIST/namibmfrc/     fileanchovy,filesardine,catch_anc_bob,catch_sar_bob,fileprobadistrib_anc,nbSizeClass_anc,sizemin_anc, &
                            fileprobadistrib_sar,nbSizeClass_sar,sizemin_sar
#endif


#include "compute_auxiliary_bounds.h"
    !!----------------------------------------------------------------------
    !! * Executable part

    ! namelists in paraibm.txt
    !--------------------------
    lstr = lenstr(debibmname)
    OPEN(50,file=debibmname(1:lstr),status='old',form='formatted',access='sequential')
    READ(50,namibmin)
    READ(50,namibmrestart)
    READ(50,namibmbio)
#ifdef IBM_SPECIES
    READ(50,namibmpop)
    READ(50,namibmdeb)
    READ(50,namibmfrc)
#endif
    CLOSE(50)

    CALL LAGRANGIAN_init(Istr,Iend,Jstr,Jend)

    ! save into simu.log
    !-------------------
    IF_MPI ( MASTER ) THEN
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) '***************************************************'
        WRITE(iscreenlog,*) '*****************   IBM_INIT.F90   ****************'
        WRITE(iscreenlog,*) '***************************************************'
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) 'fichier definissant les caracteristiques des trajectoires : ',trim(file_trajec)
    ENDIF_MPI      

    CALL tool_decompdate(tool_sectodat(time),current_day,mm_clock,current_year,hh,minu,sec)

    patch => patches % first
    DO n = 1, patches%nb      
        ! Init patch general data
        patch%yearref  = current_year - 1
        patch%t_spawn  = patch%t_beg ! clara : why ???
        patch%dt_spawn = dt_spawn*3600.0_rlg
        patch%dt_save  = patch%dt_save ! clara : interet de cette ligne ? sauf si dt_save tout court


        ! ------------------------- 
        ! --- Restart 
        IF ( ibm_restart ) THEN
            file_inp = trim(patch%file_inp)

            ! nb_part_nc = patch%nb_part_total ! denis
            CALL ionc4_openr(file_inp, .false.) ! clara
            CALL ionc4_read_dimtraj(file_inp, nb_part_nc) !clara

            ALLOCATE( flag_nc(nb_part_nc), temp_nc(nb_part_nc),  size_nc(nb_part_nc), stage_nc(nb_part_nc) )
            ALLOCATE( dens_nc(nb_part_nc), super_nc(nb_part_nc), drate_nc(nb_part_nc), dayb_nc(nb_part_nc) )
            ALLOCATE( age_nc(nb_part_nc), ageClass_nc(nb_part_nc), num_nc(nb_part_nc) )

            ! CALL ionc4_openr(trim(file_inp), .false.)
            CALL ionc4_gatt_char_read(file_inp, 'run_id', patch%run_id)
            ! Read time dimension in input file to open last time in restart file
            idimt = ionc4_read_dimt(file_inp)
            CALL ionc4_read_trajt(file_inp, "flag",      flag_nc,     1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "TEMP",      temp_nc,     1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "SIZE",      size_nc,     1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "DENSITY",   dens_nc,     1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "STAGE",     stage_nc,    1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "NUMBER",    super_nc,    1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "DRATE",     drate_nc,    1, nb_part_nc, idimt)
            CALL ionc4_read_traj (file_inp, "DAYBIRTH",  dayb_nc,     1, nb_part_nc)
            CALL ionc4_read_trajt(file_inp, "AGE",       age_nc,      1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "AGECLASS",  ageClass_nc, 1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "NUM",  num_nc,      1, nb_part_nc, idimt)

            DO m = 1,patch%nb_part_alloc
                IF (patch%nb_part_alloc == 0) CYCLE     ! To avoid an error because of a proc without any particle at restart
                IF ( .NOT. patch%particles(m)%active ) CYCLE 
                num = patch % particles(m) % num 
                ! CLARA, get index of num_nc
                ! index_num = findloc(num_nc, num, dim=1)
                index_num = -1
                do il = 1, nb_part_nc
                    if (num_nc(il) == num) then
                        index_num = il
                        exit
                    end if
                end do
                IF ( index_num == -1 ) THEN
                    PRINT *, 'ERROR not found num=', num
                    CYCLE
                END IF
                patch % particles(m) % flag      = flag_nc(index_num)
                patch % particles(m) % temp      = temp_nc(index_num)
                patch % particles(m) % size      = size_nc(index_num)
                patch % particles(m) % density   = dens_nc(index_num)
                patch % particles(m) % stage     = stage_nc(index_num)
                patch % particles(m) % super     = super_nc(index_num)
                patch % particles(m) % Drate     = drate_nc(index_num)
                patch % particles(m) % date_orig = dayb_nc(index_num)
                patch % particles(m) % age       = age_nc(index_num)
                patch % particles(m) % AgeClass  = ageClass_nc(index_num)
            END DO
            
#ifdef IBM_SPECIES
            CALL ionc4_read_traj(trim(file_inp), "DAYJUV",   dayb_nc, 1, nb_part_nc)
            CALL ionc4_read_traj(trim(file_inp), "DENSPAWN", dens_nc, 1, nb_part_nc)

            DO m = 1,patch%nb_part_alloc
                IF (patch%nb_part_alloc == 0) CYCLE
                IF ( .NOT. patch%particles(m)%active ) CYCLE

                num = patch%particles(m)%num          
                ! CLARA, get index of num_nc
                ! index_num = findloc(num_nc, num, dim=1)
                index_num = -1
                do il = 1, nb_part_nc
                    if (num_nc(il) == num) then
                        index_num = il
                        exit
                    end if
                end do
                IF ( index_num == -1 ) THEN
                    PRINT *, 'ERROR not found num=', num
                    CYCLE
                END IF
                patch%particles(m)%dayjuv   = dayb_nc(index_num)
                patch%particles(m)%denspawn = dens_nc(index_num)

                IF ( patch%particles(m)%stage >= 5 ) patch%particles(m)%itypevert = 0
            END DO

#endif /*IBM_SPECIES*/

            DEALLOCATE( flag_nc, temp_nc, size_nc, stage_nc, dens_nc, super_nc, drate_nc, dayb_nc )
            DEALLOCATE( age_nc, ageClass_nc, num_nc )
 
            ! update the date of restart, and savetraj is delayed not to have twice same time step in output
            IF ( ibm_l_time ) THEN
            !    CALL ionc4_read_time(trim(file_inp), 1, patch%t_beg)
                CALL ionc4_read_time(trim(file_inp), idimt, patch%t_beg)
                patch%t_save = patch%t_beg + patch%dt_save*3600.0_rlg
            END IF
            CALL ionc4_close(file_inp)
            
        ! ------------------------- 
        ! --- If not a restart
        ELSE     
            ! patch%particles(:)%date_orig = time
            nb_part_nc = patch%nb_part_alloc
            patch%run_id = generate_run_id()
            !--- Update FOIL.info file
            call write_run_info(patch%run_id)

            !Huret et al. 2016 (egg density =f(surf.density))
            DO m = 1, nb_part_nc
                particle => patch%particles(m)
                particle%date_orig = time - particle%age*24.0_rlg*3600.0_rlg
                IF (patch%nb_part_alloc == 0) CYCLE ! To avoid an error because of a proc without any particle at init
                IF ( .NOT. particle%active )  CYCLE
#ifdef IBM_SPECIES
                ! Init some biological parameters if not a restart
                CALL ibm_parameter_init(particle,patch%species,temp,sal,xe,Istr,Iend,Jstr,Jend)
#endif /* IBM_SPECIES */
            ENDDO      ! loop on patch%nb_part_alloc

        ENDIF  ! restart or not
        
        patch => patch%next
    ENDDO         ! loop on patches%nb
    
#ifdef IBM_SPECIES
    duration = (/duration_ibm_anc, duration_ibm_sar /) ! Store in one variable life expectancy for both species
    
    ! No need of loop to initialize DEB parameters
    IF (debuse) CALL deb_init(ibm_restart)      ! Init DEB  
    CALL fish_move_init(Istr,Iend,Jstr,Jend)    ! Init fish_move module

    yearclass = current_year + 1                ! Init yearclass to update fish's Ageclass

    ! Init reproduction by creating the egg-laying subgrid
    IF ( repro ) THEN

        ! Init repro booleans  
        newseason   = .true.
        first_spawn = .false.
        spawn       = .false.

        ! Contrainte du nombre de particules a chaque generation selon max_part

        ! Calcul nb particules a creer a chaque evenement de ponte
        nb_indv_ponte(1) = NINT(max_part*(dt_spawn/24.d0)/(138.d0)) ! 138 = nb de jours de ponte pendant une annee pour l'anchois
        nb_indv_ponte(2) = NINT(max_part*(dt_spawn/24.d0)/(211.d0)) ! 211 = nb de jours de ponte pendant une annee pour la sardine
    END IF
#endif
    
    ! save initialization
    CALL ibm_save

 END SUBROUTINE ibm_init





 !!======================================================================
 SUBROUTINE ibm_3d(xe,uz,vz,sal,temp,Istr,Iend,Jstr,Jend)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE ibm_3d  ***
    !&E
    !&E ** Purpose : Main routine with call to advection, behaviour, growth...
    !&E **          To be adapted by the user depending on case study
    !&E
    !&E 
    !&E ** Called by      : ibm_update_main
    !&E ** External calls : LAGRANGIAN_update,init_patch,indices_loc2glob,tool_latlon2i,tool_latlon2j
    !&E                     loc_h0,define_pos,ztosiggen,h0int,xeint,hc_sigint
    !&E                     ex_traj,ADD_ALL_MPI_INT,init_mpi_type_particle
    !&E                     ibm_loc_xyz,ibm_buoy,ibm_traint,ibm_proftraint,selec_dome_or_asymp,tool_julien
    !&E                     ibm_nycth_mig,death_by_fishing,eggs_grid,ibm_parameter_init,fish_move,deb_egg_init
    !&E                     deb_cycle,readtemp3d(unused anymore),readfood3d
    !&E
    !&E ** Reference      : Menu et al. (2023), Bueno-Pardo et al. (2020), Gatti et al. (2017), Huret et al. (2010), 
    !&E                     Boussouar et al. (2001), , Regner (1996), Zweifel&Lasker (1976), Rose el al. (2015)
    !&E                     Boyra (2013), Ospina (2012), Somarakis (2007)
    !&E
    !&E ** History :
    !&E       !  2008-12 (Marc Sourisseau)
    !&E       !  2010-11 (M. Huret)
    !&E       !  2011-10 (M. Huret) Introduction of MPI
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2023    (C. Menu) Densite-dependance, mortalities
    !&E       !  2024    (M. Caillaud, D. Gourves, M. Huret) Coupled with CROCO
    !&E       !  2024    (M. Huret, D. Gourves) Adding density-dependance and mortality by fishing, spatialisation of reproduction
    !&E---------------------------------------------------------------------
    !! * Modules used
#ifdef PASSIVE_TRACERS
    !USE parameters,   ONLY : nb_var
#endif /* PASSIVE_TRACERS */
    USE traject3d,    ONLY : LAGRANGIAN_update
    USE trajinitsave, ONLY : init_patch, indices_loc2glob
    USE trajectools,  ONLY : tool_latlon2i, tool_latlon2j
    USE trajectools,  ONLY : loc_h0, define_pos, ztosiggen, h0int, xeint, hc_sigint
#ifdef MPI
    USE toolmpi,      ONLY : ex_traj, ADD_ALL_MPI_INT
    USE comtraj,      ONLY : init_mpi_type_particle
    USE comtraj,      ONLY : down_give, up_give, right_give, left_give
#endif
    USE ibmtools,     ONLY : ibm_loc_xyz, ibm_buoy, ibm_traint, ibm_proftraint
    USE ibmtools,     ONLY : tool_julien
#ifdef IBM_SPECIES
    USE ibmtools,     ONLY : ibm_nycth_mig, death_by_fishing
    USE ibmtools,     ONLY : ibm_parameter_init
    USE ibmmove,      ONLY : fish_move
    USE debmodel,     ONLY : deb_egg_init, deb_cycle
    USE debmodel,     ONLY : readtemp3d,readfood3d
#endif /* IBM_SPECIES */
    USE comtraj,      ONLY : type_particle, type_patch, patches, patch_list_append, resize_patch
#ifdef IBM_SPECIES
    USE comtraj,      ONLY : file_pathout
    USE comtraj,      ONLY : jjulien, struc_ad, struc_ad_dd_DEB
    USE comtraj,      ONLY : debuse, F_Fix
    USE comtraj,      ONLY : number_tot, weight_tot, biom_tot, Wdeb_mean
    USE comtraj,      ONLY : fishing_strategy
    USE debmodel,     ONLY : Hp, Hj
    USE debmodel,     ONLY : Zaa, Zas, Zea, Zes, za, zs

    USE comtraj,      ONLY : imax, jmax
    USE comtraj,      ONLY : init_anchovy_egg, init_sardine_egg
#endif /* IBM_SPECIES */
    

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,4),           INTENT( in ) :: xe
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax),        INTENT( in ) :: sal, temp
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax,3),      INTENT( in ) :: uz, vz
    INTEGER,                                                INTENT( in ) :: Istr,Iend,Jstr,Jend


    !! * Local declarations
    INTEGER                         :: jj,mm_clock,aaaa,hh,minu,sec   ! To convert date to seconds or seconds to date
    INTEGER                         :: n,m,nn,i,j,k,ind         ! Integers for loops
    INTEGER                         :: nb_part                  ! Number of particles inside a patch for loop

    REAL(KIND=rlg)                  :: dtm                      ! Model time step
    INTEGER                         :: time_step
   
    ! Indexes for temporary particle position
    INTEGER                         :: igg, idd, jhh, jbb
    INTEGER                         :: hlb, hlt, hrb, hrt  
    INTEGER                         :: kp, km
    REAL(KIND=rsh)                  :: px, py 

    ! Indexes to calculate num in reproduction part
    INTEGER                         :: idx_s, idx_e, new_size, part_num, last_ind
    LOGICAL                         :: time_to_spawn, has_spawn                ! Logical for compliance with the reproduction time condition 
    
    CHARACTER(len=19)               :: tool_sectodat
    CHARACTER(len=8)                :: fileout_suffix

    REAL(KIND=rsh)                  :: zlag, slag

    TYPE(type_patch),    POINTER    :: patch,    child_patch
    TYPE(type_particle), POINTER    :: particle, new_particle

    ! To save a local and global position of particle for MPI and Sequential compatibility
    TYPE(type_position)             :: pos,pos_ad

    INTEGER                         :: ierr_mpi                     ! Integer returned by MPI functions
    INTEGER                         :: ind_species, child_ind_species, nb_new_particle

#ifdef IBM_SPECIES
    ! To save a local and global position of particle for MPI and Sequential compatibility in repro
    TYPE(type_position)             :: new_pos

    LOGICAL                         :: update, night
    REAL(KIND=rsh)                  :: depth_day, depth_night
    REAL(KIND=rsh)                  :: mindepth, maxdepth
    REAL(KIND=rsh)                  :: tempo, sal_part
    REAL(KIND=rsh)                  :: new_posx,new_posy            ! Eggs spawn position
    REAL(KIND=rsh)                  :: sal_surf,temp_surf,dens_surf ! Physical variables at particle's location
    REAL(KIND=rsh)                  :: dh

    REAL(KIND=rsh)                  :: harvest,rx,ry

    ! Variables for reproduction
    REAL(KIND=rsh), DIMENSION(imax+2,jmax+2,nb_species) :: mat_eggs, mat_all_eggs
    INTEGER, DIMENSION(imax+2,jmax+2)                   :: MAT_new_indv
    REAL(KIND=rsh)                                      :: ratio
    REAL(KIND=rsh)                                      :: new_super

    ! Density-dependence variables
    REAL(KIND=rsh), DIMENSION(nb_species)   :: dd_number_eggs, concentration, max_mort
    REAL(KIND=rlg), DIMENSION(nb_species)   :: Z1, Z2
    LOGICAL, SAVE                           :: first_timestep_ibm = .true.

#endif /* IBM_SPECIES */


    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! Save particle properties (before any change for getting exact initial properties) 
    CALL ibm_save
 
    ! If advection, depend de itypevert passive transport
    CALL LAGRANGIAN_update(xe,uz,vz,Istr,Iend,Jstr,Jend)
 
    ! Save old year for reproduction
    CALL tool_decompdate(tool_sectodat(time),jj,mm_clock,aaaa,hh,minu,sec)

#ifdef IBM_SPECIES
    jjulien = tool_julien(jj,mm_clock,aaaa) - tool_julien(1,1,aaaa) + 1
    IF ( debuse .AND. .NOT. F_Fix ) THEN
        CALL readfood3d(Istr,Iend,Jstr,Jend,first_timestep_ibm)
        first_timestep_ibm = .false.
        !CALL readtemp3d(Istr,Iend,Jstr,Jend,first_timestep_ibm)
    ENDIF
#endif /* IBM_SPECIES */

 !--------------------------------------------------------------
 ! Common part to all species
 !---------------------------------------------------------------
    dtm = dt   ! dt = CROCO time step
    time_step=nrhs

#ifdef IBM_SPECIES
    ! Menu et al., but script inside for density dependence from Bueno Pardo
    IF ( FIRST_TIME_STEP ) number_eggs(:)    = 0._rsh       ! For initialisation

    IF ( repro .and. density_dependent  ) THEN
        dd_number_eggs   = number_eggs
        number_eggs(:)   = 0.0_rlg
        concentration    = dd_number_eggs / 50.d9 ! 50,000 km2: area considered

        max_mort(1)      = (Zaa+(Zea-Zaa)*exp(-za*(0.35_rlg - 0.0855_rlg)))  ! For anchovy
        max_mort(2)      = (Zas+(Zes-Zas)*exp(-zs*(0.35_rlg - 0.0855_rlg)))  ! For sardine
    ENDIF

    ! Remise a 0 de la matrice des oeufs pondus avant boucle sur les patches
    IF ( repro ) mat_eggs       = 0._rsh
#endif


#ifdef MPI
   down_give  = 0
   up_give    = 0
   right_give = 0
   left_give  = 0
#endif

    has_spawn = .false.  ! Bool to know if a particle has released eggs during the time step

    ! ================================
    ! ===     Start loop on patches
    patch => patches%first
    DO n = 1,patches%nb

        IF ( (time < patch%t_beg).OR.(time > patch%t_end) ) THEN
            patch => patch%next
            CYCLE
        ENDIF

#ifdef IBM_SPECIES
        IF (patch%species == 'anchovy') ind_species = 1  
        IF (patch%species == 'sardine') ind_species = 2

        ! Implement time_to_spawn for the patch and then update patch%t_spawn
        time_to_spawn = (time >= patch%t_spawn) .AND. (time <= patch%t_end) .AND. (hh == 0) .AND. (current_day /= jj)
        IF ( time_to_spawn ) THEN 
            patch%t_spawn = patch%t_spawn + patch%dt_spawn
            has_spawn     = .true.
        ENDIF
#endif

        ! ================================
        ! ===     Start loop on particles 
        nb_part = patch%nb_part_alloc
        DO m = 1,nb_part

            particle => patch%particles(m)

            ! Skip killed or flagged fish
            IF ( .NOT. particle%active )  CYCLE
            IF (particle%flag == -valmanq) particle%super = 0._rlg
            IF ( (particle%super <= 1.0_rsh) .OR. (particle%flag == -valmanq) )  CYCLE

            ! Update age of fish
            IF (current_day /= jj) particle%age = particle%age + 1
            IF (yearclass == aaaa .and. particle%age > 101.0_rlg) particle%AgeClass = particle%AgeClass + 1

            pos%xp = particle%xpos ; pos%yp = particle%ypos
            call define_pos(pos)
            CALL ibm_loc_xyz(pos%idx_r,pos%idy_r,particle%spos,             &
                             px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,kp,km,   &
                             Istr,Iend,Jstr,Jend)

        !---------------------------------------------------------------
        ! Specific models with calls to subroutines
        !---------------------------------------------------------------

        ! ===================================
        ! ===   Fish life changes
#ifdef IBM_SPECIES
            ! size in cm
            update = .TRUE.

            ! --------------------------------------------------------------
            ! -----            STAGE 1  : eggs, buyoancy               -----
            !--------------------------------------------------------------- 
            IF (particle%stage == 1) THEN 
                particle%temp = ibm_traint(temp,xe,particle%spos,kp,km,px,py, &
                                           igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,   &
                                           Istr,Iend,Jstr,Jend)  

                IF (patch%species == 'anchovy') THEN
                    ! Egg density dependent on development (Adapted from Ospina 2012 with our data)
                    particle%density = particle%denspawn + 3.2029_rsh*particle%Drate - 7.4937_rsh*(particle%Drate)**2 &
                                       - 3.0858_rsh*(particle%Drate)**3 + 9.0376_rsh*(particle%Drate)**4 

                    ! Egg development rate (Regner, 1996)
                    particle%Drate = particle%Drate + (dtm/86400.0_rsh)*0.000559153_rsh*(particle%temp)**2.290236_rsh
                    ! Egg mortality
                    !particle % super = (1.0_rsh-0.266_rsh*dtm/86400.0_rsh)*particle % super
                ENDIF

                IF (patch%species == 'sardine') THEN
                    ! Egg density dependent on development (Garcia-Garcia)
                    particle%density = particle%denspawn + 3.58_rsh*(0.0012_rsh - 1.8528_rsh*particle%Drate &
                                       + 15.507_rsh*(particle%Drate)**2 - 41.1312_rsh*(particle%Drate)**3   &
                                       + 40.5289_rsh*(particle%Drate)**4) 

                    ! Egg developement miranda Peck (dvpt prédit pour dt en jour) dtm en s
                    particle%Drate = particle%Drate + ((dtm/86400.0_rsh)*0.17_rsh*(particle%temp)**1.9286_rsh)/100.0_rsh
                ENDIF

                ! Salinity at particle location
                sal_part = ibm_traint(sal,xe,particle%spos,kp,km,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)

                ! Buyoncy
                particle%w = ibm_buoy(particle%density, particle%size, particle%temp, sal_part)

                IF ( density_dependent ) number_eggs(ind_species) = number_eggs(ind_species) + particle%super

                IF ( particle%Drate >= 1.0_rsh ) THEN
                    particle%Drate = 0.0_rsh
                    particle%stage = 2
                    particle%size  = 0.35_rsh ! on force eclosion taille fixe
                    particle%w     = 0.0_rsh
                    update         = .FALSE.
                ENDIF  
            ENDIF ! end of stage 1

            ! --------------------------------------------------------------
            ! -----       STAGE 2  :  yolk sac larva, buoyancy         -----
            !---------------------------------------------------------------

            IF ( particle%stage == 2 .AND. update ) THEN
                particle%temp = ibm_traint(temp,xe,particle%spos,kp,km,px,py, &
                                           igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,   &
                                           Istr,Iend,Jstr,Jend) 

                ! Buoyancy
                ! sal_part = ibm_traint(sal, xe, particle % spos, kp, km, &
                !                                  px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                !                                  Istr,Iend,Jstr,Jend)
                ! particle % w = ibm_buoy(particle%density, particle%size, particle%temp, sal_part)

                IF (patch%species == 'anchovy') THEN
                    ! Yolk sac larva development rate (Boussouar et al., 2001), Emmanuelle
                    particle%Drate = particle%Drate + (dtm/3600_rsh)*0.0000311717*(particle%temp)**2.1749
                ENDIF
                IF (patch%species == 'sardine') THEN
                    ! Yolk sac larvae developement rate Zweifel and Lasker 1976 + Rose et al 2015 for sardinops sagax
                    particle%Drate = particle%Drate + (dtm/86400)*0.001_rsh*exp(6.19_rsh*(1.0_rsh - exp(-0.050_rsh*particle%temp)))&
                                                                           /log(5.15_rsh/(5.15_rsh - log(5.970_rsh/3.740_rsh)))     
                ENDIF

                ! Growth (simple as function of Drate as not handled in DEB model so far)
                particle%size = 0.35_rsh + (0.4_rsh - 0.35_rsh)*particle%Drate
                IF (particle%Drate >= 1.0_rsh) THEN
                    particle%stage = 3
                    particle%size  = 0.4_rsh ! on force ouverture bouche a taille fixe (cm) tant que yolk sac pas dans DEB, Regner 1996
                    update         = .FALSE.
                ENDIF
            ENDIF ! end of stage 2

            ! --------------------------------------------------------------
            ! -----    STAGE 3  : larva without vertical migration     -----
            !---------------------------------------------------------------
            IF (particle%stage == 3 .AND. update) THEN 
                particle%temp = ibm_traint(temp,xe,particle%spos,kp,km,px,py, &
                                           igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,   &
                                           Istr,Iend,Jstr,Jend) 

                depth_day   = 25.0_rsh
                depth_night = 25.0_rsh
                                
                !sal_part = ibm_traint(sal, xe, particle % spos, kp, km,    &
                !                      px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                !                      Istr,Iend,Jstr,Jend)
                !particle%w = ibm_buoy(particle%density, particle%size, particle%temp, sal_part)
                
                w_max      = particle%size / 10.0_rsh / 100.0_rsh ! m.s-1
                particle%w = ibm_nycth_mig(depth_day, depth_night, w_max, w_max, particle, alpha_w)
                ! to avoid deep diving at slope, to remove later on   
                  
                ! Growth
                IF ( debuse ) CALL deb_cycle(particle,dtm,aaaa,mm_clock,jj,patch%species)
                IF (particle%size >= 0.6_rsh) THEN 
                   particle%stage = 4
                   update         = .FALSE.
                ENDIF
            ENDIF ! end of stage 3
   
            ! --------------------------------------------------------------
            ! -----           STAGE 4  : nycthemeral migration         -----
            !---------------------------------------------------------------

            IF ( particle%stage == 4 .AND. update ) THEN
                particle%temp = ibm_traint(temp,xe,particle%spos,kp,km,px,py, &
                                             igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,   &
                                             Istr,Iend,Jstr,Jend) 
                depth_day   = 25.0_rsh
                depth_night = 0.0_rsh

                ! Calculate an optimal depth ?
                w_max      = particle%size/1d3 ! m.s-1
                particle%w = ibm_nycth_mig(depth_day, depth_night, w_max, w_max, particle, alpha_w)       

                ! Growth
                IF ( debuse ) CALL deb_cycle(particle,dtm,aaaa,mm_clock,jj,patch%species)
                IF (particle%H >= particle%Hj) THEN
                    particle%stage = 5
                    particle%itypevert = 0 ! Stop vertical advection/diffusion from stage 5
                    particle%zpos   = 0.0_rsh
                    particle%w      = 0.0_rsh
                    particle%dayjuv = jjulien
                    update          = .FALSE.
                ENDIF
            ENDIF ! end of stage 4

            ! -------------------------------------------------------------------------
            ! -----    STAGE 5  : Horizontal swimming, distribution from file     -----
            !--------------------------------------------------------------------------
            ! Juvenile
            IF (particle%stage == 5 .AND. update) THEN
                ! Vertical position
                night = .FALSE.
                IF(srflx(NINT(pos%idx_r),NINT(pos%idy_r))*rho0*Cp == 0.0_rsh) THEN
                   night=.TRUE.
                ENDIF
                IF (particle%H < Hp .AND. jjulien >= particle%dayjuv .AND. jjulien < 330 ) THEN
                    IF ( night ) THEN ! Juvenile Nuit
                        mindepth = 0.0_rsh
                        maxdepth = 30.0_rsh
                    ELSE! Juvenile day
                        mindepth = min(30.0_rsh, 2.27_rsh*particle%size) ! Boyra, 2013
                        maxdepth = 30.0_rsh + (120.0_rsh - 30.0_rsh)*(1.0_rsh - (330.0_rsh - REAL(jjulien,rsh)) &
                                                                               /(330.0_rsh - REAL(particle%dayjuv,rsh)))
                    ENDIF
                ELSE
                   IF (jjulien <= 90 .OR. jjulien >= 330) THEN !winter 
                      mindepth = 40.0_rsh
                      maxdepth = 120.0_rsh
                   ELSE
                      IF ( night ) THEN ! spring/summer
                         ! Adulte Nuit
                         mindepth = 0.0_rsh
                         maxdepth = 40.0_rsh
                      ELSE ! Adult day
                         mindepth = 30.0_rsh 
                         maxdepth = 120.0_rsh
                      ENDIF
                   ENDIF
                ENDIF

                !PRINT*, 'Min max Depth :', particle%num, particle%stage,mindepth,maxdepth,jjulien,particle%dayjuv,particle%H,Hp
                particle%temp = ibm_proftraint(temp,mindepth,maxdepth,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

                ! Growth
                IF ( debuse ) CALL deb_cycle(particle,dtm,aaaa,mm_clock,jj,patch%species)

                ! Time step to calculate fish movement, depending on size of the fish
                dh = NINT((100.0/(1.5*particle%size*3600.0))*        &
                     MAX(om_r(nint(pos%idx_r),nint(pos%idy_r)), on_r(nint(pos%idx_r),nint(pos%idy_r))))

                IF (hh >= particle%hmove + dh) THEN
#ifdef MPI
                        CALL fish_move(particle,ind_species)
#endif    
                    particle%hmove = hh ! update of the saved hour

                    pos_ad%xp = particle%xpos; pos_ad%yp = particle%ypos
                    CALL define_pos(pos_ad)

                    ! total depth at particle s location
                    CALL loc_h0(pos_ad%idx_r,pos_ad%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, &
                                        Istr,Iend,Jstr,Jend)
                    particle%xe = xeint(xe(:,:,time_step),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, &
                                        Istr,Iend,Jstr,Jend)
                    particle%h0 = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                    particle%d3 = particle%h0 + particle%xe
                    particle%hc = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

                ENDIF
                IF (hh == 0 .and. particle%hmove > 0) particle%hmove = particle%hmove - 24

                ! SI le super individu n'atteint pas le stade 6 et que jjulien = dayjuv, 
                ! alors mindepth = maxdepth et l'interpolation de la temperature plante.
                IF (particle%H >= particle%Hp .or. particle%age == 364) THEN
                    particle%stage = 6
                    update         = .FALSE.
                ENDIF
            ENDIF  ! end of stage 5 

            ! -------------------------------------------------------------------------
            ! -----    STAGE 6  : Horizontal swimming, distribution from file     -----
            !--------------------------------------------------------------------------
            ! Adult
            IF (particle%stage == 6) THEN
                ! Vertical position
                night = .FALSE.

                IF(srflx(NINT(pos%idx_r),NINT(pos%idy_r))*rho0*Cp <= 0.0_rsh) THEN       ! Valeur seuil de radiation pour l'alternance jour/nuit
                    night=.TRUE.
                ENDIF
                IF (jjulien <= 90 .OR. jjulien >= 330) THEN !winter 
                    mindepth = 40.0_rsh
                    maxdepth = 120.0_rsh
                ELSE
                    IF (night) THEN ! spring/summer
                        ! Adulte Nuit
                        mindepth = 0.0_rsh
                        maxdepth = 40.0_rsh
                    ELSE ! Adult day
                        mindepth = 30.0_rsh 
                        maxdepth = 120.0_rsh
                    ENDIF
                ENDIF

                particle%temp = ibm_proftraint(temp,mindepth,maxdepth,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

                IF (debuse) CALL deb_cycle(particle,dtm,aaaa,mm_clock,jj,patch%species)

                ! Time step to calculate fish movement, depending on size of the fish and size of the grid
                dh = NINT((100.0/(1.5*particle%size*3600.0))*        &
                     MAX(om_r(nint(pos%idx_r),nint(pos%idy_r)), on_r(nint(pos%idx_r),nint(pos%idy_r))))

                IF (hh >= particle%hmove + dh) THEN
#ifdef MPI
                        CALL fish_move(particle,ind_species)
#endif    
                    particle%hmove = hh ! update of the saved hour

                    pos_ad%xp = particle%xpos; pos_ad%yp = particle%ypos
                    CALL define_pos(pos_ad)

                    ! total depth at particle s location
                    CALL loc_h0(pos_ad%idx_r,pos_ad%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, &
                                        Istr,Iend,Jstr,Jend)
                    particle%xe = xeint(xe(:,:,time_step),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, &
                                        Istr,Iend,Jstr,Jend)
                    particle%h0 = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                    particle%d3 = particle%h0 + particle%xe
                    particle%hc = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

                ENDIF
                IF (hh == 0 .and. particle%hmove > 0) particle%hmove = particle%hmove - 24 
            ENDIF

            ! ===   Fin evolution stades de vie
            ! ===================================


            ! ===========================================
            ! ===    Divers mortalites sur les poissons 

            !---------------------------------------------------------------
            ! Density dependance, name struc_ad misleading, we take everybody (Menu et al. 2023)
            IF (particle%stage >= 3) struc_ad = struc_ad + (particle%WV*particle%super)  

            !---------------------------------------------------------------
            ! Mortality (Menu et al.)
            Z1(ind_species) = 0.0_rlg
            IF (patch%species == 'anchovy') THEN        ! Equivalent to (ind_species == 1)
                Z1(ind_species) =  (Zaa + (Zea - Zaa)*exp(-za*(particle%size - 0.0855_rlg))) ! size egg EN DUR
            ELSE IF (patch%species == 'sardine') THEN   ! Equivalent to (ind_species == 2)
                Z1(ind_species) =  (Zas + (Zes - Zas)*exp(-zs*(particle%size - 0.1608_rlg))) ! size egg EN DUR
            ENDIF





            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!         DENSITE - DEPENDANCE NON SAPTIALISEE        !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! --- Density dependent mortality of larvae in relation to the concentration of eggs (Somarakis, 2007) - density dependence from bueno pardo
            IF ( density_dependent ) THEN
                IF (particle%stage == 3 .OR. particle%stage == 4) THEN
                    IF (concentration(ind_species) == 0._rsh) THEN
                        Z2(ind_species) = 0
                    ELSE 
                        Z2(ind_species) =  -0.154_rlg + slope*log10(concentration(ind_species))
                    ENDIF
                    Z1(ind_species) = max(Z1(ind_species), Z1(ind_species) - max_mort(ind_species) + Z2(ind_species))
                ENDIF
            ENDIF

            ! The mortality of Somarakis is daily, so 86400 is OK
            particle%Death_NAT = particle%Death_NAT + particle%super*(1 - exp(-(Z1(ind_species)*dtm/86400_rlg))) 
            particle%super     = particle%super*exp(-Z1(ind_species)*dtm/86400_rlg)



            ! ================================================
            ! ===                                          ===
            ! ===             FISHING MORTALITY            ===
            IF ( fish_mort ) CALL death_by_fishing(particle,patch%species,aaaa,mm_clock,dtm)

            ! Some variables on whole population to calculate fishing
            IF (particle%stage >= 5 .and. particle % AgeClass >= 1) THEN
                number_tot(ind_species) = number_tot(ind_species) + particle%super!*                 &
                                                                    !selec_dome_or_asymp(particle%size,alpha_sel,beta_sel)
                weight_tot(ind_species) = weight_tot(ind_species) + (particle%Wdeb*particle%super)!* &
                                                                    !selec_dome_or_asymp(particle%size,alpha_sel,beta_sel)
            ENDIF

            IF (ind_species == 1 .and. jjulien == 135 .and. fishing_strategy == 'HCR' .and. particle%AgeClass >= 1 ) THEN
                number_tot(ind_species) = number_tot(ind_species) + particle%super
                weight_tot(ind_species) = weight_tot(ind_species) + (particle%Wdeb*particle%super)
            ENDIF 

            !---------------------------------------------------------------
            ! Add all eggs in mat_eggs, matrix for spawn, looking at the species
            IF( repro ) THEN
                ! When time to spawn comes, particle releases all accumulated eggs (Neggs) in mat_eggs where it is  
                IF (time_to_spawn .and. particle%Neggs > 0.0_rsh) THEN

                    mat_eggs(NINT(particle%xpos),NINT(particle%ypos),ind_species) = particle%Neggs +   &
                                          mat_eggs( NINT(particle%xpos), NINT(particle%ypos),ind_species)

                    ! Collect number of eggs for density-dependence on larvae
                    IF ( density_dependent ) number_eggs(ind_species) = number_eggs(ind_species) + particle%Neggs

                    particle%Neggs = 0.0_rsh
                END IF

            END IF

#endif  /* IBM_SPECIES */

            
            !--------------------------------------------------------------------------------
            ! Vertical movement and check for bottom and surface boundaries 
            zlag = -particle%zpos + particle%xe ! switch from immersion to real z
            
            zlag = zlag + particle%w*dtm                       ! update position
            zlag = MIN(MAX(zlag, -particle%h0), particle%xe)   ! check boundaries
            
            particle%zpos = -zlag + particle%xe  ! immersion
            
            !IF ( particle%stage <= 2 ) particle%zpos = min(particle%zpos, 40.0_rsh) !patch to avoid going too deep...
            CALL ztosiggen(zlag, slag, particle%xe, particle%h0, particle%hc)
            particle%spos = slag
        ENDDO ! particle

        patch => patch%next

    ENDDO ! patch

        


        ! ===================================================================================================
        ! =====                                                                                         =====
        ! =====                                     REPRODUCTION                                        =====
        ! =====                                                                                         =====
        ! ===================================================================================================
#ifdef IBM_SPECIES

        IF ( repro ) THEN

            ! Test if we have changed the year and update newseason for all species
            IF (current_year /= aaaa)  newseason = .TRUE.

            ! Init des objets avant MPI pour eviter toute valeur par defaut
            mat_all_eggs     = 0._rsh

            ! Si temps a pondu, partage des informations de reproduction obtenus par chaque proc 
            IF ( has_spawn ) THEN
#ifdef MPI
                ! Somme tous les oeufs relaches dans mat_all_eggs et partage cette matrice a tous les procs
                CALL MPI_ALLREDUCE(mat_eggs,mat_all_eggs,(imax+2)*(jmax+2)*nb_species,MPI_DOUBLE_PRECISION, &
                                   MPI_SUM,MPI_COMM_WORLD,ierr_mpi)
#else
                mat_all_eggs    = mat_eggs
#endif
            END IF  ! has_spawn

            DO ind = 1,nb_species 
                spawn(ind)       = (has_spawn .AND. SUM(mat_all_eggs(:,:,ind)) > 0.0_rlg) ! Some particles in the patch are spawning...
                first_spawn(ind) = spawn(ind) .AND. newseason(ind)                        ! ... and this is the first time . 

                ! --- Create new patch and new file to store new patch (at each first spawn of each year for each species)
                IF ( first_spawn(ind) ) THEN
                    newseason(ind) = .FALSE. 

                    ! If spawning happen for the first time, we have to create an appropriate data structure :
                    ! it is a new patch (a child one).
                    child_patch => patch_list_append(patches)
                    child_patch%t_beg       = time
                    child_patch%t_end       = child_patch%t_beg + duration(ind)*24.0_rlg*3600.0_rlg
                    child_patch%t_save      = time
                    child_patch%dt_save     = dt_save
                    child_patch%t_spawn     = time
                    child_patch%dt_spawn    = dt_spawn*3600._rsh
                    child_patch%parent_id   = patches%nb + 1                ! keep track of the childs parent
                    child_patch%yearref     = current_year
                    call read_run_info(child_patch%run_id)  ! read file FOIL.info

                    ! Inherit default values for new particles
                    IF (ind == 1) THEN
                        child_patch%init_particle = init_anchovy_egg
                        WRITE( fileout_suffix, '("_",i0,".nc")' ) current_year  ! child_patch % generation
                        child_patch%file_out      = trim(file_pathout) // "anchovy" // fileout_suffix
                        child_patch%species       = "anchovy"    ! Initialize species of the patch
                    ELSE IF (ind == 2) THEN
                        child_patch%init_particle = init_sardine_egg
                        WRITE( fileout_suffix, '("_",i0,".nc")' ) current_year  ! child_patch % generation
                        child_patch%file_out      = trim(file_pathout) // "sardine" // fileout_suffix
                        child_patch%species       = "sardine"    ! Initialize species of the patch
                    ENDIF

                    ! Allocate memory for new particles
                    CALL init_patch(child_patch, 0)               ! No allocation yet
                    child_patch%nb_part_max  = max_part           ! Initialize max part allowed in patch
                    first_spawn(ind_species) = .FALSE.

                END IF    ! first_spawn  

                ! --- Create new particles in the patch
                IF ( spawn(ind) ) THEN
                    ! Look for the child patch
                    child_patch => patches%first
                    IF(child_patch%species == 'anchovy') child_ind_species = 1
                    IF(child_patch%species == 'sardine') child_ind_species = 2

                    DO WHILE ((child_patch%yearref /= aaaa) .or. (child_ind_species /= ind))
                        child_patch => child_patch%next
                        IF(child_patch%species == 'anchovy') child_ind_species = 1
                        IF(child_patch%species == 'sardine') child_ind_species = 2
                        IF ( .NOT. ASSOCIATED(child_patch) ) THEN
                            ! This shall never happen !
                            PRINT*, " ERROR: child patch not found !!! I am ", child_patch%id
                            CALL_MPI MPI_FINALIZE(ierr_mpi)
                            STOP
                        END IF
                    ENDDO

                    ! Determination du nombre d'individus qui vont apparaitre, selon le nombre vise et le nombre d'oeufs pondus
                    ratio = nb_indv_ponte(ind)/SUM(mat_all_eggs(:,:,ind))
                    MAT_new_indv = NINT(mat_all_eggs(:,:,ind)*ratio)  

                    ! Recherche du nombre de particules a creer ( avec condition dans le domaine du proc si MPI)
                    nb_new_particle = 0
                    DO i=1,imax+2
                    DO j=1,jmax+2
                        IF(MAT_new_indv(i,j) .ne. 0.0_rsh) THEN
                            IF_MPI ( iminmpi <= i .AND. i <= imaxmpi .AND. jminmpi <= j .AND. j <= jmaxmpi ) THEN
                                nb_new_particle = nb_new_particle + MAT_new_indv(i,j)
                            ENDIF_MPI
                        ENDIF
                    ENDDO
                    ENDDO

                    !! Compute start index for numbering new particles
                    CALL indices_loc2glob(child_patch%nb_part_total, nb_new_particle, idx_s, idx_e)
                    ! Add number of new particles proc per proc if MPI, all at once otherwise
                    CALL_MPI ADD_ALL_MPI_INT(nb_new_particle)
                    
                    child_patch%nb_part_total = nb_new_particle + child_patch%nb_part_total

                    ! --- Increase size of child patch
                    new_size = child_patch%nb_part_total
                    CALL resize_patch(child_patch, new_size)

                    ! Si aucune nouvelle particule (notamment si MPI proc sans nouveau super indv), alors on passe au proc suivant
                    IF (nb_new_particle == 0) CYCLE
                    
                    part_num  = 0    ! Initialisation pour le calcul de NUM de la particule

                    ! Calcul du nombre d'indv par super, egal pour tous les nouveaux supers sur un pas de temps
                    ! Cond. IF Juste in case, ne devrait jamais se produire.
                    IF (nb_new_particle > 0) new_super = SUM(mat_all_eggs(:,:,ind))/nb_new_particle

                    ! Parcours de MAT_new_indv pour avoir les nouvelles particules
                    DO i=1,imax+2
                    DO j=1,jmax+2
                        IF (MAT_new_indv(i,j) .ne. 0.d0) THEN
                        ! Check if inside good MPI proc or if inside the domain (if not MPI)
#ifdef MPI                  
                            IF ( iminmpi <= i .AND. i <= imaxmpi .AND. &
                                 jminmpi <= j .AND. j <= jmaxmpi ) THEN
#else
                            IF ( i >= Istr .AND. i <= Iend .AND. &
                                 j >= Jstr .AND. j <= Jend ) THEN
#endif
                                ! Boucle sur le nombre d'indv a creer selon la valeur de MAT_new_indv
                                DO k=1,MAT_new_indv(i,j)
                                    ! New particle number
                                    part_num = part_num + 1
                                    last_ind = idx_s + part_num

                                    ! Check if we have not reached the limit allowed by NetCDF output (+10%)
                                    IF (child_patch%nb_part_total > child_patch%nb_part_max + child_patch%nb_part_max*0.1) THEN
                                        PRINT*, " ERROR: Not enough space in NetCDF file to store additional particles." //     &
                                                " We have", child_patch%nb_part_max, " and need", child_patch%nb_part_total +   &
                                                child_patch%nb_part_total*0.1,                                                  &
                                                " Please consider increasing nb_part_max."
                                        CALL_MPI MPI_FINALIZE(ierr_mpi)
                                        STOP
                                    ENDIF

                                    ! --- Initialize bio from patch save info
                                    new_particle => child_patch%particles(last_ind)
                                    new_particle = child_patch%init_particle
                                    new_particle % active    = .True.
                                    new_particle % num       = last_ind
                                    new_particle % date_orig = time

                                    ! Initialize spawn position of the particle randomly in the new cell
                                    CALL random_number(harvest)
                                    new_particle%xpos = min(real(i,kind=rsh) + harvest - 0.49_rsh, real(i,kind=rsh) + 0.49_rsh)
                                    CALL random_number(harvest)
                                    new_particle%ypos = min(real(j,kind=rsh) + harvest - 0.49_rsh, real(j,kind=rsh) + 0.49_rsh)

                                    ! Randomly place the particle in z betwee, 0 and 20m
                                    CALL random_number(harvest)
                                    new_particle%zpos       = 20*abs(harvest)

                                    new_particle%stage      = 1

                                    new_pos%xp = new_particle%xpos; new_pos%yp = new_particle%ypos    ! New particle global position
                                    call define_pos(new_pos)                                          ! Convert to local position inside proc if MPI, doesn't change anything in sequential

                                    ! Get bathymetry information for ztosiggen function
                                    CALL loc_h0(new_pos%idx_r,new_pos%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, &
                                                Istr,Iend,Jstr,Jend)
                                    new_particle%h0 = h0int(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                                    new_particle%xe = xeint(xe(:,:,nrhs),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,    &
                                                            Istr,Iend,Jstr,Jend)
                                    new_particle%hc = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                                    new_particle%d3 = new_particle%xe + new_particle%h0

                                    ! Init IBM parameters of the particle
                                    IF (debuse) CALL ibm_parameter_init(new_particle,child_patch%species,xe,sal,temp, &
                                                                        Istr,Iend,Jstr,Jend)

                                    new_particle%super = new_super

                                    ! Init DEB parameters
                                    IF ( debuse ) CALL deb_egg_init(new_particle, child_patch%species)

                                    ! Biological parameters
                                    IF ( debuse )  new_particle%yearspawn = aaaa + 1  ! on ne pond qu'a 1 an minimum
                                    CALL_MPI init_mpi_type_particle                   ! Init created particle for MPI

                                ENDDO
                            ENDIF   ! localisation in domain
                        ENDIF   ! Mat_new_eggs
                    ENDDO
                    ENDDO   ! End do on new_eggs_sum matrix

                END IF      ! spawn
                spawn = .false.
            ENDDO   ! Nb species
            
            ! Remise a 0 de la matrice des oeufs pondus avant fin de ce pas de temps
            mat_eggs       = 0._rsh

        END IF      ! repro
#endif /* IBM_SPECIES */

        !!! Partie commentee dans le code de Clara
        !IF (jjulien == 135 .and. fishing_strategy == 'HCR') THEN
        !           
        !    IF (weight_tot <= 24000000000.0_rlg) THEN
        !       Zfishing =  0.0_rlg
        !    ENDIF
        !            
        !    IF (weight_tot > 24000000000.0_rlg .AND. weight_tot <= 89000000000.0_rlg) THEN
        !       allowed_catch_g = (-2600000000.0_rlg + (0.40 * weight_tot)) / 365.0_rlg
        !       allowed_catch_g = allowed_catch_g * multiplier_tac
        !       mean_individual_weight = weight_tot / number_tot
        !       allowed_catch_n = allowed_catch_g / mean_individual_weight
        !       Zfishing = -log ( (number_tot - allowed_catch_n) / number_tot)
        !    ENDIF
        !            
        !    IF (weight_tot > 89000000000.0_rlg) THEN
        !       allowed_catch_g = 33000000000.0_rlg / 365.0_rlg
        !       allowed_catch_g = allowed_catch_g * multiplier_tac
        !       mean_individual_weight = weight_tot / number_tot
        !       allowed_catch_n = allowed_catch_g / mean_individual_weight
        !       Zfishing = -log ( (number_tot - allowed_catch_n) / number_tot)
        !    ENDIF
        !                                            
        !    print*, weight_tot, allowed_catch_n, number_tot, mean_individual_weight, Zfishing, situation, year
        !            
        !ENDIF

    ! Update current day and year once all particles were updated
    current_year = aaaa
    current_day  = jj
    IF (yearclass == aaaa) yearclass = yearclass + 1

#ifdef IBM_SPECIES
    ! Parametre pour mortalite et densite-dependance
    struc_ad_dd_DEB = struc_ad
    struc_ad        = 0._rsh

    ! For catches
    DO i=1,nb_species
        IF (number_tot(i) > 0._rsh) THEN
            Wdeb_mean(i) = weight_tot(i)/number_tot(i)
            biom_tot(i)  = weight_tot(i)

            weight_tot(i) = 0._rsh
            number_tot(i) = 0._rsh
        ENDIF
    ENDDO
#endif /* IBM_SPECIES*/

#ifdef MPI
    ! Exchange particles if it changed proc domain, because of fish_move
    CALL ex_traj(down_give, up_give, right_give, left_give)
#endif

  END SUBROUTINE ibm_3d





  !!======================================================================
  SUBROUTINE ibm_save
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE ibm_save  ***
    !&E
    !&E ** Purpose : Save trajectories and DEB-IBM variables 
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_init, ibm_3d
    !&E ** External calls : ionc4_createfile_traj,ionc4_createvar_traj,ionc4_write_traj
    !&E                     ionc4_write_trajt,ionc4_write_time,ionc4_sync
    !&E                     indices_loc2glob,tool_ind2lat,tool_ind2lon
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E       !  2011-01 (M. Huret) 
    !&E       !  2011-10 (M. Huret) Introduction of MPI
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves, M. Huret) Coupled with CROCO
    !&E---------------------------------------------------------------------

    !! * Modules used
    USE ionc4,          ONLY : ionc4_createfile_traj, ionc4_createvar_traj, &
                               ionc4_write_traj, ionc4_write_trajt,         &
                               ionc4_write_time, ionc4_sync, ionc4_gatt_char, &
                               ionc4_gatt_char_read, ionc4_open
    USE comtraj,        ONLY : patches, type_patch, type_particle
 
    USE trajinitsave,   ONLY : indices_loc2glob
    USE trajectools,    ONLY : tool_ind2lat,tool_ind2lon
 
    !! * Arguments
 
    !! * Local declarations
    CHARACTER(LEN=lchain)                       :: file_out
    LOGICAL                                     :: l_out_nc4par
    INTEGER                                     :: n,m,num1,num2,nb_part,nb_part_nc,p
    REAL(kind=out)                              :: fillval

    INTEGER,        ALLOCATABLE, DIMENSION(:)   :: num_out,stage_out
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: lat_out,lon_out,dateo_out
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: temp_out,flag_out,spos_out,zpos_out,xpos_out,ypos_out 
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: h0pos_out,size_out,nb_out,dens_out,Drate_out
    INTEGER,        ALLOCATABLE, DIMENSION(:)   :: age_out,ageClass_out
#ifdef IBM_SPECIES
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: food_out,f_out,Wdeb_out,Denspawn_out
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: E_out,H_out,R_out,Neggs_out,NRJ_out,Gam_out
    INTEGER,        ALLOCATABLE, DIMENSION(:)   :: dayjuv_out, dayspawn_out,yearspawn_out,season_out
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: deaddeb_out, deadfishing_out, deadnatural_out
    REAL(KIND=out), ALLOCATABLE, DIMENSION(:)   :: zoom_out
#endif /*IBM_SPECIES*/
    TYPE(type_patch),    POINTER    :: patch
    TYPE(type_particle), POINTER    :: particle
    INTEGER :: idx_s, idx_e
    LOGICAL                                     :: out_ex
    character(len=32)                           :: fileinfo_run_id
    character(len=64) :: run_id_out
 
    !!----------------------------------------------------------------------
    !! * Executable part 

    fillval = dg_valmanq_io 
#ifdef MPI
     l_out_nc4par = .true.
#else
     l_out_nc4par = .false.
#endif

    patch => patches % first

    DO n = 1,patches % nb
        IF ( (time < patch%t_save) .OR. (time > patch%t_end) ) THEN
            patch => patch%next
            CYCLE
        END IF
 
        file_out = trim(patch%file_out)

        CALL read_run_info(fileinfo_run_id)  ! read file FOIL.info

        INQUIRE(file=file_out,exist=out_ex) ! does the file exist ?
        IF (out_ex .AND. .NOT. patch%file_out_init) THEN ! file exists but not yet opened, need to check run_id for restart
            patch%file_out_init = .TRUE.
            CALL ionc4_open(file_out, .false.)
            CALL ionc4_gatt_char_read(file_out, 'run_id', run_id_out)

            ! Vérifier la cohérence du run_id
            IF (trim(run_id_out) /= trim(patch%run_id)) THEN
                print *, 'ERROR: File ', trim(file_out), ' belongs to another run: ', trim(run_id_out)
                print *, 'Current run_id is: ', trim(patch%run_id)
                out_ex = .FALSE.  ! Force re-creation of the file
            ENDIF
        ENDIF !out_ex .and. .NOT. patch%file_out_init

        IF (out_ex .AND. patch%file_out_init) THEN ! file has been already opened
            CALL ionc4_open(file_out, .false.)
        ENDIF
                

        IF (.NOT. out_ex) THEN
            ! Create output file
            nb_part_nc = patch%nb_part_max

            CALL ionc4_createfile_traj(file_out, nb_part_nc, 0, 0, l_out_nc4par=l_out_nc4par)
            CALL ionc4_gatt_char(file_out, 'run_id', trim(patch%run_id)) ! Add a global attribute run_id
            CALL ionc4_createvar_traj(file_out, "latitude", "degrees_north","latitude",                    &
                                                fill_value=REAL(dg_valmanq_io,kind=out), l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "longitude","degrees_east", "longitude",                   &
                                                fill_value=REAL(dg_valmanq_io,kind=out), l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DEPTH","m","depth as immersion",                          &
                                                fill_value=-fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "H0","m","h0",                                             &
                                                fill_value=-fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "NUM","","number of the particle",                         &
                                                fill_value=0,        l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "flag","nbr","flag",                                       &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "TEMP","degrees_Celsius","temperature",                    &
                                                fill_value=fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "SIZE","Centimeters","Size",                               &
                                                fill_value=fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "STAGE","","Stage",                                        &
                                                fill_value=-1,       l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "NUMBER","Number","Number of individuals",                 &
                                                fill_value=fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DENSITY","sigma","Density",                               &
                                                fill_value=fillval, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DRATE","","Development rate of egg or larva",             &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DAYBIRTH","","Date of birth",                             &
                                                fill_value=REAL(dg_valmanq_io,kind=out), l_out_nc4par=l_out_nc4par, ndims=1)
            CALL ionc4_createvar_traj(file_out, "AGE","","Age in days",                                    &
                                                fill_value=-1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "AGECLASS","","Age in year of fish",                       &
                                                fill_value=-1,l_out_nc4par=l_out_nc4par)

#ifdef IBM_SPECIES
            CALL ionc4_createvar_traj(file_out, "FOOD","mg/m3","food",                                     &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "F","","f",                                     &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "EDEB","Joules","Energy",                                  &
                                                     fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "HDEB","Joules","Maturity",                                &
                                                     fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "RDEB","Joules","Repro",                                   &
                                                     fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "GAM","Joules","Energy gametes",                           &
                                                     fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "WEIGHT","g","weight",                                     &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "NEGGS","Number","Number of eggs",                         &
                                                fill_value=fillval,  l_out_nc4par=l_out_nc4par)
            !CALL ionc4_createvar_traj(file_out, "NRJ","J/g","Energy Density",                              &
            !                                    fill_value=fillval,  l_out_nc4par=l_out_nc4par)

            CALL ionc4_createvar_traj(file_out, "YEARSPAWN","","Year authorised to spawn",                 &
                                                fill_value=0, ndims=1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DAYSPAWN","","Julian day start spawning",                 &
                                                fill_value=0, ndims=1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DAYJUV","","Julien day at metamorphosis",                 &
                                                fill_value=0, ndims=1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "ZOOM","","Zoom value",                                    &
                                                fill_value=0, ndims=1, l_out_nc4par=l_out_nc4par)
            !CALL ionc4_createvar_traj(file_out, "SEASON","","Wether within spawning season",               &
            !                                    fill_value=-1, ndims=1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DENSPAWN","sigma","Density of egg at spawning",           &
                                                fill_value=fillval, ndims=1, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "Death_DEB","","Number dead by starvation",                &
                                                fill_value=fillval,l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "Death_FISH","","Number dead by fishing",                  &
                                                fill_value=fillval,l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "Death_NAT","","Number dead by natural mortality",         &
                                                fill_value=fillval,l_out_nc4par=l_out_nc4par)
#endif /*IBM_SPECIES*/
        ENDIF  ! (.NOT. out_ex)
 
        nb_part = patch%nb_part_alloc
        ALLOCATE( lat_out(nb_part),   lon_out(nb_part) )
        ALLOCATE( xpos_out(nb_part),  ypos_out(nb_part) )
        ALLOCATE( spos_out(nb_part),  zpos_out(nb_part) )
        ALLOCATE( num_out(nb_part),   h0pos_out(nb_part) )
        ALLOCATE( flag_out(nb_part),  dens_out(nb_part) )
        ALLOCATE( temp_out(nb_part),  Drate_out(nb_part) )
        ALLOCATE( size_out(nb_part),  dateo_out(nb_part) )
        ALLOCATE( stage_out(nb_part), nb_out(nb_part) )
        ALLOCATE( age_out(nb_part),   ageClass_out(nb_part) )
    
        lat_out(:)  = REAL(dg_valmanq_io,kind=out) ; lon_out(:) = REAL(dg_valmanq_io,kind=out)
        dateo_out(:) = REAL(dg_valmanq_io,kind=out)
        xpos_out(:) = fillval ; ypos_out(:) = fillval ; spos_out(:) = fillval ; zpos_out(:) = -fillval
        h0pos_out(:)=-fillval ; flag_out(:) = fillval ; num_out(:)  = 0
        size_out(:) = fillval ; stage_out(:) = -1 ; dens_out(:) = fillval ; Drate_out(:) = fillval
        temp_out(:) = fillval ; nb_out(:) = fillval; age_out(:) = -1; ageClass_out(:) = -1 
 
#ifdef IBM_SPECIES
        ALLOCATE( dayjuv_out(nb_part), dayspawn_out(nb_part) )
        ALLOCATE( yearspawn_out(nb_part), season_out(nb_part) )
        ALLOCATE( zoom_out(nb_part) )
        ALLOCATE( Denspawn_out(nb_part) )
        ALLOCATE( food_out(nb_part))
        ALLOCATE( H_out(nb_part), E_out(nb_part), Gam_out(nb_part), R_out(nb_part))
        ALLOCATE( Wdeb_out(nb_part) )
        ALLOCATE( Neggs_out(nb_part) )
        ALLOCATE( f_out(nb_part) )
        ALLOCATE( deaddeb_out(nb_part) )
        ALLOCATE( deadfishing_out(nb_part) )
        ALLOCATE( deadnatural_out(nb_part) ) 
    
        dayjuv_out(:) = 0 ; dayspawn_out(:) = 0 ; yearspawn_out(:) = 0 ; season_out(:) = -1
        Denspawn_out(:) = fillval ; food_out(:) = fillval ; Wdeb_out(:) = fillval
        H_out(:) = fillval ; E_out(:) = fillval ; R_out(:) = fillval ; Gam_out(:) = fillval
        f_out(:) = fillval ; zoom_out(:) = 0 ; Neggs_out(:) = fillval
        deaddeb_out(:) = fillval ; deadfishing_out(:) = fillval ; deadnatural_out(:) = fillval
#endif /*IBM_SPECIES*/
 
        p = 0
        DO m = 1,nb_part
            particle => patch % particles(m)
            IF ( .NOT. particle % active ) CYCLE
            p = p+1
            xpos_out(p)  = REAL(particle % xpos,      kind=out)
            ypos_out(p)  = REAL(particle % ypos,      kind=out)
            lat_out(p)   = REAL(tool_ind2lat(particle%xpos, particle%ypos), kind=out)
            lon_out(p)   = REAL(tool_ind2lon(particle%xpos, particle%ypos), kind=out)
            dateo_out(p) = REAL(particle % date_orig, kind=out)
            spos_out(p)  = REAL(particle % spos,      kind=out)
            zpos_out(p)  = REAL(particle % zpos,      kind=out) ! for output as immersion
            h0pos_out(p) = REAL(particle % d3,        kind=out)
            flag_out(p)  = REAL(particle % flag,      kind=out)
            temp_out(p)  = REAL(particle % temp,      kind=out)
            size_out(p)  = REAL(particle % size,      kind=out)
            stage_out(p) = REAL(particle % stage,     kind=out)
            dens_out(p)  = REAL(particle % density,   kind=out)
            nb_out(p)    = REAL(particle % super,     kind=out)
            num_out(p)   = REAL(particle % num,       kind=out)
            Drate_out(p) = REAL(particle % Drate,     kind=out)
            age_out(p)   = REAL(particle % age,       kind=out)
            ageClass_out(p) = REAL(particle % AgeClass,  kind=out)
#ifdef IBM_SPECIES
            food_out(p)  = REAL(particle % X,         kind=out)
            f_out(p)     = REAL(particle % f,         kind=out)
            E_out(p)     = REAL(particle % E,         kind=out)
            H_out(p)     = REAL(particle % H,         kind=out)
            R_out(p)     = REAL(particle % R,         kind=out)
            Gam_out(p)   = REAL(particle % Gam,       kind=out)
            Wdeb_out(p)  = REAL(particle % Wdeb,      kind=out)
            Neggs_out(p) = REAL(particle % Neggs,     kind=out)
            yearspawn_out(p) = REAL(particle % yearspawn, kind=out)
            dayspawn_out(p)  = REAL(particle % dayspawn,  kind=out)
            dayjuv_out(p)    = REAL(particle % dayjuv,    kind=out)
            denspawn_out(p)  = REAL(particle % denspawn,  kind=out)
            zoom_out(p)      = REAL(particle % zoom,      kind=out)
            deaddeb_out(p)      = REAL(particle % Death_DEB,    kind=out)
            deadfishing_out(p)  = REAL(particle % Death_FISH,   kind=out)
            deadnatural_out(p)  = REAL(particle % Death_NAT,    kind=out)
#endif /*IBM_SPECIES*/
        END DO
 
        CALL ionc4_write_time(file_out,0,time)
    
        nb_part = count(patch%particles(:)%active)
        CALL indices_loc2glob(1, nb_part, idx_s, idx_e)
 
#ifdef MPI
        num1 = idx_s
        num2 = idx_e
#else
        ! num1 = 1
        ! num2 = nb_part
        num1 = idx_s
        num2 = idx_e
        print *, ' idx_s=', idx_s, ' idx_e=', idx_e
#endif

        print *, 'save lat', lat_out(1:nb_part)
        print *, 'save flag', flag_out(1:nb_part)
        print *, 'save NUMBER', nb_out(1:nb_part)

        CALL ionc4_write_trajt(file_out, 'latitude',  lat_out(1:nb_part),num1,num2,0,REAL(dg_valmanq_io,kind=out))
        CALL ionc4_write_trajt(file_out, 'longitude', lon_out(1:nb_part),num1,num2,0,REAL(dg_valmanq_io,kind=out))
        CALL ionc4_write_trajt(file_out, 'DEPTH',     zpos_out(1:nb_part),num1,num2,0,-fillval)
        CALL ionc4_write_trajt(file_out, 'H0',        h0pos_out(1:nb_part),num1,num2,0,-fillval)
        CALL ionc4_write_trajt(file_out, 'NUM',       num_out(1:nb_part),num1,num2,0,0)
        CALL ionc4_write_trajt(file_out, 'flag',      flag_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'TEMP',      temp_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'STAGE',     stage_out(1:nb_part),num1,num2,0,-1)
        CALL ionc4_write_trajt(file_out, 'SIZE',      size_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'NUMBER',    nb_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'DENSITY',   dens_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'DRATE',     Drate_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_traj(file_out,  'DAYBIRTH',  dateo_out(1:nb_part),num1,num2,REAL(dg_valmanq_io,kind=out))
        CALL ionc4_write_trajt(file_out, 'AGE',       age_out(1:nb_part),num1,num2,0,-1)
        CALL ionc4_write_trajt(file_out, 'AGECLASS',  ageClass_out(1:nb_part),num1,num2,0,-1)

#ifdef IBM_SPECIES
        CALL ionc4_write_trajt(file_out, 'FOOD',      food_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'F',         f_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'EDEB',      E_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'HDEB',      H_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'RDEB',      R_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'GAM',       Gam_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'WEIGHT',    Wdeb_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'NEGGS',     Neggs_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_traj(file_out,  'YEARSPAWN', yearspawn_out(1:nb_part),num1,num2,0)
        CALL ionc4_write_traj(file_out,  'DAYSPAWN',  dayspawn_out(1:nb_part),num1,num2,0)
        CALL ionc4_write_traj(file_out,  'DAYJUV',    dayjuv_out(1:nb_part),num1,num2,0)
        CALL ionc4_write_traj(file_out,  'DENSPAWN',  denspawn_out(1:nb_part),num1,num2,fillval)
        CALL ionc4_write_trajt(file_out,  'ZOOM',      zoom_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'Death_DEB', deaddeb_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'Death_FISH', deadfishing_out(1:nb_part),num1,num2,0,fillval)
        CALL ionc4_write_trajt(file_out, 'Death_NAT',  deadnatural_out(1:nb_part),num1,num2,0,fillval)
#endif /*IBM_SPECIES*/
 
        ! To write the data on the disk and not loose data in case of run crash
        CALL ionc4_sync(file_out)
        DEALLOCATE(lat_out, lon_out, xpos_out, ypos_out, spos_out, zpos_out)
        DEALLOCATE(num_out, h0pos_out, flag_out, dens_out, temp_out, Drate_out)
        DEALLOCATE(size_out, dateo_out, stage_out, nb_out,age_out,ageClass_out)
#ifdef IBM_SPECIES
        DEALLOCATE(dayjuv_out, dayspawn_out, yearspawn_out, season_out)
        DEALLOCATE(Denspawn_out, food_out, Wdeb_out, zoom_out)
        DEALLOCATE(Gam_out, H_out, E_out, R_out, Neggs_out)!, NRJ_out)
        DEALLOCATE(f_out)
        DEALLOCATE(deaddeb_out, deadfishing_out, deadnatural_out)
#endif /*IBM_SPECIES*/

        patch%t_save = time + patch%dt_save*3600.0_rlg
        patch => patch%next
       
    END DO
 
  END SUBROUTINE ibm_save





  !!======================================================================
  FUNCTION selec_dome_or_asymp(l,alpha1,beta1)
    !&E-------------------------------------------------------------------
    !&E                 ***  ROUTINE ibm_save  ***
    !&E
    !&E ** Purpose : Selectivity dome or asymptotic when computing death by fishing
    !&E
    !&E ** Called by      : ibm_3d
    !&E ** External calls : 
    !&E
    !&E ** History :
    !&E       !  2023   (C. Menu) 
    !&E-------------------------------------------------------------------
    !&E Selectivity dome_or_asymptotic

    REAL(KIND=rsh), INTENT( in ) :: l,alpha1,beta1
    REAL(KIND=rsh)               :: selec_dome_or_asymp

    selec_dome_or_asymp = 1.0_rsh/( 1.0_rsh + exp(alpha1*(beta1-l)) )
  END FUNCTION selec_dome_or_asymp

  !!======================================================================
!&E-----------------------------------------------------------------------
!&E                 ***  ROUTINE generate_run_id  ***
!&E
!&E ** Purpose :
!&E     Generates a unique run identifier (RUN_ID) for each model execution.
!&E     The identifier is based on the system date, time, and clock counter
!&E     to ensure uniqueness across multiple simulation launches.
!&E
!&E ** Called by      : ibm_init
!&E ** External calls : date_and_time, system_clock
!&E
!&E ** Output:
!&E     run_id : character string like 'RUN_20251023_070234_000123456789'
!&E
!&E ** History :
!&E     2025-10-21  (C. Menu)  First version
!&E-----------------------------------------------------------------------
FUNCTION generate_run_id() RESULT(run_id)
  implicit none
  !-----------------------------------------------------------------------
  !> Purpose:
  !>   Generate a unique run identifier based on current date, time,
  !>   and system clock count, ensuring uniqueness across runs.
  !>
  !> Output:
  !>   run_id : character string like 'RUN_20251023_070234_000123456789'
  !-----------------------------------------------------------------------
  character(len=64) :: run_id             ! generated run identifier
  character(len=8)  :: date_str           ! date string: YYYYMMDD
  character(len=10) :: time_str           ! time string: HHMMSS.SS
  character(len=12) :: count_str          ! formatted system clock count
  integer           :: count              ! system clock count (for uniqueness)
  integer           :: count_rate, count_max
  !-----------------------------------------------------------------------
  !> Get current system date and time
  call date_and_time(date=date_str, time=time_str)

  !-----------------------------------------------------------------------
  !> Get system clock counter for additional uniqueness
  call system_clock(count=count, count_rate=count_rate, count_max=count_max)

  !-----------------------------------------------------------------------
  !> Convert the clock count to string safely
  write(count_str, '(I12.12)') count

  !-----------------------------------------------------------------------
  !> Build the RUN_ID string safely
  run_id = 'RUN_' // trim(date_str) // '_' // time_str(1:6) // '_' // trim(count_str)

END FUNCTION generate_run_id
    !!======================================================================

    !&E-----------------------------------------------------------------------
!&E                 ***  ROUTINE write_run_info  ***
!&E
!&E ** Purpose :
!&E     Creates or updates a 'FOIL.info' file in the working directory.
!&E     Stores the RUN_ID and optionally the last model step written.
!&E
!&E ** Called by      : ibm_init
!&E ** External calls : none
!&E
!&E ** History :
!&E     2025-10-21  (C. Menu)  First version
!&E-----------------------------------------------------------------------
SUBROUTINE write_run_info(run_id, last_step)
  implicit none
  character(len=*), intent(in) :: run_id
  integer, intent(in), optional :: last_step
  integer :: unit

  open(newunit=unit, file='FOIL.info', status='replace', action='write')

  write(unit, '(A,1X,A)') 'RUN_ID=', trim(run_id)
  if (present(last_step)) write(unit, '(A,I10)') 'LAST_STEP=', last_step

  close(unit)
END SUBROUTINE write_run_info
!&E-----------------------------------------------------------------------


!&E-----------------------------------------------------------------------
!&E                 ***  ROUTINE read_run_info  ***
!&E
!&E ** Purpose :
!&E     Reads an existing 'FOIL.info' file if available, to retrieve
!&E     the RUN_ID and (optionally) the last written model step.
!&E
!&E ** Called by      : ibm_init, ibm_3d
!&E ** External calls : generate_run_id
!&E
!&E ** Output:
!&E     run_id : character string like 'RUN_20251023_070234_000123456789'
!&E
!&E ** History :
!&E     2025-10-21  (C. Menu)  First version
!&E-----------------------------------------------------------------------
SUBROUTINE read_run_info(run_id, last_step)
  implicit none
  !-----------------------------------------------------------------------
  !> Declarations
  character(len=*), intent(out) :: run_id
  integer, intent(out), optional :: last_step
  integer :: unit, ios
  character(len=128) :: line
  !-----------------------------------------------------------------------
  !> Try to open existing FOIL.info file
  open(newunit=unit, file='FOIL.info', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'No existing FOIL.info found.'
  endif

  !-----------------------------------------------------------------------
  !> Read RUN_ID
  read(unit, '(A)', iostat=ios) line
  if (ios == 0) then
    read(line, '(6X,A)', iostat=ios) run_id
  endif

  !-----------------------------------------------------------------------
  !> Read LAST_STEP if available
  if (present(last_step)) then
    read(unit, '(A)', iostat=ios) line
    if (ios == 0) then
      read(line, '(10X,I10)', iostat=ios) last_step
      if (ios /= 0) last_step = 0
    else
      last_step = 0
    endif
  endif

  close(unit)
END SUBROUTINE read_run_info
!&E-----------------------------------------------------------------------

#endif /* DEB_IBM */
END MODULE
