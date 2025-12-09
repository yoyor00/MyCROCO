MODULE trajinitsave
    !!======================================================================
    !!                   ***  MODULE trajinitsave  ***
    !! Initialisation and saving of particles position
    !! Read circle or rectangle boundaries to calculate initial position, 
    !! or directly real position in netcdf file
    !! Save in a netcdf file
    !&E       !  01-2012 (M. Huret) separation from traject3d
    !&E
    !&E Description :
    !&E    traj_init3d
    !&E    traj_save3d
    !&E
    !&E ** History :
    !&E       !  2011-10 (T. Odaka, M. Huret, V. Garnier) Introduction of MPI
    !&E       !  2012-01 (M. Huret) Read netcdf
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&3
    !!======================================================================
#include "cppdefs.h"
#include "toolcpp.h"
#if defined LAGRANGIAN || defined DEB_IBM

    !! * Modules used
#ifdef MPI
    USE mpi
#endif
    USE module_lagrangian !, ONLY : LLm,MMm, mynode, NNODES
    USE comtraj,    ONLY : rsh,rlg,riosh,lchain,rg_valmanq_io,dg_valmanq_io, &
                           imin,imax,jmin,jmax,kmax
    USE ionc4,      ONLY : ionc4_createfile_traj, ionc4_createvar_traj,   &
                           ionc4_write_time,ionc4_sync,ionc4_close,       &
                           ionc4_write_trajt,ionc4_read_dimtraj,          &
                           ionc4_read_trajt,ionc4_openr,ionc4_init

    IMPLICIT NONE
    PRIVATE

    !! * Accessibility
    PUBLIC   :: LAGRANGIAN_init          ! routine called by LAGRANGIAN_init_main, ibm_init
    PUBLIC   :: traj_save3d              ! routine called by LAGRANGIAN_update
    PUBLIC   :: indices_loc2glob
    PUBLIC   :: init_patch

    !! * Shared module variables
    !! * Private variables

 CONTAINS
 
   !!======================================================================

   SUBROUTINE indices_loc2glob(in_start, nb_part, out_start, out_end)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE indices_loc2glob  ***
    !&E
    !&E ** Purpose : Gather indexes of all particles in a patch, especially for MPI purpose
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : LAGRANGIAN_init, traj_save3d, ibm_save, ibm_3d
    !&E
    !&E ** External calls : exchange_vectcpu_int
    !&E
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E       !  2015-01  (M. Honnorat) IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used

#ifdef MPI
    USE toolmpi, ONLY : exchange_vectcpu_int
#endif

    !! * Arguments
    INTEGER, INTENT(in)  :: in_start
    INTEGER, INTENT(in)  :: nb_part
    INTEGER, INTENT(out) :: out_start
    INTEGER, INTENT(out) :: out_end

    !! * Local declarations
#ifdef MPI 
    INTEGER, DIMENSION(0:NNODES-1)   :: nb_part_percpu
#endif

   !!----------------------------------------------------------------------
   !! * Executable part

    out_start = in_start
#ifdef MPI 
    ! Init nb_part_percpu
    nb_part_percpu = 0 
    nb_part_percpu(mynode) = nb_part
    CALL exchange_vectcpu_int(nb_part_percpu)

    IF ( mynode /= 0 ) THEN
       out_start = out_start + SUM(nb_part_percpu(0:mynode-1))
    END IF
#endif
    out_end = out_start-1 + nb_part

   END SUBROUTINE indices_loc2glob



 !!======================================================================
 SUBROUTINE init_patch(patch, nb_part)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE init_patch  ***
    !&E
    !&E ** Purpose : Routine to init a patch object
    !&E
    !&E ** Called by : LAGRANGIAN_init, ibm_3d
    !&E
    !&E ** External calls : ADD_ALL_MPI_INT
    !&E
    !&E ** History :
    !&E       !  2015-01 (M. Honnorat) IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj, ONLY : type_patch
#ifdef MPI
    USE toolmpi, ONLY : ADD_ALL_MPI_INT
#endif
    !! * Arguments
    TYPE(type_patch), INTENT(inout)  :: patch
    INTEGER,          INTENT(in)     :: nb_part

    !! * Local declarations
    INTEGER                          :: nb_part_total
    INTEGER                          :: nb_part_alloc

    !!----------------------------------------------------------------------
    !! * Executable part
    IF ( nb_part > 0 ) THEN
        nb_part_alloc = patch%nb_part_batch*( INT(nb_part/patch%nb_part_batch) + 1)
    ELSE
        nb_part_alloc = 0
    ENDIF
    patch%nb_part_alloc = nb_part_alloc

    ! Allocation of particles of each patch
    ALLOCATE( patch%particles(patch%nb_part_alloc) )

    patch%particles(:) = patch%init_particle
    nb_part_total = nb_part

#ifdef MPI
    CALL ADD_ALL_MPI_INT(nb_part_total)
#endif
    patch % nb_part_total = nb_part_total
    patch % nb_part_max   = nb_part_total    ! Will be overwritten in ibm_init

 END SUBROUTINE init_patch



 !!======================================================================
 SUBROUTINE LAGRANGIAN_init(Istr,Iend,Jstr,Jend)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE traj_init3d  ***
    !&E
    !&E ** Purpose : Read file traject.dat or ibm.dat to initialize trajectories variables.
    !&E              There are three type of inputs patches : circle patch, rectangular patch
    !&E              and a netcdf patch. Depending on Lagrangian or Foil to fulfill patch info
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by :  LAGRANGIAN_init_main
    !&E
    !&E ** External calls : h0int,xeint,loc_h0,update_htot,update_wz,compute_dsig_dcuds
    !&E                     ztosiggen,hc_sigint,lonlat2ij,tool_latlon2i,tool_latlon2j
    !&E                     set_htot_bc,define_pos,patch_list_append,init_mpi_type_particle
    !&E                     init_patch,indices_loc2glob
    !&E                     
    !&E
    !&E ** History :
    !&E       !  12-2006 (M. Chiffle, V. Garnier)  Original code
    !&E       !  05-2010 (M. Sourisseau)  Introduction of derivate type
    !&E       !  01-2011 (M. Huret)  Adaptation to siggen
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(Istrm2:Iendp2,Jstrm2:Jendp2)
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E       !  2024    (D. Gourves) Add reading FOIL patch data for individual variability for DEB-IBM
    !&E---------------------------------------------------------------------
    !! * Modules used

    USE module_lagrangian !,  ONLY : pi,stdout,start_time,time_end, &
                         !          Eradius,h,latr,lonu,latv,sc_r,sc_w,N
    USE trajectools, ONLY : h0int,xeint,loc_h0,update_htot,update_wz,compute_dsig_dcuds,&
                            ztosiggen,hc_sigint,lonlat2ij,tool_latlon2i,tool_latlon2j,&
                            set_htot_bc, define_pos
#ifdef MPI
    USE toolmpi, ONLY : MPI_glob2loc,MPI_loc2glob
    USE comtraj, ONLY : init_mpi_type_particle
#endif
    USE comtraj, ONLY : patch_list_append,patches,type_patch,file_trajec,itypetraj,ndtz
#ifdef DEB_IBM
    USE comtraj, ONLY : ibm_restart
#endif
    USE comtraj, ONLY : dsigu,dsigw,kmax,ierrorlog,iscreenlog
    USE comtraj, ONLY : lonwest,latsouth,dlonr,dlatr,htx,hty,wz
    USE comtraj, ONLY : type_position

    !! * Arguments
    INTEGER,INTENT(in) :: Istr,Iend,Jstr,Jend 

    !! * Local declarations
    ! For reading input file
    LOGICAL                                     :: ex,l_posit
    INTEGER                                     :: lstr, lenstr
    CHARACTER(LEN=5)                            :: comment
    CHARACTER(LEN=19)                           :: dateread,tool_sectodat

    TYPE(type_patch), POINTER                   :: new_patch, patch

    ! Time info of patches
    REAL(KIND=rlg)                              :: t_traj_beg,t_traj_end,dt_traj
    REAL(KIND=rlg)                              :: tool_datosec

    ! Indexes for loops
    CHARACTER(LEN=lchain)                       :: rec
    INTEGER                                     :: eof,i,j,k,nn,npa,m1,m2,kk,l

    ! For circle patch, center of the circle patch
    REAL(KIND=rsh)                              :: dxc,dyc 
    REAL(KIND=rsh)                              :: fi,radeg,ray_patch,phimin,phimax
    REAL(KIND=rsh)                              :: xpos_patch,ypos_patch
    REAL(KIND=rsh)                              :: rx,ry,dl

    ! To read data from rectangular patch
    REAL(KIND=rlg)                              :: gmin,gmax

    ! Position indexes for particles and their number 
    INTEGER                                     :: imin_patch, imax_patch,   &
                                                   jmin_patch, jmax_patch,   &
                                                   kmin_patch, kmax_patch
    INTEGER                                     :: istep_patch, jstep_patch, kstep_patch
    INTEGER                                     :: nb_patch, nb_part_nc, nb_part_intro, nbpart_patch
    CHARACTER(LEN=19)                           :: species
    INTEGER                                     :: nb_part

    ! Variables to fix horizontal and vertical position
    REAL(KIND=rsh)                              :: xtemp,ytemp,xe_lag,h0_lag
    REAL(KIND=rsh)                              :: d3,kint,spos,hc_sig_lag

    ! DEB-IBM and SPECIES 
#ifdef DEB_IBM
    INTEGER                                     :: age,ageClass,super,number_particle,stage
    REAL(KIND=rlg)                              :: size, density
#ifdef IBM_SPECIES
    REAL(KIND=rlg)                              :: E_deb,H_deb,R_deb,Gam_deb
#endif
#endif
    ! To read data from netcdf patch
    REAL(KIND=rlg), ALLOCATABLE, DIMENSION(:)   :: lon_nc,lat_nc,depth_nc

    !  parameters for interpolation at particle location
    REAL(KIND=rsh)                              :: px,py
    INTEGER                                     :: igg,idd,jbb,jhh,hlb,hlt,hrb,hrt
    INTEGER                                     :: idx_s,idx_e,ierr
    ! To save a local and global position of particle for MPI and Sequential compatibility
    TYPE(type_position)                         :: pos, pos1, pos2

    ! MPI parameters
    INTEGER                                     :: proc_mpi, ierr_mpi

    REAL(KIND=rlg),DIMENSION(5)                 :: buff_mpi
    
    NAMELIST/namtraj/file_trajec,itypetraj,ndtz

# include "compute_auxiliary_bounds.h"
    !!----------------------------------------------------------------------
    !! * Executable part

    !translate some variables from MARS to CROCO
    !-------------------------------------------
    iscreenlog = stdout
    ierrorlog  = stdout
    imin       = 1
    jmin       = 1
    imax       = LLm
    jmax       = MMm
    kmax       = N
    
    time_start = start_time 

! start by exchange lonwest and latsouth for all procs
#ifdef MPI
    if (ii==0 .and. jj==0 .and. Istr==1.and. Jstr==1) then
        lonwest  = lonr(1,1)
        latsouth = latr(1,1)
    endif
    CALL MPI_Bcast(lonwest,  1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
    CALL MPI_Bcast(latsouth, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
#else 
    lonwest  = lonr(1,1)
    latsouth = latr(1,1)
#endif

    !calcul of lon,lat resolution
    !----------------------------
#ifdef MPI
    IF ( MASTER ) THEN
        dlonr = lonr(2,1) - lonr(1,1)
        dlatr = latr(1,2) - latr(1,1)
    ENDIF
    CALL MPI_Bcast(dlonr, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
    CALL MPI_Bcast(dlatr, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
#else
    dlonr = lonr(2,1) - lonr(1,1)
    dlatr = latr(1,2) - latr(1,1)
#endif

    !init ionc4
    !----------
    CALL ionc4_init()

    !ALLOC VAR
    !----------
    CALL ALLOC_VAR()

    !Compute sigma parameters
    !------------------------
    CALL compute_dsig_dcuds()
    
    !update water level
    !------------------
    CALL update_htot(Istr,Iend,Jstr,Jend,IstrU,JstrV)
    CALL set_htot_bc(Istr,Iend,Jstr,Jend,IstrU,JstrV)

    !update wz
    !------------------
    CALL update_wz(Istr,Iend,Jstr,Jend)
#ifdef MPI
    CALL exchange_w3d_tile (Istr,Iend,Jstr,Jend, wz(START_2D_ARRAY,0))
#endif 

#ifdef LAGRANGIAN
    ! Open paratraj.dat file, given in croco.in file if LAGRANGIAN key is defined
    ! Otherwise, file is given in ibm_init subroutine and we skip this part of the code
    lstr = lenstr(lagname)
    OPEN(50,file=lagname(1:lstr),status='old',form='formatted',access='sequential')
    READ(50,namtraj)
#endif
   
    ! save into simu.log
    !-------------------
    IF_MPI (MASTER) THEN
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) '**************************************************'
        WRITE(iscreenlog,*) '***************** TRAJ_INIT3D.F90 ****************'
        WRITE(iscreenlog,*) '**************************************************'
        WRITE(iscreenlog,*) ' '
        WRITE(iscreenlog,*) 'fichier definissant les caracteristiques des trajectoires : ',trim(file_trajec)
    ENDIF_MPI

   ! Make sure the type of trajectory is correct
    INQUIRE(file=file_trajec, exist=ex)
    IF (.NOT. ex) THEN
        PRINT*, "Trajectory file '" // trim(file_trajec) // "' does not exist."
        PRINT*, "Check in 'paraspec.txt' or 'paraibm.txt' if you use key_ibm."
        PRINT*, "Simulation stopped."
        CALL_MPI MPI_FINALIZE(ierr_mpi)
        STOP
    END IF

    IF (itypetraj/=1 .AND. itypetraj/=2 .AND. itypetraj/=3) THEN
        PRINT*, "Type of trajectory is not defined correctly."
        PRINT*, "Must be 1 (circle patch), 2 (rectangle patch) or 3 (Netcdf)"
        PRINT*, "Simulation stopped."
        CALL_MPI MPI_FINALIZE(ierr_mpi)
        STOP
    END IF

    OPEN(49, file=file_trajec, form='formatted')
    comment='debut'
    DO WHILE ( comment /= '*****' )
        READ(49,'(a)',iostat=eof) comment
    END DO

    nb_patch = 0
    DO WHILE ( eof == 0 )
        READ(49,'(a)',iostat=eof) comment
        IF ( (comment == '*****') .AND. (eof == 0) ) nb_patch = nb_patch+1
    END DO
    IF_MPI (MASTER) THEN
        WRITE(iscreenlog,*) 'Number of patches for estimation of trajectories = ',nb_patch
    ENDIF_MPI

    ! Rewind
    REWIND(49)
    comment='debut'
    DO WHILE (comment /= '*****')
        READ(49,'(a)',iostat=eof) comment
    END DO

    ! Loop over patches in trajectory file
    DO npa= 1,nb_patch

        ! Create new patch data structure
        new_patch => patch_list_append(patches)

        ! Read name of patch
        READ(49,'(a)',iostat=eof) ! name_patch

        ! Read starting date of trajectory
        READ(49,'(a)',iostat=eof) dateread
        t_traj_beg = tool_datosec(dateread)
        IF( t_traj_beg < time_start .OR. t_traj_beg > time_end ) THEN
            IF_MPI (MASTER) THEN
                WRITE(iscreenlog,*) ' '
                WRITE(iscreenlog,*) 'WARNING : PATCH NUMBER : ',npa
                WRITE(iscreenlog,*) 'WARNING : Starting date ',trim(dateread),' is uncorrect (before 3D simulation starting).'
                WRITE(iscreenlog,*) 'WARNING : Consequently trajectory starts with departure of 3D run on ' &
                                    ,tool_sectodat(time_start)
            ENDIF_MPI
            t_traj_beg = time_start
        ENDIF

        ! Read ending date of trajectory
        READ(49,'(a)',iostat=eof) dateread
        t_traj_end = tool_datosec(dateread)
        IF( t_traj_end < time_start .OR. t_traj_end > time_end ) THEN
            IF_MPI (MASTER) THEN
                WRITE(iscreenlog,*) ' '
                WRITE(iscreenlog,*) 'WARNING : PATCH NUMBER : ',npa
                WRITE(iscreenlog,*) 'WARNING : Ending date ',trim(dateread),' is uncorrect (after 3D simulation ending).'
                WRITE(iscreenlog,*) 'WARNING : Consequently trajectory ends with run on ' &
                                    , tool_sectodat(time_end)
            ENDIF_MPI
            t_traj_end = time_end
        ENDIF

        ! Read time step for outputs
        READ(49,*,iostat=eof) dt_traj

        new_patch%t_beg   = t_traj_beg
        new_patch%t_end   = t_traj_end
        new_patch%t_save  = t_traj_beg
        new_patch%dt_save = dt_traj

        IF_MPI (MASTER) THEN
            WRITE(iscreenlog,*) 'PATCH NUMBER : ',npa, new_line(''),                                    &
                                '   trajectory from ' // trim(tool_sectodat(t_traj_beg)), new_line(''), &
                                '                to ' // trim(tool_sectodat(t_traj_end)), new_line(''), &
                                '   with a ', dt_traj, 'hours time step.'
        ENDIF_MPI

        ! Depending on itypetraj in paratraj or paraibm, initialise patches with good patch
        IF (itypetraj == 1) THEN

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Initialize circle patches
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Spatial extension of initial patch in indexes or long. lati. ?
            READ(49,*,iostat=eof) l_posit

            ! Read location of initial patch
            IF ( l_posit ) THEN
                READ(49,*,iostat=eof) imin_patch, jmin_patch
            ELSE
                READ(49,*,iostat=eof) gmin,phimin
                imin_patch = NINT(tool_latlon2i(gmin,phimin))
                jmin_patch = NINT(tool_latlon2j(gmin,phimin))
                IF (imin_patch < imin .OR. imin_patch > imax) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number', new_patch%id
                    imin_patch = MIN(MAX(imin_patch,imin),imax)
                    PRINT*,' its longitude is :', imin_patch
                END IF
                IF ( jmin_patch < jmin .OR. jmin_patch > jmax ) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number', new_patch%id
                    jmin_patch = MIN(MAX(jmin_patch,jmin),jmax)
                    PRINT*,' its latitude is :', jmin_patch
                END IF
            END IF
#ifdef MPI
            IF( imin_patch >= iminmpi .and. imin_patch <= imaxmpi .and. &
                jmin_patch >= jminmpi .and. jmin_patch <= jmaxmpi ) THEN
                !indexes from global grid to local one (local is the proc's grid)
                CALL MPI_glob2loc(imin_patch,jmin_patch)
                proc_mpi=mynode  
            ELSE
                imin_patch = -1
                jmin_patch = -1
                proc_mpi   = -1
            ENDIF
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
            IF (mynode == proc_mpi) then
                !we send the number of the proc who has the center of the patch
                CALL MPI_Send(proc_mpi,1,MPI_INTEGER,0,22,MPI_COMM_WORLD,ierr_mpi) 
            ELSE IF(mynode == 0) THEN 
                !let anyone receive this information
                CALL MPI_Recv(proc_mpi,1,MPI_INTEGER,MPI_ANY_SOURCE,22,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
            END IF
            ! then broadcast to everyone
            CALL MPI_Bcast(proc_mpi,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
#endif
            ! Read ray of initial patch
            READ(49,*,iostat=eof) ray_patch

            ! Read nbpart_patch
            READ(49,*,iostat=eof) nbpart_patch

            ! Read depth of initial patch (read as immersion in meters)
            READ(49,*,iostat=eof) kmin_patch,kmax_patch

            ! Read resolution depth of initial patch
            READ(49,*,iostat=eof) kstep_patch

            ! Correct if depth not set as imersion
            kstep_patch = ABS(kstep_patch)
            kmin_patch  = ABS(kmin_patch)
            kmax_patch  = ABS(kmax_patch)

            ! Number of particles set at each initial position
            READ(49,*,iostat=eof) nb_part_intro

            ! Type of vertical behavior (integer):
             ! itypevert = 0 if constant depth
             ! itypevert < 0 if no random walk (vertical advection only)
             ! itypevert > 0 if random walk (advection + diffusion)
             ! abs(itypevert) = 1 if no vertical swimming 
             ! abs(itypevert) > 1 if vertical swimming (larval behavior):
             !                     = 2 for nycthemeral migration
             !                     = 3 for ontogenic migration (sakina), ...
            READ(49,*,iostat=eof) new_patch%init_particle%itypevert

            ! Read output file
            READ(49,'(a)',iostat=eof) rec
            kk = index(rec,',|')
            IF (kk > 0 ) THEN
                new_patch%file_out = rec(1:kk-1)
            ELSE
                new_patch%file_out = rec
            END IF
            READ(49,*,iostat=eof)
            ! == End of file reading

            ! Estimate number of particle inside the circle patch
            ! mean spatial size of the cell 
            radeg = pi/180.0_rsh
            fi = REAL(latr(imin_patch,jmin_patch),rsh)
#ifdef MPI
            IF (imin_patch + iminmpi - 1 == imax) THEN
#else
            IF (imin_patch == imax) THEN
#endif
                dxc = Eradius*COS(fi*radeg)*REAL(lonu(imin_patch,  jmin_patch)  &
                                                -lonu(imin_patch-1,jmin_patch),rsh)*radeg
            ELSE
                dxc = Eradius*COS(fi*radeg)*REAL(lonu(imin_patch+1,jmin_patch)  &
                                                -lonu(imin_patch,  jmin_patch),rsh)*radeg
            END IF
#ifdef MPI
            IF (jmin_patch + jminmpi - 1 == jmax) THEN
#else
            IF (jmin_patch == jmax) THEN
#endif
                dyc = Eradius*REAL(latv(imin_patch,jmin_patch)        &
                                  -latv(imin_patch,jmin_patch-1),rsh)*radeg
            ELSE
                dyc = Eradius*REAL(latv(imin_patch,jmin_patch+1)      &
                                  -latv(imin_patch,jmin_patch),rsh)*radeg
            END IF

            nn = INT(sqrt(4/pi*nbpart_patch)) + 1
            dl = 2*ray_patch/nn

#ifdef MPI
            IF ( iminmpi <= imin_patch .AND. imin_patch <= imaxmpi .AND. &
                 jminmpi <= jmin_patch .AND. jmin_patch <= jmaxmpi ) THEN
#endif
                IF (h(imin_patch,jmin_patch) < kmin_patch) THEN
                    PRINT*,'Patch on land or too deep'
                    PRINT*,'Check patch number', new_patch%id, 'in file "',trim(file_trajec),'"'
                    PRINT*,'Simulation stopped.'
                    STOP
                END IF
#ifdef MPI 
            ENDIF
#endif

            !send all the patch parameters to all procs
            !------------------------------------------
#ifdef MPI
            CALL MPI_loc2glob(imin_patch,jmin_patch)
            dxc=REAL(dxc,rsh)
            ! store data in a tab
            buff_mpi(1) = fi; buff_mpi(2) = dxc; buff_mpi(3) = dyc
            buff_mpi(4) = imin_patch; buff_mpi(5) = jmin_patch
            !broadcast data to all procs
            CALL MPI_Bcast(buff_mpi,5, MPI_DOUBLE_PRECISION,proc_mpi, MPI_COMM_WORLD, ierr)
            ! get data from broadcast
            fi = buff_mpi(1); dxc = buff_mpi(2); dyc = buff_mpi(3)
            imin_patch = buff_mpi(4); jmin_patch = buff_mpi(5)
#endif  
            nb_part = 0
            DO j = 1,nn
            DO i = 1,nn
                rx = -ray_patch + REAL(i-1)*dl
                ry = -ray_patch + REAL(j-1)*dl
                IF ( sqrt(rx*rx + ry*ry) <= ray_patch ) THEN 
                    xtemp = imin_patch + rx/dxc
                    ytemp = jmin_patch + ry/dyc
#ifdef MPI
                    IF ( iminmpi <= NINT(xtemp) .AND. NINT(xtemp) <= imaxmpi .AND. &
                         jminmpi <= NINT(ytemp) .AND. NINT(ytemp) <= jmaxmpi ) THEN
#endif
                        pos1%xp = xtemp; pos1%yp = ytemp
                        CALL define_pos(pos1)
                        DO k=kmin_patch,kmax_patch,kstep_patch
                            IF ( h(NINT(pos1%idx_r),NINT(pos1%idy_r)) > k ) THEN
                                CALL loc_h0(pos1%idx_r,pos1%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                                xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                                h0_lag = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                                d3 = h0_lag + xe_lag
                                IF (d3 > k) nb_part = nb_part + nb_part_intro
                            END IF
                        END DO
#ifdef MPI
                    END IF
#endif
                END IF
            END DO 
            END DO
            CALL init_patch(new_patch, nb_part)
            CALL indices_loc2glob(0, nb_part, idx_s, idx_e)
            m2 = 0
            DO j = 1,nn
            DO i = 1,nn
                rx = -ray_patch + REAL(i-1)*dl
                ry = -ray_patch + REAL(j-1)*dl
                IF ( sqrt(rx*rx+ry*ry) <= ray_patch ) THEN
                    xtemp = imin_patch + rx/dxc
                    ytemp = jmin_patch + ry/dyc
                    pos%xp = xtemp ; pos%yp = ytemp
#ifdef MPI
                    IF( iminmpi <= NINT(xtemp) .AND. NINT(xtemp) <= imaxmpi .AND. &
                        jminmpi <= NINT(ytemp) .AND. NINT(ytemp) <= jmaxmpi ) THEN
#endif
                        CALL define_pos(pos)                ! Take local and global position of particle
                        DO k=kmin_patch,kmax_patch,kstep_patch
                            IF ( h(NINT(pos%idx_r),NINT(pos%idy_r)) > k ) THEN
                                m1 = m2 + 1
                                m2 = m2 + nb_part_intro
                                xpos_patch = imin_patch + rx/dxc
                                ypos_patch = jmin_patch + ry/dyc
                                new_patch%particles(m1:m2)%xpos = xpos_patch ! i index of initial location
                                new_patch%particles(m1:m2)%ypos = ypos_patch ! j index of initial location

                                pos2%xp = xpos_patch ; pos2%yp = ypos_patch
                                CALL define_pos(pos2)                        ! Use position type for sequential/MPI compatibility

                                CALL loc_h0(pos2%idx_r,pos2%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                                xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                                h0_lag = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                                d3 = h0_lag + xe_lag
                                IF (d3 > k) THEN
                                    new_patch%particles(m1:m2)%active = .TRUE.
                                    new_patch%particles(m1:m2)%zpos   = k      ! immersion depth
                                    kint = -float(k)+xe_lag ! switch from immersion to real z
                                    hc_sig_lag = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt) 
                                    CALL ztosiggen(kint,spos,xe_lag,h0_lag,hc_sig_lag)
                                    new_patch%particles(m1:m2)%spos = spos
                                    new_patch%particles(m1:m2)%hc   = hc_sig_lag
                                    new_patch%particles(m1:m2)%d3   = d3
                                    new_patch%particles(m1:m2)%h0   = h0_lag
                                    new_patch%particles(m1:m2)%xe   = xe_lag
                                    DO l=0,nb_part_intro-1
                                       new_patch%particles(m1+l)%num = idx_s + m1 + l
                                    ENDDO
                                ELSE
                                    m2 = m2 - nb_part_intro
                                END IF
                            END IF
                        END DO
#ifdef MPI
                    END IF
#endif
                END IF
            END DO
            END DO


        ELSEIF ( itypetraj == 2 ) THEN
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Initialize rectangular patches
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Spatial extension of initial patch in indexes or long. lati. ?
            READ(49,*,iostat=eof) l_posit

            ! Read spatial extension of initial patch
            IF (l_posit) THEN
                READ(49,*,iostat=eof) imin_patch, imax_patch, jmin_patch, jmax_patch
            ELSE
                READ(49,*,iostat=eof) gmin,gmax,phimin,phimax
                imin_patch = NINT(tool_latlon2i(gmin,phimin))
                imax_patch = NINT(tool_latlon2i(gmax,phimax))
                jmin_patch = NINT(tool_latlon2j(gmin,phimin))
                jmax_patch = NINT(tool_latlon2j(gmax,phimax))

                IF ( imin_patch < imin .OR. imin_patch > imax ) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number',new_patch % id
                    imin_patch = MIN(MAX(imin_patch,imin),imax)
                    PRINT*,' its western longitude is :', imin_patch
                END IF
                IF ( jmin_patch < jmin .OR. jmin_patch > jmax ) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number',new_patch % id
                    jmin_patch = MIN(MAX(jmin_patch,jmin),jmax)
                    PRINT*,' its southern latitude is :', jmin_patch
                END IF
                IF ( imax_patch < imin .OR. imax_patch > imax ) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number',new_patch % id
                    imax_patch=MIN(MAX(imax_patch,imin),imax)
                    PRINT*,' its eastern longitude is :', imax_patch
                END IF
                IF ( jmax_patch < jmin .OR. jmax_patch > jmax ) THEN
                    PRINT*,'LOCATION OUT OF THE DOMAIN'
                    PRINT*,'Have a look at file ',trim(file_trajec),' patch number',new_patch % id
                    jmax_patch=MIN(MAX(jmax_patch,jmin),jmax)
                    PRINT*,' its northern latitude is :', jmax_patch
                END IF
            END IF
            ! Read spatial dispersion of particles inside initial patch
            READ(49,*,iostat=eof) istep_patch, jstep_patch

            ! Read depth of initial patch (read as immersion in meters)
            READ(49,*,iostat=eof) kmin_patch,kmax_patch

            ! Read resolution depth of initial patch
            READ(49,*,iostat=eof) kstep_patch

            ! Number of particles set at each initial position
            READ(49,*,iostat=eof) nb_part_intro

            ! Type of vertical behavior (integer):
                ! itypevert = 0 if constant depth
                ! itypevert < 0 if no random walk (vertical advection only)
                ! itypevert > 0 if random walk (advection + diffusion)
                ! abs(itypevert) = 1 if no vertical swimming 
                ! abs(itypevert) > 1 if vertical swimming (larval behavior):
                !                     = 2 for nycthemeral migration
                !                     = 3 for ontogenic migration (sakina), ...
            READ(49,*,iostat=eof) new_patch % init_particle % itypevert

            ! Read output file
            READ(49,'(a)',iostat=eof) rec
            kk = index(rec,',|')
            IF (kk > 0 ) THEN
                new_patch%file_out = rec(1:kk-1)
            ELSE
                new_patch%file_out = rec
            END IF
            READ(49,*,iostat=eof)
            ! == End of file reading

            ! Estimate number of particle inside the rectangular patch
            nb_part = 0
            DO j = MAX0(Jstr,jmin_patch),MIN0(Jend,jmax_patch),jstep_patch
            DO i = MAX0(Istr,imin_patch),MIN0(Iend,imax_patch),istep_patch
#ifdef MPI 
            IF ( iminmpi <= i .AND. i <= imaxmpi .AND. jminmpi <= j .AND. j <= jmaxmpi ) THEN
#endif
            pos%xp = i ; pos%yp = j
            CALL define_pos(pos)
            DO k = kmin_patch,kmax_patch,kstep_patch
                IF (h(NINT(pos%idx_r),NINT(pos%idy_r)) > k) THEN
                    CALL loc_h0(pos%idx_r,pos%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                    xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                    h0_lag = h0int(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                    d3 = h0_lag + xe_lag
                    IF (d3 > k) THEN
                        nb_part = nb_part + nb_part_intro
                    END IF
                END IF
            END DO
#ifdef MPI
            ENDIF
#endif
            END DO
            END DO

            CALL init_patch(new_patch, nb_part)
            CALL indices_loc2glob(0, nb_part, idx_s, idx_e)

            ! Place particle at their location 
            m2 = 0
            DO j = MAX0(Jstr,jmin_patch),MIN0(Jend,jmax_patch),jstep_patch
            DO i = MAX0(Istr,imin_patch),MIN0(Iend,imax_patch),istep_patch
#ifdef MPI 
            IF ( iminmpi <= i .AND. i <= imaxmpi .AND. jminmpi <= j .AND. j <= jmaxmpi ) THEN
#endif
            pos1%xp = i ; pos1%yp = j
            CALL define_pos(pos1)
            DO k = kmin_patch,kmax_patch,kstep_patch
                IF (h(NINT(pos1%idx_r),NINT(pos1%idy_r)) > k) THEN
                    m1 = m2 + 1
                    m2 = m2 + nb_part_intro
                    new_patch%particles(m1:m2)%xpos = pos1%idx_r   ! position at initial location    
                    new_patch%particles(m1:m2)%ypos = pos1%idy_r  
                    ! total depth at particle s location
                    CALL loc_h0( pos1%idx_r,pos1%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                    xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                    h0_lag = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                    d3 = h0_lag + xe_lag
                    IF (d3 > k) THEN
                        new_patch%particles(m1:m2)%active = .TRUE.
                        new_patch%particles(m1:m2)%zpos   = k      ! immersion depth
                        kint = -float(k) + xe_lag ! switch from immersion to real z
                        hc_sig_lag = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                        CALL ztosiggen(kint,spos,xe_lag,h0_lag,hc_sig_lag)
                        new_patch%particles(m1:m2)%spos = spos
                        new_patch%particles(m1:m2)%hc   = hc_sig_lag
                        new_patch%particles(m1:m2)%d3   = d3
                        new_patch%particles(m1:m2)%h0   = h0_lag
                        new_patch%particles(m1:m2)%xe   = xe_lag
                        DO l = 0,nb_part_intro-1
                            new_patch%particles(m1+l)%num = idx_s + m1 + l
                        ENDDO
                    ELSE
                        m2 = m2 - nb_part_intro
                    END IF
                END IF
            END DO
#ifdef MPI
            ENDIF
#endif
            END DO
            END DO



        ELSEIF ( itypetraj == 3 ) THEN
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Initialize netcdf patches
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Read name of input netcdf file
            READ(49,'(a)',iostat=eof) rec
            kk = index(rec,',|')
            IF (kk > 0 ) THEN
                new_patch%file_inp = rec(1:kk-1)
            ELSE
                new_patch%file_inp = rec
            END IF

            ! Read output file
            READ(49,'(a)',iostat=eof) rec
            kk = index(rec,',|')
            IF (kk > 0 ) THEN
                new_patch%file_out = rec(1:kk-1)
            ELSE
                new_patch%file_out = rec
            END IF

            ! Number of particles set at each initial position
            READ(49,*,iostat=eof) nb_part_intro

            ! Type of vertical behavior (integer):
                ! itypevert = 0 if constant depth
                ! itypevert < 0 if no random walk (vertical advection only)
                ! itypevert > 0 if random walk (advection + diffusion)
                ! abs(itypevert) = 1 if no vertical swimming 
                ! abs(itypevert) > 1 if vertical swimming (larval behavior):
                !                     = 2 for nycthemeral migration
                !                     = 3 for ontogenic migration (sakina), ...
            READ(49,*,iostat=eof) new_patch%init_particle%itypevert

#ifdef DEB_IBM
            ! Read some parameters if DEB_IBM module is used from init file
            ! Done here because starting values are given in patch file which 
            ! is read in this routine
            READ(49,'(a)',iostat=eof) species
            READ(49,*, iostat=eof) number_particle
            READ(49,*, iostat=eof) stage
            READ(49,*, iostat=eof) size
            READ(49,*, iostat=eof) super
            READ(49,*, iostat=eof) density
            READ(49,*, iostat=eof) age
            READ(49,*, iostat=eof) ageclass
#ifdef IBM_SPECIES
            READ(49,*, iostat=eof) H_deb
            READ(49,*, iostat=eof) E_deb
            READ(49,*, iostat=eof) R_deb
            READ(49,*, iostat=eof) Gam_deb
            new_patch % species = species
#endif
#endif
            READ(49,*,iostat=eof)
            ! == End of file reading

            ! Estimate/correct number of particle inside the patch

            ! Open input file and read number of particles
            CALL ionc4_openr(trim(new_patch % file_inp),.false.)
            CALL ionc4_read_dimtraj(trim(new_patch % file_inp), nb_part_nc)
         

            ALLOCATE( lon_nc(nb_part_nc), lat_nc(nb_part_nc), depth_nc(nb_part_nc) )

            ! read lat,lon,depth of particles in file
            CALL ionc4_read_trajt(trim(new_patch%file_inp), "longitude",lon_nc,  1,nb_part_nc,1)
            CALL ionc4_read_trajt(trim(new_patch%file_inp), "latitude", lat_nc,  1,nb_part_nc,1)
            CALL ionc4_read_trajt(trim(new_patch%file_inp), "DEPTH",    depth_nc,1,nb_part_nc,1)

            nb_part = 0
            DO nn = 1,nb_part_nc            

                IF ( depth_nc(nn) < 0.0_rsh ) THEN
                    WRITE(ierrorlog,*) 'Function INIT_TRAJ : depth of particle has to be > 0'
                    WRITE(ierrorlog,*) 'Check the netcdf traj file :', new_patch % file_inp
                    WRITE(ierrorlog,*) 'The simulation is stopped'
                    CALL_MPI MPI_FINALIZE(ierr_mpi)
                    STOP
                END IF

                xtemp = tool_latlon2i(lon_nc(nn),lat_nc(nn))
                ytemp = tool_latlon2j(lon_nc(nn),lat_nc(nn))
#ifdef MPI
                IF (iminmpi <= NINT(xtemp) .AND. NINT(xtemp) <= imaxmpi .AND. &
                    jminmpi <= NINT(ytemp) .AND. NINT(ytemp) <= jmaxmpi) THEN
#endif
                    pos1%xp = xtemp; pos1%yp = ytemp
                    CALL define_pos(pos1)
                    IF ( h(NINT(pos1%idx_r),NINT(pos1%idy_r)) > depth_nc(nn) ) THEN
                        CALL loc_h0(pos1%idx_r,pos1%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                        xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                        h0_lag = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                        d3     = h0_lag + xe_lag
                        IF (d3 > depth_nc(nn)) THEN   ! patch tempo pour restart (martin)
                            nb_part = nb_part + nb_part_intro
                        END IF
                   END IF
#ifdef MPI
                END IF
#endif
            END DO

            CALL init_patch(new_patch, nb_part)
            CALL indices_loc2glob(0, nb_part, idx_s, idx_e)

            m2 = 0
            DO nn = 1,nb_part_nc
                xtemp = tool_latlon2i(lon_nc(nn),lat_nc(nn))
                ytemp = tool_latlon2j(lon_nc(nn),lat_nc(nn))
                pos%xp = xtemp ; pos%yp = ytemp
#ifdef MPI
                IF (iminmpi <= NINT(xtemp) .AND. NINT(xtemp) <= imaxmpi .AND. &
                    jminmpi <= NINT(ytemp) .AND. NINT(ytemp) <= jmaxmpi) THEN
#else
                IF ( xtemp>=Istr .and. xtemp<=Iend .and. ytemp>=Jstr .and. ytemp<=Jend ) THEN
#endif
                    CALL define_pos(pos)

                    IF ( h(NINT(pos%idx_r),NINT(pos%idy_r)) > depth_nc(nn) ) THEN
                        m1 = m2 + 1
                        m2 = m2 + nb_part_intro
                        new_patch%particles(m1:m2)%xpos = tool_latlon2i(lon_nc(nn),lat_nc(nn))
                        new_patch%particles(m1:m2)%ypos = tool_latlon2j(lon_nc(nn),lat_nc(nn))

                        ! total depth at particle s location
                        CALL loc_h0( pos%idx_r, pos%idy_r,  &
                                     px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                        xe_lag = xeint(zeta(:,:,nstp),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                        h0_lag = h0int(   px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                        d3 = h0_lag + xe_lag
                        IF ( d3 < depth_nc(nn) ) depth_nc(nn) = d3 - 0.1_rsh ! patch tempo pour restart (martin)
                        IF ( d3 > depth_nc(nn) ) THEN
                            new_patch%particles(m1:m2)%active = .TRUE.
                            new_patch%particles(m1:m2)%zpos = depth_nc(nn)   ! immersion depth
                            kint = -depth_nc(nn) + xe_lag ! switch from immersion to real z
                            hc_sig_lag = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
                            CALL ztosiggen(kint,spos,xe_lag,h0_lag,hc_sig_lag)
                            new_patch%particles(m1:m2)%spos = spos
                            new_patch%particles(m1:m2)%hc   = hc_sig_lag
                            new_patch%particles(m1:m2)%d3   = d3
                            new_patch%particles(m1:m2)%h0   = h0_lag
                            new_patch%particles(m1:m2)%xe   = xe_lag
                            DO l = 0,nb_part_intro-1
                                new_patch%particles(m1+l)%num = idx_s + m1 + l
                            ENDDO
                        ELSE
                            m2 = m2 - nb_part_intro
                        END IF
                    END IF
                END IF
            END DO
            DEALLOCATE(lon_nc,lat_nc,depth_nc)

            ! close netcdf file
            CALL ionc4_close(new_patch%file_inp)
            
#ifdef DEB_IBM  
            IF (.not. ibm_restart) THEN
                DO nn = 1,new_patch%nb_part_alloc        
                    ! Init some variables from ibm.dat file for fish
                    new_patch%particles(nn)%super    = super
                    new_patch%particles(nn)%stage    = stage
                    new_patch%particles(nn)%size     = size
                    new_patch%particles(nn)%density  = density
                    new_patch%particles(nn)%age      = age
                    new_patch%particles(nn)%ageClass = ageClass
#ifdef IBM_SPECIES
                    new_patch%particles(nn)%H   = H_deb
                    new_patch%particles(nn)%E   = E_deb
                    new_patch%particles(nn)%R   = R_deb
                    new_patch%particles(nn)%Gam = Gam_deb

                ENDDO
            ENDIF
#endif
#endif
        END IF  ! end test on itypetraj

    END DO  ! loop on patches

    CLOSE(49)

    ! Complete MPI initialisation
    IF_MPI (MASTER) THEN
        INQUIRE(file=file_trajec,exist=ex)
        IF (ex) THEN
            patch => patches%first
            write(iscreenlog,*) 'VERTICAL COMPONENT OF THE DISPLACEMENT:  '
            DO npa = 1,patches%nb
                write(iscreenlog,*) ' ndtz : ', ndtz
                IF     ( patch%init_particle%itypevert == 0 ) THEN
                    WRITE(iscreenlog,*) ' Patch number',npa,': Trajectories at constant depth ' 
                ELSEIF ( patch%init_particle%itypevert < 0  ) THEN
                    WRITE(iscreenlog,*) ' Patch number',npa,': Advection only (no random walk) '
                ELSE
                    WRITE(iscreenlog,*) ' Patch number',npa,': With random walk (Advection + Diffusion) '
                ENDIF
                patch => patch%next
            END DO
        ENDIF
    ENDIF_MPI

    CALL_MPI init_mpi_type_particle


#ifdef LAGRANGIAN
    ! Save initialization only if LAGRANGIAN. 
    ! If we save here when DEB-IBM is activated, we will create a file with not 
    ! enough variables inside, which will create an error while calling ibm_save
    CALL traj_save3d
#endif
  
 END SUBROUTINE LAGRANGIAN_init 



 !!======================================================================
 SUBROUTINE traj_save3d
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE traj_save3d  ***
    !&E
    !&E ** Purpose : Save trajectories
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : LAGRANGIAN_init, LAGRANGIAN_update 
    !&E
    !&E ** External calls : tool_ind2lat,tool_ind2lon, indices_loc2glob
    !&E                     ionc4 library
    !&E
    !&E ** History :
    !&E       !  2006-12 (M. Chiffle, V. Garnier)  Original code
    !&E       !  2010-05 (M. Sourisseau, M. Huret) Coupling with environment variables
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE module_lagrangian
    USE comtraj,     ONLY : patches, type_patch, type_particle
    USE trajectools, ONLY : tool_ind2lat,tool_ind2lon
    !! * Arguments

    !! * Local declarations
    CHARACTER(LEN=lchain)    :: file_out
    LOGICAL                                     :: l_out_nc4par
    INTEGER                                     :: npa,npart,num1,num2,nb_part,nb_part_nc,p
    REAL(kind=rsh)                              :: fillval
    INTEGER,        ALLOCATABLE, DIMENSION(:)   :: num_out
    REAL(KIND=rlg), ALLOCATABLE, DIMENSION(:)   :: lat_out,lon_out
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:)   :: xpos_out,ypos_out,spos_out,zpos_out,h0pos_out,flag_out
    TYPE(type_patch),    POINTER                :: patch
    TYPE(type_particle), POINTER                :: particle
    INTEGER                                     :: idx_s, idx_e

    !!----------------------------------------------------------------------
    !! * Executable part

    IF(rsh==8) THEN 
        fillval = dg_valmanq_io
    ELSE
        fillval = rg_valmanq_io
    ENDIF

#ifdef MPI
    l_out_nc4par = .true.
#else
    l_out_nc4par = .false.
#endif

    ! Loop on patches to initialise particles position
    patch => patches%first

    DO npa = 1,patches%nb

        IF ( (time < patch%t_save) .OR. (time > patch%t_end) ) THEN
            patch => patch%next
            CYCLE
        END IF

        file_out = trim(patch%file_out)
        ! --- Create output file
        IF ( .NOT. patch%file_out_init ) THEN   
            patch%file_out_init = .TRUE.
            nb_part_nc = patch%nb_part_total
            CALL ionc4_createfile_traj(file_out, nb_part_nc, 0, 0, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "latitude", "degrees_north","latitude",                     &
                                                fill_value=dg_valmanq_io, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "longitude","degrees_east", "longitude",                    &
                                                fill_value=dg_valmanq_io, l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "DEPTH","m","depth",fill_value=-fillval,l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "H0","m","h0",fill_value=-fillval,l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "NUM","","number of the particle",fill_value=0,l_out_nc4par=l_out_nc4par)
            CALL ionc4_createvar_traj(file_out, "flag","nbr","flag",fill_value=fillval,l_out_nc4par=l_out_nc4par)
        END IF    ! (.NOT. patch%file_out_init)

        nb_part = patch % nb_part_alloc
        ALLOCATE( lat_out(nb_part),  lon_out(nb_part) )
        ALLOCATE( xpos_out(nb_part), ypos_out(nb_part) )
        ALLOCATE( spos_out(nb_part), zpos_out(nb_part) )
        ALLOCATE( num_out(nb_part),  h0pos_out(nb_part) )
        ALLOCATE( flag_out(nb_part) )

        lat_out(:)  = dg_valmanq_io ; lon_out(:) = dg_valmanq_io
        xpos_out(:) = fillval ; ypos_out(:) = fillval ; spos_out(:) = fillval ; zpos_out(:) = -fillval
        h0pos_out(:)=-fillval ; flag_out(:) = fillval ; num_out(:)  = 0

        p = 0
        DO npart = 1,nb_part
            particle => patch % particles(npart)
            IF ( .NOT. particle % active ) CYCLE
            p = p+1
            xpos_out(p)  = particle % xpos
            ypos_out(p)  = particle % ypos
            lat_out(p)   = tool_ind2lat(particle%xpos, particle%ypos)
            lon_out(p)   = tool_ind2lon(particle%xpos, particle%ypos)
            spos_out(p)  = particle % spos
            zpos_out(p)  = particle % zpos ! for output as immersion
            h0pos_out(p) = particle % d3
            flag_out(p)  = particle % flag
            num_out(p)   = particle % num
        END DO

        CALL ionc4_write_time(file_out,0,time)

        nb_part = p
        CALL indices_loc2glob(1, nb_part, idx_s, idx_e)

#ifdef MPI
        num1 = idx_s
        num2 = idx_e
#else
        num1 = 1
        num2 = nb_part
#endif
        CALL ionc4_write_trajt(file_out, 'latitude',  lat_out(1:nb_part),num1,num2,0,dg_valmanq_io)
        CALL ionc4_write_trajt(file_out, 'longitude', lon_out(1:nb_part),num1,num2,0,dg_valmanq_io)
        CALL ionc4_write_trajt(file_out, 'DEPTH',     zpos_out(1:nb_part),num1,num2,0,-fillval)
        CALL ionc4_write_trajt(file_out, 'H0',        h0pos_out(1:nb_part),num1,num2,0,-fillval)
        CALL ionc4_write_trajt(file_out, 'NUM',       num_out(1:nb_part),num1,num2,0,0)
        CALL ionc4_write_trajt(file_out, 'flag',      flag_out(1:nb_part),num1,num2,0,fillval)

        ! To write the data on the disk and not loose data in case of run crash
        CALL ionc4_sync(file_out)
        DEALLOCATE(xpos_out, ypos_out, zpos_out, spos_out)
        DEALLOCATE(lat_out, lon_out, h0pos_out, flag_out, num_out)

#ifdef LAGRANGIAN
        ! Only if LAGRANGIAN, so we are not interfering with ibm_save when using DEB_IBM key
        patch%t_save = time + patch%dt_save*3600.0_rlg
#endif
        patch => patch%next
    END DO

  END SUBROUTINE traj_save3d



  !!======================================================================
  SUBROUTINE ALLOC_VAR
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE ALLOC_VAR  ***
    !&E
    !&E ** Purpose : Allocate arrays used in MARS for LAGRANGIAN over whole domain which
    !&E              don't exist in CROCO
    !&E
    !&E ** Called by : LAGRANGIAN_init
    !&E
    !&E ** External calls : tool_ind2lat,tool_ind2lon, indices_loc2glob
    !&E                     ionc4 library
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) Added for CROCO purpose
    !&E
    !&E---------------------------------------------------------------------
    USE module_lagrangian ! for GLOBAL_2D_ARRAY
    USE comtraj, ONLY : htx,hty,hc_sig,wz,dsigu,dsigw,dcusds,dcwsds

    IMPLICIT NONE

    ALLOCATE(htx(GLOBAL_2D_ARRAY))
    ALLOCATE(hty(GLOBAL_2D_ARRAY))
    ALLOCATE(hc_sig(GLOBAL_2D_ARRAY))
    ALLOCATE(wz(GLOBAL_2D_ARRAY,0:N))
    ALLOCATE(dsigu(kmax))
    ALLOCATE(dsigw(kmax))
    ALLOCATE(dcusds(kmax))
    ALLOCATE(dcwsds(kmax))

    htx    = 0.0
    hty    = 0.0
    hc_sig = 0.0
    wz     = 0.0
    dsigu  = 0.0
    dsigw  = 0.0
    dcwsds = 0.0
    dcusds = 0.0
 
  END SUBROUTINE ALLOC_VAR



#endif

END MODULE
