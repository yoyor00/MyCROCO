MODULE comtraj


#include "cppdefs.h"
#include "toolcpp.h"

  !!======================================================================
  !!                   ***  MODULE comtraj  ***
  !!
  !!======================================================================
#if defined LAGRANGIAN || defined DEB_IBM 

    IMPLICIT NONE
    PUBLIC

    !! * Accessibility
    PUBLIC    :: patch_list_append
    PUBLIC    :: get_patch
    PUBLIC    :: enlarge_patch, resize_patch
#if defined MPI
    PUBLIC    :: init_mpi_type_particle
#endif
 
    !! General parameters (use from comsubstance later)
    ! -------------------------------------------------------------------------
    ! Definition of rsh, rlg, riosh, riolg, lchain
    ! -------------------------------------------------------------------------
    INTEGER,PARAMETER                           :: riosh = 8, riolg = 8, rlg = 8, rsh = 8
    REAL(kind=rsh), PARAMETER                   :: valmanq       = 999.0
    REAL(kind=riosh),PARAMETER                  :: rg_valmanq_io = 999.0_riosh
    REAL(kind=riolg),PARAMETER                  :: dg_valmanq_io = -1.7e+38
    REAL(kind=riolg)                            :: time_start
    INTEGER,PARAMETER                           :: lchain = 200
    INTEGER                                     :: ierrorlog, iwarnlog, iscreenlog
    INTEGER                                     :: imin,imax,jmin,jmax,kmax
    INTEGER                                     :: jjulien

    !---------------------------------------------
    ! Definition of allocatable imported from MARS
    !---------------------------------------------
    ! To compute at each time step
    REAL(KIND=rsh),DIMENSION(:,:)  ,ALLOCATABLE :: htx,hty
    REAL(KIND=rsh),DIMENSION(:,:,:),ALLOCATABLE :: wz 
    ! To compute only once in trajinit
    REAL(KIND=rsh),DIMENSION(:)  ,  ALLOCATABLE :: dsigw,dsigu,dcusds,dcwsds
    REAL(KIND=rsh),DIMENSION(:,:),  ALLOCATABLE :: hc_sig
 
    !----------------------------------------
    ! Other variables to compute at beginning
    !----------------------------------------
    REAL(kind=rlg)                             :: lonwest,latsouth,dlonr,dlatr



    ! =====================================================================
    ! =====                                                           =====
    ! =====                     TYPE type_position                    =====
    ! =====                                                           =====
    ! =====================================================================
    TYPE,PUBLIC :: type_position

        !Position in the local grid as integer (depend on MPI) 
        INTEGER             ::  idx,idy
        !Position in the local grid as real (depend on MPI) 
        REAL(kind=rlg)      :: idx_r,idy_r
        !Position in the global domain
        REAL(kind=rlg)      :: xp,yp

    END TYPE type_position



    ! =====================================================================
    ! =====                                                           =====
    ! =====                     TYPE type_particle                    =====
    ! =====                                                           =====
    ! =====================================================================
    TYPE, PUBLIC :: type_particle

        LOGICAL :: active = .False. ! .True.  if particle actually active
#if defined MPI
        ! --- MPI managing
        INTEGER                 :: limitbye = 0     ! Specify the direction of the boundary crossing
                                                    ! 1=SW, 2=S, 3=SE, 4=E, 5=NE, 6=N, 7=NW, 8=N
#endif
        INTEGER                 :: itypevert        ! type of trajectory (z=cst, random walk...)
        INTEGER                 :: num = 0          ! to keep track of particles when save/restart 

        ! --- Location
        REAL(KIND=rsh)          :: xpos, ypos, spos, zpos
        REAL(KIND=rsh)          :: d3, h0, xe, hc
        REAL(KIND=rsh)          :: flag = 0.0_rsh   ! wet-drying flag

#ifdef DEB_IBM
        ! --- Population parameters
        REAL(KIND=rlg)          :: date_orig   ! date of release
        INTEGER                 :: stage    = 0
        INTEGER                 :: age      = 0
        INTEGER                 :: AgeClass = 0
        INTEGER                 :: Nbatch   = 0
        REAL(KIND=rsh)          :: Drate    = 0.0_rsh
        REAL(KIND=rsh)          :: temp, w, size, density, denspawn

        REAL(KIND=rsh)          :: super    ! Number of individuals in particle (superindividual)

#ifdef IBM_SPECIES
        ! DEB state variables and parameters (with default value)
        INTEGER                 :: dayjuv, yearspawn
        INTEGER                 :: dayspawn     = 500
        LOGICAL                 :: season       = .FALSE.
        INTEGER                 :: hmove        = 0
        REAL(KIND=rsh)          :: Hj,Hb,pAm,pMi,EG,vc,kap,Kx,Hp,TA,K,shapeb,lfactor,E0,Rfbatch,SF
        REAL(KIND=rsh)          :: zoom         = 0.0_rsh
        REAL(KIND=rsh)          :: E,L,H,R,Wdebd,NRJd
        REAL(KIND=rsh)          :: Gam          = 0.0_rsh
        REAL(KIND=rsh)          :: f            = 0.0_rsh
        REAL(KIND=rsh)          :: X            = 0.0_rsh
        REAL(KIND=rsh)          :: Wdeb         = 0.0_rsh
        REAL(KIND=rsh)          :: Neggs        = 0.0_rsh
        REAL(KIND=rsh)          :: Neggs_tot    = 0.0_rsh
        REAL(KIND=rsh)          :: NRJ          = 0.0_rsh
        REAL(KIND=rsh)          :: WV           = 0.0_rsh
        REAL(KIND=rsh)          :: WE           = 0.0_rsh
        REAL(KIND=rsh)          :: WR           = 0.0_rsh
        REAL(KIND=rsh)          :: WG           = 0.0_rsh
        REAL(KIND=rsh)          :: NRJ_V        = 0.0_rsh
        REAL(KIND=rsh)          :: NRJ_g        = 0.0_rsh
        REAL(KIND=rlg)          :: Death_DEB    = 0.0_rsh
        REAL(KIND=rlg)          :: Death_FISH   = 0.0_rsh
        REAL(KIND=rlg)          :: Death_NAT    = 0.0_rsh      
#endif
#endif
  END TYPE type_particle



    ! =====================================================================
    ! =====                                                           =====
    ! =====                       TYPE type_patch                     =====
    ! =====                                                           =====
    ! =====================================================================
    TYPE, PUBLIC :: type_patch

        INTEGER                                         :: id = -1
#ifdef DEB_IBM                  
        INTEGER                                         :: parent_id = -2
#endif
        INTEGER                                         :: nb_part_alloc = 0    ! Size of allocated data array for particles
        INTEGER                                         :: nb_part_total = 0    ! Sum of particles over all MPI domains
        INTEGER                                         :: nb_part_batch = 10   ! Size of batch for new allocations
        INTEGER                                         :: nb_part_max   = -1   ! Maximum allowed number of particles
        REAL(KIND=rlg)                                  :: t_beg,  t_end
        REAL(KIND=rlg)                                  :: t_save, dt_save
#ifdef DEB_IBM                          
        REAL(KIND=rlg)                                  :: t_spawn, dt_spawn    !  
        INTEGER                                         :: yearref              ! 
#endif                  
        CHARACTER(LEN=lchain)                           :: file_inp             ! Input data file name
        CHARACTER(LEN=lchain)                           :: file_out             ! Output NetCDF file name
        LOGICAL                                         :: file_out_init = .FALSE.               
#ifdef IBM_SPECIES                  
        CHARACTER(LEN=lchain)                           :: species              ! name of species
#endif    
        TYPE(type_particle)                             :: init_particle        ! Init values used for new particles
        TYPE(type_particle), ALLOCATABLE, DIMENSION(:)  :: particles      
        TYPE(type_patch), POINTER                       :: next => NULL()       ! Next patch in the list

    END TYPE type_patch



    ! =====================================================================
    ! =====                                                           =====
    ! =====                    TYPE type_patch_list                   =====
    ! =====                                                           =====
    ! =====================================================================
    TYPE, PUBLIC :: type_patch_list
        INTEGER                     :: nb = 0           ! Number of elements in the list
        TYPE(type_patch), POINTER   :: first => NULL()  ! First patch in the list
        TYPE(type_patch), POINTER   :: last  => NULL()  ! Last  patch inserted in the list
    END TYPE type_patch_list
  


  !----------------------------------------
  !! * Shared module variables
  

    INTEGER,PARAMETER                       :: nb_species = 2       ! Number of species in DEB_IBM, for further developments
                                                                    ! Species with index 1 : anchovy
                                                                    ! Species with index 2 : sardine

    TYPE(type_patch_list), PUBLIC           :: patches

    ! From paraibm or paratraj file
    CHARACTER(LEN=lchain),  PUBLIC          :: file_trajec                  ! name of configuration file
    INTEGER,                PUBLIC          :: itypetraj                    ! initialisation type (circle, rectangle,netcdf)
    INTEGER,                PUBLIC          :: ndtz                         ! number of time step division for vertical subloop

#ifdef DEB_IBM
    LOGICAL,                PUBLIC          :: ibm_restart                  ! Logical for ibm restart 

    REAL(KIND=rsh),         PUBLIC          :: struc_ad = 0.0_rsh           ! Density-dependance parameter
    REAL(KIND=rsh),         PUBLIC          :: struc_ad_dd_DEB = 0.0_rsh    ! Density-dependance parameter
#ifdef IBM_SPECIES
    INTEGER, DIMENSION(nb_species), PUBLIC  :: duration                     ! Duree de vie des individus selon leur espece
    ! namibmdeb namelist parameters from paraibm
    LOGICAL,                PUBLIC          :: debuse, F_Fix, frac_deb_death
    REAL(kind=rsh),         PUBLIC          :: ffix                        
    CHARACTER(LEN=lchain),  PUBLIC          :: file_food                   ! Name of input file for food
    CHARACTER(LEN=lchain),  PUBLIC          :: file_NBSS                   

    ! namibmfrc namelist parameters from paraibm
    CHARACTER(LEN=lchain), PUBLIC           :: fileanchovy,filesardine     ! File for anchovy and sardine global parameters
    CHARACTER(LEN=lchain), PUBLIC           :: catch_anc_bob               ! File for anchois fishing if fishing_strategy = "Catch"
    CHARACTER(LEN=lchain), PUBLIC           :: catch_sar_bob               ! File for sardine fishing if fishing_strategy = "Catch"
    CHARACTER(LEN=lchain), PUBLIC           :: fileprobadistrib_anc        ! Probability map of achovy  distribution
    CHARACTER(LEN=lchain), PUBLIC           :: fileprobadistrib_sar        ! Probability map of sardine distribution
    INTEGER,               PUBLIC           :: nbSizeClass_anc,nbSizeClass_sar
    REAL(KIND=rlg),        PUBLIC           :: sizemin_anc,sizemin_sar


    CHARACTER(LEN=lchain), PUBLIC           :: fishing_strategy
    ! Global variable for fishing
    REAL(kind=rlg), DIMENSION(20,12,nb_species) :: mat_catch    ! Storing monthly catches from two input files above
    REAL(KIND=rlg), DIMENSION(nb_species)       :: number_tot, weight_tot
    REAL(KIND=rlg), DIMENSION(nb_species)       :: biom_tot  = 0._rsh
    REAL(KIND=rlg), DIMENSION(nb_species)       :: Wdeb_mean = 0._rsh

    REAL(KIND=rlg)                              :: time2spawn

    ! DEB parameters from input file deb_parameter_species
    TYPE(type_particle)                             :: init_anchovy_egg     ! Init values used for new anchovy's particles 
    TYPE(type_particle)                             :: init_sardine_egg     ! Init values used for new sardine's particles 
#endif
#endif

#if defined MPI
    ! For MPI exchange of particles betwreen procs
    INTEGER,               PUBLIC   :: type_mpi_particle
    INTEGER,               PUBLIC   :: down_give, up_give, right_give, left_give
    INTEGER                         :: type_mpi_rsh, type_mpi_rlg
#endif

CONTAINS

  !!======================================================================

#if defined MPI

  SUBROUTINE init_mpi_type_particle
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE  init_mpi_type_particle  ***
    !&E
    !&E ** Purpose : Defines the MPI type for 'type_particle'
    !&E
    !&E       !  2014-12 (M. Honnorat) Refactor MPI routines
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE mpi

    IMPLICIT NONE
    !! * Local declarations
#ifdef DEB_IBM
#ifdef IBM_SPECIES
    INTEGER, PARAMETER   :: nb = 68     ! IBM_SPECIES and DEB_IBM
#else
    INTEGER, PARAMETER   :: nb = 25     ! DEB_IBM only
#endif
#else
    INTEGER, PARAMETER   :: nb = 13     ! key_MPI_2D only
#endif

    INTEGER,                        DIMENSION(nb)    :: old_types, block_lengths
    INTEGER(kind=MPI_ADDRESS_KIND), DIMENSION(nb)    :: displacements, addresses
    INTEGER                                          :: ierr_mpi, i
    TYPE(type_particle)                              :: particle
    !!----------------------------------------------------------------------
    !! * Executable part

    ! Create MPI type for particles 
    old_types = (/                                                                  &
        MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER,                         &
        type_mpi_rlg, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,       &
        type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh                      &
#ifdef DEB_IBM
        , type_mpi_rsh, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER,         &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh                                                &
#ifdef IBM_SPECIES
        , MPI_INTEGER, MPI_INTEGER, MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER,          &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh, type_mpi_rsh,     &
          type_mpi_rsh, type_mpi_rsh, type_mpi_rsh                                  &
#endif
#endif
    /)
    block_lengths = (/                                                              &
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1                                       &
#ifdef DEB_IBM
        ,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1                                         &
#ifdef IBM_SPECIES
        ,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1                 &
        ,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1                 &
        ,1, 1, 1                                                                    &
#endif
#endif
    /)

    i = 1
    CALL MPI_GET_ADDRESS(particle % active,    addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % limitbye,  addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % itypevert, addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % num,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % xpos,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % ypos,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % spos,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % zpos,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % d3,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % h0,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % xe,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % hc,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % flag,      addresses(i), ierr_mpi) ; i = i+1
#ifdef DEB_IBM
    CALL MPI_GET_ADDRESS(particle % date_orig, addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % stage,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % age,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % AgeClass,  addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Nbatch,    addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % w,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % size,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % density,   addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Drate,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % denspawn,  addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % temp,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % super,     addresses(i), ierr_mpi) ; i = i+1
#ifdef IBM_SPECIES
    CALL MPI_GET_ADDRESS(particle % dayjuv,    addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % yearspawn, addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % dayspawn,  addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % season,    addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % hmove,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Hb,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Hj,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % pAm,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % pMi,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % EG,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % vc,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % kap,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Kx,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % TA,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % K,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % shapeb,    addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % lfactor,   addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % E0,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Rfbatch,   addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % SF,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % zoom,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % E,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % L,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % H,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % R,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Wdebd,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % NRJd,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Gam,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % f,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % X,         addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Wdeb,      addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Neggs,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % Neggs_tot, addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % NRJ,       addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % WV,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % WE,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % WR,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % WG,        addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % NRJ_V,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % NRJ_g,     addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % DEATH_DEB,  addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % DEATH_FISH, addresses(i), ierr_mpi) ; i = i+1
    CALL MPI_GET_ADDRESS(particle % DEATH_NAT,  addresses(i), ierr_mpi) ; i = i+1
#endif
#endif

    displacements(:) = addresses(:) - addresses(1)

    CALL MPI_TYPE_CREATE_STRUCT(nb, block_lengths, displacements, old_types, type_mpi_particle, ierr_mpi)
    IF (ierr_mpi /= 0) THEN
        CALL MPI_FINALIZE(ierr_mpi) ; STOP
    ENDIF
        CALL MPI_TYPE_COMMIT(type_mpi_particle, ierr_mpi)

    END SUBROUTINE init_mpi_type_particle

#endif
  
 !!======================================================================

 FUNCTION patch_list_append ( patch_list ) RESULT(new_patch)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE  patch_list_append  ***
    !&E
    !&E ** Purpose : Append a new patch to 'patch_list'. A pointer to the
    !&E              newly created patch is returned.
    !&E
    !&E       !  2014-12 (M. Honnorat) Refactor MPI routines
    !&E---------------------------------------------------------------------
    !! * Arguments
    TYPE(type_patch_list), INTENT(inout)    :: patch_list

    !! * Local declarations
    TYPE(type_patch), POINTER               :: new_patch
    INTEGER, SAVE                           :: id = 0

    !!----------------------------------------------------------------------
    !! * Executable part

    ALLOCATE( new_patch )

    id = id + 1
    new_patch%id = id

    IF ( ASSOCIATED(patch_list % last) ) then
        ! the list is not empty, append 'new_patch' after the last one.
        patch_list % last % next => new_patch
    ELSE
        ! the list has just been initialized. Let 'new_patch' be the first one.
        patch_list % first => new_patch
    END IF

    ! anyway, for next time 'new_patch' will be the last one. 
    patch_list%last => new_patch
    patch_list%nb = patch_list%nb + 1

 END FUNCTION patch_list_append



  !!======================================================================
  FUNCTION get_patch(patch_list, n) RESULT(patch)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE  get_patch  ***
   !&E
   !&E ** Purpose : Returns a pointer to the nth patch from 'patch_list'.
   !&E
   !&E       !  2014-12 (M. Honnorat) Refactor MPI routines
   !&E---------------------------------------------------------------------
   !! * Arguments
   TYPE(type_patch_list), INTENT(in)    :: patch_list
   INTEGER,               INTENT(in)    :: n

   !! * Local declarations
   TYPE(type_patch), POINTER            :: patch
   INTEGER                              :: i

   !!----------------------------------------------------------------------
   !! * Executable part

   patch => NULL()

   IF ( (patch_list%nb <= 0) .OR. (n > patch_list%nb) )  RETURN
   
   patch => patch_list % first

   DO i = 1, patch_list % nb

      IF (i == n) RETURN
      patch => patch % next
   
   END DO

  END FUNCTION get_patch



 !!======================================================================
 SUBROUTINE enlarge_patch(patch, n)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE enlarge_patch  ***
    !&E
    !&E ** Purpose : Reallocate array patch % particles so that it can store
    !&E              n additional particles. The size n may be negative.
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : ex_traj_1d
    !&E
    !&E ** External calls :
    !&E
    !&E ** History :
    !&E       !  2015-01  (M. Honnorat) IBM upgrade
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    TYPE(type_patch), INTENT(inout)                 :: patch
    INTEGER,          INTENT(in)                    :: n     ! number of additional particles

    !! * Local declarations             
    INTEGER                                         :: old_size, new_size
    INTEGER                                         :: i, j
    TYPE(type_particle), ALLOCATABLE, DIMENSION(:)  :: tmp_array

    !!----------------------------------------------------------------------
    !! * Executable part

    ! Compute new_size
    new_size = patch%nb_part_alloc + patch%nb_part_batch * CEILING(REAL(n)/patch%nb_part_batch)
    old_size = patch%nb_part_alloc

    ! Save old data in resized array
    ALLOCATE(tmp_array(new_size))
    IF ( new_size > old_size ) THEN
        tmp_array(1:old_size) = patch%particles
        tmp_array(old_size+1:new_size) = patch%init_particle
    ELSE
        print*,'FIXME: enlarge_patch : NOT TESTED new_size < old_size', new_size, old_size
        i = 1
        DO j = 1,old_size
            IF ( patch % particles(j) % active ) THEN 
                tmp_array(i) = patch%particles(j)
                i = i+1
            END IF
        END DO
    END IF

    ! Update patch data with new array
    CALL MOVE_ALLOC(TO=patch%particles, FROM=tmp_array)
    patch % nb_part_alloc = new_size

 END SUBROUTINE enlarge_patch



 !!======================================================================
 SUBROUTINE resize_patch(patch, required_size)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE resize_patch  ***
    !&E
    !&E ** Purpose : Resize array patch % particles so that it can store
    !&E              #new_size# particles.
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : ibm_3d
    !&E
    !&E ** External calls :
    !&E
    !&E ** History :
    !&E       !  2015-01  (M. Honnorat) IBM upgrade
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    TYPE(type_patch), INTENT(inout)                 :: patch
    INTEGER,          INTENT(in)                    :: required_size

    !! * Local declarations
    INTEGER                                         :: old_size, new_size
    TYPE(type_particle), ALLOCATABLE, DIMENSION(:)  :: tmp_array

    !!----------------------------------------------------------------------
    !! * Executable part

    old_size = patch % nb_part_alloc
    new_size = patch % nb_part_batch * CEILING(REAL(required_size)/patch % nb_part_batch)
    
    IF ( old_size >= new_size ) RETURN     ! Nothing to do
    
    ! Save old data in resized array
    ALLOCATE(tmp_array(new_size))
    tmp_array(1:old_size)           = patch % particles
    tmp_array(old_size+1:new_size)  = patch % init_particle

    ! Update patch data with new array
    CALL MOVE_ALLOC(TO=patch % particles, FROM=tmp_array)
    patch % nb_part_alloc           = new_size

 END SUBROUTINE resize_patch
 !!======================================================================

#endif

END MODULE
