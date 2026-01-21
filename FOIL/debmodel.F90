MODULE debmodel

!!======================================================================
!!                   ***  MODULE DEB  ***
!!             DEB anchois/sardine (larves et adultes)
!!             Martin 09/2011
!!======================================================================

#include "cppdefs.h"
#include "toolcpp.h" 
#ifdef MPI
   USE mpi
#endif

#if defined IBM_SPECIES

    USE module_ibm         ! time,sc_w,h
    USE comtraj, ONLY      : kmax,rsh,rlg,lchain,valmanq,               &
                             type_particle,type_patch,patches

    IMPLICIT NONE
    PRIVATE

    !! * Accessibility
    PUBLIC :: deb_init, deb_cycle, deb_egg_init, readfood3d, readtemp3d   ! routines called by ibm

    !! * Shared module variables
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:),PUBLIC    :: climatemp
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:)  ,PUBLIC    :: biomassezoo

    !! * Private variables
    REAL(kind=rsh), ALLOCATABLE, DIMENSION(:,:)             :: slope, intercept
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:)           :: temp1,temp2
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:)             :: food1,food2
    INTEGER                                                 :: ilecmemfood,ilecmemtemp,ilecf,ilecp
    REAL(kind=rlg)                                          :: tncf1, tncf2
    REAL(kind=rlg)                                          :: tncp1, tncp2 

    !! Set species parameters
    REAL(kind=rsh)                  :: pAm              ! J/d.cm2, spec max assim rate {J_EAm}
    REAL(kind=rsh)                  :: pMi              ! J/cm3/d, volume-spec maintenance costs
    REAL(kind=rsh)                  :: EG               ! J/cm3, volume-spec costs of structure
    REAL(kind=rsh)                  :: vc               ! cm/j, energy conductance
    REAL(kind=rsh)                  :: kap              ! -, Fraction to growth + maintenance
    REAL(kind=rsh)                  :: kj               ! 1/d, maturity maintenance rate coeff
    REAL(kind=rsh)                  :: Hb               ! J, maturity at birth (first feeding)(Nisbet 2012)
    REAL(kind=rsh), PUBLIC          :: Hj               ! J, maturity at metamorphosis
    REAL(kind=rsh), PUBLIC          :: Hp               ! J, maturity at puberty (Nisbet 2012)
    REAL(kind=rsh)                  :: Kx               ! Assimilation efficiency
    
    !Auxillary and compound parameters
    REAL(kind=rsh)                  :: TA               ! K, Arrhenius temperature ;  eggs devlpt -> 9800 = after Regner 1996; larval growth rate : 7500
    REAL(kind=rsh)                  :: T1               ! K, Reference temperature ;
    REAL(kind=rsh)                  :: TAL              ! K, Arrhenius temperature at low temp {Kooy2000}
    REAL(kind=rsh)                  :: TL               ! K, lower boundary temp range
    REAL(kind=rsh)                  :: shape            ! -, shape coefficient
    REAL(kind=rsh)                  :: d_V              ! g/cm^3, specific density of structure (Dw)
    REAL(kind=rsh)                  :: rho_V            ! J/g, = mu_V / w_V; Energy density of structure (DW)
    REAL(kind=rsh)                  :: rho_E            ! J/g, = mu_E / w_E; Energy density of reserve (DW)
    REAL(kind=rsh)                  :: rho_R            ! J/g, = mu_E / w_E; Energy density of Repro buffer (DW) 
    REAL(kind=rsh)                  :: rho_Gam          ! J/g, = mu_E / w_E; Energy density of Gametes (DW)
    REAL(kind=rsh)                  :: K                ! mgC.m3, half saturation constant
    REAL(kind=rsh)                  :: Kr               ! Fraction of the Gamete buffer fixed into eggs             

    ! Ajout Martin
    REAL(kind=rsh)                  :: shapeb           ! -cm, shape coefficient at birth (early larvae after correction for shrinking)
    REAL(kind=rsh)                  :: Sizeb            ! -cm, Length at first feeding (birth in DEB)
    REAL(kind=rsh)                  :: size              
    REAL(kind=rsh)                  :: lfactor          ! fraction capacity for larvae, weight of f
    REAL(kind=rsh)                  :: Lb   
    REAL(kind=rsh)                  :: Lj   

    ! Reproduction       
    REAL(kind=rsh)                  :: TR               ! K, temperature at start spawning
    REAL(kind=rsh)                  :: E0               ! J, Energy for one egg, after Heidi (for initialisation only, otherwise dynamic)
    REAL(kind=rsh)                  :: E_i
    REAL(kind=rsh)                  :: R_i
    REAL(kind=rsh)                  :: Rfbatch  
    REAL(kind=rsh)                  :: SF

    ! Mortality and density-dependency parameters (Menu et al. 2023)
    ! REAL(KIND=rsh)                  :: gammaA     = 0.447_rsh ! with selectivity
    ! REAL(KIND=rsh)                  :: K_biomassA = 10.670_rlg ! with selectivity
    REAL(KIND=rsh)                  :: gammaS     = 0.703343166175024681053_rsh ! with selectivity
    REAL(KIND=rsh)                  :: K_biomassS = 85.315712187958894219264_rlg ! with selectivity
    ! REAL(KIND=rsh), PUBLIC          :: Zea       = 0.063 ! with selectivity
    REAL(KIND=rsh), PUBLIC          :: Zes       = 0.072 ! with selectivity
    ! REAL(KIND=rsh), PUBLIC          :: za        = 0.154 ! with selectivity
    REAL(KIND=rsh), PUBLIC          :: zs        = 0.179 ! with selectivity

    REAL(KIND=rsh), PUBLIC          :: Zaa       = 0.0
    REAL(KIND=rsh), PUBLIC          :: Zas       = 0.0
    REAL(KIND=rsh)                  :: frac_deb  = 0.000999999_rlg

    ! without selectivity
    REAL(KIND=rsh)                  :: gammaA     = 0.384_rsh
    REAL(KIND=rsh)                  :: K_biomassA = 6.156_rlg
    REAL(KIND=rsh), PUBLIC          :: Zea       = 0.056
    REAL(KIND=rsh), PUBLIC          :: za        = 0.136


 !!===================================================================================================================================
 !!===================================================================================================================================
 !!===================================================================================================================================

 CONTAINS


  !!======================================================================
  SUBROUTINE deb_init(restart)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE deb_init  ***
    !&E
    !&E ** Purpose : initialize deb values for state variables, and read constants
    !&E              for species
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_init, ibm_3d
    !&E ** External calls : gasdev_s,tool_decompdate,
    !&E                     ionc4 library
    !&E
    !&E ** History :
    !&E       ! 2011-02  (M. Huret)
    !&E       ! 2023     (C. Menu) Individual variability
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE ionc4,      ONLY : ionc4_openr,ionc4_read_trajt,ionc4_read_traj,ionc4_close,ionc4_read_dimt, ionc4_read_dimtraj
    USE ibmtools,   ONLY : gasdev_s
    USE comtraj,    ONLY : fileanchovy, filesardine, catch_anc_bob, catch_sar_bob
    USE comtraj,    ONLY : mat_catch, fishing_strategy
    USE comtraj,    ONLY : init_anchovy_egg, init_sardine_egg
 
    !! * Arguments
    LOGICAL,intent(IN)                          :: restart
 
    !! * Local declarations
    CHARACTER(LEN=lchain)                       :: file_inp
    INTEGER                                     :: idimt
    INTEGER                                     :: m, n, num, nb_part_nc, is, ie, il, index_num !CLARA
    INTEGER                                     :: lstr,lenstr
    TYPE(type_patch), POINTER                   :: patch
    
    REAL(KIND=rsh)                              :: WV,WE,WR,WG,NRJ_V,NRJ_g,Wat,Wash,L,Wdeb,NRJ
    REAL(KIND=rsh)                              :: zoom
    INTEGER                                     :: jj,mm_clock,aaaa,hh,minu,sec
    CHARACTER(len=19)                           :: tool_sectodat

    INTEGER,        ALLOCATABLE, DIMENSION(:)   :: daysp_nc, yearsp_nc, season_nc, num_nc, zoom_nc !CLARA
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:)   :: edeb_nc, hdeb_nc, rdeb_nc, gam_nc, neggs_nc, wdeb_nc

    CHARACTER(len=lchain)                       :: file_catch
    INTEGER                                     :: id_species
    INTEGER                                     :: row,col
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! Initialisation of embryo stage
    !!----------------------------------------------------------------------
 
    !  IF (.not.F_Fix) THEN
    !     CALL read_NBSS(limin,limax,ljmin,ljmax)
    !  ENDIF

    NAMELIST/deb_para/        pAm,pMi,EG,vc,kap,Hp,Kx,TA,K,shapeb,lfactor,E0,Rfbatch,SF
    NAMELIST/deb_fixed_param/ kj,Hb,Hj,T1,TAL,TL,shape,d_V,rho_V,rho_E,rho_R,rho_Gam,Kr,Sizeb,Lj,TR,E_i,R_i
    
    CALL tool_decompdate(tool_sectodat(time),jj,mm_clock,aaaa,hh,minu,sec)

    patch => patches%first
    DO n = 1,patches%nb
        IF (patch%species == 'anchovy') THEN
            lstr = lenstr(fileanchovy)
            OPEN(52,file=fileanchovy(1:lstr),status='old',form='formatted', access='sequential')
            READ(52,deb_para)
            READ(52,deb_fixed_param)
            CLOSE(52)
        ENDIF
   
        IF (patch%species == 'sardine') THEN
            lstr = lenstr(filesardine)
            OPEN(51,file=filesardine(1:lstr),status='old',form='formatted', access='sequential')
            READ(51,deb_para)
            READ(51,deb_fixed_param)
            CLOSE(51)
        ENDIF

        ! Filtres sur valeurs paramètres
        IF (d_V*rho_V/EG > 0.9_rsh) THEN
            WRITE(*,*) 'Growth efficiency Kg is higher than 0.9, attention'
            STOP
        ENDIF
    
        IF ( restart ) THEN
            ! Initialize array with content of NetCDF input file
            file_inp = trim(patch%file_inp)
            CALL ionc4_openr(file_inp,.false.)
            CALL ionc4_read_dimtraj(file_inp, nb_part_nc) !clara
    
            ! nb_part_nc = patch%nb_part_total ! denis
            ALLOCATE( edeb_nc(nb_part_nc), hdeb_nc(nb_part_nc), rdeb_nc(nb_part_nc), gam_nc(nb_part_nc), wdeb_nc(nb_part_nc) )
            ALLOCATE( neggs_nc(nb_part_nc), daysp_nc(nb_part_nc), yearsp_nc(nb_part_nc), season_nc(nb_part_nc) )
            ALLOCATE( num_nc(nb_part_nc), zoom_nc(nb_part_nc) )
    
            ! CALL ionc4_openr(trim(file_inp), .false.)
            ! Read time dimension in input file to open last time in restart file
            idimt = ionc4_read_dimt(file_inp)

            CALL ionc4_read_trajt(file_inp, "EDEB",      edeb_nc,  1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "HDEB",      hdeb_nc,  1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "RDEB",      rdeb_nc,  1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "GAM",       gam_nc,   1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "NEGGS",     neggs_nc, 1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "WEIGHT",    wdeb_nc,  1, nb_part_nc, idimt)
            CALL ionc4_read_trajt(file_inp, "NUM",  num_nc,      1, nb_part_nc, idimt)
    
            CALL ionc4_read_traj (trim(file_inp), "DAYSPAWN",  daysp_nc,  1, nb_part_nc)
            CALL ionc4_read_traj (trim(file_inp), "YEARSPAWN", yearsp_nc, 1, nb_part_nc)
            CALL ionc4_read_traj (trim(file_inp), "ZOOM",      zoom_nc,   1, nb_part_nc)
            !CALL ionc4_read_traj (file_inp, "SEASON",    season_nc, 1, nb_part_nc)
            CALL ionc4_close(file_inp)
            
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
                patch % particles(m) % E         = edeb_nc(index_num)
                patch % particles(m) % H         = hdeb_nc(index_num)
                patch % particles(m) % R         = rdeb_nc(index_num)
                patch % particles(m) % L         = shape*(patch%particles(m)%size)
                patch % particles(m) % Gam       = gam_nc(index_num)
                patch % particles(m) % Neggs     = neggs_nc(index_num)
                patch % particles(m) % Wdeb      = wdeb_nc(index_num)
                patch % particles(m) % dayspawn  = daysp_nc(index_num)
                patch % particles(m) % yearspawn = yearsp_nc(index_num)
                patch % particles(m) % zoom      = zoom_nc(index_num)
            END DO
 
            DEALLOCATE( edeb_nc, hdeb_nc, rdeb_nc, gam_nc, wdeb_nc, neggs_nc, daysp_nc, yearsp_nc, season_nc )
            DEALLOCATE( num_nc, zoom_nc )
        ENDIF  ! Restart


        ! Les caracteristiques des oeufs d'une meme espece proviennent d'un seul et meme fichier de forcage.
        ! Sauvegarde dans un objet type_particle des informations afin de les donner aux bons individus sans avoir
        ! a deplacer toutes les variables une par une.
        IF (patch%species == 'anchovy') THEN
            init_anchovy_egg % Hj      = Hj
            init_anchovy_egg % pAm     = pAm
            init_anchovy_egg % pMi     = pMi
            init_anchovy_egg % EG      = EG
            init_anchovy_egg % vc      = vc
            init_anchovy_egg % kap     = kap
            init_anchovy_egg % Kx      = Kx
            init_anchovy_egg % Hp      = Hp
            init_anchovy_egg % TA      = TA
            init_anchovy_egg % K       = K
            init_anchovy_egg % shapeb  = shapeb
            init_anchovy_egg % lfactor = lfactor
            init_anchovy_egg % E0      = E0
            init_anchovy_egg % Rfbatch = Rfbatch
            init_anchovy_egg % SF      = SF

            ! to start for eggs, or when initialising a child path
            Lb = sizeb*Shapeb
            init_anchovy_egg % L = Lb
            init_anchovy_egg % E = E_i
            init_anchovy_egg % R = R_i
            init_anchovy_egg % H = Hb

            patch%init_particle = init_anchovy_egg ! transmission des infos au patch apparu

        ELSE IF (patch%species == 'sardine') THEN
            init_sardine_egg % Hj      = Hj
            init_sardine_egg % pAm     = pAm
            init_sardine_egg % pMi     = pMi
            init_sardine_egg % EG      = EG
            init_sardine_egg % vc      = vc
            init_sardine_egg % kap     = kap
            init_sardine_egg % Kx      = Kx
            init_sardine_egg % Hp      = Hp
            init_sardine_egg % TA      = TA
            init_sardine_egg % K       = K
            init_sardine_egg % shapeb  = shapeb
            init_sardine_egg % lfactor = lfactor
            init_sardine_egg % E0      = E0
            init_sardine_egg % Rfbatch = Rfbatch
            init_sardine_egg % SF      = SF

            ! to start for eggs, or when initialising a child patch
            Lb = sizeb*Shapeb
            init_sardine_egg % L = Lb
            init_sardine_egg % E = E_i
            init_sardine_egg % R = R_i
            init_sardine_egg % H = Hb

            patch%init_particle = init_sardine_egg
        ENDIF

        ! Weight
        WV = d_V*(patch%init_particle%L**3)
        WE = patch%init_particle%E/rho_E ! W of E
        WR = patch%init_particle%R/rho_R ! W of E_R
        WG = 0.0_rsh !  W of Gam

        ! Energy
        NRJ_V = rho_V*WV
        NRJ_g = (NRJ_V + patch%init_particle%E + patch%init_particle%R) !/ particle%Wdeb 

        IF (patch%species == 'anchovy') THEN 
           ! Conversion avec relation poids eau vs. taille * Energie
           Wat = exp(-3.19_rsh + 0.38_rsh*sizeb + 0.586_rsh*log(NRJ_g/1d3) - 0.036_rsh*sizeb*log(NRJ_g/1d3)) ! Paul
           ! Ash
           Wash = exp(-2.9326_rsh + 0.9125_rsh*log(Wat)) ! Wash=f(Water) Paul
        ENDIF !(species == 'anchovy')

        IF (patch%species == 'sardine') THEN  
           ! Conversion avec relation poids eau vs. taille * Energie
           Wat = exp((-0.841256_rsh) + 0.199038_rsh*sizeb + 0.277582_rsh*log(NRJ_g/1d3)    &
                    -0.008299_rsh*sizeb*log(NRJ_g/1000.0_rsh))
           ! Ash
           Wash = exp(-2.9325954_rsh + 0.912532_rsh*log(Wat)) ! Wash=f(water)
        ENDIF !( species == 'sardine')

        ! Dry weight   
        patch%init_particle%Wdebd = WV + WE + WR + WG + Wash                         

        IF (Wash / patch%init_particle%Wdebd  > 0.20_rsh) THEN
           Wash = (WV + WE + WR + WG)*0.20_rsh/(1.0_rsh - 0.20_rsh)
           patch%init_particle%Wdebd = WV + WE + WR + WG + Wash
        ENDIF

        ! Energie densite (dry)
        patch%init_particle%NRJd = (NRJ_g/patch%init_particle%Wdebd)/1d3

        ! Wet weight
        patch%init_particle%Wdeb = patch%init_particle%Wdebd + Wat         

        IF ((Wat/patch%init_particle%Wdeb) .gt. 0.85_rsh) THEN
           Wat = patch%init_particle%Wdebd*0.85_rsh/(1.0_rsh - 0.85_rsh)
           patch%init_particle%Wdeb = patch%init_particle%Wdebd + Wat
        ENDIF

        ! Energie densite (wet) 
        patch%init_particle%NRJ  = (NRJ_g/patch%init_particle%Wdeb)/1d3

        nb_part_nc = patch%nb_part_alloc
        DO m = 1,nb_part_nc
            IF (patch%nb_part_alloc == 0) CYCLE
            IF ( .NOT. patch%particles(m)%active ) CYCLE
            
            ! Initialisation des parametres DEB pour chaque particules
            L = patch%particles(m)%L   
            
            IF (patch%particles(m)%size < 2.0_rsh) patch%particles(m)%L = patch%particles(m)%size *shapeb
            IF (patch%particles(m)%size >= 2.0_rsh .and. patch%particles(m)%size <= 4.0_rsh) THEN
                patch%particles(m)%L = patch%particles(m)%size*((L - shapeb*2.0_rsh)*shape + (shape*4.0_rsh - L)*shapeb) &
                                       /(shape*4.0_rsh - shapeb*2.0_rsh)
            ENDIF ! ce if commenté chez clara
            IF (patch%particles(m)%size > 4.0_rsh) patch%particles(m)%L = patch%particles(m)%size*shape

            ! Weight
            WV = d_V*patch%particles(m)%L**3
            WE = patch%particles(m)%E/rho_E         ! W of E
            WR = patch%particles(m)%R/rho_R         ! W of E_R
            WG = patch%particles(m)%Gam/rho_Gam     ! W of G

            ! Energy (M Huret)
            NRJ_V = rho_V*WV
            NRJ_g = NRJ_V + patch%particles(m)%E + patch%particles(m)%R + patch%particles(m)%Gam
        
            size = patch%particles(m)%size


            ! Pour le calcul du poids a l'initialisation si pas de restart
            IF (patch%species == 'anchovy') THEN
                ! Water - Conversion avec relation poids eau vs. taille * Energie
                Wat  = exp(-3.19_rsh + 0.38_rsh*size + 0.586_rsh*log(NRJ_g/1000.0_rsh) - 0.036_rsh*size*log(NRJ_g/1000.0_rsh))    ! Gatti et al. 2017
                Wash = exp(-2.9326_rsh + 0.9125_rsh*log(Wat))               ! Wash=f(Water)
            ENDIF !( species == 'anchovy')
            
            IF (patch%species == 'sardine') THEN
                ! Water - Conversion avec relation poids eau vs. taille * Energie
            	Wat = exp((-0.841256_rsh) + 0.199038_rsh*size + 0.277582_rsh*log(NRJ_g/1000.0_rsh)-0.008299_rsh*size*log(NRJ_g/1d3))
            	Wash = exp(-2.9325954_rsh + 0.912532_rsh*log(Wat))          ! Wash=f(water)
            ENDIF !( species == 'sardine')    
       

            ! Random parameter for population variability
            IF (restart) THEN
                zoom = patch%particles(m)%zoom
            ELSE
                CALL gasdev_s(zoom)
                zoom = 1 + zoom*0.2_rsh/3.0_rsh
            ENDIF

            ! Affectation des nouveaux parametres
            ! From species parameters + individual variability
            patch%particles(m)%Hj       = Hj*(zoom**3.0_rsh) ! utile dans ibm, pour l'instant commun à tout le monde
            patch%particles(m)%Hb       = Hb
            patch%particles(m)%pAm      = pAm*zoom
            patch%particles(m)%pMi      = pMi
            patch%particles(m)%EG       = EG
            patch%particles(m)%vc       = vc
            patch%particles(m)%kap      = kap
            patch%particles(m)%Kx       = Kx
            patch%particles(m)%Hp       = Hp*(zoom**3.0_rsh)
            patch%particles(m)%TA       = TA
            patch%particles(m)%K        = K
            patch%particles(m)%shapeb   = shapeb
            patch%particles(m)%lfactor  = lfactor
            patch%particles(m)%E0       = E0
            patch%particles(m)%Rfbatch  = Rfbatch
            patch%particles(m)%SF       = SF
            patch%particles(m)%zoom     = zoom

            ! Just calculated 
            patch%particles(m)%WV       = WV
            patch%particles(m)%WE       = WE
            patch%particles(m)%WR       = WR
            patch%particles(m)%WG       = 0._rsh
            patch%particles(m)%Wdebd    = WV + WE + WR + WG + Wash
            patch%particles(m)%Wdeb     = patch%particles(m)%Wdebd + Wat
            patch%particles(m)%NRJd     = patch%init_particle%NRJd
            patch%particles(m)%NRJ_V    = NRJ_V
            patch%particles(m)%NRJ_g    = NRJ_g

            ! -- Reproduction
            patch % particles(m) % yearspawn = aaaa
            patch % particles(m) % dayspawn = 500 ! pour etre sur d etre superieur a jjulien
            patch % particles(m) % dayjuv   = 0
        END DO

        ! IF CATCH as fishing strategy
        IF (fishing_strategy == 'Catch') THEN
            IF (patch%species == 'anchovy') THEN
                file_catch = catch_anc_bob 
                id_species = 1
            ELSE IF (patch%species == 'sardine') THEN
                file_catch = catch_sar_bob 
                id_species = 2
            ENDIF

            lstr = lenstr(file_catch)
            OPEN(44, file=file_catch(1:lstr),status='old')
            DO row = 1,20
               READ(44,*) (mat_catch(row,col,id_species),col=1,12)
            END DO
            CLOSE(44)
        ENDIF

        ! Next patch
        patch => patch%next

    END DO  ! loop on patch (n)
 
  END SUBROUTINE deb_init



  !!======================================================================
  SUBROUTINE deb_egg_init(particle,species)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE deb_egg_init  ***
    !&E
    !&E ** Purpose : Initialize DEB part of a given particle with default values, adding 
    !&E              adding inidividual variability. Only for reproduction routine !
    !&E 
    !&E              Init Hj,pAm,pMi,EG,vc,kap,Kx,Hp,TA,K,shapeb,lfactor,E0,Rfbatch,SF,zoom
    !&E              L,H,E,R,WV,WE,WR,WG,NRJ_V,NRJ_G,Wdebd,NRJd,Wdeb,NRJ
    !&E              yearspawn,dayspawn,dayjuv    
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_3d
    !&E ** External calls : tool_decompdate,gasdev_s
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E       ! 2024     (D. Gourves) Add for CROCO, from C. Menu work on individual variability
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE ibmtools,       ONLY : gasdev_s

    !! * Arguments
    TYPE(type_particle),  INTENT(inout) :: particle
    CHARACTER(LEN=lchain),INTENT( in )  :: species

    !! * Local declarations
    INTEGER                             :: jj,mm_clock,aaaa,hh,minu,sec
    CHARACTER(len=19)                   :: tool_sectodat
    REAL(KIND=rsh)                      :: Wash,Wat
    REAL(KIND=rsh)                      :: zoom                     ! zoom factor for inter individual variability

    !!----------------------------------------------------------------------
    !! * Executable part
    CALL tool_decompdate(tool_sectodat(time),jj,mm_clock,aaaa,hh,minu,sec)

    ! Inter individual variability
    CALL gasdev_s(zoom)
    zoom = 1 + zoom*0.2_rsh/3.0_rsh

    !   Affectation des paramètres
    particle % Hj       = particle%Hj*(zoom**3.0_rsh) ! utile dans ibm, pour l'instant commun à tout le monde
    particle % Hb       = Hb*(zoom**3.0_rsh)
    particle % pAm      = particle%pAm*zoom
    particle % pMi      = particle % pMi
    particle % EG       = particle % EG
    particle % vc       = particle % vc
    particle % kap      = particle % kap
    particle % Kx       = particle % Kx
    particle % Hp       = particle % Hp*(zoom**3.0_rsh)
    particle % TA       = particle % TA
    particle % K        = particle % K
    particle % shapeb   = particle % shapeb
    particle % lfactor  = particle % lfactor
    particle % E0       = particle % E0
    particle % Rfbatch  = particle % Rfbatch
    particle % SF       = particle % SF
    particle % zoom     = zoom

    ! On initialise par défaut avec paramètre des parents.
    particle % L        = Sizeb*shapeb
    particle % H        = particle % H
    particle % E        = particle % E
    particle % R        = particle % R

    particle % WV       = d_V*(particle%L**3)
    particle % WE       = particle%E/rho_E
    particle % WR       = particle%R/rho_R
    particle % WG       = 0.0_rsh

    particle % NRJ_V    = rho_V*particle%WV
    particle % NRJ_g    = (particle%NRJ_V + particle%E + particle%R)

    IF (species == "anchovy") THEN    ! Anchovy
        ! Conversion avec relation poids eau vs. taille * Energie
        Wat = exp(-3.19_rsh + 0.38_rsh*sizeb + 0.586_rsh*log(particle%NRJ_g/1d3) - &
              0.036_rsh*sizeb*log(particle%NRJ_g/1d3)) ! Paul
        ! Ash
        Wash=exp(-2.9326_rsh + 0.9125_rsh*log(Wat)) ! Wash=f(Water) Paul
    ENDIF !(species == 'anchovy')

    IF (species == "sardine") THEN    ! Sardines 
        ! Conversion avec relation poids eau vs. taille * Energie
        Wat = exp((-0.841256_rsh) + 0.199038_rsh*sizeb + 0.277582_rsh*log(particle%NRJ_g/1d3)    &
                 -0.008299_rsh*sizeb*log(particle%NRJ_g/1d3))
        ! Ash
        Wash = exp(-2.9325954_rsh + 0.912532_rsh*log(Wat)) ! Wash=f(water)
    ENDIF !( species == 'sardine')
            
    ! Dry weight
    particle%Wdebd = particle%WV + particle%WE + particle%WR + particle%WG + Wash                             
            
    IF (Wash / particle%Wdebd .gt. 0.20_rsh) THEN
        Wash = (particle%WV + particle%WE + particle%WR + particle%WG)*0.20_rsh/(1.0_rsh - 0.20_rsh)
        particle%Wdebd = particle%WV + particle%WE + particle%WR + particle%WG + Wash
    ENDIF

    ! Energie densite (dry)
    particle%NRJd = particle%NRJ_g/(particle%Wdebd/1d3)

    ! Wet weight
    particle%Wdeb = particle%Wdebd + Wat           

    IF (Wat / particle%Wdeb .gt. 0.85_rsh) THEN
        Wat = particle%Wdebd*0.85_rsh/(1.0_rsh - 0.85_rsh)
        particle % Wdeb = particle%Wdebd + Wat
    ENDIF

    ! Energie densite (wet) 
    particle % NRJ  = particle%NRJ_g/(particle%Wdeb/1d3)

    ! -- Reproduction
    particle % yearspawn = aaaa+1
    particle % dayspawn  = 500 ! pour etre sur d etre superieur a jjulien
    particle % dayjuv    = 0
 
  END SUBROUTINE deb_egg_init



  !!======================================================================
  SUBROUTINE deb_cycle(particle,dt,year,month,day,species)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE deb_cycle  ***
    !&E
    !&E ** Purpose :  Modele D.E.B. generique adapte anchois
    !&E **            Calculate variation in DEB compartments with ODE
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_3d
    !&E ** External calls : get_Xdeb
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       ! 2011-02  (M. Huret) from DEB-Larve L. Pecquerie
    !&E       ! 2023     (C. Menu)  Update some variables and add density-dependance
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,    ONLY : jjulien, struc_ad, struc_ad_dd_DEB
    USE comtraj,    ONLY : type_particle
    USE comtraj,    ONLY : F_fix, ffix, frac_deb_death

    !! * Arguments
    TYPE(type_particle),   INTENT( inout )  :: particle
    REAL(KIND=rlg),        INTENT( in )     :: dt
    CHARACTER(LEN=lchain), INTENT( in )     :: species
    INTEGER,               INTENT( in )     :: year, month, day

    !! * Local declarations
    REAL(kind=rsh)                          :: cor_T, K_food,cor_l,shape_fun
    REAL(kind=rsh)                          :: pA,pC,pG,pJ,pR,pMT,pM,pM2,pM3,pR2,pGam,Gamres
    REAL(kind=rsh)                          :: kJT,pAmT,vT,f,pXmT,Eb
    REAL(kind=rsh)                          :: E,L,R,H,WV,WE,WR,WG,NRJ_V,NRJ_g,Ebatch,Ebatchlast,Gam
    REAL(KIND=rsh)                          :: size
    REAL(kind=rsh)                          :: T,X,WWoDW,Ker,kgamT
    REAL(kind=rsh)                          :: dR,dH,dE,dL,dGam
    REAL(KIND=rsh)                          :: Wat, Wash
    LOGICAL                                 :: lastbatch
    INTEGER                                 :: ierr_mpi

    !!----------------------------------------------------------------------
    !! * Executable part

    ! Remise a 0 pour flux non ecrases a chaque pas de temps
    pR2         = 0.0_rsh
    pM2         = 0.0_rsh
    pM3         = 0.0_rsh
    pGam        = 0.0_rsh
    Ebatch      = 0.0_rsh
    Ebatchlast  = 0.0_rsh
    Gamres      = 0.0_rsh

    ! Affectation
    E   = particle % E
    L   = particle % L
    R   = particle % R
    H   = particle % H
    Gam = particle % Gam 

    ! parametres DEB specifiques
    pAm     = particle % pAm
    pMi     = particle % pMi
    EG      = particle % EG
    vc      = particle % vc
    kap     = particle % kap
    Hp      = particle % Hp
    TA      = particle % TA
    K       = particle % K
    shapeb  = particle % shapeb
    lfactor = particle % lfactor
    E0      = particle % E0
    Rfbatch = particle % Rfbatch
    SF      = particle % SF
 
    ! Temperature ------------------------------------------------------------------
    T = particle%temp + 273.0_rsh     ! En Kelvin
    cor_T = exp(TA/T1 - TA/T)


    ! Food -------------------------------------------------------------------------
    X      = 0.0_rsh
    K_food = K !  + L * 10; % account for changes in preference
    f      = ffix

    IF (.not. F_Fix .and. particle%size > 0.0_rsh) THEN
        X = get_Xdeb(particle%xpos, particle%ypos, particle%spos, particle%size, 0.0_rsh, 0.0_rsh, particle%num)
     
        ! -- Densite-dependance  (Menu et al. 2023)
        IF (particle%WV > 0) THEN
                IF (species == 'anchovy') THEN
                    f = X / (X + K_food + struc_ad_dd_DEB/(K_biomassA*10d8*particle%WV**gammaA) ) !BD weight effect     
                ENDIF
                IF (species == 'sardine') THEN
                    f = X / (X + K_food + struc_ad_dd_DEB/(K_biomassS*10d8*particle%WV**gammaS) ) !BD weight effect     
                ENDIF    
        ELSE
            f = X/(X + K_food) 
        ENDIF

        IF (X < 0) THEN
            WRITE (*,*)  "w_dry", particle%Wdebd, "size", particle%size, "L", particle%L, "f", f,             &
                         "nb", particle%super, "strc", struc_ad, "Kbiom", K_biomassA, "Kfood", K_food, "X", X
            CALL_MPI MPI_FINALIZE(ierr_mpi)
            STOP
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------------
    ! shape correction function
    ! correction factor, can be related to increase of assimilation efficiency, or visual acuity
    !cor_l = MIN(1.0_rsh,MAX(((H-Hb)+(Hj-H)*lfactor) / (Hj-Hb),lfactor))
    cor_l = MIN(1.0_rsh, MAX(((L - Lb) + (Lj - L)*lfactor)/(Lj - Lb), lfactor))
    

    ! CORRECTED PARAMETERS (Mars3D, Menu et al. 2023) ------------------------------
    pAmT = pAm*cor_T
    vT   = vc*cor_T
    pMT  = pMi*cor_T
    kJT  = pMi/EG*cor_T
    pXmT = pAmT/0.8_rsh

    ! FLUXES------------------------------------------------------------------------
    pA = pAmT*f*(L**2.0_rsh)*cor_l           ! Assimilation
    pM = pMT*(L**3)                          ! J/d somatic maintenance;

    pC = (E/(L**3))*((vT*EG*(L**2)) + pM)/(EG + (kap*E/(L**3)))  ! J/d, Reserve Mobilisation, ref=Kooijman (Book,p37);            
    pG = (kap*pC) - pM                                           ! J/d, Energy allocated to somatic growth;
    pJ = kJT*H                                                   ! J/d, allocation to maturity maintenance;
    pR = (1._rsh - kap)*pC - pJ                                  ! J/d, allocation to maturity or reproduction;

    ! -- Mortality (Menu et al. 2023)
    ! pR can be negative, then energy goes from R to maturity maintenance, but there need to be enough in R
    Ker = rho_R/rho_E
    IF (pR + R*86400.0_rsh/dt/Ker < 0.0_rsh) THEN
        IF ( frac_deb_death ) THEN
            particle%DEATH_DEB = particle%DEATH_DEB + frac_deb*particle%super
            particle%super = (1._rsh - frac_deb)*particle%super     ! maintenance pas assuree then mortality
        ELSE
            particle%DEATH_DEB = particle%DEATH_DEB + particle%super
            particle%super = 0._rsh     ! maintenance pas assuree then mortality
            particle%flag = -valmanq
        ENDIF
    ENDIF

    ! Reproduction and Emergency maintenance (M Huret)
    IF (H >= Hp) THEN
        IF (species == 'anchovy') THEN
            IF (.not. particle%season .and. month > 3 .and. month < 9 .and. day > 15 .and. year >= particle%yearspawn) THEN
                particle%season = .TRUE.
            ENDIF
              
            Eb = E0*Rfbatch*(particle%Wdeb/1.1064_rsh + 0.5221_rsh) ! AZTI DEPM report, to remove gonad weight
       
            IF ( particle%season ) THEN 
                kgamT = Eb*SF*exp(TA/(273.0_rsh + 13.5_rsh) - TA/T) 
            ENDIF
        ENDIF   ! if (species == 'anchovy') 
    
        IF (species == 'sardine') THEN
            IF (.not. particle%season .and. ((month >= 2 .and. month < 7) .or. (month >= 10 .and. month < 12))) THEN ! &
                !.and. (particle%age>=624)) THEN ! printemps + automne - hiver
                particle%season = .TRUE.
            ENDIF
       
            Eb = E0*Rfbatch*(particle%Wdeb - 1.468302_rsh)/1.095110_rsh ! regression ieo
    
            IF ( particle%season ) THEN
                kgamT = Eb*SF*exp(TA/(273.0_rsh + 11.88_rsh) - TA/T) ! ponte sard   
            ENDIF
        ENDIF ! if (species == 'sardine') 
    
        IF (pG < 0.0_rsh) THEN     !Bad condition, need to remobilise energy 
            pM2 = min(-pG,R*86400.0_rsh/dt)
            IF (particle%season .and. Gam > 0.0_rsh) pM3 = min(Gam*86400.0_rsh/dt, (-pG - pM2)/Kr) ! J/d, second mobilization from gametes 
        ENDIF
    
        IF ( particle%season ) THEN
            pR2 = min(kgamT/Kr, R*86400.0_rsh/dt - pM2)
            pGam = max(0.0_rsh, Kr*pR2)
        ENDIF
    ENDIF   ! if (H>=Hp) then


    !ODE---------------------------------------------------------------------
    dE = pA - pC                            ! J/d, Dynamics of Reserve Energy;
    dL = 0.0_rsh

    IF (pG >= 0.0_rsh .and. particle%super > 0._rsh) THEN
        dL = pG/(3.0_rsh*EG*(L**2))                               ! cm/d, growth increment;
    ELSE                                                          ! Test for mortality 
        IF ((pG + pM2 + pM3*Kr)<0.0_rsh) THEN 
            IF ( frac_deb_death ) THEN
                particle%DEATH_DEB = particle%DEATH_DEB + frac_deb*particle%super 
                particle%super     = (1.0_rlg - frac_deb)*particle%super     ! maintenance pas assuree then mortality 
            ELSE
                particle%DEATH_DEB = particle%DEATH_DEB + particle%super 
                particle%super = 0.0_rlg     ! maintenance pas assuree then mortality
                particle%flag = -valmanq
            ENDIF

       ENDIF
    ENDIF

    ! -- M Huret 
    IF ( H < Hp ) THEN
       dH = pR*Ker      ! J/d, variation in level of maturity; if pR<0 (not enough for maturity maintenance pJ) then rejuvenation
       dR = 0.0_rsh		! J/d, variation in reproduction buffer energy;
    ELSE
        dH = 0.0_rsh
        H  = Hp
        IF ( particle%season ) THEN
            IF (Gam > Eb*2.0_rsh) THEN  ! enough energy to spawn 1 batch and same amount remains for vitellogenesis (Somarakis Reproduce Del.), explicit
                Ebatch = Eb
                particle%Nbatch = particle%Nbatch + 1

                IF(species=='anchovy') THEN
                    particle%Neggs = particle%Neggs + Ebatch/E0*(particle%super/2.0_rlg)           ! total egg spawned per particle over dt_spawn (see ibm) 
                ENDIF
                IF(month <= 7 .and. species == 'sardine') THEN
                    particle%Neggs = particle%Neggs + Ebatch/E0*(particle%super/2.0_rlg)           ! total egg spawned per particle over dt_spawn (see ibm) 
                ENDIF
                particle%Neggs_tot = particle%Neggs_tot + Ebatch/E0*(particle%super/2.0_rlg)       ! total egg spawned per particle over life cycle

                !test if first batch...
                IF (year >= particle%yearspawn) THEN
                    particle%yearspawn = year + 1
                    particle%dayspawn  = jjulien
                ENDIF
            ENDIF

            dGam = pGam - pM3
            Gam  = Gam + dGam*dt/86400.0_rsh - Ebatch

            ! Test if last batch
            lastbatch = .false.
            IF (species == 'anchovy') THEN
                ! last batch, implicit, more than 1 batch should remain as spawning occurs when 2 batch in the gonads...
                ! just to make sure we don't spawn again within the same year
                IF ((jjulien > particle%dayspawn + 120 .or. month > 8)  ) THEN    
                    particle%dayspawn = 367                                     
                    lastbatch           = .true.
                ENDIF
            ENDIF !if (species == 'anchovy')

            IF (species == 'sardine') THEN
                IF(((month >= 7 .and. month < 10) .or. (month == 1) .or. (month == 12))) lastbatch = .true.
            ENDIF !if (species == 'sardine')

            IF ( lastbatch ) THEN
                Ebatchlast = min(Gam, Eb)
                Gamres = Gam - Ebatchlast
                particle%season = .FALSE.
                particle%Nbatch = particle%Nbatch + 1

                IF(species == 'anchovy') THEN
                    particle%Neggs = particle%Neggs + Ebatchlast/E0*(particle%super/2.0_rlg)      ! total egg spawned per particle over dt_spawn (see ibm) 
                ENDIF
                IF(month <= 7 .and. species == 'sardine') THEN
                    particle%Neggs = particle%Neggs + Ebatchlast/E0*(particle%super/2.0_rlg)      ! total egg spawned per particle over dt_spawn (see ibm) 
                ENDIF
                particle%Neggs_tot = particle%Neggs_tot + Ebatchlast/E0*(particle%super/2.0_rlg)  ! total egg spawned per particle over life cycle
            ENDIF

            dR = pR*Ker - pR2 - pM2 ! pR can be negative, then energy goes from R to maturity maintenance (like pM2) 
            particle%Gam = Gam - Ebatchlast - Gamres

        ELSE
            dR = pR*Ker - pM2 ! only necessary fraction is mobilized for emergency, pR can be negative, then energy goes from R to maturity maintenance (like pM2) 
        ENDIF   ! if (particle%season)
    ENDIF       ! if (H<Hp)


    particle%E = E + dE*dt/86400.0_rsh
    particle%L = L + dL*dt/86400.0_rsh
    particle%H = H + dH*dt/86400.0_rsh
    particle%R = R + dR*dt/86400.0_rsh + Gamres*Kr  ! on recupere energie dernier batch dans R

    !---------------------------------------------------------
    ! Calcul du shape si variable en fonction du stade de developement
    ! transition entre shapeb et shape entre 2 et 4 cm
    IF (particle%L <  2.0_rsh*shapeb) shape_fun = shapeb
    IF (particle%L >= 2.0_rsh*shapeb .and. particle%L <= 4.0_rsh*shape)  THEN
        shape_fun = ((L - 2.0_rsh*shapeb)*shape + (4.0_rsh*shape - L)*shapeb)/(4.0_rsh*shape - 2.0_rsh*shapeb)
    ENDIF
    IF (particle%L > 4.0_rsh*shape) shape_fun = shape

    ! Physical length
    particle%size = particle%L/shape_fun !  in cm

    ! Weight
    WV = d_V*particle%L**3
    WE = particle%E/rho_E         ! W of E
    WR = particle%R/rho_R         ! W of E_R
    WG = particle%Gam/rho_Gam       ! W of Gam

    particle%WV = WV
    particle%WE = WE
    particle%WR = WR
    particle%WG = WG

    ! Energy (M Huret)
    NRJ_V = rho_V*WV
    NRJ_g = NRJ_V + particle%E + particle%R + particle%Gam

    size = particle%size

    IF (species == 'anchovy') THEN
        ! Water - Conversion avec relation poids eau vs. taille * Energie
        Wat  = exp(-3.19_rsh + 0.38_rsh*size + 0.586_rsh*log(NRJ_g/1000.0_rsh) - 0.036_rsh*size*log(NRJ_g/1000.0_rsh))    ! Gatti et al. 2017
    
        Wash = exp(-2.9326_rsh + 0.9125_rsh*log(Wat))               ! Wash=f(Water)
    ENDIF !( species == 'anchovy')
    
    
    IF (species == 'sardine') THEN
        ! Water - Conversion avec relation poids eau vs. taille * Energie
    	Wat  = exp((-0.841256_rsh) + 0.199038_rsh*size + 0.277582_rsh*log(NRJ_g/1000.0_rsh) - 0.008299_rsh*size*log(NRJ_g/1d3))
    
    	Wash = exp(-2.9325954_rsh + 0.912532_rsh*log(Wat))          ! Wash=f(water)
    ENDIF !( species == 'sardine')
    
    ! Dry weight
    particle%Wdebd = WV + WE + WR + WG + Wash
    IF (Wash/particle%Wdebd .gt. 0.20_rsh) THEN
        Wash = (WV + WE + WR + WG)*0.20_rsh/(1.0_rsh - 0.20_rsh)
        particle%Wdebd = WV + WE + WR + WG + Wash
    ENDIF
    
    IF (particle%Wdebd .gt. 0.0_rsh) THEN
        particle%NRJd = NRJ_g/particle%Wdebd/1d3    ! Energie densite (dry)
    ENDIF
    
    ! Wet weight
    particle%Wdeb = particle%Wdebd + Wat
    IF (Wat/particle%Wdeb .gt. 0.85_rsh) THEN
        Wat = particle%Wdebd*0.85_rsh/(1.0_rsh - 0.85_rsh)
        particle%Wdeb = particle%Wdebd + Wat
    ENDIF
    
    IF (particle%Wdeb .gt. 0.0_rsh) THEN
        particle%NRJ = NRJ_g/particle%Wdeb/1d3  ! Energie densite (wet)
    ENDIF
    
    IF (particle%super <= 1._rsh) THEN
        particle%Wdeb  = 0.0_rsh
        particle%flag  = -valmanq
        particle%super = 0.0_rsh
        particle%stage = 0
    ENDIF

    ! Auxillary save
    particle%f = f
    particle%X = X

 END SUBROUTINE deb_cycle



 !======================================================================
 SUBROUTINE read_NBSS(limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE read_NBSS  ***
    !&E
    !&E ** Purpose : Read NBSS (spatially explicit, but not temporally) and
    !&E              provide 2D slope and intercept
    !&E
    !&E ** Description    :
    !&E ** Called by      : deb_init
    !&E ** External calls : ionc4_openr,ionc4_read_subxy,ionc4_close
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !  2012-03 (M. Huret) Original code
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE ionc4,      ONLY : ionc4_openr,ionc4_read_subxy,ionc4_close
    USE comtraj,    ONLY : file_NBSS

    !! * Arguments
    INTEGER, INTENT ( in )   :: limin,limax,ljmin,ljmax

    !! * Local declarations
    CHARACTER(LEN=lchain)    :: name_in_slope
    CHARACTER(LEN=lchain)    :: name_in_inter

    !!----------------------------------------------------------------------
    !! * Executable part

    name_in_slope = 'slope'
    name_in_inter = 'intercept'

    CALL ionc4_openr(file_NBSS,l_in_nc4par=.true.)

    ALLOCATE(slope(GLOBAL_2D_ARRAY),intercept(GLOBAL_2D_ARRAY))

    CALL ionc4_read_subxy(file_NBSS,TRIM(name_in_slope),slope,limin,limax,ljmin,ljmax,0,0)
    CALL ionc4_read_subxy(file_NBSS,TRIM(name_in_inter),intercept,limin,limax,ljmin,ljmax,0,0)

    CALL ionc4_close(file_NBSS)

 END SUBROUTINE read_NBSS



 !!======================================================================
 SUBROUTINE readfood3d(limin,limax,ljmin,ljmax,timestep_ibm)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE readfood3d  ***
    !&E
    !&E ** Purpose : read model variables from a NetCDF file
    !&E
    !&E ** Description    : read micro and zooplankton from a MARS 3D file, and sum those
    !&E ** Called by      : ibm_3d
    !&E ** External calls : ionc4_openr,ionc4_read_subzxyt, ionc4_read_dimt, 
    !&E                     ionc4_read_time
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !          (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------

    !! * Modules used
    USE comtraj, ONLY : ierrorlog
    USE ionc4,   ONLY : ionc4_openr,ionc4_read_subzxyt, ionc4_read_dimt, ionc4_read_time, &
                        ionc4_read_subxyt
    USE comtraj, ONLY : file_food
    
    !! * Arguments
    INTEGER,INTENT( in )                                 :: limin,limax,ljmin,ljmax 

    !! * Local declarations
    CHARACTER(LEN=lchain)                                :: name1_in_food,name2_in_food
    CHARACTER(LEN=19)                                    :: tool_sectodat,tdate
    REAL(KIND=rsh),ALLOCATABLE,DIMENSION(:,:)            :: food1_1,food2_1
    REAL(kind=rlg)                                       :: dt1,dt2,torigin, tfood 
    LOGICAL                                              :: l_pb, timestep_ibm
    INTEGER                                              :: valimin,valimax,valjmin,valjmax
    INTEGER                                              :: imin,jmin,imax,jmax
    INTEGER                                              :: IERR_MPI
    INTEGER                                              :: i,j,k,idimt


    !!----------------------------------------------------------------------
    !! * Executable part
    tfood=time
    torigin=0.0_rlg

    name1_in_food = 'zooc'

    CALL ionc4_openr(file_food,l_in_nc4par=.true.)

    ! Definit les indices de lecture en fonction du proc mpi dans le fichier de forcage
    ! Lit sur tout le domaine en sequentiel sinon
    imin = 0 ; jmin = 0
    valimin = 1 ; valjmin = 1 ! version initiale Denis
    ! valimin = 0 ; valjmin = 0 ! version modifiée Clara

#ifdef MPI
    if (ii .gt. 0) then
        valimin = 1 - imin + iminmpi
        imin    = 1
    endif
    if (ii .eq. NP_XI-1) then
        imax = Lmmpi + 1
    else
        imax = Lmmpi
    endif
    if (jj .gt. 0) then
        valjmin = 1 - jmin + jminmpi
        jmin    = 1
    endif
    if (jj .eq. NP_ETA-1) then
        jmax = Mmmpi+1
    else
        jmax = Mmmpi
    endif

#else
    imax = Lm+1
    jmax = Mm+1
#endif

    valimax = imax - imin + valimin
    valjmax = jmax - jmin + valjmin

    ! Read food at first time step from file
    ! IF( FIRST_TIME_STEP ) THEN
    IF ( timestep_ibm ) THEN 

        ALLOCATE(biomassezoo(GLOBAL_2D_ARRAY))
        ALLOCATE(food1(GLOBAL_2D_ARRAY), food2(GLOBAL_2D_ARRAY))
        
        ALLOCATE(food1_1(valimin:valimax,valjmin:valjmax))
        ALLOCATE(food2_1(valimin:valimax,valjmin:valjmax))

        food1 = 0._rsh ; food2 = 0._rsh

        ! Time index in the food file
        ! -------------------------------------
        idimt=ionc4_read_dimt(file_food)
        ilecf=0
        tncf1=0.0_rlg
        tncf2=0.0_rlg
        l_pb = .false.

        DO WHILE((tfood-(tncf1-torigin))*(tfood-(tncf2-torigin)) > 0.0_rlg .AND. ilecf < idimt)
            ilecf=ilecf+1
            IF (ilecf == idimt) THEN
                WRITE(ierrorlog,*) ' '
                WRITE(ierrorlog,*) 'ERROR: routine readfood3d.F90'
                WRITE(ierrorlog,*) 'All records are read, dates are not compatible with simulation time.'
                WRITE(ierrorlog,*) 'Make sure your food file is correct: ',TRIM(file_food)
                WRITE(ierrorlog,*) 'Make sure the origin time is 01/01/1900 00:00:00 in the food file.'
                WRITE(ierrorlog,*) 'If not, modify variable "torigin" in the routine readfood3d.F90.'
                WRITE(ierrorlog,*) ''
                tdate = tool_sectodat(tfood)
                WRITE(ierrorlog,*) 't in MARS model =',tdate
                tdate = tool_sectodat(tncf1-torigin)
                WRITE(ierrorlog,*) 't1 in food file =',tdate
                tdate = tool_sectodat(tncf2-torigin)
                WRITE(ierrorlog,*) 't2 in food file =',tdate
                WRITE(ierrorlog,*) 'Simulation stopped.'
                WRITE(ierrorlog,*) ''
                l_pb = .true.
            END IF
            IF (.NOT. l_pb) THEN
                CALL ionc4_read_time(file_food,ilecf,tncf1)
                CALL ionc4_read_time(file_food,ilecf+1,tncf2)
            END IF
        ENDDO

        CALL ionc4_read_subxyt(file_food,TRIM(name1_in_food),food1_1,valimin,valimax,valjmin,valjmax,ilecf  ,1,1)
        CALL ionc4_read_subxyt(file_food,TRIM(name1_in_food),food2_1,valimin,valimax,valjmin,valjmax,ilecf+1,1,1)
        
        ! Copie des lectures dans les bons indices pour CROCO, ie entre 0 et imax-1
        ! food1(1:valimax-valimin+1,1:valjmax-valjmin+1) = food1_1 ! version initiale Denis
        ! food2(1:valimax-valimin+1,1:valjmax-valjmin+1) = food2_1 ! version initiale Denis
        food1(imin:imax,jmin:jmax) = food1_1(valimin:valimax,valjmin:valjmax) ! version modifiée Clara
        food2(imin:imax,jmin:jmax) = food2_1(valimin:valimax,valjmin:valjmax) ! version modifiée Clara

        
! #ifdef MPI
!             food1(1:valimax-valimin+1,1:valjmax-valjmin+1) = food1_1 ! version initiale Denis
!             food2(1:valimax-valimin+1,1:valjmax-valjmin+1) = food2_1 ! version initiale Denis
! #else
!             food1(LBOUND(food1,1) : MIN(LBOUND(food1,1)+valimax-valimin, UBOUND(food1,1)), LBOUND(food1,2) : MIN(LBOUND(food1,2)+valjmax-valjmin, UBOUND(food1,2)) ) = food1_1 ! version modifiée Clara
!             food2(LBOUND(food2,1) : MIN(LBOUND(food2,1)+valimax-valimin, UBOUND(food2,1)), LBOUND(food2,2) : MIN(LBOUND(food2,2)+valjmax-valjmin, UBOUND(food2,2)) ) = food2_1 ! version modifiée Clara
! #endif

        DEALLOCATE(food1_1,food2_1)
        ilecmemfood=ilecf

    ELSE           ! -------------------------------------
        ! ---   Other read of food from file

        ! Time index in the food file
        idimt=ionc4_read_dimt(file_food)

        ! Decide if read new food data looking at time in file and in model
        DO WHILE((tfood-(tncf1-torigin))*(tfood-(tncf2-torigin)) > 0.0_rlg .AND. ilecf < idimt)
            ilecf=ilecf+1
            IF (ilecf == idimt) THEN
                WRITE(ierrorlog,*) ' '
                WRITE(ierrorlog,*) 'ERROR: routine readfood3d.F90'
                WRITE(ierrorlog,*) 'All records are already read.'
                WRITE(ierrorlog,*) 'The food file is not long enough'
                WRITE(ierrorlog,*) 'Simulation stopped.'
                WRITE(ierrorlog,*) ''
                !CALL_MPI MPI_FINALIZE(MPI_COMM_WORLD, ierr_mpi)
                STOP
            END IF

            CALL ionc4_read_time(file_food,ilecf,tncf1)
            CALL ionc4_read_time(file_food,ilecf+1,tncf2)
        ENDDO

        ! Read and save food data from file
        IF(ilecmemfood/=ilecf) THEN
            ALLOCATE(food2_1(valimin:valimax,valjmin:valjmax))
            CALL ionc4_read_subxyt(file_food,TRIM(name1_in_food),food2_1,valimin,valimax,valjmin,valjmax,ilecf+1,1,1)

            food1 = food2
            ! food2(1:valimax-valimin+1,1:valjmax-valjmin+1) = food2_1 ! version initiale Denis
            food2(imin:imax,jmin:jmax) = food2_1(valimin:valimax,valjmin:valjmax) ! version modifiée Clara
            
! #ifdef MPI
!             food2(1:valimax-valimin+1,1:valjmax-valjmin+1) = food2_1 ! version initiale Denis
! #else
!             food2(LBOUND(food2,1) : MIN(LBOUND(food2,1)+valimax-valimin, UBOUND(food2,1)), LBOUND(food2,2) : MIN(LBOUND(food2,2)+valjmax-valjmin, UBOUND(food2,2)) ) = food2_1 ! version modifiée Clara
! #endif
            
            

            DEALLOCATE(food2_1)
            ilecmemfood=ilecf
        ENDIF

    ENDIF   ! FIRST_TIME_STEP

    dt1=(tncf2-torigin-tfood)/(tncf2-tncf1)
    dt2=1.0_rlg-dt1

    biomassezoo(:,:) = food1(:,:)*dt1 + food2(:,:)*dt2

    ! convert in mgC/m3
    biomassezoo(:,:) = biomassezoo(:,:) !*5.45_rsh*12.0_rsh

 END SUBROUTINE readfood3d


 
 !!======================================================================
 ! --- Unused anymore --- !
 SUBROUTINE readtemp3d(limin,limax,ljmin,ljmax,timestep_ibm)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE readtemp3d  ***
    !&E
    !&E ** Purpose : read temperature data from a NetCDF file
    !&E   !!! USELESS SINCE ONLINE READ FROM HYDRODYNAMIC MODEL !!!
    !&E
    !&E ** Description    : read temp from a climatology file
    !&E ** Called by      : ibm_3d
    !&E ** External calls : ionc4_openr,ionc4_read_subzxyt
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !          (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
 
    !! * Modules used
    USE comtraj,   ONLY : ierrorlog
    USE ionc4,     ONLY : ionc4_openr,ionc4_read_subzxyt,ionc4_read_dimt,ionc4_read_time
 
    !! * Arguments
    INTEGER, INTENT( in )                                :: limin,limax,ljmin,ljmax 
 
    !! * Local declarations
    CHARACTER(LEN=lchain)                                :: name_in_temp
    CHARACTER(LEN=19)                                    :: tool_sectodat,tdatep
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:)        :: temp_1,temp_2
    REAL(kind=rlg)                                       :: dt1,dt2,ttemp,toriginp
    CHARACTER(LEN=lchain)                                :: file_temp
    LOGICAL                                              :: l_pb, timestep_ibm
    INTEGER                                              :: valimin,valimax,valjmin,valjmax
    INTEGER                                              :: imin,imax,jmin,jmax
    INTEGER                                              :: idimt,i,j,k, IERR_MPI
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    file_temp = 'PREPROC/ibm_temp.nc'
    name_in_temp = 'TEMP'
 
    ttemp=time
    toriginp=0.0_rlg
 
    CALL ionc4_openr(file_temp,l_in_nc4par=.true.)

    ! Definit les indices de lecture en fonction du proc mpi dans le fichier de forcage
    ! Lit sur tout le domaine en sequentiel sinon
    imin = 0 ; jmin = 0
    valimin = 1 ; valjmin = 1

#ifdef MPI
    if (ii .gt. 0) then
        valimin = 1 - imin + iminmpi
        imin    = 1
    endif
    if (ii .eq. NP_XI-1) then
        imax = Lmmpi + 1
    else
        imax = Lmmpi
    endif
    if (jj .gt. 0) then
        valjmin = 1 - jmin + jminmpi
        jmin    = 1
    endif
    if (jj .eq. NP_ETA-1) then
        jmax = Mmmpi+1
    else
        jmax = Mmmpi
    endif

#else
    imax = Lm+1
    jmax = Mm+1
#endif

    valimax = imax - imin + valimin
    valjmax = jmax - jmin + valjmin
 
    ! IF( FIRST_TIME_STEP ) THEN
    IF ( timestep_ibm ) THEN

        ALLOCATE(climatemp(GLOBAL_2D_ARRAY,kmax),   &
                 temp1(GLOBAL_2D_ARRAY,kmax),       &
                 temp2(GLOBAL_2D_ARRAY,kmax))
        ALLOCATE(temp_1(valimin:valimax,valjmin:valjmax,kmax),      &
                 temp_2(valimin:valimax,valjmin:valjmax,kmax))

        idimt=ionc4_read_dimt(file_temp)
        ilecp = 0
        tncp1 = 0.0_rlg
        tncp2 = 0.0_rlg
        l_pb  = .false.


        DO WHILE((ttemp-(tncp1-toriginp))*(ttemp-(tncp2-toriginp)) > 0.0_rlg .AND. ilecp < idimt)
            ilecp=ilecp+1
            IF (ilecp == idimt) THEN
              WRITE(ierrorlog,*) ' '
              WRITE(ierrorlog,*) 'ERROR: routine readtemp3d.F90'
              WRITE(ierrorlog,*) 'All records are read, dates are not compatible with simulation time.'
              WRITE(ierrorlog,*) 'Make sure your food file is correct: ',TRIM(file_temp)
              WRITE(ierrorlog,*) 'Make sure the origin time is 01/01/1900 00:00:00 in the temp file.'
              WRITE(ierrorlog,*) 'If not, modify variable "torigin" in the routine readtemp3d.F90.'
              WRITE(ierrorlog,*) ''
              tdatep = tool_sectodat(ttemp)
              WRITE(ierrorlog,*) 't in MARS model =',tdatep
              tdatep = tool_sectodat(tncp1-toriginp)
              WRITE(ierrorlog,*) 't1 in temp file =',tdatep
              tdatep = tool_sectodat(tncp2-toriginp)
              WRITE(ierrorlog,*) 't2 in temp file =',tdatep
              WRITE(ierrorlog,*) 'Simulation stopped.'
              WRITE(ierrorlog,*) ''
              l_pb = .true.
            END IF
            IF (.NOT. l_pb) THEN
              CALL ionc4_read_time(file_temp,ilecp,tncp1)
              CALL ionc4_read_time(file_temp,ilecp+1,tncp2)
            END IF
        ENDDO
        
        CALL ionc4_read_subzxyt(file_temp,TRIM(name_in_temp),temp_1,valimin,valimax,valjmin,valjmax,1,kmax,ilecp  ,1,1,1)
        CALL ionc4_read_subzxyt(file_temp,TRIM(name_in_temp),temp_2,valimin,valimax,valjmin,valjmax,1,kmax,ilecp+1,1,1,1)
        
        temp1(1:valimax-valimin+1,1:valjmax-valjmin+1,:) = temp_1(:,:,:)
        temp2(1:valimax-valimin+1,1:valjmax-valjmin+1,:) = temp_2(:,:,:)
        DEALLOCATE(temp_1,temp_2)
        ilecmemtemp=ilecp
 
    ELSE

        idimt=ionc4_read_dimt(file_temp)
        DO WHILE((ttemp-(tncp1-toriginp))*(ttemp-(tncp2-toriginp)) > 0.0_rlg .AND. ilecp < idimt)
          ilecp=ilecp+1
          IF (ilecp == idimt) THEN
            WRITE(ierrorlog,*) ' '
            WRITE(ierrorlog,*) 'ERROR: routine readtemp3d.F90'
            WRITE(ierrorlog,*) 'All records are already read.'
            WRITE(ierrorlog,*) 'The temp file is not long enough'
            WRITE(ierrorlog,*) 'Simulation stopped.'
            WRITE(ierrorlog,*) ''
            CALL_MPI MPI_FINALIZE(IERR_MPI)
            STOP
          END IF
            CALL ionc4_read_time(file_temp,ilecp,tncp1)
            CALL ionc4_read_time(file_temp,ilecp+1,tncp2)
        ENDDO
 
       IF(ilecmemtemp/=ilecp) THEN
          ALLOCATE(temp_2(valimin:valimax,valjmin:valjmax,kmax))
          CALL ionc4_read_subzxyt(file_temp,TRIM(name_in_temp),temp_2,valimin,valimax,valjmin,valjmax,1,kmax,ilecp+1,1,1,1)

          temp1(1:valimax-valimin+1,1:valjmax-valjmin+1,:) = temp2(:,:,:)
          temp2(1:valimax-valimin+1,1:valjmax-valjmin+1,:) = temp_2(:,:,:)
          DEALLOCATE(temp_2)
          ilecmemtemp=ilecp
       ENDIF
    ENDIF

    dt1=(tncp2-toriginp-ttemp)/(tncp2-tncp1)
    dt2=1.0_rlg-dt1

    climatemp(:,:,:) = temp1(:,:,:)*dt1 + temp2(:,:,:)*dt2
 
  END SUBROUTINE readtemp3d



 !!======================================================================
 FUNCTION get_Xdeb(xpos,ypos,spos,size,minprof,maxprof,debug)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE  get_Xdeb ***
    !&E
    !&E ** Purpose : calculate f (standardized food response) for the DEB
    !&E
    !&E ** Description    :
    !&E ** Called by      : deb_cycle
    !&E ** External calls : ibm_profmean,ksupkinf
    !&E 
    !&E ** History :
    !&E       ! 2012-07  (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,     ONLY : hc_sig,jjulien, type_position
    USE ibmtools,    ONLY : ibm_profmean
    USE trajectools, ONLY : ksupkinf, define_pos

    !! * Arguments
    REAL(kind=rsh), INTENT (in)  :: xpos,ypos,spos,size,minprof,maxprof
    INTEGER,        INTENT(in)   :: debug
    TYPE(type_position)          :: pos 

    REAL(kind=rsh)               :: get_Xdeb

    !! * Local declarations
    REAL(kind=rsh)               ::  slope_max_preywidth,slope_min_preywidth,maxwidth,esdmax,esdmin,minwidth
    REAL(kind=rsh)               ::  biom,sl,inter
    INTEGER                      ::  kwm,kw,i,j

   !!----------------------------------------------------------------------
   !! * Executable part

   !slope_max_preywidth = 0.035_rsh    ! Delivrable FACTS
   !slope_min_preywidth = 0.005        ! Test
   !maxwidth = slope_max_preywidth*size*10.0_rsh ! cm -> mm
   !maxwidth = min(maxwidth,3.0_rsh) ! arbitraire ! pour coherence a maxzoowidth issu du modele ci-dessous
   !esdmax = size2ESD(maxwidth)
   !minwidth = 0.04_rsh+slope_min_preywidth*(size*10.0_rsh-4.0_rsh)
   !minwidth = min(minwidth,0.2_rsh) !toute gamme mesozoo
   !esdmin = size2ESD(minwidth)

   !sl = slope(NINT(xpos),NINT(ypos)) ! slope climato LOPC
   !inter = intercept(NINT(xpos),NINT(ypos))

   !IF ( maxprof == 0.0_rsh )  THEN  ! larva
   !   kwm = 0
   !   kw  = kmax
   !   CALL ksupkinf(spos,sc_w,kmax+1,kw,kwm,12)
   !   biom = biomassezoo(NINT(xpos),NINT(ypos),kw)
   !ELSE
   !   i = NINT(xpos)
   !   j = NINT(ypos)
   !   IF ( h0_g(i,j) > 0.0_rsh ) THEN
   !      biom = ibm_profmean(biomassezoo(i,j,:),minprof,maxprof,h(i,j),zeta(i,j),hc_sig(i,j))
   !   ELSE
   !      biom = 2.0_rsh
   !   ENDIF
   !ENDIF

   !IF( biom <= 0.0_rsh .or. biom >= 999.0_rsh) THEN ! not exactly same grid, so lacking cells along the coast
   !   sl   = -1.0_rsh   ! slope=-1, default
   !   biom = 2.0_rsh  ! 2mg/m3 = mini Poulet1996
   !ENDIF

   !biom = max(2.0_rsh,biom) ! pour eviter valeurs trop basses
   !IF ( jjulien > 90 .and. jjulien < 300 ) THEN
   !   biom = max(10.0_rsh,biom) ! pour eviter valeurs trop basses, mini observe ete
   !ENDIF
   !inter = biom_sl2inter(sl,biom) ! on recalcule l'intercept du NBSS a partir de biomasse modele
   !biom  = NBSS_2biom(sl,inter,esdmin,esdmax) ! on recalcule l'abondance dispo en fonction de la taille
    
    pos%xp = xpos ; pos%yp = ypos
    call define_pos(pos)
    
    i = NINT(pos%idx_r)
    j = NINT(pos%idy_r)

    biom = biomassezoo(i,j)

    !! -- Ajout M Huret 
    !IF (nbss) THEN
    !    slope_max_preywidth=0.035_rsh   ! Delivrable FACTS
    !    slope_min_preywidth=0.005       ! Test 
    !
    !    maxwidth=slope_max_preywidth*particle%size*10.0_rsh ! cm -> mm
    !    maxwidth=min(maxwidth,3.0_rsh) ! ordre de grandeur du max proie observe, ~1cm en longueur, ~ 3mm en largeur. Coherent avec representation  taille modele NPZD, cf routine ci-dessous
    !    
    !    IF (species == 'anchovy') THEN
    !        minwidth=0.04_rsh+slope_min_preywidth*(particle%size*10.0_rsh-4.0_rsh)
    !        minwidth=min(minwidth,0.2_rsh) !toute gamme mesozoo 
    !    ENDIF !(species == 'anchovy')
    !    
    !    IF (species == 'sardine') THEN
    !        IF (particle%size<3.66_rsh) THEN        ! fish length in cm, caldeira 2014
    !            minwidth=(5.27_rsh+0.007_rsh*particle%size*10.0_rsh*1d3)/1d3/3.0_rsh ! conversion: cm-->mm--> mu
    !        ENDIF
    !       
    !        IF (particle%size>=3.66_rsh) THEN       ! bachiller 2012
    !            minwidth=exp(log(particle%size*10.0_rsh)*0.0974_rsh-1.524_rsh)/3.0_rsh ! size to width, hyp : copepod form L=3l
    !        ENDIF
    !    ENDIF !(species == 'sardine')
    !    
    !    esdmax=size2ESD(maxwidth)
    !    esdmin=size2ESD(minwidth)
    !    
    !    newinter=biom_sl2inter(slope(n,dayfood),biom)                ! on recalcule l'intercept du NBSS à partir de biomasse modele
    !    biom=NBSS_2biom(slope(n,dayfood),newinter,esdmin,esdmax)     ! on recalcule l'abondance dispo en fonction de la taille  
    !
    !ENDIF

    biom = max(2.0_rsh,biom) ! pour s'assurer un minimum de zoo quand même, attention en mg/m3, a adapter si autre unité
    
    get_Xdeb = biom

 END FUNCTION get_Xdeb



 ! Appeles dans la partie commentee de get_Xdeb
 !!======================================================================
 FUNCTION biom_sl2inter(sl,biom)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE  biom_sl2inter ***
    !&E
    !&E ** Purpose : calculate the intercept from NBSS slope and model biomass
    !&E
    !&E ** Description :
    !&E ** Called by   : get_Xdeb
    !&E
    !&E ** History :
    !&E       ! 2012-07  (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
 
    !! * Arguments
    REAL(kind=rsh), INTENT (in)  :: sl, biom
    REAL(kind=rsh)               :: biom_sl2inter
 
    !! * Local declarations
    REAL(kind=rsh)   :: minesd,maxesd,minwidth,maxwidth,slin,xx1,xx2
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    minwidth = 0.02_rsh  !  min width zoo dans modele en mm
    maxwidth = 3.0_rsh
    minesd = size2ESD(minwidth) ! transfo en ESD
    maxesd = size2ESD(maxwidth)
    slin = sl
    IF (slin==-1) slin= -1.01
 
    xx1 = ESD2weight(minesd)
    xx2 = ESD2weight(maxesd)
 
    biom_sl2inter = log((biom*(slin+1))/(xx2**(slin+1)-xx1**(slin+1)))
 
 END FUNCTION biom_sl2inter
 


 !!======================================================================
 FUNCTION NBSS_2biom(sl,int,size1,size2)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE  NBSS_2biom ***
    !&E
    !&E ** Purpose : calculate the biomasse within a size range given a NBSS (slope + intercept)
    !&E
    !&E ** Description :
    !&E ** Called by   : get_Xdeb
    !&E
    !&E ** History :
    !&E       ! 2012-07  (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
 
    !! * Arguments
    REAL(kind=rsh), INTENT(in)   :: sl,int,size1,size2
    REAL(kind=rsh)               :: NBSS_2biom
 
    !! * Local declarations
    REAL(kind=rsh)               :: xmid,xx1,xx2
    REAL(kind=rsh)               :: slin, interin,biom
 
    !!----------------------------------------------------------------------
    !! * Executable part
    slin    = sl
    interin = int
    IF (slin==-1) slin = -1.01
 
    interin = exp(int)  ! log(y) = sl * log(x) + int     eq.    y = exp(int) * x^sl
 
    ! convert size1 and size2 from mm esd to mgC:
    xx1 = ESD2weight(size1)
    xx2 = ESD2weight(size2)
 
    ! -- compute xmid in mugC (mean so as to have biom[x1,xmid]=biom[xmid,x2]
    !xmid = ((xx1_2[2]^(sl+1) + xx1_2[1]^(sl+1))/2)^(1/(sl+1));
 
    ! -- compute biom in mgC.m^{-3} and ab in #.m^{-3}
    biom = (interin*xx2**(slin+1))/(slin+1) - (interin*xx1**(slin+1))/(slin+1)
    !ab = biom / (xmid/1000)
 
    NBSS_2biom = biom
 
 END FUNCTION NBSS_2biom
 


 !=========================================================================
 FUNCTION size2ESD(width)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE size2esd ***
    !&E
    !&E ** Purpose : Calculate esd from width
    !&E
    !&E ** Description :
    !&E!      ESD=2*sqrt(ab) ou ESD=sqrt(WL), Herman92 avec w et l largeur et longueur. Pour a=3b, ESD=sqrt(3*w^2)=w*sqrt(3)
    !&E ** Called by   : get_Xdeb
    !&E ** References  : Herman, 1992
    !&E
    !&E ** History :
    !&E       ! 2012-07  (M. Huret)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
 
    !! * Arguments
    REAL(kind=rsh), INTENT(in)   :: width
    REAL(kind=rsh)               :: size2ESD
 
    !! * Local declarations
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    size2ESD = width*3.0_rsh**0.5_rsh
 
 END FUNCTION size2ESD

    

 !=========================================================================
 FUNCTION ESD2weight(esd)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE size2esd ***
    !&E
    !&E ** Purpose : calculate Carbon weight from ESD
    !&E
    !&E ** Description :
    !&E ** Called by   : biomsl2inter, NBSS_2biom
    !&E ** References  : ! Lehette & Hernandez-Leon, 2009 (DW general mesozooplankton)
    !&E                  ! Mauchline, 1998                (average carbon content of copepods)
    !&E
    !&E ** History :
    !&E       !  M. Huret (07-2012)
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
 
    !! * Arguments
    REAL(kind=rsh), INTENT(in)   :: esd
    REAL(kind=rsh)               :: ESD2weight
 
    !! * Local declarations
    REAL(kind=rsh)               :: S,DW
 
    !!----------------------------------------------------------------------
    !! * Executable part
  
    S  = pi*(esd/2.0_rsh)**2.0_rsh
    DW = 43.38_rsh*S**1.5_rsh   !from Lehette & Hernandez-Leon 2009 (general mesozooplankton)
    ESD2weight = 0.447_rsh*DW   !from Mauchline 1998 (average carbon content of copepods)
 
 END FUNCTION ESD2weight


#endif /* IBM_SPECIES */

END MODULE
