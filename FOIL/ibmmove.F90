MODULE ibmmove

!!======================================================================
!!                   ***  MODULE movement   ***
!! 
!!         M. Huret 03/2012
!!         Movement for fish IBM
!!         Several movement rules are implemented
!!        
!&E
!&E ** History :
!&E       !  2024     (M. Huret, D. Gourves) Coupled with CROCO
!&E
!!======================================================================
#include "cppdefs.h"
#include "toolcpp.h" 

#if defined DEB_IBM && defined IBM_SPECIES

    USE module_ibm         ! time,h,om_r,on_r
    USE comtraj,        ONLY : imin,imax,jmin,jmax,kmax,      &
                               rsh,rlg,lchain,valmanq,type_particle
    USE comtraj,        ONLY : 
    use typeSizes

    IMPLICIT NONE
    PRIVATE
  
    !! * Accessibility
    PUBLIC :: fish_move_init         ! routine called by ibm_init
    PUBLIC :: fish_move              ! routine called by ibm

    !! * Shared module variables
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: fish_anc
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: fish_sar

 
    !! * Private variables
    CHARACTER(LEN=lchain)        :: file_fish1,file_fish2
    CHARACTER(LEN=lchain)        :: name_in_fish
 
 



 !!===================================================================================================================================
 !!===================================================================================================================================
 !!===================================================================================================================================

 CONTAINS

  !!===================================================================================================================================
  SUBROUTINE fish_move_init(limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE fish_move_init  ***
    !&E
    !&E ** Purpose : Read fish distribution file
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_init
    !&E ** External calls : ionc4_openr, ionc4_read_time,ionc4_read_subxyt, ionc4_read_subzxyt,
    !&E                     ionc4_read_dimt,ionc4_close
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E       !  2012-03  (M. Huret) Original code
    !&E       !  2024     (M. Huret, D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    
    !! * Modules used

    USE ionc4,   ONLY : ionc4_openr,ionc4_read_dimt,ionc4_read_time,        &
                        ionc4_read_subxyt,ionc4_read_subzxyt,ionc4_close
    USE comtraj, ONLY : fileprobadistrib_anc,nbSizeClass_anc,fileprobadistrib_sar,nbSizeClass_sar
  
    !! * Arguments
    INTEGER, INTENT(in)                             :: limin,limax,ljmin,ljmax
 
    !! * Local declarations
    REAL(KIND=rsh), ALLOCATABLE, DIMENSION(:,:,:)   :: fish1
    INTEGER                                         :: idimt,t
    INTEGER                                         :: valimin,valimax,valjmin,valjmax
    INTEGER                                         :: imin,jmin,imax,jmax
    INTEGER                                         :: lstr, lenstr
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! Open probadistrib for anchovies
    file_fish1   = fileprobadistrib_anc
    write(*,*) 'file proba distrib = ',file_fish1
    name_in_fish = 'Probability'
    
    CALL ionc4_openr(file_fish1,l_in_nc4par=.true.)
    
    ! read number of population distribution along the year
    idimt = ionc4_read_dimt(file_fish1)  
    ALLOCATE(fish_anc(GLOBAL_2D_ARRAY,nbSizeClass_anc,idimt)) 

    
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
    ALLOCATE(fish1   (valimin:valimax,valjmin:valjmax,nbSizeClass_anc)) 

    ! ! Vérifie que la taille du bloc à lire ne dépasse pas fish_anc
    ! ! Modif Clara, pas sur pourquoi les indices de fish_anc commencent à 0 et fish1 à 1 dans 3D-1DV
    ! ! Attention : si fish1 plus petit que fish_anc, pas de message d'erreur
    ! if ((valimax-valimin+1 > UBOUND(fish_anc,1)-LBOUND(fish_anc,1)+1) .or. &
    !     (valjmax-valjmin+1 > UBOUND(fish_anc,2)-LBOUND(fish_anc,2)+1)) then
    !     WRITE(*,*) 'Erreur : le bloc à lire dépasse les dimensions de fish_anc'
    !     STOP
    ! endif !Fin modif Clara

    DO t = 1, idimt
        CALL ionc4_read_subzxyt(file_fish1,TRIM(name_in_fish),fish1,valimin,valimax,valjmin,valjmax,1,nbSizeClass_anc,t,1,1,1)
        ! fish_anc(1:valimax-valimin+1,1:valjmax-valjmin+1,:,t) = fish1 ! version initiale Denis
        fish_anc(0:valimax-valimin,0:valjmax-valjmin,:,t) = fish1 ! version modifiée Clara
! #ifdef MPI
            ! fish_anc(1:valimax-valimin+1,1:valjmax-valjmin+1,:,t) = fish1 ! version initiale Denis
! #else
            ! fish_anc(LBOUND(fish_anc,1) : MIN(LBOUND(fish_anc,1)+valimax-valimin, UBOUND(fish_anc,1)), LBOUND(fish_anc,2) : MIN(LBOUND(fish_anc,2)+valjmax-valjmin, UBOUND(fish_anc,2)), :, t ) = fish1 ! version modifiée Clara
! #endif 

    END DO
    CALL ionc4_close(file_fish1)
    DEALLOCATE(fish1)
    
    ! Open probadistrib for sardines
    file_fish2   = fileprobadistrib_anc
    CALL ionc4_openr(file_fish2,l_in_nc4par=.true.)

    ! read number of population distribution along the year
    idimt = ionc4_read_dimt(file_fish2)              
 
    ALLOCATE(fish_sar(GLOBAL_2D_ARRAY,nbSizeClass_sar,idimt))
    ALLOCATE(fish1   (valimin:valimax,valjmin:valjmax,nbSizeClass_sar))

    ! ! Vérifie que la taille du bloc à lire ne dépasse pas fish_sar
    ! ! Modif Clara, pas sur pourquoi les indices de fish_sar commencent à 0 et fish1 à 1 dans 3D-1DV
    ! ! Attention : si fish1 plus petit que fish_sar, pas de message d'erreur
    ! if ((valimax-valimin+1 > UBOUND(fish_sar,1)-LBOUND(fish_sar,1)+1) .or. &
    !     (valjmax-valjmin+1 > UBOUND(fish_sar,2)-LBOUND(fish_sar,2)+1)) then
    !     WRITE(*,*) 'Erreur : le bloc à lire dépasse les dimensions de fish_sar'
    !     STOP
    ! endif !Fin modif Clara
    
    DO t = 1, idimt
        CALL ionc4_read_subzxyt(file_fish2,TRIM(name_in_fish),fish1,valimin,valimax,valjmin,valjmax,1,nbSizeClass_sar,t,1,1,1)
        ! fish_sar(1:valimax-valimin+1,1:valjmax-valjmin+1,:,t) = fish1 ! version initiale Denis
        fish_sar(0:valimax-valimin,0:valjmax-valjmin,:,t) = fish1 ! version modifiée Clara
! #ifdef MPI
!             fish_sar(1:valimax-valimin+1,1:valjmax-valjmin+1,:,t) = fish1 ! version initiale Denis
! #else
!             fish_sar(LBOUND(fish_sar,1) : MIN(LBOUND(fish_sar,1)+valimax-valimin, UBOUND(fish_sar,1)), LBOUND(fish_sar,2) : MIN(LBOUND(fish_sar,2)+valjmax-valjmin, UBOUND(fish_sar,2)), :, t ) = fish1 ! version modifiée Clara
! #endif

    END DO

    CALL ionc4_close(file_fish2)

    DEALLOCATE(fish1)
 
  END SUBROUTINE fish_move_init



  !!======================================================================
  SUBROUTINE fish_move(particle,ind_species)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE fish_move  ***
    !&E
    !&E ** Purpose : routine for random walk within a given distribution area
    !&E              following a Metropolis algorithm  
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_3d
    !&E ** External calls : tool_decompdate, define_pos
    !&E ** Reference      :
    !&E
    !&E ** History        :
    !&E       !  2012-03  (M. Huret) Original code    
    !&E       !  2024     (M. Huret, D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------    
    !! * Modules used
    USE trajectools, ONLY : define_pos
    USE comtraj,     ONLY : type_position
    USE comtraj,     ONLY : sizemin_anc,sizemin_sar,nbSizeClass_anc,nbSizeClass_sar
#ifdef MPI
    USE comtraj,     ONLY : down_give, up_give, right_give, left_give
#endif
    
    !! * Arguments
    TYPE(type_particle), INTENT(inout)   :: particle
    INTEGER, INTENT(in)                  :: ind_species
    
    !! * Local declarations
    LOGICAL              :: move
    CHARACTER(len=19)    :: tool_sectodat
    REAL(rsh)            :: v1,v2,Pi,Pj,harvest,xpos_n,ypos_n,cell
    REAL(KIND=rsh)       :: fpos_x, fpos_y
    INTEGER              :: index,icell,jcell,ipos,jpos,saison,icells,jcells
    INTEGER              :: jj,mm_clock,aaaa,hh,minu,sec
    INTEGER              :: ierr_mpi
    TYPE(type_position)  :: pos,pos_n 

    REAL(KIND=rsh)       :: sizemin
    
    INTEGER, PARAMETER   :: n = 1     ! FIXME : a l'origine, c'etait le numero du patch... (mh-2015)
     
    !!----------------------------------------------------------------------
    !! * Executable part

    IF (ind_species == 1) THEN
        sizemin = sizemin_anc
    ELSE IF (ind_species == 2) THEN
        sizemin = sizemin_sar
    ENDIF

    pos%xp = particle%xpos ; pos%yp = particle%ypos
    CALL define_pos(pos)
    
    !Get the index of size in the probability distribution matrix
    index  = min(max(1, NINT(particle%size*2) - INT(sizemin*2) + 1), 34) ! 34=20cm, au dela rien a l automne

    icells = NINT(pos%idx_r)
    jcells = NINT(pos%idy_r)
    
    !Get the seasonal distribution
    CALL tool_decompdate(tool_sectodat(time),jj,mm_clock,aaaa,hh,minu,sec)
    
    ! when transition is let free, it takes 2 to 3 months to come to a reasonable R2 for comparison
    ! between obs and model distribution by size (March-April and July-August)
    IF (n == 1) THEN ! normal migration
        saison = 1
        IF(mm_clock <= 2 .OR. mm_clock >= 7) THEN
            saison = 2
        END IF
        IF (mm_clock <= 2 .OR. mm_clock >= 10) THEN 
            IF (particle%age < 365) THEN ! to avoid remaining offshore for juveniles (still small) in winter...
                saison = 1
            END IF
        END IF
        IF (particle%size < 7.25) saison = 2 ! pour cas ou juvenile avant aout, pas de size_class correspondante saison 1
    END IF
    
    !IF (n == 2) THEN ! remain in the south
    !   saison = 1
    !   IF (particle%size < 7.25) saison = 2 ! pour cas ou juvenile avant aout, pas de size_class correspondante saison 1
    !END IF
    !
    !IF (n == 3) THEN ! remain in the north
    !    saison=2
    !    IF(mm_clock <= 2  .OR. mm_clock >= 10) THEN 
    !        IF ((time - particle%date_orig) < 86400.0_rlg*365.0_rlg) THEN ! to avoid remaining offshore for juveniles in winter...
    !            saison=1
    !        END IF
    !    END IF
    !    IF (particle%size < 7.25) saison = 2 ! pour cas ou juvenile avant aout, pas de size_class correspondante saison 1
    !END IF
    IF (ind_species == 1) THEN
        Pi = fish_anc(icells,jcells,index,saison)
    ELSE IF (ind_species == 2) THEN
        Pi = fish_sar(icells,jcells,index,saison)
    ENDIF

    !Select a neighbouring cell j at random among the 4 neighbouring cells
    CALL random_number(v1) 
    !CALL random_number(v2)
    !icell = icells+(2*NINT(v1)-1) ! -1 or 1
    !jcell = jcells+(2*NINT(v2)-1)
    
    icell = icells
    jcell = jcells
    cell = NINT(v1*4.0_rsh-0.5_rsh)

    ! f_posx and f_posy allow to keep global position for movement in case of MPI
    fpos_x = NINT(particle%xpos) ; fpos_y = NINT(particle%ypos) 
    IF (cell == 0 .or. cell == -1) THEN
        icell  = icells - 1
        fpos_x = fpos_x - 1
    ELSE IF (cell == 1 ) THEN
        jcell  = jcells - 1 
        fpos_y = fpos_y - 1
    ELSE IF (cell == 2 ) THEN
        icell  = icells + 1 
        fpos_x = fpos_x + 1
    ELSE IF (cell == 3 .or. cell == 4) THEN
        jcell  = jcells + 1 
        fpos_y = fpos_y + 1
    ENDIF

    icell = MIN(MAX(icell,imin),imax) ! check IF in boundaries
    jcell = MIN(MAX(jcell,jmin),jmax)
    
    ! Get the probability at this new cell
    IF (ind_species == 1) THEN
        Pj = fish_anc(icell,jcell,index,saison)
    ELSE IF (ind_species == 2) THEN
        Pj = fish_sar(icell,jcell,index,saison)
    ENDIF
    
    ! Metropolis algorithm to decide on whether to move or not
    CALL random_number(harvest) 
    move = .FALSE.
    IF (Pi == 0.0_rsh) move = .TRUE. ! adapted
    IF (Pi >  0.0_rsh) move = ((Pj >= Pi) .or. (harvest < Pj/Pi))
    
    ! Calculate displacement (dt forced at 1h, see ibm.f90)
    IF ( move ) THEN
    
        ! To modulate movement with fish velocity, but not really compatible with Metropolis...
        !IF (Pi==0.0_rsh .and. h(icells,jcells)>1000.0_rsh .and. particle % size>8.0_rsh) THEN ! to have the vagrants above a given size (8) offshore back on shelf
        !    xpos_n=particle % xpos+2.0_rsh*particle % size*0.01_rsh*3600.0_rsh/om_r(icells,jcells)
        !    ypos_n=particle % ypos
        !ELSE IF (jcell==jcells) THEN
        !    xpos_n=particle % xpos+2.0_rsh*particle % size*0.01_rsh*3600.0_rsh/om_r(icells,jcells)*(icell-icells)
        !    ypos_n=particle % ypos
        !ELSE IF (icell==icells) THEN
        !    ypos_n=particle % ypos+2.0_rsh*particle % size*0.01_rsh*3600.0_rsh/on_r(icells,jcells)*(jcell-jcells)
        !    xpos_n=particle % xpos
        !ELSE
        !    xpos_n=particle % xpos+2.0_rsh*particle % size*0.01_rsh*3600.0_rsh/om_r(icells,jcells)*(icell-icells)/sqrt(2.0_rsh)
        !    ypos_n=particle % ypos+2.0_rsh*particle % size*0.01_rsh*3600.0_rsh/on_r(icells,jcells)*(jcell-jcells)/sqrt(2.0_rsh)
        !END IF
        !IF (h(NINT(xpos_n),NINT(ypos_n))>0.0_rsh) THEN ! test si on reste en mer
        !    IF (NINT(xpos_n) .gt. (NINT(particle % xpos)+1) .or. NINT(xpos_n) .lt. (NINT(particle % xpos)-1) .or. NINT(ypos_n) .gt. (NINT(particle % ypos)+1) .or. NINT(ypos_n) .lt. (NINT(particle % ypos)-1)) THEN
        !        WRITE(*,*)  particle % xpos,xpos_n,particle % ypos,ypos_n
        !        WRITE(*,*) 'saut de CPU'
        !        STOP
        !    END IF
        !    particle % xpos=xpos_n
        !    particle % ypos=ypos_n
        !END IF
    
        IF (Pi == 0.0_rsh) THEN ! to have the vagrants above a given size offshore back on shelf
            ypos_n = particle%ypos
            xpos_n = particle%xpos+1
        ELSE
            CALL random_number(harvest)
            xpos_n = min(real(fpos_x,kind=rsh) + harvest - 0.49_rsh, real(fpos_x,kind=rsh) + 0.49_rsh) ! to have it randomly within new cell
            CALL random_number(harvest)
            ypos_n = min(real(fpos_y,kind=rsh) + harvest - 0.49_rsh, real(fpos_y,kind=rsh) + 0.49_rsh)
        END IF
    
       ! detection d un franchissement limite de domaine
        IF (ypos_n <= 1.0_rsh.or.ypos_n >= real(jmax)) THEN
            particle%flag  = -valmanq ! on prefere flagger la particule et la laisser ou elle est (ne bougera plus) 
            particle%super = 0
        END IF
        IF (xpos_n <= 1.0_rsh.or.xpos_n >= real(imax)) THEN
            particle%flag  = -valmanq ! on prefere flagger la particule et la laisser ou elle est (ne bougera plus) 
            particle%super = 0
        END IF

        pos_n%xp = xpos_n ; pos_n%yp = ypos_n
        CALL define_pos(pos_n)
    
        IF (h(NINT(pos_n%idx_r),NINT(pos_n%idy_r)) > 0.0_rsh .and. particle%flag /= -valmanq) THEN
            ! test si on reste en mer (possibilite si incompatibilite de grille)
            IF ( (NINT(pos_n%idx_r) > icells+1) .or. (NINT(pos_n%idx_r) < icells-1) .or. &
                 (NINT(pos_n%idy_r) > jcells+1) .or. (NINT(pos_n%idy_r) < jcells-1) ) THEN 

                PRINT*, ' ERREUR MPI : ', pos_n%idx_r,pos_n%idy_r,icells+1,icells-1,jcells+1,jcells-1
                WRITE(*,*) 'STOP jump MPI'
                CALL_MPI MPI_FINALIZE(ierr_mpi)
                STOP
            END IF
            particle%xpos = xpos_n
            particle%ypos = ypos_n
        END IF
    END IF    ! move

#ifdef MPI     
    ! check IF we need to exchange or not between CPUs
    ! -----------------------------------
    ipos = NINT(particle%xpos)
    jpos = NINT(particle%ypos)
    
    ! here we verify IF the particle stays in your zone
    ! although it gets out of its zone; it will only move max to 1+-1 and j+-1 
    ! thus we can compute the particle s new value in this time step without pb
    !
    ! IF the particle went out the zone; we count on which side it went out and
    ! note it in the down give up give r give l give 
    ! because this value is needed when all the particle finished calculated for this time step
    ! to CALL ex_traj
    !
    ! The neighbor domain where the particle went is coded in variable 'limitbye' according
    ! to this scheme :
    !         7  |  6  |  5
    !       -----+-----+-----
    !         8  |  0  |  4
    !       -----+-----+-----
    !         1  |  2  |  3wo

    IF ( jpos < jminmpi ) THEN
        particle%limitbye = 2
        down_give = down_give+1
        IF ( ipos < iminmpi ) particle%limitbye = 1
        IF ( ipos > imaxmpi ) particle%limitbye = 3

    ELSE IF ( jpos > jmaxmpi ) THEN
        particle%limitbye = 6
        up_give = up_give+1
        IF ( ipos < iminmpi ) particle%limitbye = 7
        IF ( ipos > imaxmpi ) particle%limitbye = 5

    ELSE IF ( ipos < iminmpi ) THEN
        particle%limitbye = 8
        left_give = left_give + 1

    ELSE IF ( ipos > imaxmpi ) THEN
        particle%limitbye = 4
        right_give = right_give + 1
    END IF

#endif
     
  END SUBROUTINE fish_move   
  !!======================================================================

#endif /* DEB_IBM */

END MODULE


