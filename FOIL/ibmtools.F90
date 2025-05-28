MODULE ibmtools

!==================================================================================   
!!                                   *** MODULE ibmtools ***
!!
!&E ** Purpose : Generic routines and functions for IBM 
!&E
!&E ** Description : 
!&E
!&E ** Called by : ibm_init, ibm_3d
!&E
!&E ** External calls : trajectools (interpvit,h_int_siggen,ksupkinf,siggentoz) 
!&E
!&E ** History :
!&E       !  2011-01 (M. Huret) Original code
!&E       !  2013-04 (M. Huret) Separation from IBM
!&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
!&E       !  2024    (M. Caillaud, M. Huret, D. Gourves) Coupled with CROCO
!
!==================================================================================   

#include "cppdefs.h"
#include "toolcpp.h"
#ifdef MPI
   USE mpi
#endif

#if defined DEB_IBM
    !! * Module used
    USE module_ibm          ! time,sc_r,sc_w,Cs_r,h,hc,g,srflx,zeta
    USE comtraj, ONLY       : imin,imax,jmin,jmax,kmax,rsh,rlg,type_particle,lchain
   
    IMPLICIT NONE
    PRIVATE 

   !! * Accessibility
    PUBLIC w_dens, ibm_buoy, ibm_traint, ibm_nycth_mig, ibm_proftraint,   &
#ifdef IBM_SPECIES
           ibm_parameter_init,death_by_fishing,                           &
#endif
          ibm_profmean, ibm_loc_xyz, gasdev_s,tool_julien               
          !ibm_profuint, ibm_profvint                                  ! non utilise

    !! * Shared module variables
    REAL(kind=rsh),PUBLIC, dimension(GLOBAL_2D_ARRAY)   :: opt_depth_part
    !! * Private variables
    REAL(kind=rsh),dimension(GLOBAL_2D_ARRAY)           :: depth_strat  !-1 no stratif; >0 stratif

    INTEGER, PARAMETER                                  :: track = 1  

   



 !!===================================================================================================================================
 !!===================================================================================================================================
 !!===================================================================================================================================

 CONTAINS



 FUNCTION tool_julien(jou,moi,ia)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_julien  ***
    !&E
    !&E ** Purpose : compute the Julian day number from the date entered as 
    !&E              function arguments in the form of day, month and year
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm3d
    !&E ** External calls : 
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !  2004-08-18
    !&E       !  2024    (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    IMPLICIT NONE

    !! * Declarations function
    INTEGER              :: tool_julien

    !! * Arguments
    INTEGER, INTENT(in)  :: ia,jou,moi

    !! * Local declarations
    INTEGER,DIMENSION(12,2),PARAMETER :: jo=                             &
        RESHAPE((/0,31,59,90,120,151,181,212,243,273,304,334             &
                ,0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
    INTEGER                           :: i,ibs,m,iy,j,mois,jour
    REAL(kind=rsh)                    :: a,b

   !!----------------------------------------------------------------------
   !! * Executable part

    mois=moi
    jour=jou
    IF(mois == 1) THEN
        ibs=1
        IF (ia < 1582 .OR. (ia == 1582 .AND. jou < 283)) GOTO 10
        ibs=MOD(ia,4)+2
        IF(ibs /= 2) ibs=1

10 CONTINUE
        DO i=1,12
            mois=12-i+1
            jour=jou-jo(mois,ibs)
            IF(jour > 0) GOTO 20
        ENDDO
    ENDIF

20 CONTINUE
    a=REAL(ia,rsh)+REAL(mois,rsh)/100.0_rsh+REAL(jour,rsh)/10000.0_rsh
    b=0.0_rsh
    iy=ia
    m=mois
    IF(mois <= 2) THEN
        iy=iy-1
        m=m+12
    ENDIF

    IF(a >= 1582.1015_rsh) THEN
        j=iy/100
        b=2.0_rsh-REAL(j,rsh)+REAL(j,rsh)/4.0_rsh
    ENDIF

    tool_julien=INT(365.25_rsh*iy)+INT(30.6001_rsh*(m+1))+jour+1720995+INT(b)

 END FUNCTION tool_julien



#ifdef IBM_SPECIES

  !!======================================================================
  SUBROUTINE ibm_parameter_init(particle,species,xe,sal,temp,Istr,Iend,Jstr,Jend)
    !&E---------------------------------------------------------------------
    !&E                *** SUBROUTINE death_by_fishing ***
    !&E ** Purpose : Init IBM parameters, at start simulation and for particles created
    !&E              by reproduction
    !&E              Drate, w, temp, density, denspawn and size (size only in case of repro)
    !&E ** 
    !&E ** Called by      : ibm_3d, ibm_init
    !&E ** External calls : define_pos,ztosiggen,ibm_loc_xyz,gasdev_s
    !&E ** Reference      : 
    !&E
    !&E ** History        : 
    !&E       ! 2024     (D. Gourves) From MARS3D gathered into a function
    !&E---------------------------------------------------------------------
    !! Modules used
    USE trajectools,  ONLY : define_pos,ztosiggen
    USE comtraj,      ONLY : type_position

    !! Arguments
    TYPE(type_particle),                               INTENT( inout )  :: particle
    CHARACTER(LEN=lchain),                             INTENT( in )     :: species
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,4),      INTENT( in )     :: xe
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax),   INTENT( in )     :: sal, temp
    INTEGER,                                           INTENT( in )     :: Istr,Iend,Jstr,Jend 

    ! Local declaration
    ! To save a local and global position of particle for MPI and Sequential compatibility
    TYPE(type_position)                             :: pos
    
    ! Temporary indexes of particle's location
    INTEGER                                         :: igg,idd,jhh,jbb,kp,km
    INTEGER                                         :: hlb,hlt,hrb,hrt
    REAL(KIND=rsh)                                  :: px,py
    
    REAL(KIND=rsh)                                  :: sal_surf,temp_surf,dens_surf ! Physical variables particle's location
    REAL(KIND=rsh)                                  :: prof,sprof                   ! Z and Sigma temporary depths
    
    REAL(KIND=rsh)                                  :: tir                          ! Random values between in [-1,1]
    REAL(KIND=rsh)                                  :: ec_type_size, ec_type        ! Egg density and size variability

    
    ! Egg density ----------------
    prof = min(4.0_rsh, particle%d3) ! surface
    prof = -prof + particle%xe ! switch from immersion to real z
    
    CALL ztosiggen(prof, particle%spos, particle%xe, particle%h0, particle%hc)

    ! -- Locate the particle
    pos%xp = particle%xpos ; pos%yp = particle%ypos
    call define_pos(pos)
    CALL ibm_loc_xyz(pos%idx_r, pos%idy_r, particle%spos,           &
                     px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,kp,km,   &
                     Istr,Iend,Jstr,Jend)

    ! Initialize particle's stage, Drate and w
    particle%Drate     = 0.0_rsh 
    particle%w         = 0.0_rsh

    ! Initialize particle's size, checking at species
    !Huret et al. 2016
    IF (species == 'anchovy') THEN
        ! Egg density variability ----------------
        ec_type      = 0.714_rsh
        ! Egg size and variability ----------------------
        IF (particle%stage == 1) particle%size = 0.855_rsh*0.1_rsh ! egg diameter (mm ->cm)    
        ec_type_size = 0.036_rsh*0.1_rsh
    ENDIF
    IF (species == 'sardine') THEN
        ! Egg density variability ----------------
        ec_type      = 0.832_rsh
        ! Egg size and variability ----------------------
        IF (particle%stage == 1) particle%size = 1.608_rsh*0.1_rsh ! egg diameter (mm ->cm)
        ec_type_size = 0.095_rsh*0.1_rsh
    ENDIF

    ! Initialize particle's density and temperature at its location
    ! Egg density = f(surf. density)
    sal_surf  = ibm_traint(sal,xe,particle%spos,kp,km,px,py,       &
                           igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,        &
                           Istr,Iend,Jstr,Jend)
    temp_surf = ibm_traint(temp,xe,particle%spos,kp,km,px,py,      &
                           igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,        &
                           Istr,Iend,Jstr,Jend)
    dens_surf = w_dens(temp_surf, sal_surf)

    ! Si temperature realiste a l'initialisation, on prend sa valeur, sinon on prend 0
    particle%temp    = 0._rsh
    particle%density = dens_surf

    ! Initialize particle's denspawn
    IF (species == "anchovy") THEN
       particle%denspawn = -3.434_rsh + 1.1_rsh*(dens_surf - 1000.0_rsh) + 1000.0_rsh ! with thermal expansion correction
       !particle % denspawn = particle % denspawn + 0.5_rsh ! cf Paul G.
    ENDIF
    IF (species == "sardine") THEN
       particle%denspawn = 5.343_rsh + 0.748_rsh*(dens_surf - 1000.0_rsh) + 1000.0_rsh ! sardine
       !particle % denspawn = particle % denspawn + 0.78_rh ! sardine
    ENDIF

    ! -- Randomly modify denspawn and size of particles
    CALL gasdev_s(tir)
    particle%denspawn = particle%denspawn + tir*ec_type
    particle%density  = particle%denspawn
    particle%size     = particle%size     + tir*ec_type_size

  END SUBROUTINE ibm_parameter_init



  !!======================================================================
  SUBROUTINE death_by_fishing(particle,species,year,month,dtm)
    !&E---------------------------------------------------------------------
    !&E                *** SUBROUTINE death_by_fishing ***
    !&E ** Purpose : Calculate death by fishing of a patch, deoending on chosen scenario
    !&E ** 
    !&E ** Called by      : ibm_3d
    !&E ** External calls : 
    !&E ** Reference      : Clara Menu's thesis, 2023
    !&E
    !&E ** History        : 
    !&E       ! 2024     (D. Gourves) From Menu's 1D code to CROCO
    !&E---------------------------------------------------------------------
    !&E Calculate death by fishing of a patch

    USE comtraj,            ONLY : catch_anc_bob, catch_sar_bob, fishing_strategy
    USE comtraj,            ONLY : mat_catch, number_tot, weight_tot, biom_tot, Wdeb_mean

    !! Arguments
    TYPE(type_particle),    INTENT( inout ) :: particle
    REAL(KIND=rsh),         INTENT( in )    :: dtm
    INTEGER,                INTENT( in )    :: year,month 
    CHARACTER(len=lchain),  INTENT( in )    :: species

    ! Local arguments
    INTEGER                                 :: id_species
    REAL(KIND=rsh)                          :: Zfishing
    REAL(KIND=rsh)                          :: c_F     = 1.0_rsh
    REAL(KIND=rlg)                          :: percent = 1.0_rsh  
    REAL(KIND=rsh)                          :: coeff1

    !!----------------------------------------------------------------------
    !! * Executable part
    Zfishing = 0.0_rlg

    ! ============================================================================
    ! =====                                                                  =====
    ! =====                     Fishing strategy : F_eval                    =====
    ! =====                                                                  =====
    ! Fishing of First Age Class
    IF (particle%AgeClass >= 1 .and. fishing_strategy == 'F_eval') THEN !ageclass>=1 to be sure not to fish newborns
    
        IF (particle%stage == 6 .and. month < 7 .and. species == 'anchovy') THEN
            !IF (year < 1990) Zfishing = f_spin
            !IF (year == 1990) Zfishing = 0.00344_rlg
            !IF (year == 1991) Zfishing = 0.00321_rlg
            !IF (year == 1992) Zfishing = 0.00335_rlg 
            !IF (year == 1993) Zfishing = 0.00242_rlg 
            !IF (year == 1994) Zfishing = 0.00327_rlg
            !IF (year == 1995) Zfishing = 0.00436_rlg
            !IF (year == 1996) Zfishing = 0.00349_rlg
            !IF (year == 1997) Zfishing = 0.00185_rlg
            !IF (year == 1998) Zfishing = 0.00131_rlg
            !IF (year == 1999) Zfishing = 0.00150_rlg
            IF (year <= 2000) Zfishing = 0.00200_rlg
            IF (year == 2001) Zfishing = 0.00185_rlg
            IF (year == 2002) Zfishing = 0.00155_rlg
            IF (year == 2003) Zfishing = 0.00108_rlg
            IF (year == 2004) Zfishing = 0.00242_rlg
            IF (year == 2005) Zfishing = 0.00044_rlg
            IF (year == 2006) Zfishing = 0.00068_rlg
            IF (year == 2007) Zfishing = 0.000038_rlg
            IF (year == 2008) Zfishing = 0.0_rlg
            IF (year == 2009) Zfishing = 0.0_rlg
            IF (year == 2010) Zfishing = 0.00115_rlg
            IF (year == 2011) Zfishing = 0.00083_rlg
            IF (year == 2012) Zfishing = 0.00055_rlg
            IF (year == 2013) Zfishing = 0.00101_rlg
            IF (year == 2014) Zfishing = 0.00128_rlg
            IF (year == 2015) Zfishing = 0.00119_rlg
            IF (year == 2016) Zfishing = 0.00098_rlg
            IF (year == 2017) Zfishing = 0.00178_rlg
            IF (year == 2018) Zfishing = 0.00160_rlg
            IF (year == 2019) Zfishing = 0.00135_rlg

            IF (particle%AgeClass == 1) Zfishing = Zfishing*0.46_rlg
        ENDIF

        IF (particle%stage >= 5 .and. month >= 7 .and. species == 'anchovy') THEN
            !   IF (year < 1990) Zfishing = f_spin
            !   IF (year == 1990) Zfishing = 0.00208_rlg
            !   IF (year == 1991) Zfishing = 0.00081_rlg
            !   IF (year == 1992) Zfishing = 0.00108_rlg 
            !   IF (year == 1993) Zfishing = 0.00163_rlg 
            !   IF (year == 1994) Zfishing = 0.00180_rlg
            !   IF (year == 1995) Zfishing = 0.00107_rlg
            !   IF (year == 1996) Zfishing = 0.00204_rlg
            !   IF (year == 1997) Zfishing = 0.00170_rlg
            !   IF (year == 1998) Zfishing = 0.00146_rlg
            !   IF (year == 1999) Zfishing = 0.00119_rlg
            IF (year <= 2000) Zfishing = 0.00110_rlg
            IF (year == 2001) Zfishing = 0.00144_rlg
            IF (year == 2002) Zfishing = 0.00144_rlg
            IF (year == 2003) Zfishing = 0.00184_rlg
            IF (year == 2004) Zfishing = 0.00177_rlg
            IF (year == 2005) Zfishing = 0.0_rlg
            IF (year == 2006) Zfishing = 0.000033_rlg
            IF (year == 2007) Zfishing = 0.0_rlg
            IF (year == 2008) Zfishing = 0.0_rlg
            IF (year == 2009) Zfishing = 0.0_rlg
            IF (year == 2010) Zfishing = 0.00054_rlg
            IF (year == 2011) Zfishing = 0.00019_rlg
            IF (year == 2012) Zfishing = 0.00043_rlg
            IF (year == 2013) Zfishing = 0.00033_rlg
            IF (year == 2014) Zfishing = 0.00041_rlg
            IF (year == 2015) Zfishing = 0.00046_rlg
            IF (year == 2016) Zfishing = 0.00029_rlg
            IF (year == 2017) Zfishing = 0.00025_rlg
            IF (year == 2018) Zfishing = 0.00028_rlg
            IF (year == 2019) Zfishing = 0.00026_rlg

            IF (particle%AgeClass == 1) Zfishing = Zfishing*1.035_rlg
        ENDIF

        IF (particle%stage == 6 .and. particle%AgeClass >= 2 .and. species == 'sardine') THEN 
           !IF (year < 2000) Zfishing = f_spin
           IF (year <= 2000) Zfishing = 0.0004247_rlg
           IF (year == 2001) Zfishing = 0.0004384_rlg
           IF (year == 2002) Zfishing = 0.0005151_rlg
           IF (year == 2003) Zfishing = 0.0004164_rlg
           IF (year == 2004) Zfishing = 0.0003973_rlg
           IF (year == 2005) Zfishing = 0.0003890_rlg
           IF (year == 2006) Zfishing = 0.0004274_rlg
           IF (year == 2007) Zfishing = 0.0004521_rlg
           IF (year == 2008) Zfishing = 0.0006301_rlg
           IF (year == 2009) Zfishing = 0.0005151_rlg
           IF (year == 2010) Zfishing = 0.0005096_rlg
           IF (year == 2011) Zfishing = 0.0006849_rlg
           IF (year == 2012) Zfishing = 0.0012055_rlg
           IF (year == 2013) Zfishing = 0.0013151_rlg
           IF (year == 2014) Zfishing = 0.0016164_rlg
           IF (year == 2015) Zfishing = 0.0013973_rlg
           IF (year == 2016) Zfishing = 0.0016986_rlg
           IF (year == 2017) Zfishing = 0.0016712_rlg
           IF (year == 2018) Zfishing = 0.0019178_rlg
           IF (year == 2019) Zfishing = 0.0014521_rlg
        ENDIF

        particle % Death_FISH = particle%Death_FISH + particle%super*(1 - exp(-(Zfishing*c_F)*dtm/86400_rlg))               
        particle % super      = particle%super*exp(- (Zfishing*c_F)*dtm/86400_rlg)

    ENDIF   ! fishing_strategy == F_eval



    ! ============================================================================
    ! =====                                                                  =====
    ! =====                     Fishing strategy : Catch                     =====
    ! =====                                                                  =====
    IF (particle%stage >= 5 .and. fishing_strategy == 'Catch' .and. particle%AgeClass >= 1 ) THEN !ageclass>=1 to be sure not to fish newborns THEN
        
        ! mat_catch is read in ibm_init routine
        IF (species == 'anchovy') id_species = 1
        IF (species == 'sardine') id_species = 2

        coeff1 = dtm/86400_rlg*c_F

        IF (biom_tot(id_species) .gt. 0.0_rlg) THEN
            !date_start : to avoid 1st t where biom and wdebmean =0
            IF(year > 1970 .and. year < 2000) THEN
                particle%Death_FISH = particle%Death_FISH + particle%Wdeb*particle%super*mat_catch(1,month,id_species)         &
                                        *coeff1/(Wdeb_mean(id_species)*biom_tot(id_species)) 
                particle%super = particle%super - particle%Wdeb*particle%super*mat_catch(1,month,id_species)*coeff1   &
                                /(Wdeb_mean(id_species)*biom_tot(id_species)) 
            ELSE IF (year >= 2000) THEN
                particle%Death_FISH = particle%Death_FISH + particle%Wdeb*particle%super*mat_catch(year-1999,month,id_species) &
                                        *coeff1/(Wdeb_mean(id_species)*biom_tot(id_species))
                particle%super = particle%super - particle%Wdeb*particle%super*mat_catch(year-1999,month,id_species)             &
                                 *coeff1/(Wdeb_mean(id_species)*biom_tot(id_species))
            ENDIF
        ENDIF
    ENDIF   ! fishing_strategy == Catch


    ! ============================================================================
    ! =====                                                                  =====
    ! =====                      Fishing strategy : HCR                      =====
    ! =====                                                                  =====
    IF (fishing_strategy == 'HCR') THEN ! Attention: time step is second!!!! - Bueno Pardo
        IF (year == 1996 .and. month < 7 .and. particle%stage >= 5)  Zfishing = 0.0033_rlg
        IF (year == 1997 .and. month < 7 .and. particle%stage >= 5)  Zfishing = 0.0017_rlg
        IF (year == 1998 .and. month < 7 .and. particle%stage >= 5)  Zfishing = 0.0013_rlg
        IF (year == 1999 .and. month < 7 .and. particle%stage >= 5)  Zfishing = 0.0015_rlg

        IF (year == 1996 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0017_rlg
        IF (year == 1997 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0015_rlg
        IF (year == 1998 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0012_rlg
        IF (year == 1999 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0010_rlg

        !! FUTURE SCENARIOS:
        IF (year >= 2036 .and. year < 2040 .and. month < 7  .and. particle%stage >= 5) Zfishing = 0.0008_rlg
        IF (year >= 2036 .and. year < 2040 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0008_rlg
        IF (year >= 2076 .and. year < 2080 .and. month < 7  .and. particle%stage >= 5) Zfishing = 0.0008_rlg
        IF (year >= 2076 .and. year < 2080 .and. month >= 7 .and. particle%stage == 6) Zfishing = 0.0008_rlg

        particle%Death_FISH = particle%Death_FISH + particle%super*(1 - exp(-(Zfishing*percent)*dtm/86400_rlg))               
        particle%super      = particle%super*exp(-(Zfishing*percent)*dtm/86400_rlg)                   
    ENDIF ! fishing_strategy == HCR

    ! ============================================================================
    ! =====                                                                  =====
    ! =====                   Fishing strategy : historical                  =====
    ! =====                                                                  =====

    IF (fishing_strategy == 'historical') THEN
                
        !! Average fishing mortality first semester (between 1988 and 2015, without 0's): 0.001953 (half of this = 0.0009765)
        !! Average fishing mortality second semester (between 1988 and 2015, without 0's): 0.00099 (half of this = 0.000495)
        !! Average fishing mortality first semester (between 2000 and 2015, without 0's): 0.001207
        !! Average fishing mortality second semester (between 2000 and 2015, without 0's): 0.000825
                
        !! Commenté Par Clara pour test à 2100
        !IF (year > 2030) THEN
        !   IF (month < 7) THEN
        !      IF (particle % stage == 6 .OR. particle % stage == 5) THEN
        !         particle % number = particle % number * exp(- (0.001953_rlg*percent) * dtm/86400_rlg) !see /IBMDEB/TimeSeriesAnchovy/CBBMdata.ods (juan)
         !     ENDIF
         !  ENDIF
        !ENDIF
                
        !IF (year > 2030) THEN
        !   IF (month >= 7) THEN
        !      IF (particle % stage == 6) THEN
        !         particle % number = particle % number * exp(- (0.00099_rlg*percent) * dtm/86400_rlg) !see /IBMDEB/TimeSeriesAnchovy/CBBMdata.ods (juan)
        !      ENDIF
        !   ENDIF
        !ENDIF
                      
        IF (month < 7 .and. particle%stage >= 5) THEN
            ! Zfishing = 0.001207_rlg
            IF (year == 1980)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1981)  Zfishing = 0.0021_rlg !average of the first 10 years  
            IF (year == 1982)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1983)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1984)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1985)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1986)  Zfishing = 0.0021_rlg !average of the first 10 years
            IF (year == 1987)  Zfishing = 0.0033_rlg
            IF (year == 1988)  Zfishing = 0.0027_rlg
            IF (year == 1989)  Zfishing = 0.0025_rlg
            IF (year == 1990)  Zfishing = 0.0033_rlg
            IF (year == 1991)  Zfishing = 0.0031_rlg
            IF (year == 1992)  Zfishing = 0.0032_rlg 
            IF (year == 1993)  Zfishing = 0.0023_rlg 
            IF (year == 1994)  Zfishing = 0.0031_rlg
            IF (year == 1995)  Zfishing = 0.0041_rlg
            IF (year == 1996)  Zfishing = 0.0033_rlg
            IF (year == 1997)  Zfishing = 0.0017_rlg
            IF (year == 1998)  Zfishing = 0.0013_rlg
            IF (year == 1999)  Zfishing = 0.0015_rlg
            IF (year == 2000)  Zfishing = 0.0019_rlg
            IF (year == 2001)  Zfishing = 0.0018_rlg
            IF (year == 2002)  Zfishing = 0.0015_rlg
            IF (year == 2003)  Zfishing = 0.0010_rlg
            IF (year == 2004)  Zfishing = 0.0023_rlg
            IF (year == 2005)  Zfishing = 0.0004_rlg
            IF (year == 2006)  Zfishing = 0.0007_rlg
            IF (year == 2007)  Zfishing = 0.0000_rlg
            IF (year == 2008)  Zfishing = 0.0000_rlg
            IF (year == 2009)  Zfishing = 0.0000_rlg
            IF (year == 2010)  Zfishing = 0.0011_rlg
            IF (year == 2011)  Zfishing = 0.0008_rlg
            IF (year == 2012)  Zfishing = 0.0005_rlg
            IF (year == 2013)  Zfishing = 0.0010_rlg
            IF (year == 2014)  Zfishing = 0.0013_rlg
            IF (year == 2015)  Zfishing = 0.0011_rlg

            IF (year > 2015) Zfishing = 0.001207_rlg !! Fishing scenario: CONSTANT-FISHING (test Clara à 2100)
                       
            particle%Death_FISH = particle%Death_FISH + particle%super*(1 - exp(-(Zfishing*percent)*dtm/86400_rlg))               
            particle%super      = particle%super*exp(-(Zfishing*percent)*dtm/86400_rlg)

        ELSE IF (month >= 7 .and. particle%stage == 6) THEN
            ! Zfishing = 0.000825_rlg
            IF (year == 1980)  Zfishing = 0.0009_rlg !average of the first 10 years
            IF (year == 1981)  Zfishing = 0.0009_rlg !average of the first 10 years  
            IF (year == 1982)  Zfishing = 0.0009_rlg !average of the first 10 years                    
            IF (year == 1983)  Zfishing = 0.0009_rlg !average of the first 10 years
            IF (year == 1984)  Zfishing = 0.0009_rlg !average of the first 10 years
            IF (year == 1985)  Zfishing = 0.0009_rlg !average of the first 10 years
            IF (year == 1986)  Zfishing = 0.0009_rlg !average of the first 10 years
            IF (year == 1987) Zfishing = 0.0009_rlg                
            IF (year == 1988) Zfishing = 0.0009_rlg        
            IF (year == 1989) Zfishing = 0.0005_rlg
            IF (year == 1990) Zfishing = 0.0018_rlg
            IF (year == 1991) Zfishing = 0.0007_rlg 
            IF (year == 1992) Zfishing = 0.0009_rlg 
            IF (year == 1993) Zfishing = 0.0014_rlg 
            IF (year == 1994) Zfishing = 0.0015_rlg
            IF (year == 1995) Zfishing = 0.0009_rlg
            IF (year == 1996) Zfishing = 0.0017_rlg  
            IF (year == 1997) Zfishing = 0.0015_rlg
            IF (year == 1998) Zfishing = 0.0012_rlg
            IF (year == 1999) Zfishing = 0.0010_rlg
            IF (year == 2000) Zfishing = 0.0010_rlg
            IF (year == 2001) Zfishing = 0.0013_rlg                
            IF (year == 2002) Zfishing = 0.0013_rlg                
            IF (year == 2003) Zfishing = 0.0016_rlg                
            IF (year == 2004) Zfishing = 0.0016_rlg
            IF (year == 2005) Zfishing = 0.0000_rlg                
            IF (year == 2006) Zfishing = 0.0000_rlg
            IF (year == 2007) Zfishing = 0.0000_rlg
            IF (year == 2008) Zfishing = 0.0000_rlg
            IF (year == 2009) Zfishing = 0.0000_rlg
            IF (year == 2010) Zfishing = 0.0005_rlg
            IF (year == 2011) Zfishing = 0.0002_rlg
            IF (year == 2012) Zfishing = 0.0004_rlg
            IF (year == 2013) Zfishing = 0.0003_rlg
            IF (year == 2014) Zfishing = 0.0003_rlg                
            IF (year == 2015) Zfishing = 0.0004_rlg

            IF (year > 2015) Zfishing = 0.000825_rlg !! Fishing scenario: CONSTANT-FISHING (test Clara à 2100)

            particle%Death_FISH = particle%Death_FISH + particle%super*(1 - exp(-(Zfishing*percent)*dtm/86400_rlg))                
            particle%super      = particle%super*exp(-(Zfishing*percent)*dtm/86400_rlg)
        ENDIF
    ENDIF ! fishing_strategy == historical

  END SUBROUTINE death_by_fishing
#endif  /* IBM_SPECIES */



  !!===========================================================================
  FUNCTION w_dens(tempw,salw)
    !&E---------------------------------------------------------------------
    !&E ** Purpose : Estimate the water density
    !&E ** 
    !&E ** Description    : 
    !&E ** Called by      : ibm_buoy, ibm_init
    !&E ** External calls : 
    !&E ** Reference      : Mellor, 1991
    !&E
    !&E ** History        : 
    !&E       ! 2010-11  (M. Huret)
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh)   , INTENT( in )  :: salw,tempw
    REAL(KIND=rsh)                   :: w_dens

    !! * Local declarations
    REAL(KIND=rsh)                  :: TP,S,P
    REAL(KIND=rlg)                  :: row,rost,c,pcc

    !!----------------------------------------------------------------------
    !! * Executable part
    ! water density (Mellor, 1991)
    P = 10.0_rsh                  ! pressure in decibars = 1000mbar reference atmospherical pressure

    TP = MAX(tempw,0.001_rsh)
    S  = MAX(salw,0.0001_rsh)

    ! density of pure water (UNESCO)
    row = 999.842594_rlg + (6.793952d-02 + (-9.095290d-03 + (1.001685d-04 - 1.120083d-06*TP + 6.536332d-09*TP*TP)*TP)*TP)*TP

    ! density of sea water at the surface (UNESCO)
    rost = row + S*(0.824493_rlg + (-4.0899d-03 + (7.6438d-05 + (-8.2467d-07 + 5.3875d-09*TP)*TP)*TP)*TP)    &
           + S**1.5_rlg*( -5.72466d-03 + (1.0227d-04 - 1.6546d-06*TP)*TP) + S*S*4.8314d-04

    ! density of sea water at pressure p (decibars)
    c      = 1449.2_rlg + 1.34_rlg*(S-35.0_rlg) + 4.55_rlg*TP - 0.045_rlg*TP*TP + 0.00821_rlg*p + 15.0d-09*p*p
    pcc    = p/(c*c)
    w_dens = rost + 1.0d4*(1.0_rlg-0.2_rlg*pcc)*pcc
  END FUNCTION w_dens



  !!===========================================================================
  FUNCTION ibm_buoy(density,size,tempw,salw)

    !&E---------------------------------------------------------------------
    !&E ** Purpose : Estimate the terminal velocity from buoyancy scheme
    !&E ** 
    !&E ** Description    : 
    !&E ** Called by      : ibm_3d
    !&E ** External calls : w_dens
    !&E ** Reference      : Adlandsvik, 2000
    !&E
    !&E ** History        : 
    !&E       ! 2010-11  (M. Huret)
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh)   , INTENT( in ) :: salw,tempw,density,size
    REAL(KIND=rsh)                  :: ibm_buoy

    !! * Local declarations
    REAL(KIND=rsh)                  :: w_buoy
    REAL(KIND=rsh)                  :: wsto,wdal,cor,rayoeuf,diam,difdens
    REAL(KIND=rsh)                  :: TP,S,wdens,dynvisc,reyn
    
    !!----------------------------------------------------------------------
    !! * Executable part
    cor     = 1.0_rsh
    rayoeuf = size*0.01_rsh/2.0_rsh    ! diameter to radius

    TP = MAX(tempw,0.001_rsh)
    S  = MAX(salw,0.001_rsh)
    
    wdens = w_dens(tempw,salw)
    
    ! Dynamic viscosity (kg.m-1.s-1) (approx. of table values, vertegg v1.0 toolbox, Adlandsvik, 2000)
    dynvisc = (1.7915_rsh - 0.0538_rsh*TP + 0.0007_rsh*TP*TP + 0.0023_rsh*S)*0.001_rsh
    
    !  calcul de la vitesse de remontee des oeufs de poissons (Stokes)             
    !  unites: m s-1
    difdens=density-REAL(wdens,rsh)
    wsto   = 2.0_rsh*g*rayoeuf**2.0_rsh*difdens/(9.0_rsh*dynvisc)
    wsto   = wsto*cor
    w_buoy = -wsto

    ! calcul du nombre de Reynolds (Re)
    reyn = (wdens*abs(wsto)*2.0_rsh*rayoeuf)/dynvisc
    
    ! choix de la vitesse de Dallavalle si Re > 0.5     
    IF (reyn.gt.0.5_rsh) THEN
        diam   = (9.0_rsh*(dynvisc**2)/(wdens*g*abs(difdens)))**(1.0_rsh/3.0_rsh)      
        wdal   = 0.0875_rsh*abs(difdens)**(2.0_rsh/3.0_rsh)*((2.0_rsh*rayoeuf) - 0.4_rsh*diam)/dynvisc**(1.0_rsh/3.0_rsh)
        w_buoy = wdal*cor
        IF(difdens .ge. 0.0_rsh) THEN
            w_buoy = -wdal*cor
        ENDIF
    ENDIF
    
    ibm_buoy = w_buoy
    
  END FUNCTION ibm_buoy



  !!======================================================================
  FUNCTION ibm_nycth_mig(depth1,depth2,w1,w2,p,alpha)
    !&E---------------------------------------------------------------------
    !&E
    !&E ** Purpose : Nycthemeral migration between depth1 (day) and depth2 (night)
    !&E **           Migration from depth1 to depth2 starts when rad=0 with w=w2
    !&E **           Migration from depth2 to depth1 starts when rad>0 with w=w1
    !&E **
    !&E ** Description    :
    !&E ** Called by      : ibm_3d
    !&E ** External calls : define_pos
    !&E ** Reference      :
    !&E ** 
    !&E ** History        : 
    !&E       ! 2010-11  (M. Huret)
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE trajectools,        ONLY : define_pos
    USE comtraj,            ONLY : type_position
 
    !! * Arguments
    REAL(KIND=rsh), INTENT( in )          :: depth1,depth2,w1,w2,alpha
    TYPE(type_particle), INTENT(inout)    :: p
    REAL(KIND=rsh)                        :: ibm_nycth_mig
 
    !! * Local declarations
    REAL(kind=rsh)                        :: depth,w
    character(len=19)                     :: tool_sectodat,date
    integer                               :: jj,mm,aaaa,hh,minu,sec
    TYPE(type_position)                   :: pos  
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    depth = depth1
    w     = w1

    ! Switch to local coordinates inside each MPI domain
    pos%xp = p%xpos ; pos%yp = p%ypos 
    CALL define_pos(pos)

    IF(srflx(NINT(pos%idx_r),NINT(pos%idy_r)) == 0.0_rsh) THEN
        depth = depth2
        w     = w2
    ENDIF
    IF (depth .lt. p%d3) THEN  
        ! p%zpos is immersion of particle
        ibm_nycth_mig = w*tanh(alpha*(p%zpos - depth)) 
    ELSE
        ibm_nycth_mig = 0.
    ENDIF
   
  END FUNCTION ibm_nycth_mig



  !!===============================================================================
  FUNCTION ibm_traint(varint,xe,spos,khi,klo,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,&
                     limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION traint  ***
    !&E
    !&E ** Purpose : Linear interpolation of tracer values 
    !&E              at the location (xpos,ypos,zpos)
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_3d,ibm_parameter_init
    !&E ** External calls : h_int_siggen,f_lag_sigz_tra,ibm_interp,varint
    !&E ** Reference      :
    !&E **
    !&E ** History :
    !&E       !  2002-08 (F. Dumas, P. Lazure) 
    !&E       !  2010-11 (M. Huret) 
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,     ONLY : dsigw,ierrorlog
    USE trajectools, ONLY : h_int_siggen


    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )  :: varint
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),      INTENT( in )  :: xe

    REAL(KIND=rsh), INTENT( in )                :: spos,px,py
    INTEGER, INTENT( in )                       :: khi,klo
    INTEGER, INTENT( in )                       :: i1,i2,j1,j2
    INTEGER, INTENT( in )                       :: rh11,rh21,rh12,rh22
    INTEGER, INTENT( in )                       :: limin,limax,ljmin,ljmax

    REAL(KIND=rsh)                              :: ibm_traint

    !! * Local declarations
    REAL(KIND=rsh)                              :: traintm,traintp
    REAL(KIND=rsh)                              :: xe_lag,hc_sig_lag,h0_lag
    REAL(KIND=rsh)                              :: f_lag_hzp,f_lag_hzm,f_lag_hz
    REAL(KIND=rsh), DIMENSION(2,2)              :: val4
    INTEGER                                     :: IERR_MPI



    !!----------------------------------------------------------------------
    !! * Executable part
    IF ((khi>kmax+1) .or. (klo>kmax+1) .or. (khi<0) .or. (klo<0)) THEN
        WRITE(ierrorlog,*) 'Function IBM_TRAINT'
        WRITE(ierrorlog,*) 'bad layer from sigma'
        WRITE(ierrorlog,*) 'The simulation is stopped'
        CALL_MPI MPI_FINALIZE(ierr_mpi)
        STOP
    END IF

    ibm_traint = 0.0_rsh
    traintm    = 0.0_rsh   
    traintp    = 0.0_rsh

    !-- interpolation proprement dite sur la nappe sigma au-dessus
    IF (khi /= kmax+1) THEN
        val4(:,:)  = varint(i1:i2,j1:j2,khi)
        traintp    = ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)
    ELSE
        val4(:,:)  = varint(i1:i2,j1:j2,klo)
        ibm_traint = ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)
    END IF
    !-- interpolation proprement dite sur la nappe sigma au-dessous
    IF (klo /= 0) THEN
        val4(:,:)  = varint(i1:i2,j1:j2,klo)
        traintm    = ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)
    ELSE
        val4(:,:)  = varint(i1:i2,j1:j2,khi)
        ibm_traint = ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)
    END IF
    IF (khi /= kmax+1 .and. klo /= 0) THEN
    !-- interpolation sur la verticale entre la valeur au-dessus (traintp)
    !-- et la valeur au-dessous (traintm).
        f_lag_hzp = 1.0_rsh
        f_lag_hzm = 1.0_rsh
        f_lag_hz  = 1.0_rsh
        CALL h_int_siggen(xe_lag,h0_lag,hc_sig_lag,xe,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,&
                          limin,limax,ljmin,ljmax)
        CALL f_lag_sigz_tra(spos,klo,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)
        ibm_traint=(f_lag_hzp*(sc_r(khi) - spos)*traintm + f_lag_hzm*(spos - sc_r(klo))*traintp) &
                   /(f_lag_hz*dsigw(klo)) 
    ENDIF

  END FUNCTION ibm_traint



  !!======================================================================
  FUNCTION ibm_proftraint(var,minprof,maxprof,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION proftraint  ***
    !&E
    !&E ** Purpose :Interpolate on the horizontal and average over minprof - maxprof
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_3d
    !&E ** External calls : ibm_profmean
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       ! 2013-04  (M. Huret) 
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,    ONLY : hc_sig


    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )  :: var
    REAL(KIND=rsh),                                  INTENT( in )  :: px,py,minprof, maxprof
    INTEGER,                                         INTENT( in )  :: i1,i2,j1,j2
    INTEGER,                                         INTENT( in )  :: rh11,rh21,rh12,rh22

    REAL(KIND=rsh)                                                 :: ibm_proftraint
    !! * Local declarations
    REAL(KIND=rsh), DIMENSION(2,2)                                 :: val4
    
    !!----------------------------------------------------------------------
    !! * Executable part
    val4(:,:) = 0.0_rsh

    ! Processing prof mean on 4 neighbooring locations
    IF (rh11 == 1) THEN
        val4(1,1) = ibm_profmean(var(i1,j1,:),minprof,maxprof,h(i1,j1),zeta(i1,j1,kstp),hc_sig(i1,j1))
    END IF

    IF (rh21 == 1) THEN
        val4(2,1) = ibm_profmean(var(i2,j1,:),minprof,maxprof,h(i2,j1),zeta(i2,j1,kstp),hc_sig(i2,j1))
    END IF

    IF (rh22 == 1) THEN
        val4(2,2) = ibm_profmean(var(i2,j2,:),minprof,maxprof,h(i2,j2),zeta(i2,j2,kstp),hc_sig(i2,j2))
    END IF

    IF (rh12 == 1) THEN
        val4(1,2) = ibm_profmean(var(i1,j2,:),minprof,maxprof,h(i1,j2),zeta(i1,j2,kstp),hc_sig(i1,j2))
    END IF

    ! Interpolation
    ibm_proftraint = ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)

  END FUNCTION



  !!======================================================================
  FUNCTION ibm_profmean(var,minprof,maxprof,h,xe,hc_sig)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION ibm_profmean  ***
    !&E
    !&E ** Purpose : Average of values of var over a given prof range (minprof, maxprof)
    !&E               minprof and maxprof are provided as immersion values (positive)
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_proftraint,ibm_profuint,ibm_profvint,
    !&E ** External calls : siggentoz, ksupkinf 
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       ! 2013-04  (M. Huret) 
    !&E       ! 2024     (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE trajectools, ONLY : ksupkinf,siggentoz
 
    !! * Arguments
    REAL(KIND=rsh), DIMENSION(kmax), INTENT( in )  :: var
    REAL(KIND=rsh), INTENT( in )                   :: minprof, maxprof,h,xe
    REAL(KIND=rsh), INTENT( in )                   :: hc_sig

    REAL(KIND=rsh)                                 :: ibm_profmean
 
    !! * Local declarations
    REAL(KIND=rsh)                                 :: zmin, zmax, valint,d3
    REAL(KIND=rsh), DIMENSION(0:kmax)              :: zpos
    INTEGER                                        :: kw,kwm,l,ksub,ksup
    
    !---------------------------------------------------------
    !! * Executable part
    valint = 0.0_rsh 
 
    ! switch from immersion to real z
    d3   = h+xe
    zmin = -minprof+xe  
    zmax = -maxprof+xe
 
   ! To make sure minprof is smaller than maxprof
    IF (minprof >= maxprof) THEN 
       write(*,*) 'attention minprof > maxprof', minprof, maxprof, d3
    end IF
 
    ! Convert to z because we have minprof and maxprof in z
    DO l=0,kmax
        CALL siggentoz(zpos(l),sc_w(l),xe,h,hc_sig)
    END DO
 
    !  Calcul de la prof sigma min et max
    IF (minprof >= d3) THEN
        ksup = 1
        zmin = zpos(1)
    ELSE
        kw=kmax
        kwm=0
        CALL ksupkinf(zmin,zpos,kmax+1,kw,kwm,10)
        ksup=kw
    END IF
    
    IF (maxprof >= d3) THEN
        ksub = 0
        zmax = zpos(0)
    ELSE 
        kw  = kmax
        kwm = 0
        CALL ksupkinf(zmax,zpos,kmax+1,kw,kwm,11)
        ksub = kwm
    END IF
 
    ! Cumul des valeurs entre minprof et maxprof
    IF (ksup == ksub+1) THEN
       valint = var(ksup)*(zmin-zmax)
    ELSE IF (ksup == ksub+2) THEN
       valint = valint + var(ksup)*(zmin - zpos(ksup-1))
       valint = valint + var(ksup-1)*(zpos(ksup-1) - zmax)
    ELSE
       valint = valint + var(ksup)*(zmin - zpos(ksup-1))
       DO l = ksub+2,ksup-1 
            valint = valint + var(l)*(zpos(l) - zpos(l-1))
       END DO
       valint = valint + var(ksub+1)*(zpos(ksub+1) - zmax)
    ENDIF
 
    ! Average over minprof - maxprof
    ibm_profmean = valint/(zmin-zmax)
 
  END FUNCTION ibm_profmean



  !!======================================================================
  SUBROUTINE ibm_loc_xyz(xpos,ypos,spos,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,khi,klo,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE location  ***
    !&E
    !&E ** Purpose : search index of the cells around particle and location of the particle 
    !&E              inside the cell for later interpolation 
    !&E              Valid relative to h0,xe,wkz... 
    !&E
    !&E ** Description    : 
    !&E ** Called by      : ibm_3d,ibm_parameter_init
    !&E ** External calls : ksupkinf
    !&E ** Reference      : 
    !&E
    !&E ** History :
    !&E       !  2012-07  (M. Huret)
    !&E       !  2011-11  (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       ! 2024      (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,     ONLY : valmanq
    USE trajectools, ONLY : ksupkinf

    ! * Arguments
    REAL(KIND=rsh),                 INTENT( in )   :: xpos,ypos,spos
    REAL(KIND=rsh),                 INTENT( out )  :: px,py 
    INTEGER,                        INTENT( out )  :: i1,i2,j1,j2
    INTEGER,                        INTENT( out )  :: rh11,rh21,rh12,rh22
    INTEGER,                        INTENT( out )  :: khi,klo
    INTEGER,                        INTENT( in )   :: limin,limax,ljmin,ljmax

    !! * Local declarations
    REAL(KIND=rsh)                                 :: xe_lag,h0_lag
    REAL(KIND=rsh), DIMENSION(0:kmax+1)            :: spos_w  
   
    !!----------------------------------------------------------------------
    !! * Executable part
    
    ! indexes of the cell where the particle is
    ! reperage de la maille XE ou se situe la particule
    i1 = INT(xpos) 
    i2 = INT(xpos) + 1   
    j1 = INT(ypos) 
    j2 = INT(ypos) + 1
     
    i1 = MIN(MAX(i1,limin-1),limax+1)
    i2 = MIN(MAX(i2,limin-1),limax+1)
    j1 = MIN(MAX(j1,ljmin-1),ljmax+1)
    j2 = MIN(MAX(j2,ljmin-1),ljmax+1)


    ! location of the particle inside the cell
    ! calcul des positions px et py dans la maille XE
    px = xpos - INT(xpos)
    py = ypos - INT(ypos)
    
    ! check for next cells in water
    rh11 = 1
    rh12 = 1
    rh21 = 1
    rh22 = 1

    ! If rmask = 0 : means the cell is on land. If rmask = 1, means cell is on sea
#ifdef MASKING
    IF (rmask(i1,j1)==0)  rh11 = 0
    IF (rmask(i1,j2)==0)  rh12 = 0
    IF (rmask(i2,j1)==0)  rh21 = 0
    IF (rmask(i2,j2)==0)  rh22 = 0
#endif
    
    ! get neighbooring indices on the vertical
    spos_w(0)      = -1.0_rsh
    spos_w(kmax+1) = 0.0_rsh
    spos_w(1:kmax) = sc_r(1:kmax)
    
    klo = 0
    khi = kmax+1
    CALL ksupkinf(spos,spos_w,kmax+2,khi,klo,6)

  END SUBROUTINE ibm_loc_xyz



  !!====================================================================
  SUBROUTINE gasdev_s(harvest)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE gasdev_s  ***
    !&E
    !&E ** Purpose : Returns in harvest a normally distributed deviate with zero mean and unit variance, 
    !&E              using ran1 as the source of uniform deviates.
    !&E
    !&E ** Description    : 
    !&E ** Called by      : ibm_parameter_init,deb_egg_init,deb_init
    !&E ** External calls : 
    !&E ** Reference      : 
    !&E
    !&E ** History :
    !&E
    !&E---------------------------------------------------------------------

    REAL(rsh), INTENT(OUT) :: harvest

    REAL(rsh)              :: rsq,v1,v2
    REAL(rsh), SAVE        :: g
    LOGICAL,   SAVE        :: gaus_stored = .false.
    
    IF (gaus_stored) THEN       ! We have an extra deviate handy,
       harvest     = g          ! so return it,
       gaus_stored = .false.    ! and unset the flag.
    
    ELSE        !We dont have an extra deviate handy, so
        DO
            CALL random_number(v1)      ! pick two uniform numbers in the square extending from 
            CALL random_number(v2)      !-1 to +1 in each direction,
            v1  = 2.0_rsh*v1 - 1.0_rsh
            v2  = 2.0_rsh*v2 - 1.0_rsh
            rsq = v1**2 + v2**2         !see IF they are in the unit circle,
            IF (rsq > 0.0 .and. rsq < 1.0) EXIT
        END Do       !otherwise try again.

        !Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
        rsq     = sqrt(-2.0_rsh*log(rsq)/rsq)  
        harvest = v1*rsq
        g       = v2*rsq
        gaus_stored = .true.               ! Set flag.
    END IF

  END SUBROUTINE gasdev_s




  ! ######################################################################
  ! #####
  ! #####                   Interne au module ibmtools

  !!======================================================================
  FUNCTION ibm_interp(val4,px,py,rh11,rh21,rh12,rh22)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE ibm_interp  ***
    !&E
    !&E ** Purpose : Horizontal interpolation with verification of boundary limits
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_traint,ibm_proftraint
    !&E ** External calls :
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !  2011-01  (M. Huret)
    !&E       ! 2024      (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
 
    !! * Arguments
    REAL(KIND=rsh), DIMENSION(2,2), INTENT( in )  :: val4
    REAL(KIND=rsh),                 INTENT( in )  :: px,py
    INTEGER,                        INTENT( in )  :: rh11,rh21,rh12,rh22

    REAL(KIND=rsh)                                :: ibm_interp
 
    !! * Local declarations
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! check IF boundary condition
    IF (rh11 > 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        ibm_interp = px*(1-py)*val4(2,1) + px*py*val4(2,2) + (1-px)*py*val4(1,2) + (1-px)*(1-py)*val4(1,1)
 
    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        ibm_interp = 0
 
    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        ibm_interp = val4(2,2)
 
    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        ibm_interp = val4(1,2)
 
    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        ibm_interp = px*val4(2,2) + (1-px)*val4(1,2)
       
    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        ibm_interp = val4(2,1)
       
    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        IF (px > 0.5_rsh) THEN
            ibm_interp = val4(2,1) 
        ELSE
            ibm_interp = val4(1,2)
        END IF
       
    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        ibm_interp = (1-py)*val4(2,1) + py*val4(2,2)
       
    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        ibm_interp = px*(val4(2,2) - val4(1,2)) + py*(val4(2,2)-val4(2,1)) - val4(2,2) + val4(1,2) + val4(2,1)
       
    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        ibm_interp = val4(1,1)
       
    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        IF (px > 0.5_rsh) THEN
            ibm_interp = val4(2,2)
        ELSE
            ibm_interp = val4(1,1)
        END IF
       
    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        ibm_interp = px*val4(2,1)+(1-px)*val4(1,1)
       
    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        ibm_interp = px*(val4(2,1)-val4(1,1)) + py*(val4(2,2)-val4(2,1)) + val4(1,1)
       
    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        ibm_interp = py*val4(1,2) + (1-py)*val4(1,1)
       
    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        ibm_interp = px*(val4(2,1)-val4(1,1)) + py*(val4(1,2)-val4(1,1)) + val4(1,1)
       
    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        ibm_interp = px*(val4(2,2)-val4(1,2)) + py*(val4(1,2)-val4(1,1)) + val4(1,1)
    END IF
    
  END FUNCTION ibm_interp



  !!======================================================================
  SUBROUTINE f_lag_sigz_tra(spos,k,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE f_lag  ***
    !&E
    !&E ** Purpose : Scale factor from generalized sigma to z
    !&E              for vertical interpolation at particle location
    !&E
    !&E ** Description    :
    !&E ** Called by      : ibm_traint
    !&E ** External calls :
    !&E ** Reference      :
    !&E
    !&E ** History : 
    !&E       ! 2010-12   (M. Huret)
    !&E       ! 2024      (D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,     ONLY : dcwsds
    USE trajectools, ONLY : CSF
 
    !! * Arguments
    REAL(KIND=rsh), INTENT( in ) :: spos,xe_lag,h0_lag,hc_sig_lag
    REAL(KIND=rsh), INTENT(out ) :: f_lag_hzp,f_lag_hzm,f_lag_hz
    INTEGER,        INTENT( in ) :: k
 
    !! * Local declarations
    REAL(KIND=rsh)               :: Csusig,dcwsdsp,dcwsdsm,hinv
    
    !!----------------------------------------------------------------------
    !! * Executable part
    Csusig=CSF(spos,theta_s,theta_b)
    hinv=1./(h0_lag+hc_sig_lag)
    dcwsdsp = 0.0_rsh
    dcwsdsm = 0.0_rsh
 
    IF (spos /= sc_r(k+1)) dcwsdsp = (Cs_r(k+1) - Csusig)/(sc_r(k+1) - spos)
    IF (spos /= sc_r(k))   dcwsdsm = (Csusig - Cs_r(k))/(spos - sc_r(k))
 
    f_lag_hzp = hinv*(hc_sig_lag + (h0_lag-hc_sig_lag)*dcwsdsp)
    f_lag_hzm = hinv*(hc_sig_lag + (h0_lag-hc_sig_lag)*dcwsdsm)
    f_lag_hz  = hinv*(hc_sig_lag + (h0_lag-hc_sig_lag)*dcwsds(k))
 
  END SUBROUTINE f_lag_sigz_tra





 !!======================================================================================================
 !!======================================================================================================
 !!======================================================================================================

#ifdef key_ibm_unused

 ! #####################################################################
 ! #####     Reprendre les variables hx_g, hcx_sig de MARS et      #####
 ! #####    les adapter a CROCO pour faire fonctionner profuint    #####
 ! #####################################################################


 !!======================================================================
 FUNCTION ibm_profuint(var,xpos,ypos,minprof,maxprof)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION ibm_profuint  ***
    !&E
    !&E ** Purpose : Interpolate vertically averaged current 
    !&E              from the 4 neighbooring location 
    !&E
    !&E ** Description    : 
    !&E ** Called by      :
    !&E ** External calls : interpvit,ibm_profmean
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !  2013-04  (M. Huret) 
    !&E
    !&E---------------------------------------------------------------------
   !! * Modules used
   USE comtraj,  ONLY : htx
   USE trajectools, ONLY : interpvit
   !USE comvars2d, ONLY : hx_g
   !USE comsiggen, ONLY : hcx_sig

   !! * Arguments
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )  :: var
   REAL(KIND=rsh), INTENT( in )                                   :: xpos,ypos
   REAL(KIND=rsh), INTENT( in )                                   :: minprof, maxprof
   REAL(KIND=rsh)                                                 :: ibm_profuint

   !! * Local declarations
   REAL(KIND=rsh), DIMENSION(2,2)                                 :: val4
   INTEGER                                                        :: igg,idd,jhh,jbb
   INTEGER                                                        :: hlb,hlt,hrb,hrt  
   REAL(KIND=rsh)                                                 :: xe,px,py


   !!----------------------------------------------------------------------
   !! * Executable part

   !---- reperage de la maille ou se situe la particule pour les U
   ! location of ssh : xpos=REAL(i), ypos=REAL(j)
   igg=nint(xpos)
   idd=igg+1
   jbb=int(ypos)
   jhh=jbb+1

   idd=min(max(idd,imin),imax)
   igg=min(max(igg,imin),imax)
   jbb=min(max(jbb,jmin),jmax)
   jhh=min(max(jhh,jmin),jmax)

   ! calcul des positions px et py dans la maille U
   px=xpos-real(igg)-0.5_rsh
   py=ypos-real(jbb)


   val4(:,:)=0.0_rsh

   ! Processing prof mean on 4 neighbooring locations
   IF (htx(igg,jbb)>0.0_rsh) THEN
      xe=htx(igg,jbb)-hx_g(igg,jbb)
      val4(1,1)=ibm_profmean(var(:,igg,jbb),minprof,maxprof,hx_g(igg,jbb),xe,hcx_sig(igg,jbb))
   END IF

   IF (htx(idd,jbb)>0.0_rsh) THEN
      xe=htx(idd,jbb)-hx_g(idd,jbb)
      val4(2,1)=ibm_profmean(var(:,idd,jbb),minprof,maxprof,hx_g(idd,jbb),xe,hcx_sig(idd,jbb))
   END IF

   IF (htx(idd,jhh)>0.0_rsh) THEN
      xe=htx(idd,jhh)-hx_g(idd,jhh)
      val4(2,2)=ibm_profmean(var(:,idd,jhh),minprof,maxprof,hx_g(idd,jhh),xe,hcx_sig(idd,jhh))
   END IF

   IF (htx(igg,jhh)>0.0_rsh) THEN
      xe=htx(igg,jhh)-hx_g(igg,jhh)
      val4(1,2)=ibm_profmean(var(:,igg,jhh),minprof,maxprof,hx_g(igg,jhh),xe,hcx_sig(igg,jhh))
   END IF

   ! Interpolation
   ibm_profuint=interpvit(px,py,val4(1,1),val4(2,1),val4(1,2),val4(2,2), &
                          htx(igg,jbb),htx(idd,jbb),htx(igg,jhh),htx(idd,jhh))

 END FUNCTION



 ! #####################################################################
 ! #####     Reprendre les variables hy_g, hcy_sig de MARS et      #####
 ! #####    les adapter a CROCO pour faire fonctionner profvint    #####
 ! #####################################################################

 !!======================================================================
 FUNCTION ibm_profvint(var,xpos,ypos,minprof,maxprof)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION ibm_profvint  ***
    !&E
    !&E ** Purpose :Interpolate vertically averaged current 
    !&E             from the 4 neighbooring location 
    !&E
    !&E ** Description    : 
    !&E ** Called by      :
    !&E ** External calls : interpvit,ibm_profmean
    !&E ** Reference      :
    !&E
    !&E ** History :
    !&E       !  2013-04  (M. Huret) 
    !&E
    !&E---------------------------------------------------------------------
   !! * Modules used
   USE comtraj,     ONLY : hty
   USE trajectools, ONLY : interpvit
   !USE comvars2d,  ONLY : hy_g
   !USE comsiggen,  ONLY : hcy_sig

   !! * Arguments
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )  :: var
   REAL(KIND=rsh), INTENT( in )                                   :: xpos,ypos
   REAL(KIND=rsh), INTENT( in )                                   :: minprof, maxprof
   REAL(KIND=rsh)                                                 :: ibm_profvint

   !! * Local declarations
   REAL(KIND=rsh), DIMENSION(2,2)            :: val4
   INTEGER                                   :: igg,idd,jhh,jbb
   INTEGER                                   :: hlb,hlt,hrb,hrt  
   REAL(KIND=rsh)                            :: xe,px,py


   !!----------------------------------------------------------------------
   !! * Executable part

   !---- reperage de la maille ou se situe la particule pour les V
   ! location of ssh : xpos=REAL(i), ypos=REAL(j)
   igg=int(xpos)
   idd=igg+1
   jbb=nint(ypos)
   jhh=jbb+1

   idd=min(max(idd,imin),imax)
   igg=min(max(igg,imin),imax)
   jbb=min(max(jbb,jmin),jmax)
   jhh=min(max(jhh,jmin),jmax)

   ! calcul des positions px et py dans la maille V
   px=xpos-real(igg)
   py=ypos-real(jbb)-0.5_rsh


   val4(:,:)=0.0_rsh
   ! Processing prof mean on 4 neighbooring locations
   IF (hty(igg,jbb)>0.0_rsh) THEN
      xe=hty(igg,jbb)-hy_g(igg,jbb)
      val4(1,1)=ibm_profmean(var(:,igg,jbb),minprof,maxprof,hy_g(igg,jbb),xe,hcy_sig(igg,jbb))
   ENDIF

   IF (hty(idd,jbb)>0.0_rsh) THEN
      xe=hty(idd,jbb)-hy_g(idd,jbb)
      val4(2,1)=ibm_profmean(var(:,idd,jbb),minprof,maxprof,hy_g(idd,jbb),xe,hcy_sig(idd,jbb))
   ENDIF

   IF (hty(idd,jhh)>0.0_rsh) THEN
      xe=hty(idd,jhh)-hy_g(idd,jhh)
      val4(2,2)=ibm_profmean(var(:,idd,jhh),minprof,maxprof,hy_g(idd,jhh),xe,hcy_sig(idd,jhh))
   ENDIF

   IF (hty(igg,jhh)>0.0_rsh) THEN
      xe=hty(igg,jhh)-hy_g(igg,jhh)
      val4(1,2)=ibm_profmean(var(:,igg,jhh),minprof,maxprof,hy_g(igg,jhh),xe,hcy_sig(igg,jhh))
   ENDIF

   ! Interpolation
   ibm_profvint=interpvit(px,py,val4(1,1),val4(2,1),val4(1,2),val4(2,2), &
                          hty(igg,jbb),hty(idd,jbb),hty(igg,jhh),hty(idd,jhh))


 END FUNCTION



 !!======================================================================
 FUNCTION ibm_tidal_mig(xe,depth1,depth2,w,p,alpha)

    !&E---------------------------------------------------------------------
    !&E
    !&E ** Purpose : Tidal migration between depth1 (flood) and depth2 (ebb) 
    !&E ** 
    !&E ** Description    : 
    !&E ** Called by      : 
    !&E ** External calls :
    !&E ** Reference      : 
    !&E
    !&E ** History        : 
    !&E                ! 2010-11  (M. Huret) (A ameliorer)
    !&E                ! 2011-11  (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E---------------------------------------------------------------------

    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(liminm2:limaxp2,ljminm2:ljmaxp2),INTENT( in )  :: xe
    REAL(KIND=rsh), INTENT( in )                                             :: depth1,depth2,alpha,w
    REAL(KIND=rsh)                                                           :: ibm_tidal_mig
    ! INTEGER, INTENT( in )                                                   :: n,m
    TYPE(type_particle)                                                      :: p

    !! * Local declarations
    INTEGER                 :: jj,mm,aaaa,hh,minu,sec
    REAL(kind=rsh)          :: depth
    
    REAL(KIND=rsh),  ALLOCATABLE, DIMENSION(:,:)  :: xep
    REAL(KIND=rsh),  ALLOCATABLE, DIMENSION(:,:)  :: dxe
    REAL(KIND=rsh),  ALLOCATABLE, DIMENSION(:,:)  :: dxep
    REAL(KIND=rlg),  ALLOCATABLE, DIMENSION(:,:)  :: deb_descend
    REAL(KIND=rlg),  ALLOCATABLE, DIMENSION(:,:)  :: deb_monte  

    !!----------------------------------------------------------------------
    !! * Executable part

    dxe(n,m) = xe(NINT(p%xpos),NINT(p%ypos))-xep(n,m) 
    
    ! Si chgmt d elevation decroissant apres au moins 6h
    IF(dxep(n,m) > 0.0 .AND. dxe(n,m) <= 0.0 .AND. time >= deb_monte(n,m) + 6.0*3600.0) THEN
       deb_monte(n,m) = time + 4.0*3600.0 + 15.0*60.0 ! chgt de courant dans 4heures 15
    ENDIF

    ! Si chgmt d elevation croissant apres au moins 6h
    IF(dxep(n,m) <= 0.0 .AND. dxe(n,m) > 0.0 .AND. time >= deb_descend(n,m) + 6.0*3600.0) THEN
       deb_descend(n,m) = time + 4.0*3600.0 + 15.0*60.0 ! chgt de courant dans 4heures 15
    ENDIF

    depth=depth1 ! en surface qd maree monte
    IF(time < deb_monte(n,m) .AND. time >= deb_descend(n,m)) THEN 
       depth=depth2 ! au fond qd maree descend
    ENDIF
  
    dxep(n,m) = dxe(n,m)
    xep(n,m)=xe(NINT(p%xpos),NINT(p%ypos))
    
    ibm_tidal_mig= w*tanh(alpha*(p%zpos-depth)) 
  
 END FUNCTION ibm_tidal_mig



 ! #####################################################################
 ! #####       Reprendre les variables hm,ig et id de MARS et      #####
 ! #####  les adapter a CROCO pour faire fonctionner detect_pycno  #####
 ! #####################################################################

 !!======================================================================
 SUBROUTINE ibm_detect_pycno(xe,bz,sal,temp,limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E ** Purpose : Detect pycnocl for part behaviour
    !&E 
    !&E ** Description    :
    !&E ** Called by      : 
    !&E ** External calls :
    !&E ** Reference      :
    !&E
    !&E ** History : 
    !&E       ! 2008     (M. Sourisseau) non terminée mais pas encore utilisée
    !&E       ! 2011-11  (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E---------------------------------------------------------------------
   !! * Modules used
   !USE comvars2d, ONLY : hm,ig,id

   !! * Arguments
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),INTENT( in )            :: xe
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )      :: sal,temp
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )      :: bz
   INTEGER,                                                           :: limin,limax,ljmin,ljmax  

   !! * Local declarations
   INTEGER                                          :: k_strat_max,i,j,k
   REAL(kind=rsh)                                   :: d,roref,dzu
   REAL(kind=rsh),dimension(GLOBAL_2D_ARRAY,kmax-1) :: gradv_bz
   INTEGER,dimension(GLOBAL_2D_ARRAY)               :: k_stratif

   !!----------------------------------------------------------------------
   !! * Executable part
   !
   ! calcul du gradient vertical de densite
   !-----------------------
   !       sig au centre des cellules et sigw position verticale entre des limites (calcul des w)

   !cval Make sure you use the same value of roref than in tsstateeqblumberg.F90
   roref=1040.0_rsh

   DO  j=ljmin,ljmax
      !do  i=ig(j),id(j)
      DO i=MAX(limin,ig(j)+1),MIN(limax,id(j)-1)


       depth_strat(i,j)=0.0_rsh
       d=h(i,j)+xe(i,j)

       IF (d .gt. hm) THEN
        k_strat_max=1
        gradv_bz(1,i,j)=0.0_rsh            !grad initial ; grad > 0
         DO k=2,kmax

            dzu=REAL(hz(h(i,j),xe(i,j),k,i,j),rsh)   !dsigu=sigw(k)-sigw(k-1) => dzu>0
            gradv_bz(k,i,j)=((roref-bz(k-1,i,j)/g*roref) - (roref-bz(k,i,j)/g*roref)) &
                              /dzu !grad between k-1 and k ; grad > 0 en m-1


            IF (gradv_bz(k,i,j) .gt. gradv_bz(k_strat_max,i,j)) THEN
               k_strat_max=k           ! numero de la couche ou le gradv de flottab. (bz) est max.sur sa frontiere inf.
            END IF
         ENDDO !endo on k

        IF ((abs(gradv_bz(k_strat_max,i,j)) .gt. 0.02_rsh) ) THEN 
                                 ! .and. (d*(1.0_rsh-(sigw(k_strat_max-1))) .gt. 60.0_rsh)
                                 ! 0.02 m-1 uniquement pour les panaches / stratif estivale non visible
                                 ! 0.000005 m-1 stratif obs partout de mars à sept

            k_stratif(i,j)  = k_strat_max-1  !indic of sigma levels with the highest upper gradient
            depth_strat(i,j)= d*(1.0_rsh-(sc_w(k_stratif(i,j))))   !depth in m of the highest buoanc. gradient
        ELSE
            depth_strat(i,j)=-1.0_rsh
        END IF


       END IF
      ENDDO
   ENDDO
 
 END SUBROUTINE ibm_detect_pycno



 ! ##################################################################### 
 ! #####       Reprendre les variables hm,ig et id de MARS et      #####
 ! #####   les adapter a CROCO pour faire fonctionner opt_depth    #####
 ! #####################################################################

 !!======================================================================
 SUBROUTINE ibm_opt_depth(xe,limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E
    !&E ** Purpose : Estimate the optimal depth for part. behaviour.
    !&E ** Description :
    !&E ** Called by : w_behav_part
    !&E ** External calls :
    !&E ** Reference :
    !&E ** History : 
    !&E       ! 2008-11  (M. Sourisseau) pour l instant opt_depth fixe
    !&E       ! 2011-11  (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E---------------------------------------------------------------------
   !! * Modules used
   !USE comvars2d, ONLY : hm,ig,id

   !! * Arguments
   REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY),INTENT( in )   :: xe
   INTEGER,                                                 :: limin,limax,ljmin,ljmax

   !! * Local declarations
   INTEGER           :: k_strat_max,i,j,k
   REAL(kind=rsh)    :: d
 
   !!----------------------------------------------------------------------
   !! * Executable part
   !
   ! estimation of the optimal vertical depth
   !-----------------------
   !
   DO  j=ljmin,ljmax
      DO  i=MAX(limin,ig(j)+1),MIN(limax,id(j)-1)

       d= h(i,j)+xe(i,j)
       IF (d .gt. hm) THEN
        !  depth_strat(i,j)=15.0_rsh;

         IF (depth_strat(i,j) .gt. 0.0_rsh) THEN
          opt_depth_part(i,j)=depth_strat(i,j)
         ELSE
          opt_depth_part(i,j)=0.0_rsh ! dpeth in m
         END IF

       END IF
      ENDDO
   ENDDO

 END SUBROUTINE ibm_opt_depth
#endif /* key_ibm__unused */


#endif  /* DEB_IBM */


END MODULE