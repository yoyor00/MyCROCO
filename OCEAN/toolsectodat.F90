  FUNCTION tool_sectodat(time)

   !&E---------------------------------------------------------------------
   !&E                 ***  FUNCTION tool_sectodat  ***
   !&E
   !&E ** Purpose : convert a time expressed in seconds (relative to calendar
   !&E              epoch ~year 1900) into a date string "yyyy-mm-dd hh:mm:ss".
   !&E              Supports: gregorian (default), 360_day, 365_day/no_leap.
   !&E
   !&E---------------------------------------------------------------------
   USE croco_namelist, ONLY: calendar_type
   IMPLICIT NONE
   INTEGER,PARAMETER     :: rlg=8
   ! Gregorian epoch offset
   REAL(kind=rlg),PARAMETER :: tref     = 59958230400_rlg
   ! Non-Gregorian epoch offsets (must match tooldatosec.F90)
   REAL(kind=rlg),PARAMETER :: tref_360 = 1900_rlg * 360.0_rlg * 86400.0_rlg
   REAL(kind=rlg),PARAMETER :: tref_365 = 1900_rlg * 365.0_rlg * 86400.0_rlg

   !! * Declarations function
   CHARACTER(len=19)          :: tool_sectodat

   !! * Arguments
   real :: time

   !! * Local declarations
   LOGICAL               :: bissext
   CHARACTER(len=19)     :: date
   INTEGER               :: annee,mois,jour,heure,minute,seconde,tot_jours
   INTEGER               :: nb_siecles,nb_annees,nb_4siecles,nb_4annees
   INTEGER,DIMENSION(12) :: jours_avt_mois
   INTEGER,DIMENSION(365):: mois_from_jour
   REAL(kind=rlg),PARAMETER :: secs_in_minute=60.0_rlg,secs_in_heure=3600.0_rlg,secs_in_jour=86400.0_rlg
   REAL(kind=rlg),PARAMETER :: secs_in_annee=secs_in_jour*365.0_rlg
   REAL(kind=rlg),PARAMETER :: secs_in_siecle=secs_in_jour*(76.0_rlg*365.0_rlg+24.0_rlg*366.0_rlg)
   REAL(kind=rlg),PARAMETER :: secs_in_4annees=secs_in_jour*(3.0_rlg*365.0_rlg+366.0_rlg)
   REAL(kind=rlg),PARAMETER :: secs_in_4siecles=4.0_rlg*secs_in_siecle+secs_in_jour
   REAL(kind=rlg)           :: total_secs

   DATA mois_from_jour /31*01,28*02,31*03,30*04,31*05,30*06,31*07,31*08,30*09,31*10,30*11,31*12/
   DATA jours_avt_mois /0,31,59,90,120,151,181,212,243,273,304,334/

   !!----------------------------------------------------------------------
   !! * Executable part

   IF (TRIM(calendar_type) == '360_day') THEN
     ! 360-day calendar: 12 months of 30 days, no leap years
     total_secs = REAL(time,rlg) + tref_360

     annee      = INT(total_secs / (360.0_rlg * secs_in_jour))
     total_secs = total_secs - REAL(annee,rlg) * 360.0_rlg * secs_in_jour

     mois       = INT(total_secs / (30.0_rlg * secs_in_jour)) + 1
     total_secs = total_secs - REAL(mois-1,rlg) * 30.0_rlg * secs_in_jour

     jour       = INT(total_secs / secs_in_jour) + 1
     total_secs = MOD(total_secs, secs_in_jour)

     heure      = INT(total_secs / secs_in_heure)
     total_secs = MOD(total_secs, secs_in_heure)
     minute     = INT(total_secs / secs_in_minute)
     seconde    = INT(MOD(total_secs, secs_in_minute))

   ELSE IF (TRIM(calendar_type) == '365_day' .OR. &
            TRIM(calendar_type) == 'no_leap') THEN
     ! 365-day calendar: standard month lengths, no leap years
     total_secs = REAL(time,rlg) + tref_365

     annee      = INT(total_secs / secs_in_annee)
     total_secs = total_secs - REAL(annee,rlg) * secs_in_annee

     tot_jours  = INT(total_secs / secs_in_jour)
     total_secs = MOD(total_secs, secs_in_jour)

     ! Find month (no leap year, use standard month lengths)
     mois = mois_from_jour(MIN(tot_jours+1, 365))
     jour = tot_jours - jours_avt_mois(mois) + 1

     heure      = INT(total_secs / secs_in_heure)
     total_secs = MOD(total_secs, secs_in_heure)
     minute     = INT(total_secs / secs_in_minute)
     seconde    = INT(MOD(total_secs, secs_in_minute))

   ELSE
     ! Default: proleptic Gregorian calendar (original implementation)
     total_secs = time+tref
     annee      = 0
     bissext = (MOD(3,2) == 0)

     IF(total_secs  >=  (secs_in_annee+secs_in_jour)) THEN
       total_secs = total_secs - (secs_in_annee+secs_in_jour)
       nb_4siecles= INT(total_secs/secs_in_4siecles)
       total_secs = MOD(total_secs,secs_in_4siecles)
       nb_siecles = INT(total_secs/secs_in_siecle)

       IF(nb_siecles==4 .AND. total_secs >= secs_in_4siecles-secs_in_jour) nb_siecles = 3
       total_secs = total_secs - nb_siecles*secs_in_siecle
       annee      = 400*nb_4siecles + 100*nb_siecles
       nb_4annees = INT(total_secs/secs_in_4annees)
       total_secs = MOD(total_secs,secs_in_4annees)
       nb_annees  = INT(total_secs/secs_in_annee)

       IF(nb_annees==4 .AND. total_secs >= secs_in_4annees-secs_in_jour) nb_annees=3
       total_secs = total_secs - nb_annees*secs_in_annee
       annee      = annee + 4*nb_4annees + nb_annees + 1

     ENDIF

     bissext =(MOD(annee,400) == 0 .OR. (MOD(annee,4) == 0 .AND. MOD(annee,100) /= 0))

     tot_jours  = INT(total_secs/secs_in_jour)
     total_secs = MOD(total_secs,secs_in_jour)

     IF (bissext .AND. tot_jours >= 59) THEN
       mois = mois_from_jour(tot_jours)
     ELSE
       mois = mois_from_jour(tot_jours+1)
     ENDIF

     IF (bissext .AND. mois>2) THEN
       jour = tot_jours - jours_avt_mois(mois)
     ELSE
       jour = tot_jours - jours_avt_mois(mois) + 1
     ENDIF

     heure      = INT(total_secs/secs_in_heure)
     total_secs = MOD(total_secs,secs_in_heure)
     minute     = INT(total_secs/secs_in_minute)
     total_secs = MOD(total_secs,secs_in_minute)
     seconde = total_secs

   END IF

   WRITE (date,800) annee,mois,jour,heure,minute,seconde
   tool_sectodat = date

800 FORMAT(i4.4,'-',i2.2,'-',i2.2,' ',2(i2.2,':'),i2.2)

  END FUNCTION tool_sectodat
