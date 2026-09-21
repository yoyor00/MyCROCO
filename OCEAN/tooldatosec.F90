  FUNCTION tool_datosec(date)

   !&E---------------------------------------------------------------------
   !&E                 ***  FUNCTION tool_datosec  ***
   !&E
   !&E ** Purpose : return the number of seconds elapsed since a calendar
   !&E              epoch (~year 1900) for the current calendar_type.
   !&E              Supports: gregorian (default), 360_day, 365_day/no_leap.
   !&E
   !&E---------------------------------------------------------------------
   USE croco_namelist, ONLY: calendar_type
   IMPLICIT NONE
   INTEGER,PARAMETER     :: rlg=8
   ! Gregorian epoch offset (~year 1900 in Gregorian seconds from year 0)
   REAL(kind=rlg),PARAMETER :: tref     = 59958230400_rlg
   ! Non-Gregorian epoch offsets centred near year 1900 for single-precision accuracy
   REAL(kind=rlg),PARAMETER :: tref_360 = 1900_rlg * 360.0_rlg * 86400.0_rlg
   REAL(kind=rlg),PARAMETER :: tref_365 = 1900_rlg * 365.0_rlg * 86400.0_rlg

   !! * Declaration function
   REAL(kind=rlg)             :: tool_datosec

   !! * Arguments
   CHARACTER(len=19),INTENT(in) :: date

   !! * Local declarations
   INTEGER              :: annee, mois, jour, heure, minute, seconde
   INTEGER,DIMENSION(12),PARAMETER :: jours_avt_mois=(/0,31,59,90,120,151,181,212,243,273,304,334/)
   REAL(kind=rlg),PARAMETER   :: secs_in_minute=60.0_rlg,secs_in_heure=3600.0_rlg
   REAL(kind=rlg),PARAMETER   :: secs_in_jour=86400.0_rlg,secs_in_annee=secs_in_jour*365.0_rlg
   REAL(kind=rlg),PARAMETER   :: secs_in_siecle=secs_in_jour*36524.0_rlg
   REAL(kind=rlg)             :: total_secs

   !!----------------------------------------------------------------------
   !! * Executable part

   CALL tool_decompdate(date,jour,mois,annee,heure,minute,seconde)

   IF (TRIM(calendar_type) == '360_day') THEN
     ! 360-day calendar: 12 months of 30 days, no leap years
     total_secs = REAL(annee,rlg) * 360.0_rlg * secs_in_jour
     total_secs = total_secs + REAL(mois-1,rlg) * 30.0_rlg * secs_in_jour
     total_secs = total_secs + REAL(jour-1,rlg) * secs_in_jour
     total_secs = total_secs + REAL(heure,rlg)  * secs_in_heure
     total_secs = total_secs + REAL(minute,rlg) * secs_in_minute
     total_secs = total_secs + REAL(seconde,rlg)
     tool_datosec = total_secs - tref_360

   ELSE IF (TRIM(calendar_type) == '365_day' .OR. &
            TRIM(calendar_type) == 'no_leap') THEN
     ! 365-day calendar: standard month lengths, no leap years ever
     total_secs = REAL(annee,rlg) * secs_in_annee
     total_secs = total_secs + REAL(jours_avt_mois(mois),rlg) * secs_in_jour
     total_secs = total_secs + REAL(jour-1,rlg) * secs_in_jour
     total_secs = total_secs + REAL(heure,rlg)  * secs_in_heure
     total_secs = total_secs + REAL(minute,rlg) * secs_in_minute
     total_secs = total_secs + REAL(seconde,rlg)
     tool_datosec = total_secs - tref_365

   ELSE
     ! Default: proleptic Gregorian calendar (original implementation)

     ! on ajoute le nombre de secondes pour chaque siecle depuis l origine
     total_secs = secs_in_siecle * INT(annee/100)

     ! on ajoute un jour tous les 400 ans (siecle bissextile)
     total_secs = total_secs + secs_in_jour*INT(DBLE(annee)/400.0_rlg+0.9975_rlg)

     ! on ajoute chaque annee depuis le debut du dernier siecle
     total_secs = total_secs + secs_in_annee*MOD(annee,100)

     ! on ajoute un jour pour chaque annee bissextile depuis le dernier siecle
     total_secs = total_secs + secs_in_jour*INT((MOD(annee,100)-1)/4)

     ! on ajoute chaque mois depuis le debut de l annee
     total_secs = total_secs + jours_avt_mois(mois)*secs_in_jour

     ! on ajoute un jour si on est apres fevrier et que l annee est bissextile
     IF(mois > 2) THEN
       IF(MOD(annee,400) == 0) THEN
         total_secs = total_secs + secs_in_jour
       ELSE
         IF ((MOD(annee,4) == 0) .AND. (MOD(annee,100) /= 0)) &
           total_secs = total_secs + secs_in_jour
       ENDIF
     ENDIF

     ! on ajoute le nombre de jours entres en parametre
     total_secs = total_secs + secs_in_jour*(jour-1)

     ! on ajoute le nombre d heures entrees en parametre
     total_secs = total_secs + secs_in_heure*(heure)

     ! on ajoute le nombre de minutes entrees en parametre
     total_secs = total_secs + secs_in_minute*(minute)

     ! on ajoute le nombre de secondes entrees en parametre
     total_secs = total_secs + seconde
     tool_datosec = total_secs - tref

   END IF

  END FUNCTION tool_datosec
