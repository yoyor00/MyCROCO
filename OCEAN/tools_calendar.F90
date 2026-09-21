# include "cppdefs.h"
MODULE tools_calendar
!!---------------------------------------------------------------------
!!                 ***  MODULE tools_calendar  ***
!!
!! ** Purpose : calendar conversion utilities for CROCO.
!!              Supports gregorian (default), 360_day, 365_day/no_leap.
!!
!!---------------------------------------------------------------------

   USE croco_namelist, ONLY: calendar_type, start_date
   IMPLICIT NONE
   PRIVATE

   ! Kind parameter for double-precision calendar arithmetic
   INTEGER, PARAMETER :: rlg = 8

   ! Epoch offsets: seconds from year 0 to year 1900 for each calendar.
   ! Centering on 1900 keeps CROCO's single-precision 'time' (seconds since
   ! start_date) small enough (~1e9 for multi-decade runs) to avoid precision loss.
   REAL(kind=rlg), PARAMETER :: tref = 59958230400_rlg
   REAL(kind=rlg), PARAMETER :: tref_360 = 1900_rlg*360.0_rlg*86400.0_rlg
   REAL(kind=rlg), PARAMETER :: tref_365 = 1900_rlg*365.0_rlg*86400.0_rlg

   ! Time unit constants
   REAL(kind=rlg), PARAMETER :: secs_in_minute = 60.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_hour = 3600.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_day = 86400.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_year = secs_in_day*365.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_century = secs_in_day*36524.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_4years = secs_in_day*(3.0_rlg*365.0_rlg + 366.0_rlg)
   REAL(kind=rlg), PARAMETER :: secs_in_4centuries = 4.0_rlg*secs_in_century + secs_in_day

   ! Days elapsed before each month in a non-leap year
   INTEGER, DIMENSION(12), PARAMETER :: days_before_month = &
                                        (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

   PUBLIC :: tool_sectodat, tool_datosec, tool_datetosec, &
             tool_decompdate, tool_origindate

CONTAINS

   !====================================================================
   CHARACTER(len=19) FUNCTION tool_sectodat(time)
   !! ** Purpose : convert a time in seconds since the year-1900 epoch
   !!              to a date string "yyyy-mm-dd hh:mm:ss".

      IMPLICIT NONE
      REAL, INTENT(in) :: time

      LOGICAL               :: is_leap
      CHARACTER(len=19)     :: date
      INTEGER               :: year, month, day, hour, minute, second, tot_days
      INTEGER               :: nb_centuries, nb_years, nb_4centuries, nb_4years
      INTEGER, DIMENSION(365):: month_from_day
      REAL(kind=rlg)        :: total_secs

      DATA month_from_day/31*01, 28*02, 31*03, 30*04, 31*05, 30*06, 31*07, 31*08, 30*09, 31*10, 30*11, 31*12/

      IF (TRIM(calendar_type) == '360_day') THEN
         ! 360-day calendar: 12 months of 30 days, no leap years
         total_secs = REAL(time, rlg) + tref_360

         year = INT(total_secs/(360.0_rlg*secs_in_day))
         total_secs = total_secs - REAL(year, rlg)*360.0_rlg*secs_in_day

         month = INT(total_secs/(30.0_rlg*secs_in_day)) + 1
         total_secs = total_secs - REAL(month - 1, rlg)*30.0_rlg*secs_in_day

         day = INT(total_secs/secs_in_day) + 1
         total_secs = MOD(total_secs, secs_in_day)

         hour = INT(total_secs/secs_in_hour)
         total_secs = MOD(total_secs, secs_in_hour)
         minute = INT(total_secs/secs_in_minute)
         second = INT(MOD(total_secs, secs_in_minute))

      ELSE IF (TRIM(calendar_type) == '365_day' .OR. &
               TRIM(calendar_type) == 'no_leap') THEN
         ! 365-day calendar: standard month lengths, no leap years
         total_secs = REAL(time, rlg) + tref_365

         year = INT(total_secs/secs_in_year)
         total_secs = total_secs - REAL(year, rlg)*secs_in_year

         tot_days = INT(total_secs/secs_in_day)
         total_secs = MOD(total_secs, secs_in_day)

         month = month_from_day(MIN(tot_days + 1, 365))
         day = tot_days - days_before_month(month) + 1

         hour = INT(total_secs/secs_in_hour)
         total_secs = MOD(total_secs, secs_in_hour)
         minute = INT(total_secs/secs_in_minute)
         second = INT(MOD(total_secs, secs_in_minute))

      ELSE
         ! Default: proleptic Gregorian calendar
         total_secs = time + tref
         year = 0
         is_leap = (MOD(3, 2) == 0)

         IF (total_secs >= (secs_in_year + secs_in_day)) THEN
            total_secs = total_secs - (secs_in_year + secs_in_day)
            nb_4centuries = INT(total_secs/secs_in_4centuries)
            total_secs = MOD(total_secs, secs_in_4centuries)
            nb_centuries = INT(total_secs/secs_in_century)

            IF (nb_centuries == 4 .AND. total_secs >= secs_in_4centuries - secs_in_day) nb_centuries = 3
            total_secs = total_secs - nb_centuries*secs_in_century
            year = 400*nb_4centuries + 100*nb_centuries
            nb_4years = INT(total_secs/secs_in_4years)
            total_secs = MOD(total_secs, secs_in_4years)
            nb_years = INT(total_secs/secs_in_year)

            IF (nb_years == 4 .AND. total_secs >= secs_in_4years - secs_in_day) nb_years = 3
            total_secs = total_secs - nb_years*secs_in_year
            year = year + 4*nb_4years + nb_years + 1

         END IF

         is_leap = (MOD(year, 400) == 0 .OR. (MOD(year, 4) == 0 .AND. MOD(year, 100) /= 0))

         tot_days = INT(total_secs/secs_in_day)
         total_secs = MOD(total_secs, secs_in_day)

         IF (is_leap .AND. tot_days >= 59) THEN
            month = month_from_day(tot_days)
         ELSE
            month = month_from_day(tot_days + 1)
         END IF

         IF (is_leap .AND. month > 2) THEN
            day = tot_days - days_before_month(month)
         ELSE
            day = tot_days - days_before_month(month) + 1
         END IF

         hour = INT(total_secs/secs_in_hour)
         total_secs = MOD(total_secs, secs_in_hour)
         minute = INT(total_secs/secs_in_minute)
         total_secs = MOD(total_secs, secs_in_minute)
         second = total_secs

      END IF

      WRITE (date, 800) year, month, day, hour, minute, second
      tool_sectodat = date

800   FORMAT(i4.4, '-', i2.2, '-', i2.2, ' ', 2(i2.2, ':'), i2.2)

   END FUNCTION tool_sectodat

   !====================================================================
   REAL(kind=rlg) FUNCTION tool_datosec(date)
   !! ** Purpose : return seconds elapsed since the year-1900 epoch
   !!              for a date string "yyyy-mm-dd hh:mm:ss".

      IMPLICIT NONE
      CHARACTER(len=19), INTENT(in) :: date

      INTEGER        :: year, month, day, hour, minute, second
      REAL(kind=rlg) :: total_secs

      CALL tool_decompdate(date, day, month, year, hour, minute, second)

      IF (TRIM(calendar_type) == '360_day') THEN
         ! 360-day calendar: 12 months of 30 days, no leap years
         total_secs = REAL(year, rlg)*360.0_rlg*secs_in_day
         total_secs = total_secs + REAL(month - 1, rlg)*30.0_rlg*secs_in_day
         total_secs = total_secs + REAL(day - 1, rlg)*secs_in_day
         total_secs = total_secs + REAL(hour, rlg)*secs_in_hour
         total_secs = total_secs + REAL(minute, rlg)*secs_in_minute
         total_secs = total_secs + REAL(second, rlg)
         tool_datosec = total_secs - tref_360

      ELSE IF (TRIM(calendar_type) == '365_day' .OR. &
               TRIM(calendar_type) == 'no_leap') THEN
         ! 365-day calendar: standard month lengths, no leap years ever
         total_secs = REAL(year, rlg)*secs_in_year
         total_secs = total_secs + REAL(days_before_month(month), rlg)*secs_in_day
         total_secs = total_secs + REAL(day - 1, rlg)*secs_in_day
         total_secs = total_secs + REAL(hour, rlg)*secs_in_hour
         total_secs = total_secs + REAL(minute, rlg)*secs_in_minute
         total_secs = total_secs + REAL(second, rlg)
         tool_datosec = total_secs - tref_365

      ELSE
         ! Default: proleptic Gregorian calendar
         tool_datosec = gregorian_to_sec(day, month, year, hour, minute, second)

      END IF

   END FUNCTION tool_datosec

   !====================================================================
   SUBROUTINE tool_decompdate(date, dd, mm, yyyy, hh, minu, sec)
   !! ** Purpose : decompose a 19-character date string "yyyy-mm-dd hh:mm:ss"
   !!              into integer components.

      IMPLICIT NONE
      CHARACTER(len=19), INTENT(in)   :: date
      INTEGER, INTENT(out)  :: dd, mm, yyyy, hh, minu, sec

      LOGICAL                        :: cont
      CHARACTER(len=2)               :: month_str
      CHARACTER(len=2), DIMENSION(12) :: tab_months = (/'01', '02', '03', '04', '05', '06', &
                                                        '07', '08', '09', '10', '11', '12'/)
      INTEGER                        :: i

      READ (date, 800) yyyy, month_str, dd, hh, minu, sec
      cont = .TRUE.
      i = 1

      DO WHILE (cont .AND. i <= 12)
         IF (tab_months(i) == month_str) THEN
            cont = .FALSE.
         ELSE
            i = i + 1
         END IF
      END DO

      IF (i <= 12) THEN
         mm = i
      ELSE
         PRINT *, 'tool_decompdate error: cannot determine month'
         STOP
      END IF

800   FORMAT(i4, 1x, a2, 1x, i2, 1x, 2(i2, 1x), i2)

   END SUBROUTINE tool_decompdate

   !====================================================================
   SUBROUTINE tool_datetosec(day, month, year, hour, minute, second, datetosec)
   !! ** Purpose : return seconds elapsed since epoch for a date given as
   !!              integer components. Delegates to tool_datosec so that all
   !!              calendar types (gregorian, 360_day, 365_day/no_leap) are
   !!              supported automatically.
   !! ** Called by : init_oa subroutine in module_oa_interface

      IMPLICIT NONE
      INTEGER, INTENT(in)  :: year, month, day, hour, minute, second
      REAL(kind=rlg), INTENT(out) :: datetosec

      CHARACTER(len=19) :: date_str
      WRITE (date_str, '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') &
         year, month, day, hour, minute, second
      datetosec = tool_datosec(date_str)

   END SUBROUTINE tool_datetosec

   !====================================================================
   REAL(kind=rlg) FUNCTION gregorian_to_sec(day, month, year, hour, minute, second)
      ! Private helper: proleptic Gregorian date components → seconds since tref.
      ! Called by tool_datosec (Gregorian branch) and tool_datetosec.

      IMPLICIT NONE
      INTEGER, INTENT(in) :: day, month, year, hour, minute, second

      REAL(kind=rlg) :: total_secs

      total_secs = secs_in_century*INT(year/100)
      total_secs = total_secs + secs_in_day*INT(DBLE(year)/400.0_rlg + 0.9975_rlg)
      total_secs = total_secs + secs_in_year*MOD(year, 100)
      total_secs = total_secs + secs_in_day*INT((MOD(year, 100) - 1)/4)
      total_secs = total_secs + days_before_month(month)*secs_in_day
      IF (month > 2) THEN
         IF (MOD(year, 400) == 0) THEN
            total_secs = total_secs + secs_in_day
         ELSE
            IF ((MOD(year, 4) == 0) .AND. (MOD(year, 100) /= 0)) &
               total_secs = total_secs + secs_in_day
         END IF
      END IF
      total_secs = total_secs + secs_in_day*(day - 1)
      total_secs = total_secs + secs_in_hour*(hour)
      total_secs = total_secs + secs_in_minute*(minute)
      total_secs = total_secs + second
      gregorian_to_sec = total_secs - tref

   END FUNCTION gregorian_to_sec

   !====================================================================
   SUBROUTINE tool_origindate(netcdfid, varid, date_in_sec)
   !! ** Purpose : read origin date from a NetCDF time variable 'units'
   !!              attribute and return it as seconds since epoch.

#ifdef MPI
      USE scalars, ONLY: mynode
#endif
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(in)  :: netcdfid, varid
      REAL(kind=rlg), INTENT(out) :: date_in_sec

      CHARACTER*180 :: units
      CHARACTER*40  :: file_calendar
      CHARACTER*19  :: date_str
      INTEGER       :: lenstr, luni, indst, ierr, ierr2

      ierr = nf90_get_att(netcdfid, varid, 'units', units)
      if (ierr .eq. nf90_noerr) then
         luni = lenstr(units)
         if (units(1:6) .eq. 'second') then
            indst = 15
         elseif (units(1:3) .eq. 'day') then
            indst = 12
         else
            MPI_master_only write (*, '(/1x,2A/6x,A/10x,A)') &
               'TOOL_ORIGINDATE ERROR: ', &
               'unknown units for time variable', &
               'Time variable should follow Netcdf CF format: ', &
               '''seconds(days) since YYYY-MM-DD hh:mm:ss'''
            STOP
         end if
      else
         MPI_master_only write (*, '(/1x,A/6x,2A/)') &
            'TOOL_ORIGINDATE WARNING: no units attribute in forcing file.', &
            'Assuming time axis is relative to start_date: ', TRIM(start_date)
         date_in_sec = tool_datosec(start_date)
         RETURN
      end if

      if (luni < indst) then
         MPI_master_only write (*, '(/1x,A/6x,A/10x,A/6x,2A/)') &
            'TOOL_ORIGINDATE WARNING: no date found in time var units.', &
            'Time variable should follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'Assuming time axis is relative to start_date: ', TRIM(start_date)
         date_in_sec = tool_datosec(start_date)
         RETURN
      elseif (luni - indst .eq. 3) then
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_ORIGINDATE: ', &
            'Only origin year is specified, suppose it is ', units(indst:luni), &
            '-01-01 00:00:00'
         date_str = units(indst:luni)//'-01-01 00:00:00'
      elseif (luni - indst .eq. 6) then
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_ORIGINDATE: ', &
            'Only origin year and month are specified, suppose it is ', units(indst:luni), &
            '/01 00:00:00'
         date_str = units(indst:luni)//'-01 00:00:00'
      elseif (luni - indst .eq. 9) then
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_ORIGINDATE: ', &
            'Only origin year,month,day are specified, suppose it is ', units(indst:luni), &
            ' 00:00:00'
         date_str = units(indst:luni)//' 00:00:00'
      elseif (luni - indst .eq. 12) then
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_ORIGINDATE: ', &
            'Only origin year,month,day,hour are specified, suppose it is ', units(indst:luni), &
            ':00:00'
         date_str = units(indst:luni)//':00:00'
      elseif (luni - indst .eq. 15) then
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_ORIGINDATE: ', &
            'Only origin year,month,day,hour,minute are specified, suppose it is ', units(indst:luni), &
            ':00'
         date_str = units(indst:luni)//':00'
      else
         date_str = units(indst:luni)
      end if

      file_calendar = ' '
      ierr2 = nf90_get_att(netcdfid, varid, 'calendar', file_calendar)
      if (ierr2 .eq. nf90_noerr) then
         if (TRIM(file_calendar) /= TRIM(calendar_type) .AND. &
             .NOT. (TRIM(file_calendar) == 'gregorian' .AND. TRIM(calendar_type) == 'gregorian') .AND. &
             .NOT. (TRIM(file_calendar) == 'standard' .AND. TRIM(calendar_type) == 'gregorian')) then
            MPI_master_only write (*, '(/1x,A/6x,A,A/6x,A,A/)') &
               'TOOL_ORIGINDATE WARNING: calendar mismatch between file and model:', &
               '  file calendar   = ', TRIM(file_calendar), &
               '  model calendar_type = ', TRIM(calendar_type)
         end if
      else
         if (TRIM(calendar_type) /= 'gregorian') then
            MPI_master_only write (*, '(/1x,2A,A/)') &
               'TOOL_ORIGINDATE WARNING: no calendar attribute in file,', &
               ' assuming model calendar_type: ', TRIM(calendar_type)
         end if
      end if

      date_in_sec = tool_datosec(date_str)

   END SUBROUTINE tool_origindate

END MODULE tools_calendar
