# include "cppdefs.h"
MODULE tools_calendar
!!---------------------------------------------------------------------
!!                 ***  MODULE tools_calendar  ***
!!
!! ** Purpose : calendar conversion utilities for CROCO.
!!              Supports gregorian/standard (default), 360_day,
!!              365_day/no_leap.
!!
!! ** History : creation of a module from standalone functions and
!!              subroutines, 2026, Croco Team
!!
!!---------------------------------------------------------------------
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
             tool_decompdate, tool_origindate, init_tools_calendar

   ! Module-level copies of namelist calendar settings.
   ! Set once via init_tools_calendar (called from init_calendar).
   ! Stored here so this module does not depend on AGRIF-managed croco_namelist.
   ! AGRIF conv can't handle character(len=N)
!$AGRIF_DO_NOT_TREAT
   character(len=20) :: calendar_type = 'gregorian'
   character(len=19) :: start_date    = '                   '
!$AGRIF_END_DO_NOT_TREAT

   ! (netCDF file id, variable id, warning kind) triples for which a
   ! tool_origindate warning has already been issued, so that each kind of
   ! warning (partial date, calendar mismatch, missing calendar attribute)
   ! is only printed once per time variable (repeat reads/cycling of the
   ! same variable are silenced, but distinct time variables sharing a
   ! file, and distinct warning kinds for the same variable, each still
   ! warn once). Grown dynamically (doubling) in first_origindate_warning:
   ! no fixed cap on the number of distinct entries tracked.
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: warned_ncid
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: warned_varid
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: warned_kind
   INTEGER, SAVE :: n_warned_partial = 0

   ! Warning-kind discriminators for first_origindate_warning, so that e.g.
   ! a partial-date warning for a variable does not suppress a later,
   ! genuinely first-time calendar warning for that same variable.
   INTEGER, PARAMETER :: warn_partial_date = 1
   INTEGER, PARAMETER :: warn_calendar     = 2

CONTAINS

   !====================================================================
   LOGICAL FUNCTION first_origindate_warning(ncdfid, vid, kind)
   !! ** Purpose : return .TRUE. only the first time it is called for a
   !!              given (netCDF file id, variable id, warning kind)
   !!              triple, so that a per-variable, per-warning-kind
   !!              tool_origindate warning can be printed once instead of
   !!              once per call (e.g. once per data-cycling re-read).

      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncdfid, vid, kind
      INTEGER :: i
      INTEGER, DIMENSION(:), ALLOCATABLE :: grown

      DO i = 1, n_warned_partial
         IF (warned_ncid(i) == ncdfid .AND. &
             warned_varid(i) == vid .AND. &
             warned_kind(i) == kind) THEN
            first_origindate_warning = .FALSE.
            RETURN
         END IF
      END DO

      first_origindate_warning = .TRUE.

      IF (.NOT. ALLOCATED(warned_ncid)) THEN
         ALLOCATE (warned_ncid(16))
         ALLOCATE (warned_varid(16))
         ALLOCATE (warned_kind(16))
      ELSE IF (n_warned_partial == SIZE(warned_ncid)) THEN
         ALLOCATE (grown(2*n_warned_partial))
         grown(1:n_warned_partial) = warned_ncid(1:n_warned_partial)
         CALL MOVE_ALLOC(grown, warned_ncid)

         ALLOCATE (grown(2*n_warned_partial))
         grown(1:n_warned_partial) = warned_varid(1:n_warned_partial)
         CALL MOVE_ALLOC(grown, warned_varid)

         ALLOCATE (grown(2*n_warned_partial))
         grown(1:n_warned_partial) = warned_kind(1:n_warned_partial)
         CALL MOVE_ALLOC(grown, warned_kind)
      END IF

      n_warned_partial = n_warned_partial + 1
      warned_ncid(n_warned_partial) = ncdfid
      warned_varid(n_warned_partial) = vid
      warned_kind(n_warned_partial) = kind

   END FUNCTION first_origindate_warning

   !====================================================================
   SUBROUTINE warn_if_partial_date(ncdfid, varid, date_str, context)
   !! ** Purpose : called by tool_origindate: if date_str is a partial
   !!              'yyyy[-mm[-dd[ hh[:mm[:ss]]]]]'
   !!              string, warn once (per ncdfid/varid, via the same
   !!              first_origindate_warning tracker/warn_partial_date
   !!              kind) which fields are being defaulted. Does not
   !!              rewrite date_str: tool_datosec/tool_decompdate already
   !!              default missing fields on their own. A length that is
   !!              not one of 4/7/10/13/16/19 is left for tool_decompdate
   !!              to reject as malformed.

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncdfid, varid
      CHARACTER(len=*), INTENT(in) :: date_str, context

      CHARACTER(len=16) :: filled
      INTEGER :: dlen

      dlen = LEN_TRIM(date_str)

      SELECT CASE (dlen)
      CASE (4)
         filled = '-01-01 00:00:00'
      CASE (7)
         filled = '-01 00:00:00'
      CASE (10)
         filled = ' 00:00:00'
      CASE (13)
         filled = ':00:00'
      CASE (16)
         filled = ':00'
      CASE DEFAULT
         RETURN
      END SELECT

      IF (first_origindate_warning(ncdfid, varid, warn_partial_date)) THEN
         MPI_master_only write (*, '(/1x,4A/1x)') &
            'TOOL_CALENDAR WARNING: ', TRIM(context), &
            ': partial date '//TRIM(date_str), &
            ', assuming '//TRIM(date_str)//TRIM(filled)
      END IF

   END SUBROUTINE warn_if_partial_date

   !====================================================================
   SUBROUTINE tool_fatal_stop()
   !! ** Purpose : abort the run on a fatal calendar error. Under MPI,
   !!              aborts all ranks (mpi_abort) instead of relying on a
   !!              plain STOP, which only terminates the calling process
   !!              and can leave the other ranks hanging.

      IMPLICIT NONE
#if defined MPI
      include 'mpif.h'

      call mpi_abort (MPI_COMM_WORLD, 1)
#endif
      STOP

   END SUBROUTINE tool_fatal_stop

   !====================================================================
   CHARACTER(len=19) FUNCTION tool_sectodat(time)
   !! ** Purpose : convert a time in seconds since the year-1900 epoch
   !!              to a date string "yyyy-mm-dd hh:mm:ss".

      IMPLICIT NONE
      REAL(kind=rlg), INTENT(in) :: time

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
   !!              for a date string. Accepts partial strings:
   !!              "yyyy", "yyyy-mm", "yyyy-mm-dd", "yyyy-mm-dd hh:mm",
   !!              "yyyy-mm-dd hh:mm:ss". Missing fields default to 1 (day/month)
   !!              or 0 (hour/minute/second).

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: date

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
   !! ** Purpose : decompose a date string into integer components.
   !!              Accepts partial strings; missing fields default to
   !!              1 (month/day) or 0 (hour/minute/second).
   !!              Supported formats (dlen = LEN_TRIM):
   !!                dlen >= 4  : "yyyy"
   !!                dlen >= 7  : "yyyy-mm"
   !!                dlen >= 10 : "yyyy-mm-dd"
   !!                dlen >= 16 : "yyyy-mm-dd hh:mm"
   !!                dlen >= 19 : "yyyy-mm-dd hh:mm:ss"

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)  :: date
      INTEGER, INTENT(out) :: dd, mm, yyyy, hh, minu, sec

      INTEGER :: dlen

      ! Defaults for optional fields
      mm = 1; dd = 1; hh = 0; minu = 0; sec = 0

      dlen = LEN_TRIM(date)

      IF (dlen < 4) THEN
         MPI_master_only write (*, '(/1x,3A/)') &
            'TOOL_DECOMPDATE ERROR: date string too short: "', TRIM(date), '"'
         call tool_fatal_stop()
      END IF

      READ (date(1:4),   '(i4)') yyyy
      IF (dlen >= 7)  READ (date(6:7),   '(i2)') mm
      IF (dlen >= 10) READ (date(9:10),  '(i2)') dd
      IF (dlen >= 13) READ (date(12:13), '(i2)') hh
      IF (dlen >= 16) READ (date(15:16), '(i2)') minu
      IF (dlen >= 19) READ (date(18:19), '(i2)') sec

      IF (mm < 1 .OR. mm > 12) THEN
         MPI_master_only write (*, '(/1x,3A/)') &
            'TOOL_DECOMPDATE ERROR: invalid month in date: "', TRIM(date), '"'
         call tool_fatal_stop()
      END IF

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

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(in)  :: netcdfid, varid
      REAL(kind=rlg), INTENT(out) :: date_in_sec

      CHARACTER*180 :: units
      CHARACTER*40  :: file_calendar
      CHARACTER*19  :: date_str
      CHARACTER*40  :: varname
      CHARACTER*250 :: ncfile
      CHARACTER*120 :: fix_hint
      CHARACTER*300 :: loc_info
      INTEGER       :: lenstr, luni, indst, ierr, ierr2, ierr3, ierr4
      INTEGER :: pathlen

      ! netcdfid/varid alone are meaningless to a user reading an error or
      ! warning message below: look up the actual file path and variable
      ! name so every message can name what it is complaining about.
      ! Defaults are kept if either lookup fails, so messages still make
      ! sense (just less specific) rather than referencing an undefined
      ! string.
      varname = 'time'
      ierr3 = nf90_inquire_variable(netcdfid, varid, name=varname)

      ncfile = '(unknown file)'
      ierr4 = nf90_inq_path(netcdfid, path=ncfile, pathlen=pathlen)

      loc_info = 'netCDF file: '//TRIM(ncfile)//', variable: '//TRIM(varname)

      ! Suggested fix for an existing but malformed 'units' attribute
      fix_hint = 'ncatted -a units,'//TRIM(varname)// &
                 ',o,c,''seconds since YYYY-MM-DD hh:mm:ss'' '//TRIM(ncfile)

      ierr = nf90_get_att(netcdfid, varid, 'units', units)
      if (ierr .ne. nf90_noerr) then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'no ''units'' attribute found for the time variable.', &
            TRIM(loc_info), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'You can add it with, e.g.: ', &
            'ncatted -a units,'//TRIM(varname)// &
            ',c,c,''seconds since YYYY-MM-DD hh:mm:ss'' '//TRIM(ncfile)
         call tool_fatal_stop()
      end if

      luni = lenstr(units)
      if (index(units(1:luni), 'since') == 0) then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'no ''since'' keyword found in time variable units attribute.', &
            TRIM(loc_info), &
            'units attribute value: '//TRIM(units), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end if

      if (units(1:6) .eq. 'second') then
         indst = 15
      elseif (units(1:3) .eq. 'day') then
         indst = 12
      else
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'unknown units for time variable.', &
            TRIM(loc_info), &
            'units attribute value: '//TRIM(units), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end if

      if (luni < indst) then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'no date found after ''since'' in time variable units attribute.', &
            TRIM(loc_info), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      else
         date_str = units(indst:luni)
         call warn_if_partial_date(netcdfid, varid, date_str, TRIM(loc_info))
      end if

      file_calendar = ' '
      ierr2 = nf90_get_att(netcdfid, varid, 'calendar', file_calendar)
      if (ierr2 .eq. nf90_noerr) then
         if (TRIM(file_calendar) /= TRIM(calendar_type) .AND. &
             .NOT. (TRIM(file_calendar) == 'gregorian' .AND. TRIM(calendar_type) == 'gregorian') .AND. &
             .NOT. (TRIM(file_calendar) == 'standard' .AND. TRIM(calendar_type) == 'gregorian')) then
            if (first_origindate_warning(netcdfid, varid, warn_calendar)) then
               MPI_master_only write (*, '(/1x,A/6x,A/6x,A,A/6x,A,A/)') &
                  'TOOL_ORIGINDATE WARNING: calendar mismatch between file and model:', &
                  TRIM(loc_info), &
                  '  file calendar   = ', TRIM(file_calendar), &
                  '  model calendar_type = ', TRIM(calendar_type)
            end if
         end if
      else
         if (TRIM(calendar_type) /= 'gregorian') then
            if (first_origindate_warning(netcdfid, varid, warn_calendar)) then
               MPI_master_only write (*, '(/1x,2A,A/6x,A/)') &
                  'TOOL_ORIGINDATE WARNING: no calendar attribute in file,', &
                  ' assuming model calendar_type: ', TRIM(calendar_type), &
                  TRIM(loc_info)
            end if
         end if
      end if

      ! Reject non-standard date strings (must start with 4-digit year YYYY-)
      if (date_str(1:1) < '0' .or. date_str(1:1) > '9' .or. &
          date_str(2:2) < '0' .or. date_str(2:2) > '9' .or. &
          date_str(3:3) < '0' .or. date_str(3:3) > '9' .or. &
          date_str(4:4) < '0' .or. date_str(4:4) > '9') then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'non-standard date format in time variable units.', &
            TRIM(loc_info), &
            'origin date read: '//TRIM(date_str), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end if

      date_in_sec = tool_datosec(date_str)

   END SUBROUTINE tool_origindate

   !====================================================================
   SUBROUTINE init_tools_calendar(cal_type, s_date)
   !! ** Purpose : initialise module-level calendar settings from croco_namelist.
   !!              Must be called once during setup (from init_calendar).
   !!              calendar_type is mandatory and must be one of the
   !!              calendars this module implements: an unrecognized value
   !!              is a fatal error rather than a silent fallback to
   !!              gregorian.

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: cal_type, s_date

      SELECT CASE (TRIM(cal_type))
      CASE ('gregorian', 'standard', '360_day', '365_day', 'no_leap')
         calendar_type = cal_type
      CASE DEFAULT
         MPI_master_only write (*, '(/1x,2A/6x,A/)') &
            'INIT_TOOLS_CALENDAR ERROR: ', &
            'unknown ''calendar_type'' in croco_calendar namelist: '''//TRIM(cal_type)//'''.', &
            'Supported values: ''gregorian'', ''standard'', ''360_day'', ''365_day'', ''no_leap''.'
         call tool_fatal_stop()
      END SELECT

      start_date = s_date

   END SUBROUTINE init_tools_calendar

END MODULE tools_calendar
