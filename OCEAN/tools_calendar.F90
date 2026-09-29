# include "cppdefs.h"
MODULE tools_calendar
!!---------------------------------------------------------------------
!!                 ***  MODULE tools_calendar  ***
!!
!! ** Purpose : calendar conversion utilities for CROCO.
!!              Supports the CF calendars (case-insensitive):
!!                - gregorian / standard / proleptic_gregorian (default,
!!                  implemented as proleptic Gregorian: identical to
!!                  CF 'standard' for all dates from 1582-10-15 on; a
!!                  'standard'/'gregorian' start date or file origin date
!!                  before 1582-10-15 is rejected, see
!!                  before_gregorian_reform)
!!                - 360_day
!!                - noleap / 365_day (legacy non-CF alias: no_leap)
!!                - all_leap / 366_day
!!                - julian
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
   REAL(kind=rlg), PARAMETER :: tref_366 = 1900_rlg*366.0_rlg*86400.0_rlg
   ! Julian: 1900 years of 365 days + 475 leap years (0, 4, ..., 1896)
   REAL(kind=rlg), PARAMETER :: tref_jul = (1900.0_rlg*365.0_rlg + 475.0_rlg)*86400.0_rlg

   ! Time unit constants
   REAL(kind=rlg), PARAMETER :: secs_in_minute = 60.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_hour = 3600.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_day = 86400.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_year = secs_in_day*365.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_century = secs_in_day*36524.0_rlg
   REAL(kind=rlg), PARAMETER :: secs_in_4years = secs_in_day*(3.0_rlg*365.0_rlg + 366.0_rlg)
   REAL(kind=rlg), PARAMETER :: secs_in_4centuries = 4.0_rlg*secs_in_century + secs_in_day

   ! Days elapsed before each month in a non-leap / leap year
   INTEGER, DIMENSION(12), PARAMETER :: days_before_month = &
                                        (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)
   INTEGER, DIMENSION(12), PARAMETER :: days_before_month_leap = &
                                        (/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/)

   ! Canonical calendar identifiers: every CF name/alias maps to one of
   ! these (see cf_calendar_id), so that equivalent names never trigger
   ! a mismatch and branches test an integer instead of strings.
   INTEGER, PARAMETER :: cal_unknown   = 0
   INTEGER, PARAMETER :: cal_gregorian = 1
   INTEGER, PARAMETER :: cal_360day    = 2
   INTEGER, PARAMETER :: cal_noleap    = 3
   INTEGER, PARAMETER :: cal_allleap   = 4
   INTEGER, PARAMETER :: cal_julian    = 5

   PUBLIC :: tool_sectodat, tool_datosec, tool_datetosec, &
             tool_decompdate, tool_origindate, init_tools_calendar

   ! Module-level copies of namelist calendar settings.
   ! Set once via init_tools_calendar (called from init_calendar).
   ! Stored here so this module does not depend on AGRIF-managed croco_namelist.
   ! AGRIF conv can't handle character(len=N)
   ! model_cal is kept here too: the calendar is common to all grids.
!$AGRIF_DO_NOT_TREAT
   character(len=20) :: calendar_type = 'gregorian'
   character(len=26) :: start_date    = '                          '
   integer           :: model_cal     = cal_gregorian
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
   CHARACTER(len=40) FUNCTION normalize_calendar(str)
   !! ** Purpose : return a calendar name lower-cased, left-adjusted and
   !!              cut at the first NUL character (some tools write
   !!              NUL-terminated netCDF string attributes, and TRIM
   !!              does not remove NULs).

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: str

      CHARACTER(len=40) :: s
      INTEGER :: i, c

      s = str
      DO i = 1, LEN(s)
         c = IACHAR(s(i:i))
         IF (c == 0) THEN
            s(i:) = ' '
            EXIT
         ELSE IF (c >= IACHAR('A') .AND. c <= IACHAR('Z')) THEN
            s(i:i) = ACHAR(c + 32)
         END IF
      END DO
      normalize_calendar = ADJUSTL(s)

   END FUNCTION normalize_calendar

   !====================================================================
   INTEGER FUNCTION cf_calendar_id(str)
   !! ** Purpose : map a CF calendar name (or alias) to its canonical
   !!              identifier; cal_unknown if not supported.

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: str

      SELECT CASE (TRIM(normalize_calendar(str)))
      CASE ('gregorian', 'standard', 'proleptic_gregorian')
         cf_calendar_id = cal_gregorian
      CASE ('360_day')
         cf_calendar_id = cal_360day
      CASE ('noleap', '365_day', 'no_leap')
         cf_calendar_id = cal_noleap
      CASE ('all_leap', '366_day')
         cf_calendar_id = cal_allleap
      CASE ('julian')
         cf_calendar_id = cal_julian
      CASE DEFAULT
         cf_calendar_id = cal_unknown
      END SELECT

   END FUNCTION cf_calendar_id

   !====================================================================
   LOGICAL FUNCTION is_mixed_gregorian(str)
   !! ** Purpose : .TRUE. if str names the CF mixed Julian/Gregorian
   !!              calendar ('standard' or its deprecated alias
   !!              'gregorian'), which this module implements as
   !!              proleptic Gregorian, i.e. exactly only from 1582-10-15.

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: str

      SELECT CASE (TRIM(normalize_calendar(str)))
      CASE ('standard', 'gregorian')
         is_mixed_gregorian = .TRUE.
      CASE DEFAULT
         is_mixed_gregorian = .FALSE.
      END SELECT

   END FUNCTION is_mixed_gregorian

   !====================================================================
   LOGICAL FUNCTION before_gregorian_reform(date_str)
   !! ** Purpose : .TRUE. if date_str is before 1582-10-15, the first
   !!              day of the Gregorian calendar. Before that date, CF
   !!              'standard' is Julian and differs from the proleptic
   !!              Gregorian arithmetic used here (by 10 days in 1582,
   !!              2 days at year 1).

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: date_str

      INTEGER :: yyyy, mm, dd, hh, minu, sec

      CALL tool_decompdate(date_str, dd, mm, yyyy, hh, minu, sec)
      before_gregorian_reform = yyyy < 1582 .OR. &
         (yyyy == 1582 .AND. (mm < 10 .OR. (mm == 10 .AND. dd < 15)))

   END FUNCTION before_gregorian_reform

   !====================================================================
   SUBROUTINE doy_to_month_day(doy, is_leap, month, day)
   !! ** Purpose : convert a 0-based day of year into month and day of
   !!              month, for a 365- or 366-day year.

      IMPLICIT NONE
      INTEGER, INTENT(in)  :: doy
      LOGICAL, INTENT(in)  :: is_leap
      INTEGER, INTENT(out) :: month, day

      INTEGER, DIMENSION(12) :: dbm

      IF (is_leap) THEN
         dbm = days_before_month_leap
      ELSE
         dbm = days_before_month
      END IF

      month = 12
      DO WHILE (dbm(month) > doy)
         month = month - 1
      END DO
      day = doy - dbm(month) + 1

   END SUBROUTINE doy_to_month_day

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

      SELECT CASE (model_cal)

      CASE (cal_360day)
         ! 360-day calendar: 12 months of 30 days, no leap years
         total_secs = REAL(time, rlg) + tref_360

         year = INT(total_secs/(360.0_rlg*secs_in_day))
         total_secs = total_secs - REAL(year, rlg)*360.0_rlg*secs_in_day

         month = INT(total_secs/(30.0_rlg*secs_in_day)) + 1
         total_secs = total_secs - REAL(month - 1, rlg)*30.0_rlg*secs_in_day

         day = INT(total_secs/secs_in_day) + 1
         total_secs = MOD(total_secs, secs_in_day)

      CASE (cal_noleap)
         ! 365-day calendar: standard month lengths, no leap years
         total_secs = REAL(time, rlg) + tref_365

         year = INT(total_secs/secs_in_year)
         total_secs = total_secs - REAL(year, rlg)*secs_in_year

         tot_days = INT(total_secs/secs_in_day)
         total_secs = MOD(total_secs, secs_in_day)

         month = month_from_day(MIN(tot_days + 1, 365))
         day = tot_days - days_before_month(month) + 1

      CASE (cal_allleap)
         ! 366-day calendar: every year is a leap year
         total_secs = REAL(time, rlg) + tref_366

         year = INT(total_secs/(366.0_rlg*secs_in_day))
         total_secs = total_secs - REAL(year, rlg)*366.0_rlg*secs_in_day

         tot_days = INT(total_secs/secs_in_day)
         total_secs = MOD(total_secs, secs_in_day)

         CALL doy_to_month_day(MIN(tot_days, 365), .TRUE., month, day)

      CASE (cal_julian)
         ! Julian calendar: leap year every 4 years, no century rule.
         ! Year 0 is a leap year, so each 4-year cycle is 366+3*365 days.
         total_secs = REAL(time, rlg) + tref_jul

         tot_days = INT(total_secs/secs_in_day)
         total_secs = total_secs - REAL(tot_days, rlg)*secs_in_day

         nb_4years = tot_days/1461
         tot_days = MOD(tot_days, 1461)
         IF (tot_days < 366) THEN
            nb_years = 0
         ELSE
            nb_years = 1 + (tot_days - 366)/365
            tot_days = tot_days - 366 - 365*(nb_years - 1)
         END IF
         year = 4*nb_4years + nb_years

         CALL doy_to_month_day(tot_days, nb_years == 0, month, day)

      CASE DEFAULT
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

      END SELECT

      ! total_secs now holds the seconds elapsed within the day
      hour = INT(total_secs/secs_in_hour)
      total_secs = MOD(total_secs, secs_in_hour)
      minute = INT(total_secs/secs_in_minute)
      second = INT(MOD(total_secs, secs_in_minute))

      WRITE (date, 800) year, month, day, hour, minute, second
      tool_sectodat = date

800   FORMAT(i4.4, '-', i2.2, '-', i2.2, ' ', 2(i2.2, ':'), i2.2)

   END FUNCTION tool_sectodat

   !====================================================================
   REAL(kind=rlg) FUNCTION tool_datosec(date)
   !! ** Purpose : return seconds elapsed since the year-1900 epoch
   !!              for a date string. Accepts partial strings:
   !!              "yyyy", "yyyy-mm", "yyyy-mm-dd", "yyyy-mm-dd hh:mm",
   !!              "yyyy-mm-dd hh:mm:ss", "yyyy-mm-dd hh:mm:ss.ffff...".
   !!              Missing fields default to 1 (day/month) or 0
   !!              (hour/minute/second); a trailing ".ffff..." after the
   !!              whole seconds is an optional fractional-second part.

      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: date

      INTEGER        :: year, month, day, hour, minute, second, tot_days
      REAL(kind=rlg) :: total_secs, secs_of_day, frac_sec

      CALL tool_decompdate(date, day, month, year, hour, minute, second, frac_sec)

      secs_of_day = REAL(hour, rlg)*secs_in_hour &
                    + REAL(minute, rlg)*secs_in_minute &
                    + REAL(second, rlg)

      SELECT CASE (model_cal)

      CASE (cal_360day)
         ! 360-day calendar: 12 months of 30 days, no leap years
         total_secs = REAL(year, rlg)*360.0_rlg*secs_in_day
         total_secs = total_secs + REAL(month - 1, rlg)*30.0_rlg*secs_in_day
         total_secs = total_secs + REAL(day - 1, rlg)*secs_in_day
         tool_datosec = total_secs + secs_of_day - tref_360

      CASE (cal_noleap)
         ! 365-day calendar: standard month lengths, no leap years ever
         total_secs = REAL(year, rlg)*secs_in_year
         total_secs = total_secs + REAL(days_before_month(month), rlg)*secs_in_day
         total_secs = total_secs + REAL(day - 1, rlg)*secs_in_day
         tool_datosec = total_secs + secs_of_day - tref_365

      CASE (cal_allleap)
         ! 366-day calendar: every year is a leap year
         total_secs = REAL(year, rlg)*366.0_rlg*secs_in_day
         total_secs = total_secs + REAL(days_before_month_leap(month), rlg)*secs_in_day
         total_secs = total_secs + REAL(day - 1, rlg)*secs_in_day
         tool_datosec = total_secs + secs_of_day - tref_366

      CASE (cal_julian)
         ! Julian calendar: (year+3)/4 leap years (0, 4, ...) before 'year'
         tot_days = 365*year + (year + 3)/4 + days_before_month(month) + day - 1
         IF (month > 2 .AND. MOD(year, 4) == 0) tot_days = tot_days + 1
         tool_datosec = REAL(tot_days, rlg)*secs_in_day + secs_of_day - tref_jul

      CASE DEFAULT
         ! Default: proleptic Gregorian calendar
         tool_datosec = gregorian_to_sec(day, month, year, hour, minute, second)

      END SELECT

      tool_datosec = tool_datosec + frac_sec

   END FUNCTION tool_datosec

   !====================================================================
   SUBROUTINE tool_decompdate(date, dd, mm, yyyy, hh, minu, sec, frac_sec)
   !! ** Purpose : decompose a date string into integer components.
   !!              Accepts partial strings; missing fields default to
   !!              1 (month/day) or 0 (hour/minute/second).
   !!              Supported formats (dlen = LEN_TRIM):
   !!                dlen >= 4  : "yyyy"
   !!                dlen >= 7  : "yyyy-mm"
   !!                dlen >= 10 : "yyyy-mm-dd"
   !!                dlen >= 16 : "yyyy-mm-dd hh:mm"
   !!                dlen >= 19 : "yyyy-mm-dd hh:mm:ss"
   !!                dlen >  19 : "yyyy-mm-dd hh:mm:ss.ffff..." -- a '.' right
   !!                             after the whole seconds, followed by any
   !!                             number of fractional-second digits, returned
   !!                             via the optional frac_sec (0 if absent/not
   !!                             requested). Needed for calendar_type
   !!                             sub-second precision (e.g. NBQ/acoustic test
   !!                             cases with dt well under 1 second).

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)  :: date
      INTEGER, INTENT(out) :: dd, mm, yyyy, hh, minu, sec
      REAL(kind=rlg), INTENT(out), OPTIONAL :: frac_sec

      INTEGER :: dlen
      REAL(kind=rlg) :: frac_local

      ! Defaults for optional fields
      mm = 1; dd = 1; hh = 0; minu = 0; sec = 0
      frac_local = 0.0_rlg

      dlen = LEN_TRIM(date)

      ! Accept a trailing ISO-8601 UTC 'Z' designator (e.g. a NetCDF units
      ! attribute such as "seconds since 2000-01-01T00:00:00Z"): CROCO has
      ! no other timezone handling, so it is equivalent to no suffix at all.
      IF (dlen > 0) THEN
         IF (date(dlen:dlen) == 'Z' .OR. date(dlen:dlen) == 'z') dlen = dlen - 1
      END IF

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
      IF (dlen > 19) THEN
         IF (date(20:20) /= '.') THEN
            MPI_master_only write (*, '(/1x,3A/)') &
               'TOOL_DECOMPDATE ERROR: invalid trailing characters after seconds in date: "', TRIM(date), '"'
            call tool_fatal_stop()
         END IF
         READ (date(20:dlen), *) frac_local
      END IF
      IF (PRESENT(frac_sec)) frac_sec = frac_local

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
   !!              supported calendar types are handled automatically.
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
      ! Called by tool_datosec (Gregorian branch).

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
   SUBROUTINE tool_origindate(netcdfid, varid, date_in_sec, secs_per_unit)
   !! ** Purpose : read origin date from a NetCDF time variable 'units'
   !!              attribute ('<unit> since <date>') and return it as
   !!              seconds since epoch.
   !!              Supported <unit> (case-insensitive): seconds, minutes,
   !!              hours, days, with their UDUNITS singular/short forms.
   !!              Optional secs_per_unit returns the number of seconds
   !!              in one <unit> (1, 60, 3600 or 86400), so that callers
   !!              convert the time values themselves consistently:
   !!                 time_in_sec = date_in_sec + value*secs_per_unit

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      USE netcdf
      IMPLICIT NONE

      INTEGER, INTENT(in)  :: netcdfid, varid
      REAL(kind=rlg), INTENT(out) :: date_in_sec
      REAL(kind=rlg), INTENT(out), OPTIONAL :: secs_per_unit

      CHARACTER*180 :: units
      CHARACTER*180 :: units_lc
      CHARACTER*20  :: unit_word
      REAL(kind=rlg):: unit_factor
      INTEGER       :: isince, i, ic
      CHARACTER*40  :: file_calendar
      CHARACTER*26  :: date_str
      CHARACTER*40  :: varname
      CHARACTER*250 :: ncfile
      CHARACTER*120 :: fix_hint
      CHARACTER*300 :: loc_info
      INTEGER       :: lenstr, luni, ierr, ierr2, ierr3, ierr4
      INTEGER :: pathlen
      INTEGER :: file_cal
      LOGICAL :: file_is_mixed
      CHARACTER*60 :: cal_desc

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
            '''seconds|minutes|hours|days since YYYY-MM-DD hh:mm:ss''', &
            'You can add it with, e.g.: ', &
            'ncatted -a units,'//TRIM(varname)// &
            ',c,c,''seconds since YYYY-MM-DD hh:mm:ss'' '//TRIM(ncfile)
         call tool_fatal_stop()
      end if

      luni = lenstr(units)

      ! Lower-cased copy for parsing ('Hours Since ...' is accepted);
      ! the original is kept for messages.
      units_lc = units
      do i = 1, luni
         ic = IACHAR(units_lc(i:i))
         if (ic >= IACHAR('A') .AND. ic <= IACHAR('Z')) units_lc(i:i) = ACHAR(ic + 32)
      end do

      isince = index(units_lc(1:luni)//' ', ' since ')
      if (isince == 0) then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'no ''since'' keyword found in time variable units attribute.', &
            TRIM(loc_info), &
            'units attribute value: '//TRIM(units), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds|minutes|hours|days since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end if

      ! Unit word: everything before ' since ', e.g. 'hours'
      unit_word = ADJUSTL(units_lc(1:isince))
      select case (TRIM(unit_word))
      case ('seconds', 'second', 'secs', 'sec', 's')
         unit_factor = 1.0_rlg
      case ('minutes', 'minute', 'mins', 'min')
         unit_factor = secs_in_minute
      case ('hours', 'hour', 'hrs', 'hr', 'h')
         unit_factor = secs_in_hour
      case ('days', 'day', 'd')
         unit_factor = secs_in_day
      case default
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'unknown units for time variable.', &
            TRIM(loc_info), &
            'units attribute value: '//TRIM(units), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds|minutes|hours|days since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end select

      ! Date: everything after ' since ' (7 characters), left-adjusted
      if (LEN_TRIM(units(isince+7:luni)) == 0 .OR. isince + 7 > luni) then
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/10x,A/6x,A/10x,A/)') &
            'TOOL_ORIGINDATE ERROR: ', &
            'no date found after ''since'' in time variable units attribute.', &
            TRIM(loc_info), &
            'Time variable must follow Netcdf CF format: ', &
            '''seconds|minutes|hours|days since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      else
         date_str = ADJUSTL(units(isince+7:luni))
         call warn_if_partial_date(netcdfid, varid, date_str, TRIM(loc_info))
      end if

      ! Calendar check: names are compared through their canonical id, so
      ! CF-equivalent names (e.g. 'standard' / 'proleptic_gregorian' /
      ! 'gregorian', or 'noleap' / '365_day') never trigger a mismatch.
      ! A missing attribute means 'standard' (CF default).
      file_calendar = ' '
      ierr2 = nf90_get_att(netcdfid, varid, 'calendar', file_calendar)
      if (ierr2 .eq. nf90_noerr) then
         file_calendar = normalize_calendar(file_calendar)
         file_is_mixed = is_mixed_gregorian(file_calendar)
         cal_desc = 'file calendar = '//TRIM(file_calendar)
      else
         ! CF default when the attribute is absent is 'standard'. Only
         ! relevant if the model calendar is Gregorian-like: otherwise the
         ! model calendar is assumed for the file (warning below).
         file_is_mixed = (model_cal == cal_gregorian)
         cal_desc = 'no calendar attribute (CF default: standard)'
      end if

      if (ierr2 .eq. nf90_noerr) then
         file_cal = cf_calendar_id(file_calendar)
         if (file_cal == cal_unknown) then
            if (first_origindate_warning(netcdfid, varid, warn_calendar)) then
               MPI_master_only write (*, '(/1x,A/6x,A/6x,A,A/6x,A,A/)') &
                  'TOOL_ORIGINDATE WARNING: unrecognized calendar in file:', &
                  TRIM(loc_info), &
                  '  file calendar   = ', TRIM(file_calendar), &
                  '  assuming model calendar_type = ', TRIM(calendar_type)
            end if
         else if (file_cal /= model_cal) then
            if (first_origindate_warning(netcdfid, varid, warn_calendar)) then
               MPI_master_only write (*, '(/1x,A/6x,A/6x,A,A/6x,A,A/)') &
                  'TOOL_ORIGINDATE WARNING: calendar mismatch between file and model:', &
                  TRIM(loc_info), &
                  '  file calendar   = ', TRIM(file_calendar), &
                  '  model calendar_type = ', TRIM(calendar_type)
            end if
         end if
      else
         if (model_cal /= cal_gregorian) then
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
            '''seconds|minutes|hours|days since YYYY-MM-DD hh:mm:ss''', &
            'You can fix it with, e.g.: ', &
            TRIM(fix_hint)
         call tool_fatal_stop()
      end if

      ! A 'standard' origin date before 1582-10-15 is a Julian date: the
      ! proleptic Gregorian arithmetic of tool_datosec would shift it
      ! (and every time value of the file) by up to ~10 days.
      if (file_is_mixed) then
         if (before_gregorian_reform(date_str)) then
            MPI_master_only write (*, '(/1x,A/(6x,A))') &
               'TOOL_ORIGINDATE ERROR: origin date before 1582-10-15 in a mixed Julian/Gregorian calendar.', &
               TRIM(loc_info), &
               'origin date read: '//TRIM(date_str)//', '//TRIM(cal_desc), &
               'CROCO implements ''standard''/''gregorian'' as proleptic Gregorian,', &
               'which differs from the CF mixed calendar before 1582-10-15.', &
               'If the file times are really proleptic Gregorian, declare it with, e.g.:', &
               '    ncatted -a calendar,'//TRIM(varname)//',o,c,proleptic_gregorian '//TRIM(ncfile), &
               'Otherwise, rebase the time axis on an origin from 1582-10-15 on, e.g.:', &
               '    cdo setreftime,1900-01-01,00:00:00 '//TRIM(ncfile)//' out.nc'
            call tool_fatal_stop()
         end if
      end if

      date_in_sec = tool_datosec(date_str)
      IF (PRESENT(secs_per_unit)) secs_per_unit = unit_factor

   END SUBROUTINE tool_origindate

   !====================================================================
   SUBROUTINE init_tools_calendar(cal_type, s_date)
   !! ** Purpose : initialise module-level calendar settings from croco_namelist.
   !!              Must be called once during setup (from init_calendar).
   !!              calendar_type is mandatory and must be one of the
   !!              calendars this module implements (CF names, case-
   !!              insensitive): an unrecognized value is a fatal error
   !!              rather than a silent fallback to gregorian.

#if defined MPI
      USE scalars, ONLY: mynode
#endif
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: cal_type, s_date

      CHARACTER(len=40) :: cal_norm

      model_cal = cf_calendar_id(cal_type)

      IF (model_cal == cal_unknown) THEN
         MPI_master_only write (*, '(/1x,2A/6x,A/6x,A/)') &
            'INIT_TOOLS_CALENDAR ERROR: ', &
            'unknown ''calendar_type'' in croco_calendar namelist: '''//TRIM(cal_type)//'''.', &
            'Supported values: ''gregorian'', ''standard'', ''proleptic_gregorian'', ''360_day'',', &
            '''noleap'', ''365_day'', ''all_leap'', ''366_day'', ''julian''.'
         call tool_fatal_stop()
      END IF

      cal_norm = normalize_calendar(cal_type)
      calendar_type = cal_norm(1:20)
      start_date = s_date

      ! Model time axis starts at start_date and only moves forward: a
      ! start date from 1582-10-15 on is enough for 'standard'/'gregorian'
      ! to be computed exactly by the proleptic Gregorian arithmetic.
      IF (is_mixed_gregorian(calendar_type) .AND. LEN_TRIM(s_date) >= 4) THEN
         IF (before_gregorian_reform(s_date)) THEN
            MPI_master_only write (*, '(/1x,A/(6x,A))') &
               'INIT_TOOLS_CALENDAR ERROR: start_date before 1582-10-15 with calendar_type '''// &
               TRIM(calendar_type)//'''.', &
               'start_date = '//TRIM(s_date), &
               'CROCO implements ''standard''/''gregorian'' as proleptic Gregorian,', &
               'which differs from the CF mixed Julian/Gregorian calendar before 1582-10-15.', &
               'Use calendar_type = ''proleptic_gregorian'' for such dates.'
            call tool_fatal_stop()
         END IF
      END IF

   END SUBROUTINE init_tools_calendar

END MODULE tools_calendar