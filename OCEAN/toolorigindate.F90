# include "cppdefs.h"
  SUBROUTINE tool_origindate(netcdfid,varid,date_in_sec)
#ifdef MPI
    use scalars, ONLY : mynode
#endif
    USE croco_namelist, ONLY: calendar_type, start_date
    IMPLICIT NONE
# include "netcdf.inc"

    INTEGER       :: netcdfid,varid
    CHARACTER*180 :: units
    CHARACTER*40  :: file_calendar
    CHARACTER*19  :: date_str
    INTEGER       :: lenstr, luni, indst, ierr, ierr2
    REAL*8        :: tool_datosec,date_in_sec

    ierr=nf_get_att_text(netcdfid, varid, 'units', units)
    if (ierr .eq. nf_noerr) then
      luni = lenstr(units)
      if (units(1:6).eq.'second')then
        ! Case units as the format 'seconds since YYYY-MM-DD hh:mm:ss'
        indst=15 ! position to start reading date
      elseif (units(1:3).eq.'day') then
        ! Case units as the format 'days since YYYY-MM-DD hh:mm:ss'
        indst=12 ! position to start reading date
      else
        MPI_master_only write(*,'(/1x,2A/6x,A/10x,A)') &
               'TOOL_ORIGINDATE ERROR: ', &
               'unknown units for time variable' ,&
               'Time variable should follow Netcdf CF format: ',&
               '''seconds(days) since YYYY-MM-DD hh:mm:ss'''
        STOP
      endif
    else
      ! No units attribute: assume the file's time axis is relative to
      ! model start_date (old-style climatology files without USE_CALENDAR).
      ! Use start_date as the origin so the result is consistent with model time.
      MPI_master_only write(*,'(/1x,A/6x,2A/)') &
        'TOOL_ORIGINDATE WARNING: no units attribute in forcing file.', &
        'Assuming time axis is relative to start_date: ', TRIM(start_date)
      date_in_sec = tool_datosec(start_date)
      RETURN
    endif

    if (luni<indst) then
      MPI_master_only write(*,'(/1x,A/6x,A/10x,A/6x,2A/)') &
              'TOOL_ORIGINDATE WARNING: no date found in time var units.', &
              'Time variable should follow Netcdf CF format: ', &
              '''seconds(days) since YYYY-MM-DD hh:mm:ss''', &
              'Assuming time axis is relative to start_date: ', TRIM(start_date)
      date_in_sec = tool_datosec(start_date)
      RETURN
    elseif (luni-indst.eq.3) then
      MPI_master_only write(*,'(/1x,4A/1x)') &
           'TOOL_ORIGINDATE: ',&
           'Only origin year is specified, suppose it is ', units(indst:luni),&
           '-01-01 00:00:00'
      date_str=units(indst:luni)//'-01-01 00:00:00'
    elseif (luni-indst.eq.6) then
      MPI_master_only write(*,'(/1x,4A/1x)') &
           'TOOL_ORIGINDATE: ',&
           'Only origin year and month are specified, suppose it is ', units(indst:luni),&
           '/01 00:00:00'
      date_str=units(indst:luni)//'-01 00:00:00'
    elseif (luni-indst.eq.9) then
      MPI_master_only write(*,'(/1x,4A/1x)') &
           'TOOL_ORIGINDATE: ',&
           'Only origin year,month,day are specified, suppose it is ', units(indst:luni),&
           ' 00:00:00'
      date_str=units(indst:luni)//' 00:00:00'
    elseif (luni-indst.eq.12) then
      MPI_master_only write(*,'(/1x,4A/1x)') &
           'TOOL_ORIGINDATE: ', &
           'Only origin year,month,day,hour are specified, suppose it is ', units(indst:luni),&
           ':00:00'
      date_str=units(indst:luni)//':00:00'
    elseif (luni-indst.eq.15) then
       MPI_master_only write(*,'(/1x,4A/1x)') &
           'TOOL_ORIGINDATE: ', &
           'Only origin year,month,day,hour,minute are specified, suppose it is ', units(indst:luni),&
           ':00'
      date_str=units(indst:luni)//':00'
    else
      date_str=units(indst:luni)
    endif

    ! Check calendar attribute of the file against model calendar_type
    file_calendar = ' '
    ierr2=nf_get_att_text(netcdfid, varid, 'calendar', file_calendar)
    if (ierr2 .eq. nf_noerr) then
      ! File has a calendar attribute — check it matches the model
      if (TRIM(file_calendar) /= TRIM(calendar_type) .AND. &
          .NOT. (TRIM(file_calendar)=='gregorian' .AND. TRIM(calendar_type)=='gregorian') .AND. &
          .NOT. (TRIM(file_calendar)=='standard'  .AND. TRIM(calendar_type)=='gregorian')) then
        MPI_master_only write(*,'(/1x,A/6x,A,A/6x,A,A/)') &
          'TOOL_ORIGINDATE WARNING: calendar mismatch between file and model:', &
          '  file calendar   = ', TRIM(file_calendar), &
          '  model calendar_type = ', TRIM(calendar_type)
      endif
    else
      ! No calendar attribute in file: assume it matches the model, warn if non-gregorian
      if (TRIM(calendar_type) /= 'gregorian') then
        MPI_master_only write(*,'(/1x,2A,A/)') &
          'TOOL_ORIGINDATE WARNING: no calendar attribute in file,', &
          ' assuming model calendar_type: ', TRIM(calendar_type)
      endif
    endif

    date_in_sec=tool_datosec(date_str)

  END SUBROUTINE TOOL_ORIGINDATE
