#!/bin/bash
# Source this file after defining SCRIPTDIR, start date, job duration and NBJOB.
# Calendar helpers julday.sh and caldat.sh are unchanged.

export CALTYPE=greg

# Preserve the existing helpers (and their argument order: month day year).
function valid_date() {
    local jd
    jd=$("${SCRIPTDIR}/julday.sh" "$1" "$2" "$3" "$CALTYPE") || return
    "${SCRIPTDIR}/caldat.sh" "$jd" "$CALTYPE"
}

function makedate() {
    local m d y h
    read -r m d y <<< "$(valid_date "$1" "$2" "$3")"
    h=$4
    printf '%04d%02d%02dh%02d\n' "$y" "$m" "$d" "$h"
}

function sec2hour() {
    local secs=$1
    printf '%02d:%02d:%02d\n' \
        "$((secs / 3600))" "$(((secs / 60) % 60))" "$((secs % 60))"
}

# Convert a date/time to an absolute number of hours in CALTYPE.
function date2hour() {
    local jd
    jd=$("${SCRIPTDIR}/julday.sh" "$1" "$2" "$3" "$CALTYPE") || return
    echo "$((24 * jd + 10#$4))"
}

# Convert absolute hours back to: MONTH DAY YEAR HOUR.
function hour2date() {
    local t=$1 jd hour mdy
    jd=$((t / 24))
    hour=$((t % 24))
    if ((hour < 0)); then
        jd=$((jd - 1))
        hour=$((hour + 24))
    fi
    mdy=$("${SCRIPTDIR}/caldat.sh" "$jd" "$CALTYPE") || return
    printf '%s %d\n' "$mdy" "$hour"
}

# Next job start (= current job's EXCLUSIVE end).
# Arguments: MONTH DAY YEAR HOUR. Output: MONTH DAY YEAR HOUR.
function next_job_date() {
    local m=$1 d=$2 y=$3 h=$4 mdy t
    if ((JOB_DUR_MTH > 0)); then
        mdy=$(valid_date "$((m + JOB_DUR_MTH))" "$d" "$y") || return
        read -r m d y <<< "$mdy"
    elif ((JOB_DUR_DAY > 0)); then
        mdy=$(valid_date "$m" "$((d + JOB_DUR_DAY))" "$y") || return
        read -r m d y <<< "$mdy"
    else
        t=$(date2hour "$m" "$d" "$y" "$h") || return
        mdy=$(hour2date "$((t + JOB_DUR_HOUR))") || return
        read -r m d y h <<< "$mdy"
    fi
    printf '%s %s %s %s\n' "$m" "$d" "$y" "$h"
}

# Missing HOUR_BEGIN_JOB means midnight; missing JOB_DUR_HOUR means zero.
export HOUR_BEGIN_JOB=${HOUR_BEGIN_JOB:-0}
export JOB_DUR_HOUR=${JOB_DUR_HOUR:-0}
export JOB_DUR_DAY=${JOB_DUR_DAY:-0}
export JOB_DUR_MTH=${JOB_DUR_MTH:-0}

# Avoid shell arithmetic surprises (leading zeros / negative / missing values).
for value in "$YEAR_BEGIN_JOB" "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" \
             "$HOUR_BEGIN_JOB" "$JOB_DUR_MTH" "$JOB_DUR_DAY" \
             "$JOB_DUR_HOUR" "$NBJOB"; do
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Calendar parameters must be non-negative integers" >&2
        return 1 2>/dev/null || exit 1
    fi
done
# Convert hour to decimal, including inputs like 08 and 09.
HOUR_BEGIN_JOB=$((10#$HOUR_BEGIN_JOB))
JOB_DUR_MTH=$((10#$JOB_DUR_MTH))
JOB_DUR_DAY=$((10#$JOB_DUR_DAY))
JOB_DUR_HOUR=$((10#$JOB_DUR_HOUR))
NBJOB=$((10#$NBJOB))
if ((HOUR_BEGIN_JOB > 23 || NBJOB < 1 ||
      (JOB_DUR_MTH > 0) + (JOB_DUR_DAY > 0) + (JOB_DUR_HOUR > 0) != 1)); then
    echo "ERROR: Invalid start hour, NBJOB or job duration (choose exactly one unit)" >&2
    return 1 2>/dev/null || exit 1
fi


# Job start: existing output names.
export DATE_BEGIN_JOB=$(makedate "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" "$YEAR_BEGIN_JOB" "$HOUR_BEGIN_JOB")
export JDAY_BEGIN_JOB=$("${SCRIPTDIR}/julday.sh" \
    "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" "$YEAR_BEGIN_JOB" "$CALTYPE")
export JDAY_BEGIN_EXP=$JDAY_BEGIN_JOB
export THOUR_BEGIN_JOB=$(date2hour "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" \
    "$YEAR_BEGIN_JOB" "$HOUR_BEGIN_JOB")

# First job: exclusive end is also start of job N+1.
read -r MONTH_BEGIN_JOBp1 DAY_BEGIN_JOBp1 YEAR_BEGIN_JOBp1 HOUR_BEGIN_JOBp1 \
    <<< "$(next_job_date "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" "$YEAR_BEGIN_JOB" "$HOUR_BEGIN_JOB")"
export MONTH_BEGIN_JOBp1 DAY_BEGIN_JOBp1 YEAR_BEGIN_JOBp1 HOUR_BEGIN_JOBp1
export JDAY_BEGIN_JOBp1=$("${SCRIPTDIR}/julday.sh" \
    "$MONTH_BEGIN_JOBp1" "$DAY_BEGIN_JOBp1" "$YEAR_BEGIN_JOBp1" "$CALTYPE")
export THOUR_END_JOB=$(date2hour "$MONTH_BEGIN_JOBp1" "$DAY_BEGIN_JOBp1" \
    "$YEAR_BEGIN_JOBp1" "$HOUR_BEGIN_JOBp1")
export TOTAL_JOB_DUR_HOUR=$((THOUR_END_JOB - THOUR_BEGIN_JOB))
if ((TOTAL_JOB_DUR_HOUR <= 0)); then
    echo "ERROR: Job duration must be positive" >&2
    return 1 2>/dev/null || exit 1
fi
export TOTAL_JOB_DUR_SEC=$((TOTAL_JOB_DUR_HOUR * 3600 ))
# Legacy integer duration in days, rounded up for sub-day jobs.
export TOTAL_JOB_DUR=$(((TOTAL_JOB_DUR_HOUR + 23) / 24))

# Legacy END_* variables: calendar date of the LAST INCLUDED HOUR.
# Thus a midnight-to-midnight job still ends on the previous calendar day.
read -r MONTH_END_JOB DAY_END_JOB YEAR_END_JOB HOUR_END_JOB \
    <<< "$(hour2date "$((THOUR_END_JOB - 1))")"
export MONTH_END_JOB DAY_END_JOB YEAR_END_JOB HOUR_END_JOB
export DATE_END_JOB=$(makedate "$MONTH_END_JOB" "$DAY_END_JOB" "$YEAR_END_JOB" "$HOUR_END_JOB")
export JDAY_END_JOB=$("${SCRIPTDIR}/julday.sh" \
    "$MONTH_END_JOB" "$DAY_END_JOB" "$YEAR_END_JOB" "$CALTYPE")

# Dates for namelists
printf -v START_DATE_RUN '%04d-%02d-%02d %02d:00:00' \
    "$YEAR_BEGIN_JOB" "$MONTH_BEGIN_JOB" "$DAY_BEGIN_JOB" "$HOUR_BEGIN_JOB"

printf -v END_DATE_RUN '%04d-%02d-%02d %02d:00:00' \
    "$YEAR_BEGIN_JOBp1" "$MONTH_BEGIN_JOBp1" "$DAY_BEGIN_JOBp1" "$HOUR_BEGIN_JOBp1"

export START_DATE_RUN END_DATE_RUN

# Experiment end: advance NBJOB times, not NBJOB*MONTHS (month lengths vary).
CHECK_M=$MONTH_BEGIN_JOB
CHECK_D=$DAY_BEGIN_JOB
CHECK_Y=$YEAR_BEGIN_JOB
CHECK_H=$HOUR_BEGIN_JOB
for ((i=0; i<NBJOB; i++)); do
    read -r CHECK_M CHECK_D CHECK_Y CHECK_H \
        <<< "$(next_job_date "$CHECK_M" "$CHECK_D" "$CHECK_Y" "$CHECK_H")"
done
export CHECK_M CHECK_D CHECK_Y CHECK_H
export THOUR_END_EXP=$(date2hour "$CHECK_M" "$CHECK_D" "$CHECK_Y" "$CHECK_H")
read -r MONTH_END_EXP DAY_END_EXP YEAR_END_EXP HOUR_END_EXP \
    <<< "$(hour2date "$((THOUR_END_EXP - 1))")"
export MONTH_END_EXP DAY_END_EXP YEAR_END_EXP HOUR_END_EXP
export DATE_END_EXP=$(makedate "$MONTH_END_EXP" "$DAY_END_EXP" "$YEAR_END_EXP" "$HOUR_END_EXP")
export JDAY_END_EXP=$("${SCRIPTDIR}/julday.sh" \
    "$MONTH_END_EXP" "$DAY_END_EXP" "$YEAR_END_EXP" "$CALTYPE")

if ((THOUR_END_JOB > THOUR_END_EXP)); then
    echo "ERROR: Job ends after experiment end" >&2
    return 1 2>/dev/null || exit 1
fi

