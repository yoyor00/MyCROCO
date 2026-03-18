!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
                                        ! Read/report model input
      subroutine read_inp (ierr)        ! parameters  from  startup
                                        ! script file using keywords
                                        ! to recognize variables.
#ifdef MPI
      use mpi
#endif
!     implicit none

#include "param.h"
#include "scalars.h"
#include "ncscrum.h"
#ifdef ABL1D
# include "abl1d.h"
#endif
#ifdef SEDIMENT
# include "sediment.h"
#endif
#ifdef STATIONS
# include "nc_sta.h"
#endif
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
# include "sources.h"
#endif
#ifdef WKB_WWAVE
# include "wkb_wwave.h"
#endif
#ifdef WAVE_MAKER
# include "forces.h"
#endif
#ifdef ONLINE
# include "online.h"
#endif
#ifdef NBQ
# include "nbq.h"
#endif
#include "mpi_cpl.h"
      integer kwsize, testunit, input
      parameter (kwsize=32, testunit=40, input=15)
      character end_signal*3, keyword*32, fname*180
      parameter (end_signal='end')
      integer ierr, iargc, is,ie, kwlen, lstr, lenstr
#ifdef SOLVE3D
     &                                       , itrc
#endif
#ifdef DIAGNOSTICS_BIO
     &                                       , iflux
#endif
      logical dumboolean
#ifdef USE_CALENDAR
      real  dt_his,dt_avg,dt_rst
#endif
#ifdef XIOS
      REAL*8        :: tool_datosec
#endif

!
! Use pre-set default startup filename for known applications,
! or get it as an argument from command line via iargc-getarg

! (override default). NOTE: The usage of the executable should
! be either
!           croco
!        or
!           croco startup_file_name
!
! WITHOUT the UNIX redirection (<): croco<startup_file like it
! used to be.
!
#if defined BASIN
      fname='TEST_CASES/croco.in.Basin'
#elif defined CANYON
      fname='TEST_CASES/croco.in.Canyon'
#elif defined EQUATOR
      fname='TEST_CASES/croco.in.Equator'
#elif defined GRAV_ADJ
# ifdef NBQ
      fname='TEST_CASES/croco.in.Grav_adj_nbq'
# else
      fname='TEST_CASES/croco.in.Grav_adj'
# endif
#elif defined ISOLITON
      fname='TEST_CASES/croco.in.Isoliton'
#elif defined INNERSHELF
      fname='TEST_CASES/croco.in.Innershelf'
#elif defined INTERNAL
      fname='TEST_CASES/croco.in.Internal'
#elif defined OVERFLOW
      fname='TEST_CASES/croco.in.Overflow'
#elif defined RIVER
      fname='TEST_CASES/croco.in.River'
#elif defined SEAMOUNT
      fname='TEST_CASES/croco.in.Seamount'
#elif defined SHELFRONT
      fname='TEST_CASES/croco.in.Shelfront'
#elif defined SOLITON
      fname='TEST_CASES/croco.in.Soliton'
#elif defined UPWELLING
      fname='TEST_CASES/croco.in.Upwelling'
#elif defined VORTEX
      fname='TEST_CASES/croco.in.Vortex'
#elif defined JET
      fname='TEST_CASES/croco.in.Jet'
#elif defined SHOREFACE
      fname='TEST_CASES/croco.in.ShoreFace'
#elif defined RIP
      fname='TEST_CASES/croco.in.Rip'
#elif defined FLASH_RIP
      fname='TEST_CASES/croco.in.Flashrip'
#elif defined SANDBAR
      fname='TEST_CASES/croco.in.Sandbar'
#elif defined THACKER
      fname='TEST_CASES/croco.in.Thacker'
#elif defined SINGLE_COLUMN && defined KATO_PHILIPS
      fname='TEST_CASES/croco.in.KATO_PHILLIPS'
#elif defined SINGLE_COLUMN && defined WILLIS_DEARDORFF
      fname='TEST_CASES/croco.in.WILLIS_DEARDORFF'
#elif defined SINGLE_COLUMN && defined DIURNAL_CYCLE
      fname='TEST_CASES/croco.in.DIURNALCYCLE'
#elif defined SINGLE_COLUMN && defined FORCED_EKBBL
      fname='TEST_CASES/croco.in.FORCED_EKBBL'
#elif defined SINGLE_COLUMN && defined FORCED_DBLEEK
      fname='TEST_CASES/croco.in.FORCED_DBLEEK'
#elif defined SINGLE_COLUMN && defined FORCED_NONROTBBL
      fname='TEST_CASES/croco.in.FORCED_NONROTBBL'
#elif defined SINGLE_COLUMN && defined FORCED_OSCNONROTBBL
      fname='TEST_CASES/croco.in.FORCED_OSCNONROTBBL'
#elif defined TANK
      fname='TEST_CASES/croco.in.Tank'
#elif defined MOVING_BATHY
      fname='TEST_CASES/croco.in.MovBat'
#elif defined ACOUSTIC
      fname='TEST_CASES/croco.in.Acoustic'
#elif defined SWASH
      fname='TEST_CASES/croco.in.Swash'
#elif defined KH_INST
      fname='TEST_CASES/croco.in.KH_INST'
#elif defined IGW
# ifdef NBQ
      fname='TEST_CASES/croco.in.IGW_nbq'
# else
      fname='TEST_CASES/croco.in.IGW'
# endif
#elif defined TS_HADV_TEST
# ifdef SOLID_BODY_ROT
      fname='TEST_CASES/croco.in.TS_adv_body'
# elif defined DIAGONAL_ADV
      fname='TEST_CASES/croco.in.TS_adv_diag'
# elif defined SOLID_BODY_PER
      fname='TEST_CASES/croco.in.TS_adv_rep'
# endif
#elif defined DUNE
# ifdef DUNE3D
      fname='TEST_CASES/croco.in.Dune3d'
# elif defined ANA_DUNE
      fname='TEST_CASES/croco.in.Ana_dune'
# else
      fname='TEST_CASES/croco.in.Dune'
# endif
#elif defined SED_TOY
# ifdef SED_TOY_RESUSP
      fname='TEST_CASES/croco.in.Sed_toy_resusp'
# elif defined SED_TOY_CONSOLID
      fname='TEST_CASES/croco.in.Sed_toy_consolid'
# elif defined SED_TOY_FLOC
      fname='TEST_CASES/croco.in.Sed_toy_floc'
# elif defined SED_TOY_ROUSE
      fname='TEST_CASES/croco.in.Sed_toy_rouse'
# endif
#elif defined TIDAL_FLAT
      fname='TEST_CASES/croco.in.Tidal_flat'
#elif defined KILPATRICK
      fname='TEST_CASES/croco.in.Kilpatrick'
#elif defined ESTUARY
      fname='TEST_CASES/croco.in.Estuary'
#elif defined SEAGRASS
      fname='TEST_CASES/croco.in.Seagrass'
#else
      fname='croco.in'
#endif
#ifdef MPI
      if (mynode.eq.0 .and. iargc().GT.0) call getarg(1,fname)
      call MPI_Bcast(fname,64,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#else
      if (iargc().eq.1) call getarg(1,fname)
#endif
!
! if child grid, use an input name: fname.1, .2, .3, ...
!
#ifdef AGRIF
      if (.Not.Agrif_Root()) then
        lstr=lenstr(fname)
# ifdef AGRIF_ADAPTIVE
        fname=fname(1:lstr) / / '.1'
# else
        fname=fname(1:lstr) / / '.' / / Agrif_Cfixed()
# endif
        lstr=lenstr(fname)
      endif
#endif
!
      wrthis(indxTime)=.false.
#ifdef AVERAGES
      wrtavg(indxTime)=.false.
#endif
#ifdef ABL1D
      wrtabl(indxTime)=.false.
# ifdef AVERAGES
      wrtabl_avg(indxTime)=.false.
# endif
#endif

#ifdef TRACERS
# ifdef DIAGNOSTICS_TS
      wrtdia3D(NT+1)=.false.
#  ifdef AVERAGES
      wrtdia3D_avg(NT+1)=.false.
#  endif
#  ifdef DIAGNOSTICS_TS_MLD
      wrtdia2D(NT+1)=.false.
#   ifdef AVERAGES
      wrtdia2D_avg(NT+1)=.false.
#   endif
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_UV
      wrtdiaM(3)=.false.
# ifdef AVERAGES
      wrtdiaM_avg(3)=.false.
# endif
#endif
#ifdef DIAGNOSTICS_VRT
      wrtdiags_vrt(3)=.false.
# ifdef AVERAGES
      wrtdiags_vrt_avg(3)=.false.
# endif
#endif
#ifdef DIAGNOSTICS_EK
      wrtdiags_ek(3)=.false.
# ifdef AVERAGES
      wrtdiags_ek_avg(3)=.false.
# endif
#endif
#if defined DIAGNOSTICS_PV && defined TRACERS
      wrtdiags_pv(NT+1)=.false.
# ifdef AVERAGES
      wrtdiags_pv_avg(NT+1)=.false.
# endif
#endif
#if defined DIAGNOSTICS_EDDY && ! defined XIOS
      wrtdiags_eddy(3)=.false.
# ifdef AVERAGES
      wrtdiags_eddy_avg(3)=.false.
# endif
#endif
#if defined OUTPUTS_SURFACE && ! defined XIOS
      wrtsurf(3)=.false.
# ifdef AVERAGES
      wrtsurf_avg(3)=.false.
# endif
#endif
#ifdef DIAGNOSTICS_BIO
      wrtdiabioFlux(NumFluxTerms+1)=.false.
      wrtdiabioVSink(NumVSinkTerms+1)=.false.
      wrtdiabioGasExc(NumGasExcTerms+1)=.false.
# ifdef AVERAGES
      wrtdiabioFlux_avg(NumFluxTerms+1)=.false.
      wrtdiabioVSink_avg(NumVSinkTerms+1)=.false.
      wrtdiabioGasExc_avg(NumGasExcTerms+1)=.false.
# endif
#endif
!
! Read in keyword: keep trying, until keyword is found.
! ==== == ======== ==== ======= ===== ======= == ======
!
      ierr=0    ! <-- reset error counter
      call setup_kwds (ierr)
      open (input,file=fname,status='old',form='formatted',err=97)
   1  keyword='                                '
      read(input,'(A)',err=1,end=99) keyword
      if (ichar(keyword(1:1)).eq.33) goto 1
      is=1
   2  if (is.eq.kwsize) then
        goto 1
      elseif (keyword(is:is).eq.' ') then
        is=is+1
        goto 2
      endif
      ie=is
   3  if (keyword(ie:ie).eq.':') then
        keyword(ie:ie)=' '
        goto 4           !--> recognized keyword.
      elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
        ie=ie+1
        goto 3
      endif
      goto 1
   4  kwlen=ie-is
      if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
!
! Read input parameters according to the keyword:
! ==== ===== ========== ========= == === ========
!
! Title
!
      if (keyword(1:kwlen).eq.'title') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) title

#ifdef LOGFILE
!
! LOGFILE file name. Check its availability and open it.
!
       elseif (keyword(1:kwlen).eq.'logfile') then
         call cancel_kwd (keyword(1:kwlen), ierr)
         read(input,'(A)',err=95) fname
         lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
         call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
         logname=fname(1:lstr)
# ifndef AGRIF
         MPI_master_only stdout=80
# else
         stdout=Agrif_Get_Unit()
         if (.Not.Agrif_Root()) then
           if(mynode.eq.0) then
             open (stdout,file=logname//'.'//Agrif_Cfixed(),
     &         status='REPLACE',form='formatted',err=98)
           else
             stdout=6
           endif
          else
# endif
            if(mynode.eq.0) then
# ifdef ENSEMBLE
              open (stdout,file=cmember//logname,
     &          status='REPLACE',form='formatted',err=98)
# else
              open (stdout,file=logname,
     &          status='REPLACE',form='formatted',err=98)
# endif
            else
              stdout=6
            endif
# ifdef AGRIF
          endif
# endif
#endif /*LOGFILE*/
        

#ifdef USE_CALENDAR
!
! start date
!
      elseif (keyword(1:kwlen).eq.'start_date') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) start_date
# ifdef ANA_INITIAL
        if (nrrec.eq.0) then
          origin_date=start_date
          origin_date_in_sec=tool_datosec(origin_date)
          READ(origin_date(1:4),fmt='(i4)') origin_year
          READ(origin_date(6:7),fmt='(i2)') origin_month
          READ(origin_date(9:10),fmt='(i2)') origin_day
          READ(origin_date(12:13),fmt='(i2)') origin_hour
          READ(origin_date(15:16),fmt='(i2)') origin_minute
          READ(origin_date(18:19),fmt='(i2)') origin_second
        endif
# endif
!
! end date
!
      elseif (keyword(1:kwlen).eq.'end_date') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) end_date
! output time steps
      elseif (keyword(1:kwlen).eq.'output_time_steps') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) dt_his,dt_avg,dt_rst
#endif /*USE_CALENDAR*/

#ifdef XIOS
!
! XIOS origin date
!
      elseif (keyword(1:kwlen).eq.'xios_origin_date') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) xios_origin_date
        xios_origin_date_in_sec=tool_datosec(xios_origin_date)
#endif
!
! Time-stepping parameters
!
      elseif (keyword(1:kwlen).eq.'time_stepping') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntimes,dt,ndtfast, ninfo

        dtfast=dt/float(ndtfast)     ! set barotropic time step.

        if (NWEIGHT.lt.(2*ndtfast-1)) then
          write(stdout,'(a,i3)')
     &    ' Error - Number of 2D timesteps (2*ndtfast-1): ',
     &    2*ndtfast-1
          write(stdout,'(a,i3)')
     &    '           exceeds barotopic weight dimension: ',NWEIGHT
          goto 95
        endif

#ifdef NBQ
!
! NBQ Time-stepping parameters
!
      elseif (keyword(1:kwlen).eq.'time_stepping_nbq') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ndtnbq, csound_nbq, visc2_nbq

        dtnbq=dtfast
        ndtnbq=1
#endif

#ifdef SOLVE3D
!
! Vertical S-coordinates parameters.
!
      elseif (keyword(1:kwlen).eq.'S-coord') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) theta_s, theta_b, Tcline

#endif
!
! Initial conditions file name. Check its availability (in the case
! of analytical initial conditions and nrrec=0 initial conditions are
! created internally and no file is needed).
!
      elseif (keyword(1:kwlen).eq.'initial') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrrec
#ifdef ANA_INITIAL
        if (nrrec.gt.0) then
#endif
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
#if defined MPI && defined PARALLEL_FILES
          call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#endif
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          ininame=fname(1:lstr)

#ifdef ANA_INITIAL
        endif
#endif
#ifndef ANA_GRID
!
! Grid file name. Check its availability.
!
      elseif (keyword(1:kwlen).eq.'grid') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        open(testunit,file=fname(1:lstr), status='old', err=97)
        close(testunit)
        grdname=fname(1:lstr)

#endif


! forcing file for tides
#if !defined NO_FRCFILE  && \
        ( defined TIDES \
          || (defined MRL_WCI     && !defined ANA_WWAVE && \
                                     !defined WKB_WWAVE && \
                                     !defined OW_COUPLING) \
          || (defined SFLX_CORR   && !defined ANA_SSS   && \
                                   !defined BULK_FLUX   && \
                                   !defined OA_COUPLING) \
          || (defined QCORRECTION && !defined ANA_SST   && \
                                   !defined OA_COUPLING) \
          || (defined BBL         && !defined ANA_BSEDIM  && \
                                   !defined SEDIMENT) \
          || (!defined ANA_STFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                   defined SOLVE3D) \
          || (defined SALINITY    && !defined ANA_SSFLUX && \
                                   !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                    defined SOLVE3D) \
          || (!defined ANA_SRFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING && \
                                   defined SOLVE3D) \
          || (!defined ANA_SMFLUX && !defined BULK_FLUX && \
                                   !defined OA_COUPLING) \
          || ( defined BULK_FLUX    && !defined SALINITY ) \
          || ( defined WAVE_OFFLINE && !defined MUSTANG ) )
!
!
! Forcing file name. Check its availability.
!
      elseif (keyword(1:kwlen).eq.'forcing') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        frcname=fname(1:lstr)

#endif  /* !defined NO_FRCFILE */

!

#if defined WAVE_OFFLINE && defined MUSTANG
! Forcing file name. Check its availability.
!
      elseif (keyword(1:kwlen).eq.'wave_offline') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode, NNODES, ierr)
# endif
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        wave_file=fname(1:lstr)

#endif  /* !defined WAVE_OFFLINE */

!
! Biology forcing: iron dust deposition.
!
#if defined BIOLOGY && defined PISCES
      elseif (keyword(1:kwlen).eq.'biology') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
c        open (testunit, file=fname(1:lstr), status='old', err=97)
c        close(testunit)
        bioname=fname(1:lstr)

#endif
!
! Bulk file name. Check its availability.
!
#if defined BULK_FLUX && !defined ANA_ABL_LSDATA
# ifndef ONLINE /* ONLINE FORCING */
        elseif (keyword(1:kwlen).eq.'bulk_forcing') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
          call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          bulkname=fname(1:lstr)

# endif  /* ONLINE FORCING */
#endif


#if (defined TCLIMATOLOGY  && !defined ANA_TCLIMA) || \
      (defined ZCLIMATOLOGY  && !defined ANA_SSH) || \
      (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
      (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
!
! Climatology file name. Check availability.
!
      elseif (keyword(1:kwlen).eq.'climatology') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef AGRIF
       if (Agrif_Root()) then
# endif
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        clmname=fname(1:lstr)

# ifdef AGRIF
       endif
# endif
#endif


#if !defined ANA_BRY && defined FRC_BRY
!
! Boundary file name. Check availability.
!
        elseif (keyword(1:kwlen).eq.'boundary') then
          call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef AGRIF
       if (Agrif_Root()) then
# endif
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
          call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          bry_file=fname(1:lstr)
# ifdef AGRIF
       endif
# endif
#endif



#if defined WKB_WWAVE && !defined ANA_BRY_WKB
!
! WKB boundary file name. Check availability.
!
        elseif (keyword(1:kwlen).eq.'wkb_boundary') then
          call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef AGRIF
       if (Agrif_Root()) then
# endif
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
          call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          brywkb_file=fname(1:lstr)
# ifdef AGRIF
       endif
# endif
#endif
!
! Restart file name.
!
      elseif (keyword(1:kwlen).eq.'restart') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrst, nrpfrst
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
#if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#endif
#ifdef USE_CALENDAR
        nrst=ceiling((dt_rst*3600.)/dt)
#endif
        rstname=fname(1:lstr)
!
! History file name.
!
      elseif (keyword(1:kwlen).eq.'history') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefhis, nwrt, nrpfhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#endif
#ifdef USE_CALENDAR
        nwrt=ceiling((dt_his*3600.)/dt)
#endif
        hisname=fname(1:lstr)

#ifdef AVERAGES
!
! Averages file name.
!
      elseif (keyword(1:kwlen).eq.'averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntsavg, navg, nrpfavg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
# ifdef USE_CALENDAR
        navg=ceiling((dt_avg*3600.)/dt)
# endif
        avgname=fname(1:lstr)

#endif

#if defined ABL1D && ! defined XIOS
!
! ABL file name.
!
      elseif (keyword(1:kwlen).eq.'abl') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefablhis, nwrtablhis, nrpfablhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        ablname=fname(1:lstr)

# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'abl_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefablavg, ntsablavg, nwrtablavg,
     &                       nrpfablavg
        if (nwrtablavg.eq.0) nwrtablavg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        ablname_avg=fname(1:lstr)

# endif
#endif /* ABL1D */

#if defined DIAGNOSTICS_TS
!
! Diagnostics file name.
!
      elseif (keyword(1:kwlen).eq.'diagnostics') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdia, nwrtdia, nrpfdia
        if (nwrtdia.eq.0) nwrtdia = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        dianame=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diag_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdia_avg, ntsdia_avg, nwrtdia_avg,
     &                                                nrpfdia_avg
        if (nwrtdia_avg.eq.0) nwrtdia_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        dianame_avg=fname(1:lstr)
# endif
#endif /* DIAGOSTICS_TS */

#if defined DIAGNOSTICS_UV
!
! Diagnostics Momentum file name.
!
      elseif (keyword(1:kwlen).eq.'diagnosticsM') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiaM, nwrtdiaM, nrpfdiaM
        if (nwrtdiaM.eq.0) nwrtdiaM = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        dianameM=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diagM_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiaM_avg, ntsdiaM_avg, nwrtdiaM_avg,
     &                                                  nrpfdiaM_avg
        if (nwrtdiaM_avg.eq.0) nwrtdiaM_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        dianameM_avg=fname(1:lstr)

# endif
#endif /*DIAGNOSTICS_UV */


#if defined DIAGNOSTICS_VRT
!
! Diagnostics vrt file name.
!
      elseif (keyword(1:kwlen).eq.'diags_vrt') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_vrt, nwrtdiags_vrt, nrpfdiags_vrt
        if (nwrtdiags_vrt.eq.0) nwrtdiags_vrt = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        diags_vrtname=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diags_vrt_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_vrt_avg, ntsdiags_vrt_avg,
     &                               nwrtdiags_vrt_avg, nrpfdiags_vrt_avg
        if (nwrtdiags_vrt_avg.eq.0) nwrtdiags_vrt_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        diags_vrtname_avg=fname(1:lstr)

# endif
#endif /*DIAGNOSTICS_VRT */


#if defined DIAGNOSTICS_EK
!
! Diagnostics ek file name.
!
      elseif (keyword(1:kwlen).eq.'diags_ek') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_ek, nwrtdiags_ek, nrpfdiags_ek
        if (nwrtdiags_ek.eq.0) nwrtdiags_ek = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        diags_ekname=fname(1:lstr)

# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diags_ek_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_ek_avg, ntsdiags_ek_avg,
     &                      nwrtdiags_ek_avg,  nrpfdiags_ek_avg
        if (nwrtdiags_ek_avg.eq.0) nwrtdiags_ek_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        diags_ekname_avg=fname(1:lstr)

# endif
#endif /*DIAGNOSTICS_EK */

#if defined DIAGNOSTICS_PV
!
! Diagnostics pv file name.
!
      elseif (keyword(1:kwlen).eq.'diags_pv') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_pv, nwrtdiags_pv, nrpfdiags_pv
        if (nwrtdiags_pv.eq.0) nwrtdiags_pv = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        diags_pvname=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diags_pv_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_pv_avg, ntsdiags_pv_avg,
     &                      nwrtdiags_pv_avg,  nrpfdiags_pv_avg
        if (nwrtdiags_pv_avg.eq.0) nwrtdiags_pv_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        diags_pvname_avg=fname(1:lstr)

# endif
#endif /*DIAGNOSTICS_PV */

#if defined DIAGNOSTICS_EDDY && ! defined XIOS
!
! Diagnostics eddy file name.
!
      elseif (keyword(1:kwlen).eq.'diags_eddy') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_eddy, nwrtdiags_eddy, nrpfdiags_eddy
        if (nwrtdiags_eddy.eq.0) nwrtdiags_eddy = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        diags_eddyname=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diags_eddy_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_eddy_avg, ntsdiags_eddy_avg,
     &                      nwrtdiags_eddy_avg,  nrpfdiags_eddy_avg
        if (nwrtdiags_eddy_avg.eq.0) nwrtdiags_eddy_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        diags_eddyname_avg=fname(1:lstr)

# endif
#endif /*DIAGNOSTICS_EDDY */

#if defined OUTPUTS_SURFACE && ! defined XIOS
!
! surface outputs file name.
!
      elseif (keyword(1:kwlen).eq.'surf') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefsurf, nwrtsurf, nrpfsurf
        if (nwrtsurf.eq.0) nwrtsurf = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        surfname=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'surf_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefsurf_avg, ntssurf_avg,
     &                      nwrtsurf_avg,  nrpfsurf_avg
        if (nwrtsurf_avg.eq.0) nwrtsurf_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        surfname_avg=fname(1:lstr)

# endif
#endif /*OUTPUTS_SURFACE */


#ifdef DIAGNOSTICS_BIO
!
! Diagnostics Biology file name.
!
      elseif (keyword(1:kwlen).eq.'diagnostics_bio') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiabio, nwrtdiabio, nrpfdiabio
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
        dianamebio=fname(1:lstr)
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'diagbio_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiabio_avg, ntsdiabio_avg,
     &                       nwrtdiabio_avg, nrpfdiabio_avg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
#  if defined MPI && defined PARALLEL_FILES
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
#  endif
        dianamebio_avg=fname(1:lstr)
# endif
#endif /* DIAGNOSTICS_BIO */

#ifdef STATIONS
!
! Stations file name.
!
      elseif (keyword(1:kwlen).eq.'stations') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef AGRIF
        if (Agrif_Root()) then
# endif
          read(input,*,err=95) ldefsta, nsta, nrpfsta
          read(input,'(A)',err=95) staposname
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES
          call insert_node (staposname, lstr, mynode2, NNODES2, ierr)
# endif
          staname=fname(1:lstr)
# ifdef AGRIF
        else
          ldefsta=Agrif_Parent(ldefsta)
          nsta=Agrif_Parent(nsta)
          nrpfsta=Agrif_Parent(nrpfsta)
        endif
# endif
#endif /* STATIONS */

#ifdef ASSIMILATION
!
! Assimilation input/output file names.
!
      elseif (keyword(1:kwlen).eq.'assimilation') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) aparnam
        read(input,'(A)',err=95) assname
        fname=aparnam
        lstr=lenstr(aparnam)
        open (testunit,file=aparnam(1:lstr),status='old',err=97)
        close(testunit)
        fname=assname
        lstr=lenstr(assname)
        open (testunit,file=assname(1:lstr),status='old',err=97)
        close(testunit)
#endif
!
! Switches for fields to be saved into history file.
!
      elseif (keyword(1:kwlen).eq.'primary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                       ,  wrthis(indxVb)
#ifdef SOLVE3D
     &                    ,  wrthis(indxU),  wrthis(indxV)
# ifdef TRACERS
     &                    , (wrthis(itrc), itrc=indxV+1,indxV+NT)
        do itrc=1,NT
          if (wrthis(indxV+itrc)) wrthis(indxTime)=.true.
# ifdef DIAGNOSTICS_TS_MLD_CRIT
          if (itrc.eq.iCRT2 .or. itrc.eq.iCRT3
     &   .or. itrc.eq.iCRT4) then
            wrthis(indxV+itrc)=.false.
          endif

# endif
        enddo
# endif
#endif
        if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
#ifdef SOLVE3D
     &                        .or. wrthis(indxU) .or. wrthis(indxV)
#endif
     &     ) wrthis(indxTime)=.true.

#if !defined SOLVE3D && defined RIP
      elseif (keyword(1:kwlen).eq.'auxiliary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                                             wrthis(indxBostr)
     &                                          ,  wrthis(indxBustr)
     &                                          ,  wrthis(indxBvstr)
     &                                          ,  wrthis(indxWstr)
     &                                          ,  wrthis(indxUWstr)
     &                                          ,  wrthis(indxVWstr)
#elif defined SOLVE3D
      elseif (keyword(1:kwlen).eq.'auxiliary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                                             wrthis(indxR)
     &                                          ,  wrthis(indxO)
     &                                          ,  wrthis(indxW)
     &                                          ,  wrthis(indxAkv)
# ifdef TEMPERATURE
     &                                          ,  wrthis(indxAkt)
# else
     &                                          ,  dumboolean
# endif
# ifdef SALINITY
     &                                          ,  wrthis(indxAks)
# else
     &                                          ,  dumboolean
# endif
# if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
     &                                          ,  wrthis(indxbvf)
# else
     &                                          ,  dumboolean
# endif
# ifdef VIS_COEF_3D
     &                                          ,  wrthis(indxVisc)
# else
     &                                          ,  dumboolean
# endif
# ifdef DIF_COEF_3D
     &                                          ,  wrthis(indxDiff)
# else
     &                                          ,  dumboolean
# endif
# if defined LMD_SKPP || defined GLS_MIXING
     &                                          ,  wrthis(indxHbl)
# else
     &                                          ,  dumboolean
# endif
# ifdef LMD_BKPP
     &                                          ,  wrthis(indxHbbl)
# else
     &                                          ,  dumboolean
# endif
     &                                          ,  wrthis(indxBostr)
     &                                          ,  wrthis(indxBustr)
     &                                          ,  wrthis(indxBvstr)
     &                                          ,  wrthis(indxWstr)
     &                                          ,  wrthis(indxUWstr)
     &                                          ,  wrthis(indxVWstr)
# ifdef TEMPERATURE
     &                                          ,  wrthis(indxShflx)
# else
     &                                          ,  dumboolean
# endif
# ifdef SALINITY
     &                                          ,  wrthis(indxSwflx)
# else
     &                                          ,  dumboolean
# endif
# ifdef TEMPERATURE
     &                                          ,  wrthis(indxShflx_rsw)
# else
     &                                          ,  dumboolean
# endif
# ifdef BULK_FLUX
     &                                          ,  wrthis(indxShflx_rlw)
     &                                          ,  wrthis(indxShflx_lat)
     &                                          ,  wrthis(indxShflx_sen)
# else
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
# endif
# ifdef BHFLUX
     &                                          ,  wrthis(indxBhflx)
# else
     &                                          ,  dumboolean
# endif
# if defined BWFLUX && defined SALINITY
     &                                          ,  wrthis(indxBwflx)
# else
     &                                          ,  dumboolean
# endif
# ifdef MORPHODYN
     &                                          ,  wrthis(indxHm)
# else
     &                                          ,  dumboolean
# endif
# if defined BIOLOGY && !defined PISCES
     &                                          ,  wrthis(indxHel)
#  ifdef BIO_NChlPZD
     &                                          ,  wrthis(indxChC)
#   ifdef OXYGEN
     &                                          ,  wrthis(indxU10)
     &                                          ,  wrthis(indxKvO2)
     &                                          ,  wrthis(indxO2sat)
#   endif
#  elif defined BIO_BioEBUS
     &                                          ,  wrthis(indxAOU)
     &                                          ,  wrthis(indxWIND10)
#  endif
# endif

        if ( wrthis(indxR)
     &                                        .or. wrthis(indxO)
     &                                        .or. wrthis(indxW)
     &                                        .or. wrthis(indxAkv)
# ifdef TEMPERATURE
     &                                        .or. wrthis(indxAkt)
# endif
# ifdef SALINITY
     &                                        .or. wrthis(indxAks)
# endif
# if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
     &                                        .or. wrthis(indxbvf)
# endif
# ifdef VIS_COEF_3D
     &                                        .or. wrthis(indxVisc)
# endif
# ifdef DIF_COEF_3D
     &                                        .or. wrthis(indxDiff)
# endif
# if defined LMD_SKPP || defined GLS_MIXING
     &                                        .or. wrthis(indxHbl)
# endif
# ifdef LMD_BKPP
     &                                        .or. wrthis(indxHbbl)
# endif
     &                                        .or. wrthis(indxBostr)
     &                                        .or. wrthis(indxBustr)
     &                                        .or. wrthis(indxBvstr)
     &                                        .or. wrthis(indxWstr)
     &                                        .or. wrthis(indxUWstr)
     &                                        .or. wrthis(indxVWstr)
# ifdef TEMPERATURE
     &                                        .or. wrthis(indxShflx)
# endif
# ifdef SALINITY
     &                                        .or. wrthis(indxSwflx)
# endif
# ifdef TEMPERATURE
     &                                        .or. wrthis(indxShflx_rsw)
# endif
# if defined BULK_FLUX
     &                                        .or. wrthis(indxShflx_rlw)
     &                                        .or. wrthis(indxShflx_lat)
     &                                        .or. wrthis(indxShflx_sen)
# endif
# if defined BHFLUX
     &                                        .or. wrthis(indxBhflx)
# endif
# if defined BWFLUX && defined SALINITY
     &                                        .or. wrthis(indxBwflx)
# endif
# if defined BIOLOGY && !defined PISCES
     &                                        .or. wrthis(indxHel)
#  ifdef BIO_NChlPZD
     &                                        .or. wrthis(indxChC)
#   ifdef OXYGEN
     &                                        .or. wrthis(indxU10)
     &                                        .or. wrthis(indxKvO2)
     &                                        .or. wrthis(indxO2sat)
#   endif
#  elif defined BIO_BioEBUS
     &                                        .or. wrthis(indxAOU)
     &                                        .or. wrthis(indxWIND10)
#  endif
# endif
     &     ) wrthis(indxTime)=.true.

!
! Switches for GLS fields to be saved into history file.
!
# ifdef GLS_MIXING
      elseif (keyword(1:kwlen).eq.'gls_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxTke),  wrthis(indxGls)
     &                    ,  wrthis(indxLsc)
        if (wrthis(indxTke) .or. wrthis(indxGls) .or. wrthis(indxLsc)
     &     ) wrthis(indxTime)=.true.

# endif
#endif /* SOLVE3D */

#ifdef AVERAGES
!
! Switches for fields to be saved into averages file.
!
      elseif (keyword(1:kwlen).eq.'primary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxZ),  wrtavg(indxUb)
     &                                    ,  wrtavg(indxVb)
# ifdef SOLVE3D
     &                    ,  wrtavg(indxU),  wrtavg(indxV)
#  ifdef TRACERS
     &                    , (wrtavg(itrc), itrc=indxV+1,indxV+NT)
#  endif
# endif
        if ( wrtavg(indxZ) .or. wrtavg(indxUb) .or. wrtavg(indxVb)
# ifdef SOLVE3D
     &                     .or. wrtavg(indxU)  .or. wrtavg(indxV)
# endif
     &     ) wrtavg(indxTime)=.true.

# ifdef SOLVE3D
#  ifdef TRACERS
         do itrc=1,NT
          if (wrtavg(indxV+itrc)) wrtavg(indxTime)=.true.
        enddo
#  endif

      elseif (keyword(1:kwlen).eq.'auxiliary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxR), wrtavg(indxO)
     &        ,  wrtavg(indxW),  wrtavg(indxAkv)
#  ifdef TEMPERATURE
     &                                          ,  wrtavg(indxAkt)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef SALINITY
     &                                          ,  wrtavg(indxAks)
#  else
     &                                          ,  dumboolean
#  endif
#  if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
     &                                          ,  wrtavg(indxbvf)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef VIS_COEF_3D
     &                                          ,  wrtavg(indxVisc)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef DIF_COEF_3D
     &                                          ,  wrtavg(indxDiff)
#  else
     &                                          ,  dumboolean
#  endif
#  if defined LMD_SKPP || defined GLS_MIXING
     &                                          ,  wrtavg(indxHbl)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef LMD_BKPP
     &                                          ,  wrtavg(indxHbbl)
#  else
     &                                          ,  dumboolean
#  endif
     &                                          ,  wrtavg(indxBostr)
     &                                          ,  wrtavg(indxBustr)
     &                                          ,  wrtavg(indxBvstr)
     &                                          ,  wrtavg(indxWstr)
     &                                          ,  wrtavg(indxUWstr)
     &                                          ,  wrtavg(indxVWstr)
#  ifdef TEMPERATURE
     &                                          ,  wrtavg(indxShflx)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef SALINITY
     &                                          ,  wrtavg(indxSwflx)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef TEMPERATURE
     &                                          ,  wrtavg(indxShflx_rsw)
#  else
     &                                          ,  dumboolean
#  endif
#  ifdef BULK_FLUX
     &                                          ,  wrtavg(indxShflx_rlw)
     &                                          ,  wrtavg(indxShflx_lat)
     &                                          ,  wrtavg(indxShflx_sen)
#  else
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
#  endif
#  ifdef BHFLUX
     &                                          , wrtavg(indxBhflx)
#  else
     &                                          , dumboolean
#  endif
#  if defined BWFLUX && defined SALINITY
     &                                          , wrtavg(indxBwflx)
#  else
     &                                          , dumboolean
#  endif
#  ifdef MORPHODYN
     &                                          ,  wrtavg(indxHm)
#  else
     &                                          ,  dumboolean
#  endif
#  if defined BIOLOGY && !defined PISCES
     &                                          ,  wrtavg(indxHel)
#   ifdef BIO_NChlPZD
     &                                          ,  wrtavg(indxChC)
#    ifdef OXYGEN
     &                                          ,  wrtavg(indxU10)
     &                                          ,  wrtavg(indxKvO2)
     &                                          ,  wrtavg(indxO2sat)
#    endif
#   elif defined BIO_BioEBUS
     &                                          ,  wrtavg(indxAOU)
     &                                          ,  wrtavg(indxWIND10)
#   endif
#  endif

        if ( wrtavg(indxR) .or. wrtavg(indxO) .or. wrtavg(indxW)
     &                   .or. wrtavg(indxAkv)
#  ifdef TEMPERATURE
     &                                        .or. wrtavg(indxAkt)
#  endif
#  ifdef SALINITY
     &                                        .or. wrtavg(indxAks)
#  endif
#  if defined ANA_VMIX || defined LMD_MIXING \
  || defined LMD_SKPP || defined LMD_BKPP \
  || defined GLS_MIXING
     &                                        .or. wrtavg(indxbvf)
#  endif
#  ifdef VIS_COEF_3D
     &                                        .or. wrtavg(indxVisc)
#  endif
#  ifdef DIF_COEF_3D
     &                                        .or. wrtavg(indxDiff)
#  endif
#  if defined LMD_SKPP || defined GLS_MIXING
     &                                        .or. wrtavg(indxHbl)
#  endif
#  ifdef LMD_BKPP
     &                                        .or. wrtavg(indxHbbl)
#  endif
     &                                        .or. wrtavg(indxBostr)
     &                                        .or. wrtavg(indxBustr)
     &                                        .or. wrtavg(indxBvstr)
     &                                        .or. wrtavg(indxWstr)
     &                                        .or. wrtavg(indxUWstr)
     &                                        .or. wrtavg(indxVWstr)
#  ifdef TEMPERATURE
     &                                        .or. wrtavg(indxShflx)
#  endif
#  ifdef SALINITY
     &                                        .or. wrtavg(indxSwflx)
#  endif
#  ifdef TEMPERATURE
     &                                        .or. wrtavg(indxShflx_rsw)
#  endif
     &
#  ifdef BULK_FLUX
     &                                        .or. wrtavg(indxShflx_rlw)
     &                                        .or. wrtavg(indxShflx_lat)
     &                                        .or. wrtavg(indxShflx_sen)
#  endif
#  if defined BHFLUX
     &                                        .or. wrtavg(indxBhflx)
#  endif
#  if defined BWFLUX && defined SALINITY
     &                                        .or. wrtavg(indxBwflx)
#  endif
#  if defined BIOLOGY && !defined PISCES
     &                                        .or. wrtavg(indxHel)
#   ifdef BIO_NChlPZD
     &                                        .or. wrtavg(indxChC)
#    ifdef OXYGEN
     &                                        .or. wrtavg(indxU10)
     &                                        .or. wrtavg(indxKvO2)
     &                                        .or. wrtavg(indxO2sat)
#    endif
#   elif defined BIO_BioEBUS
     &                                     .or. wrtavg(indxAOU)
     &                                     .or. wrtavg(indxWIND10)
#   endif
#  endif  /* BIOLOGY */
     &     ) wrtavg(indxTime)=.true.

!
! Switches for GLS fields to be saved into history file.
!
#  ifdef GLS_MIXING
      elseif (keyword(1:kwlen).eq.'gls_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxTke),  wrtavg(indxGls)
     &                    ,  wrtavg(indxLsc)
        if ( wrtavg(indxAkk) .or. wrtavg(indxAkp) .or. wrtavg(indxTke)
     &                       .or. wrtavg(indxGls) .or. wrtavg(indxLsc)
     &     ) wrtavg(indxTime)=.true.

#  endif
# endif /* SOLVE3D */
#endif /* AVERAGES */

#if defined ABL1D && ! defined XIOS
!
! Switches for fields to be saved into ABL output file.
!
      elseif (keyword(1:kwlen).eq.'abl_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &       wrtabl(indxabl_pu_dta)  , wrtabl(indxabl_pv_dta)  ,
     &       wrtabl(indxabl_pt_dta)  , wrtabl(indxabl_pq_dta)  ,
     &       wrtabl(indxabl_pgu_dta) , wrtabl(indxabl_pgv_dta) ,
     &       wrtabl(indxabl_u_abl)   , wrtabl(indxabl_v_abl)   ,
     &       wrtabl(indxabl_t_abl)   , wrtabl(indxabl_q_abl)   ,
     &       wrtabl(indxabl_tke_abl) , wrtabl(indxabl_mxlm_abl),
     &       wrtabl(indxabl_mxld_abl), wrtabl(indxabl_avm_abl) ,
     &       wrtabl(indxabl_avt_abl) , wrtabl(indxabl_ablh_abl),
     &       wrtabl(indxabl_zr_abl)  , wrtabl(indxabl_zw_abl)  ,
     &       wrtabl(indxabl_Hzr_abl) , wrtabl(indxabl_Hzw_abl)
!
        if ( wrtabl(indxabl_pu_dta)   .or. wrtabl(indxabl_pv_dta)   .or.
     &       wrtabl(indxabl_pt_dta)   .or. wrtabl(indxabl_pq_dta)   .or.
     &       wrtabl(indxabl_pgu_dta)  .or. wrtabl(indxabl_pgv_dta)  .or.
     &       wrtabl(indxabl_u_abl)    .or. wrtabl(indxabl_v_abl)    .or.
     &       wrtabl(indxabl_t_abl)    .or. wrtabl(indxabl_q_abl)    .or.
     &       wrtabl(indxabl_tke_abl)  .or. wrtabl(indxabl_mxlm_abl) .or.
     &       wrtabl(indxabl_mxld_abl) .or. wrtabl(indxabl_avm_abl)  .or.
     &       wrtabl(indxabl_avt_abl)  .or. wrtabl(indxabl_ablh_abl) .or.
     &       wrtabl(indxabl_zr_abl)   .or. wrtabl(indxabl_zw_abl)   .or.
     &       wrtabl(indxabl_Hzr_abl)  .or. wrtabl(indxabl_Hzw_abl)
     &     ) wrtabl(indxTime)=.true.
!
# ifdef AVERAGES
      elseif (keyword(1:kwlen).eq.'abl_averages_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &       wrtabl_avg(indxabl_pu_dta)  , wrtabl_avg(indxabl_pv_dta)  ,
     &       wrtabl_avg(indxabl_pt_dta)  , wrtabl_avg(indxabl_pq_dta)  ,
     &       wrtabl_avg(indxabl_pgu_dta) , wrtabl_avg(indxabl_pgv_dta) ,
     &       wrtabl_avg(indxabl_u_abl)   , wrtabl_avg(indxabl_v_abl)   ,
     &       wrtabl_avg(indxabl_t_abl)   , wrtabl_avg(indxabl_q_abl)   ,
     &       wrtabl_avg(indxabl_tke_abl) , wrtabl_avg(indxabl_mxlm_abl),
     &       wrtabl_avg(indxabl_mxld_abl), wrtabl_avg(indxabl_avm_abl) ,
     &       wrtabl_avg(indxabl_avt_abl) , wrtabl_avg(indxabl_ablh_abl),
     &       wrtabl_avg(indxabl_zr_abl)  , wrtabl_avg(indxabl_zw_abl)  ,
     &       wrtabl_avg(indxabl_Hzr_abl) , wrtabl_avg(indxabl_Hzw_abl)
!
        if ( wrtabl_avg(indxabl_pu_dta)   .or. wrtabl_avg(indxabl_pv_dta)   .or.
     &       wrtabl_avg(indxabl_pt_dta)   .or. wrtabl_avg(indxabl_pq_dta)   .or.
     &       wrtabl_avg(indxabl_pgu_dta)  .or. wrtabl_avg(indxabl_pgv_dta)  .or.
     &       wrtabl_avg(indxabl_u_abl)    .or. wrtabl_avg(indxabl_v_abl)    .or.
     &       wrtabl_avg(indxabl_t_abl)    .or. wrtabl_avg(indxabl_q_abl)    .or.
     &       wrtabl_avg(indxabl_tke_abl)  .or. wrtabl_avg(indxabl_mxlm_abl) .or.
     &       wrtabl_avg(indxabl_mxld_abl) .or. wrtabl_avg(indxabl_avm_abl)  .or.
     &       wrtabl_avg(indxabl_avt_abl)  .or. wrtabl_avg(indxabl_ablh_abl) .or.
     &       wrtabl_avg(indxabl_zr_abl)   .or. wrtabl_avg(indxabl_zw_abl)   .or.
     &       wrtabl_avg(indxabl_Hzr_abl)  .or. wrtabl_avg(indxabl_Hzw_abl)
     &     ) wrtabl_avg(indxTime)=.true.
!
# endif
#endif /* ABL1D && ! defined XIOS */

#if defined DIAGNOSTICS_TS
# if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
      elseif (keyword(1:kwlen).eq.'diag_mld_crit_threshold') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) mld_crit_D, mld_crit_T, mld_crit_T2
      elseif (keyword(1:kwlen).eq.'diag_mld_depth_ref') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) mld_depth_ref
# endif  /* DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT */
!
! Switches for fields to be saved into tracer diagnostics file.
!
      elseif (keyword(1:kwlen).eq.'diag3D_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef TRACERS
        read(input,*,err=95)  (wrtdia3D(itrc), itrc=1, NT)
# ifdef DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
        if (itrc.eq.iCRT2 .or. itrc.eq.iCRT3
     &      .or. itrc.eq.iCRT4) then
          wrtdia3D(itrc)=.false.
        endif
# endif  /* DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT */

        do itrc=1,NT
          if (wrtdia3D(itrc)) wrtdia3D(NT+1)=.true.
        enddo
# endif  /* TRACERS */
# ifdef DIAGNOSTICS_TS_MLD
      elseif (keyword(1:kwlen).eq.'diag2D_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
        read(input,*,err=95)  (wrtdia2D(itrc), itrc=1, NT)
        do itrc=1,NT
          if (wrtdia2D(itrc)) wrtdia2D(NT+1)=.true.
        enddo
#  endif  /* TRACERS */
# endif  /* DIAGNOSTICS_TS_MLD */

# ifdef AVERAGES
!
! Switches for fields to be saved into tracer diagnostics average file.
!
      elseif (keyword(1:kwlen).eq.'diag3D_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
        read(input,*,err=95) (wrtdia3D_avg(itrc),itrc=1,NT)
# ifdef DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
        if (itrc.eq.iCRT2 .or. itrc.eq.iCRT3
     &     .or. itrc.eq.iCRT4) then
            wrtdia3D_avg(itrc)=.false.
        endif
# endif  /* DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT */


        do itrc=1,NT
          if (wrtdia3D_avg(itrc)) wrtdia3D_avg(NT+1)=.true.
        enddo
#  endif  /* TRACERS */
!
#  ifdef DIAGNOSTICS_TS_MLD
      elseif (keyword(1:kwlen).eq.'diag2D_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#   ifdef TRACERS
        read(input,*,err=95) (wrtdia2D_avg(itrc),itrc=1,NT)
        do itrc=1,NT
          if (wrtdia2D_avg(itrc)) wrtdia2D_avg(NT+1)=.true.
        enddo
#   endif  /* TRACERS */
#  endif  /* DIAGNOSTICS_TS_MLD */
# endif  /* AVERAGES */

#endif  /* DIAGNOSTICS_TS */


#if defined DIAGNOSTICS_UV
!
! Switches for fields to be saved into momentum diagnostics file.
!
!!          wrtdiaM=.true.
      elseif (keyword(1:kwlen).eq.'diagM_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  (wrtdiaM(itrc), itrc=1,2)
         do itrc=1,2
           if (wrtdiaM(itrc)) wrtdiaM(3)=.true.
        enddo

# ifdef AVERAGES
!!         wrtdiaM_avg=.true.
      elseif (keyword(1:kwlen).eq.'diagM_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiaM_avg(itrc),itrc=1,2)
        do itrc=1,2
          if (wrtdiaM_avg(itrc)) wrtdiaM_avg(3)=.true.
       enddo
# endif /* AVERAGES  */

#endif /*DIAGNOSTICS_UV */

#if defined DIAGNOSTICS_VRT
!
! Switches for fields to be saved into vorticity diagnostics file.
!
!!          wrtdiags_vrt=.true.
      elseif (keyword(1:kwlen).eq.'diags_vrt_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  (wrtdiags_vrt(itrc), itrc=1,1)
         do itrc=1,1
           if (wrtdiags_vrt(itrc)) wrtdiags_vrt(3)=.true.
        enddo

# ifdef AVERAGES
!!         wrtdiags_vrt_avg=.true.
      elseif (keyword(1:kwlen).eq.'diags_vrt_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiags_vrt_avg(itrc),itrc=1,1)
        do itrc=1,1
          if (wrtdiags_vrt_avg(itrc)) wrtdiags_vrt_avg(3)=.true.
       enddo
# endif /* AVERAGES  */

#endif /*DIAGNOSTICS_VRT */

#if defined DIAGNOSTICS_EK
!
! Switches for fields to be saved into kinetic energy diagnostics file.
!
!!          wrtdiags_ek=.true.
      elseif (keyword(1:kwlen).eq.'diags_ek_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  (wrtdiags_ek(itrc), itrc=1,1)
         do itrc=1,1
           if (wrtdiags_ek(itrc)) wrtdiags_ek(3)=.true.
        enddo

# ifdef AVERAGES
!!         wrtdiags_ek_avg=.true.
      elseif (keyword(1:kwlen).eq.'diags_ek_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiags_ek_avg(itrc),itrc=1,1)
        do itrc=1,1
          if (wrtdiags_ek_avg(itrc)) wrtdiags_ek_avg(3)=.true.
       enddo
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_EK */


#if defined DIAGNOSTICS_PV
!
! Switches for fields to be saved into potential vorticity diagnostics file.
!
!!          wrtdiags_pv=.true.
      elseif (keyword(1:kwlen).eq.'diags_pv_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef TRACERS
        read(input,*,err=95)  (wrtdiags_pv(itrc), itrc=1,NT)
         do itrc=1,NT
           if (wrtdiags_pv(itrc)) wrtdiags_pv(3)=.true.
        enddo
# endif /* TRACERS  */

# ifdef AVERAGES
!!         wrtdiags_pv_avg=.true.
      elseif (keyword(1:kwlen).eq.'diags_pv_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
        read(input,*,err=95) (wrtdiags_pv_avg(itrc),itrc=1,NT)
        do itrc=1,NT
          if (wrtdiags_pv_avg(itrc)) wrtdiags_pv_avg(3)=.true.
       enddo
#  endif /* TRACERS  */
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_PV */


#if defined DIAGNOSTICS_EDDY && ! defined XIOS
!
! Switches for fields to be saved into eddy diagnostics file.
!
!!          wrtdiags_eddy=.true.
      elseif (keyword(1:kwlen).eq.'diags_eddy_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  wrtdiags_eddy(1)
         do itrc=1,1
           if (wrtdiags_eddy(itrc)) wrtdiags_eddy(3)=.true.
        enddo

# ifdef AVERAGES
!!         wrtdiags_eddy_avg=.true.
      elseif (keyword(1:kwlen).eq.'diags_eddy_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtdiags_eddy_avg(1)
        do itrc=1,1
          if (wrtdiags_eddy_avg(itrc)) wrtdiags_eddy_avg(3)=.true.
       enddo
# endif /* AVERAGES  */
#endif /*DIAGNOSTICS_EDDY */


#if defined OUTPUTS_SURFACE && ! defined XIOS
!
! Switches for fields to be saved into surface outputs file.
!
!!          wrtsurf=.true.
      elseif (keyword(1:kwlen).eq.'surf_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  (wrtsurf(itrc), itrc=1,1)
         do itrc=1,1
           if (wrtsurf(itrc)) wrtsurf(3)=.true.
        enddo

# ifdef AVERAGES
!!         wrtsurf_avg=.true.
      elseif (keyword(1:kwlen).eq.'surf_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtsurf_avg(itrc),itrc=1,1)
        do itrc=1,1
          if (wrtsurf_avg(itrc)) wrtsurf_avg(3)=.true.
       enddo
# endif /* AVERAGES  */

#endif /*OUTPUTS_SURFACE*/

#ifdef DIAGNOSTICS_BIO
!
! Switches for fields to be saved into biology diagnostics file.
!
!          wrtdiabio=.true.
      elseif (keyword(1:kwlen).eq.'diagbioFlux_history_fields') then
!
        call cancel_kwd (keyword(1:kwlen), ierr)
!     ==
        read(input,*,err=95) (wrtdiabioFlux(iflux),iflux=1,NumFluxTerms)
        do iflux=1,NumFluxTerms
          if (wrtdiabioFlux(iflux))
     &         wrtdiabioFlux(NumFluxTerms+1)=.true.
        enddo
!     ==
      elseif (keyword(1:kwlen).eq.'diagbioVSink_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiabioVSink(iflux),
     &       iflux=1,NumVSinkTerms)
        do iflux=1,NumVSinkTerms
          if (wrtdiabioVSink(iflux))
     &         wrtdiabioVSink(NumVSinkTerms+1)=.true.
        enddo
!     ==
# if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      elseif (keyword(1:kwlen).eq.'diagbioGasExc_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiabioGasExc(iflux),
     &                                       iflux=1,NumGasExcTerms)
        do iflux=1,NumGasExcTerms
          if (wrtdiabioGasExc(iflux))
     &         wrtdiabioGasExc(NumGasExcTerms+1)=.true.
        enddo
# endif
! === === ===
# ifdef AVERAGES
!         wrtdiabio_avg=.true.
!        write(*,*)'=============='
      elseif (keyword(1:kwlen).eq.'diagbioFlux_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
!     ==
        read(input,*,err=95) (wrtdiabioFlux_avg(iflux),
     &       iflux=1,NumFluxTerms)
        do iflux=1,NumFluxTerms
          if (wrtdiabioFlux_avg(iflux))
     &         wrtdiabioFlux_avg(NumFluxTerms+1)=.true.
        enddo
!     ==
      elseif (keyword(1:kwlen).eq.'diagbioVSink_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiabioVSink_avg(iflux),
     &                              iflux=1,NumVSinkTerms)
        do iflux=1,NumVSinkTerms
          if (wrtdiabioVSink_avg(iflux))
     &         wrtdiabioVSink_avg(NumVSinkTerms+1)=.true.
        enddo
!     ==
#  if (defined BIO_NChlPZD && defined OXYGEN) || defined BIO_BioEBUS
      elseif (keyword(1:kwlen).eq.'diagbioGasExc_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtdiabioGasExc_avg(iflux),
     &                                  iflux=1,NumGasExcTerms)
        do iflux=1,NumGasExcTerms
          if (wrtdiabioGasExc_avg(iflux))
     &         wrtdiabioGasExc_avg(NumGasExcTerms+1)=.true.
        enddo
#  endif
# endif  /* AVERAGES */
#endif /* DIAGNOSTICS_BIO */

#ifdef STATIONS
!
! Switches for fields to be saved into stations output file.
!
      elseif (keyword(1:kwlen).eq.'station_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# ifdef AGRIF
        if (Agrif_Root()) then
# endif
          read(input,*,err=95) wrtsta(indxstaGrd),
     &       wrtsta(indxstaTemp), wrtsta(indxstaSalt),
     &       wrtsta(indxstaRho), wrtsta(indxstaVel)
# ifdef AGRIF
        endif
# endif
#endif /* STATIONS */

!
! Boussinesq Approximation mean density.
!
      elseif (keyword(1:kwlen).eq.'rho0') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rho0
#if defined UV_VIS2 || defined UV_VIS4 || defined SPONGE_VIS2
!
! Horizontal viscosity coefficients.
!
      elseif (keyword(1:kwlen).eq.'lateral_visc') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) visc2, visc4
#endif

!
! Bottom drag coefficients.
!
      elseif (keyword(1:kwlen).eq.'bottom_drag') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rdrg, rdrg2, Zobt, Cdb_min, Cdb_max
! Lateral boundary slipperness.
!
      elseif (keyword(1:kwlen).eq.'gamma2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) gamma2
#ifdef SOLVE3D
# if defined TS_DIF2 || defined SPONGE_DIF2
!
! Horizontal Laplacian mixing coefficients for tracers.
!
      elseif (keyword(1:kwlen).eq.'tracer_diff2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
        read(input,*,err=95) (tnu2(itrc),itrc=1,NT)
#  endif /* TRACERS */
# endif
# ifdef TS_DIF4
!
! Horizontal biharmonic mixing coefficients for tracer.
!
      elseif (keyword(1:kwlen).eq.'tracer_diff4') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
        read(input,*,err=95) (tnu4(itrc),itrc=1,NT)
#  endif /* TRACERS */

# endif
# if !defined LMD_MIXING && !defined GLS_MIXING
!
! Background vertical viscosity and mixing coefficients for tracers.
!
      elseif (keyword(1:kwlen).eq.'vertical_mixing') then
        call cancel_kwd (keyword(1:kwlen), ierr)
#  ifdef TRACERS
         read(input,*,err=95) Akv_bak,(Akt_bak(itrc),itrc=1,NT)
#  endif /* TRACERS */
# endif
# ifdef BODYFORCE
      elseif (keyword(1:kwlen).eq.'bodyforce') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) levsfrc,levbfrc
        if (levsfrc.lt.1 .or. levsfrc.gt.N) then
          MPI_master_only write(stdout,19) 'LEVSFRC = ',levsfrc
  19      format(' READ_INP - Illegal bodyforce level, ',A,i4)
          ierr=ierr+1
        endif
        if (levbfrc.lt.1 .or. levbfrc.gt.N) then
          MPI_master_only write(stdout,19) 'LEVBFRC = ',levbfrc
          ierr=ierr+1
        endif

# endif
#endif /*SOLVE3D*/

#if (defined SPONGE && !defined SPONGE_GRID)
!
! Parameters for sponge layers
!
!  if SPONGE_GRID is defined, they are set generically
!  in set_nudgcof routine
!
      elseif (keyword(1:kwlen).eq.'sponge') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) x_sponge
!
#elif (defined SPONGE && defined SPONGE_GRID)
!
!  if SPONGE_GRID is defined, they are set generically
!  in set_nudgcof routine
!
      elseif (keyword(1:kwlen).eq.'sponge') then
         call cancel_kwd (keyword(1:kwlen), ierr)
#endif
!
! Nudging parameters for restoring toward large-scale data in ABL
!
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_TRA
      elseif (keyword(1:kwlen).eq.'abl_nudg_tra_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ltra_min,ltra_max
        ltra_min = 1./ltra_min
        ltra_max = 1./ltra_max
#endif
#if defined ABL1D && defined ABL_NUDGING && defined ABL_NUDGING_DYN
      elseif (keyword(1:kwlen).eq.'abl_nudg_dyn_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldyn_min,ldyn_max
        ldyn_min = 1./ldyn_min
        ldyn_max = 1./ldyn_max
#endif

#if defined T_FRC_BRY     || defined M2_FRC_BRY    || \
     defined M3_FRC_BRY    || defined Z_FRC_BRY     || \
     defined TCLIMATOLOGY  || defined M2CLIMATOLOGY || \
     defined M3CLIMATOLOGY || defined ZCLIMATOLOGY
!
! Nudging parameters for OBC and nudging layers
! (converted from [days] to [sec^-1]
!
      elseif (keyword(1:kwlen).eq.'nudg_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
# if defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
     !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
        if (Agrif_Root()) then
# endif
          read(input,*,err=95) tauT_in,tauT_out,tauM_in,tauM_out
          tauT_in =1./(tauT_in *86400.)
          tauT_out=1./(tauT_out*86400.)
          tauM_in =1./(tauM_in *86400.)
          tauM_out=1./(tauM_out*86400.)
# if defined AGRIF && !defined AGRIF_OBC_M2ORLANSKI && \
     !defined AGRIF_OBC_M3ORLANSKI && !defined AGRIF_OBC_TORLANSKI
        endif
# endif
#endif

#ifdef SOLVE3D
# if !defined NONLIN_EOS
!
! Parameters for linear equations of state.
!
      elseif (keyword(1:kwlen).eq.'lin_EOS_cff') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) R0, T0, S0, Tcoef, Scoef
# endif
# ifdef SEDIMENT
!
! Sediments input file name for USGS sediment model.
!
      elseif (keyword(1:kwlen).eq.'sediments') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) sedname

      elseif (keyword(1:kwlen).eq.'sediment_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                   (wrthis(itrc), itrc=indxSed-1,indxSed+NST+1
#  ifdef SUSPLOAD
     &                                       +2*NST
#  endif
#  ifdef BEDLOAD
     &                                       +2*NST
#  endif
#  if defined MIXED_BED || defined COHESIVE_BED
     &                                       + 1
#  endif
     &  )
#  ifdef AVERAGES
        do itrc=indxSed,indxSed+NST+1
#   ifdef SUSPLOAD
     &                                       +2*NST
#   endif
#   ifdef BEDLOAD
     &                                       +2*NST
#   endif
         wrtavg(itrc)=wrthis(itrc)
        enddo
#  endif

# elif defined MUSTANG
!
! Sediments input file name for Mustang.
!
      elseif (keyword(1:kwlen).eq.'sediments_mustang') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) sedname_must

# endif /* SEDIMENT */
#endif /* SOLVE3D */

#if defined SUBSTANCE
      elseif (keyword(1:kwlen).eq.'substance') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) subsfilename
#endif
#ifdef BBL
      elseif (keyword(1:kwlen).eq.'bbl_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                    (wrthis(itrc), itrc=indxBBL,indxBBL+5)

#endif /* BBL */

#if defined OBSTRUCTION
      elseif (keyword(1:kwlen).eq.'obstruction') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) obstname
        lstr=lenstr(obstname)
#endif /* OBSTRUCTION */

#ifdef MRL_WCI
!
! MRL wave-current interaction variables history fields
!
        elseif (keyword(1:kwlen).eq.'wci_history_fields') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrthis(indxSUP),  wrthis(indxUST2D)
     &                                        ,  wrthis(indxVST2D)
# ifdef SOLVE3D
     &         , wrthis(indxUST), wrthis(indxVST), wrthis(indxwST)
     &         , wrthis(indxAKB), wrthis(indxAKW), wrthis(indxKVF)
     &                        , wrthis(indxCALP), wrthis(indxKAPS)
# endif

# ifdef AVERAGES
!
! MRL wave-current interaction variables average fields
!
        elseif (keyword(1:kwlen).eq.'wci_average_fields') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrtavg(indxSUP),  wrtavg(indxUST2D)
     &                                        ,  wrtavg(indxVST2D)
#  ifdef SOLVE3D
     &         , wrtavg(indxUST), wrtavg(indxVST), wrtavg(indxwST)
     &         , wrtavg(indxAKB), wrtavg(indxAKW), wrtavg(indxKVF)
     &                        , wrtavg(indxCALP), wrtavg(indxKAPS)
#  endif

# endif /* AVERAGES */
#endif /* MRL_WCI */

#if defined MRL_WCI || defined OW_COUPLING
!
! Wave history fields
!
      elseif (keyword(1:kwlen).eq.'wave_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrthis(i), i=indxHRM,indxHRM+8)
# ifdef AVERAGES
!
! Wave average fields
!
      elseif (keyword(1:kwlen).eq.'wave_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) (wrtavg(i), i=indxHRM,indxHRM+8)
# endif
#endif /* WKB_WWAVE || OW_COUPLING || WAVE_OFFLINE */

#if defined WKB_WWAVE
!
! WKB primary waves and empirical breaking model parameters
!
      elseif (keyword(1:kwlen).eq.'wkb_wwave') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wkb_amp,wkb_ang,wkb_prd,wkb_tide,wkb_btg,
     &                       wkb_gam

# ifdef WAVE_ROLLER
!
! Svendsen's surface roller model parameters
!
      elseif (keyword(1:kwlen).eq.'wkb_roller') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wkb_rsb, wkb_roller
# endif
#endif /* WKB_WWAVE */

#ifdef WAVE_MAKER
!
! Offshore wavemaker parameters for wave-resolving simulations
!
      elseif (keyword(1:kwlen).eq.'wave_maker') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wmaker_amp,wmaker_prd,wmaker_dir,
     &                       wmaker_dsp,wmaker_fsp
#endif /* WAVE_MAKER */

#ifdef ANA_PSOURCE
!
! Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
! J-grid locations (Isrc,Jsrc), and logical switch for type of tracer
! to apply (Lsrc). Currently, the direction can be along XI-direction
! (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The mass sources are
! located at U- or V-points so the grid locations should range from
! 1 =< Isrc =< L  and  1 =< Jsrc =< M.
!
# ifdef PSOURCE_NCFILE
      elseif (keyword(1:kwlen).eq.'psource_ncfile') then
        call cancel_kwd (keyword(1:kwlen), ierr)
! river runoff file name. Check its availability.
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        qbarname=fname(1:lstr)
        read(input,*,err=95) Nsrc
        do is=1,Nsrc
#  ifdef PSOURCE_NCFILE_TS
          read(input,*,err=95) Isrc(is),Jsrc(is),Dsrc(is),qbardir(is),
     &                                    (Lsrc(is,itrc), itrc=1,NT)
#  else
          read(input,*,err=95) Isrc(is),Jsrc(is),Dsrc(is),qbardir(is),
     &                                     (Lsrc(is,itrc), itrc=1,NT),
#   ifdef TRACERS
     &                                     (Tsrc0(is,itrc), itrc=1,NT)
#   endif /* TRACERS */
#  endif /* PSOURCE_NCFILE_TS  */

#  ifdef MPI
          Isrc_mpi(is,mynode)=Isrc(is)-iminmpi+1
          Jsrc_mpi(is,mynode)=Jsrc(is)-jminmpi+1
#  endif

      enddo
# else
!<- PSOURCE_NCFILE
      elseif (keyword(1:kwlen).eq.'psource') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) Nsrc
        do is=1,Nsrc
          read(input,*,err=95) Isrc(is), Jsrc(is), Dsrc(is), Qbar(is)
#  if defined SOLVE3D && defined TRACERS
     &                                    ,(Lsrc(is,itrc), itrc=1,NT)
     &                                    ,(Tsrc0(is,itrc), itrc=1,NT)
#  endif

#  ifdef MPI
          Isrc_mpi(is,mynode)=Isrc(is)-iminmpi+1
          Jsrc_mpi(is,mynode)=Jsrc(is)-jminmpi+1
#  endif

        enddo
# endif   /*  define PSOURCE_NCFILE  */
#endif /* ANA_PSOURCE */
!
!
! Hydrothermal forcing file name. Check its availability.
!
#if defined BHFLUX || ( defined BWFLUX && defined SALINITY)
      elseif (keyword(1:kwlen).eq.'bottom_forcing') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
# if defined MPI && defined PARALLEL_FILES && defined MPI_NOLAND
        call insert_node (fname, lstr, mynode2, NNODES2, ierr)
# endif
# if defined MPI && defined PARALLEL_FILES && !defined MPI_NOLAND
        call insert_node (fname, lstr, mynode, NNODES, ierr)
# endif
        open (testunit, file=fname(1:lstr), status='old', err=97)
        close(testunit)
        btfname=fname(1:lstr)
#endif  /* BHFLUX || (SALINITY && BWFLUX )  */
!
! Online forcing
!
#ifdef ONLINE   /* JG Online */
      elseif (keyword(1:kwlen).eq.'online') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) yearnum, monthnum, recordsperday, yearend,
     &                       monthend
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
        pathbulk=fname(1:lstr)
#endif /* ONLINE */

      else
        MPI_master_only write(stdout,'(/3(1x,A)/)')
     &                  'WARNING: Unrecognized keyword:',
     &                   keyword(1:kwlen),' --> DISREGARDED.'
      endif
      if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
!
! Error while reading input parameters.
!
  95  write(stdout,'(/1x,4A/)') 'READ_INP ERROR while reading block',
     &                    ' with keyword ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  write(stdout,'(/1x,4A/)') 'READ_INP ERROR: Cannot find input ',
     &                                'file ''', fname(1:lstr), '''.'
#ifdef LOGFILE
  98  write(6,'(/1x,4A/)') 'READ_INP ERROR: Cannot open log file ',
     &                                'file   , croco.log', '''.'
#endif
      ierr=ierr+1
  99  close (input)
!
! Check that all keywords were canceled
! Complain if some of them are left
!
      if (ierr.eq.0) then
        call check_kwds (ierr)
!
! Check CPP-switches for consistency. This operation is split into
! two stages because the first subroutine, "check_switches1", is
! generated by special program cppcheck (file cppcheck.F) by
! examining and documentation all available switches in cppdefs.h.
! This subroutine creates log of all switches defined in "cppdefs.h",
! as well as traps multiply defined global configurations (project
! switches, such as REGIONAL, etc).
! The second routine, "check_switches2" is hand written and it
! contains traps for mutually exclussive definition of all other
! CPP-switches (i.e. those which are NOT project selection switches,
! for example, it traps multiply defined vertical mixing schemes or
! lateral boundary conditions).
!
! Both codes are written in transparent mode: they assumed that error
! variable (ierr) is initialized at entry and they add 1 for each
! error discovered.
!
        call check_srcs
        call check_switches1 (ierr)
        call check_switches2 (ierr)

      return
      end
