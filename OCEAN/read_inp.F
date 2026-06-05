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

# ifdef ENSEMBLE
      use ensmpi , only : ln_ens_rst_in ! do we start from an ensemble of restart
# endif
# ifdef DIAGNOSTICS_TS_MLD_DENS
      use croco_namelist, only: mld_crit_D, mld_crit_T
# endif

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
#ifdef MPI
      include 'mpif.h'
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
#elif defined SINGLE_COLUMN && defined KATO_PHILLIPS
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
# elif defined SED_TOY_FLOC_0D
      fname='TEST_CASES/croco.in.Sed_toy_floc_0d'
# elif defined SED_TOY_FLOC_1D
      fname='TEST_CASES/croco.in.Sed_toy_floc_1d'
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
      if (iargc().ge.1) call getarg(1,fname)
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
      if (1 == 0) then ! temporary, just to keep if first, TODO remove it when all is in nml
          ! never pass here, it is just to keep first if
          ! (to be remove when all is in nml)


!
! Switches for fields to be saved into history file.
!
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

#endif /* SOLVE3D */

#ifdef AVERAGES
# ifdef SOLVE3D
!
! Switches for fields to be saved into averages file.
!
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

# endif /* SOLVE3D */
#endif /* AVERAGES */




!









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
          read(input,*,err=95) Isrc(is),Jsrc(is),Dsrc(is),qbardir(is)
#   ifdef TRACERS
     &                                     ,(Lsrc(is,itrc), itrc=1,NT)
     &                                     ,(Tsrc0(is,itrc), itrc=1,NT)
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
!
! Online forcing
!

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
      endif

      ! start of summary writing
      if (ierr.eq.0) then




          MPI_master_only write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &    , wrthis(indxZ),  'write zeta ', 'free-surface.'
     &    , wrthis(indxUb), 'write UBAR ', '2D U-momentum component.'
     &    , wrthis(indxVb), 'write VBAR ', '2D V-momentum component.'
#ifdef SOLVE3D
     &    , wrthis(indxU),  'write U    ', '3D U-momentum component.'
     &    , wrthis(indxV),  'write V    ', '3D V-momentum component.'
# ifdef TRACERS
          do itrc=1,NT
          MPI_master_only write(stdout, '(6x,L1,2x,A,I2,A,I2,A)')
     &                     wrthis(indxV+itrc), 'write T(', itrc,
     &                              ')  Tracer of index', itrc,'.'
          enddo
# endif
#endif  /* SOLVE3D */

#if defined SOLVE3D
          MPI_master_only write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrthis(indxR),    'write RHO  ', 'Density anomaly.'
     &  , wrthis(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrthis(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrthis(indxAkv),  'write Akv  ', 'Vertical viscosity.'
# ifdef TEMPERATURE
     &  , wrthis(indxAkt),  'write Akt  ',
     &                      'Vertical diffusivity for temperature.'
# endif
# ifdef SALINITY
     &  , wrthis(indxAks),  'write Aks  ',
     &                      'Vertical diffusivity for salinity.'
# endif
# if defined ANA_VMIX || defined LMD_MIXING \
     || defined LMD_SKPP || defined LMD_BKPP \
     || defined GLS_MIXING
     &  , wrthis(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
# endif
# ifdef VIS_COEF_3D
     &  , wrthis(indxVisc),  'write Visc3d', 'Horizontal viscosity.'
# endif
# ifdef DIF_COEF_3D
     &  , wrthis(indxDiff),  'write Visc3d', 'Horizontal diffusivity.'
# endif
# if defined LMD_SKPP || defined GLS_MIXING
     &  , wrthis(indxHbl),  'write Hbl  ',
     &                      'Depth of model boundary layer.'
# endif
# ifdef LMD_BKPP
     &  , wrthis(indxHbbl), 'write Hbbl  ',
     &                      'Depth of bottom planetary boundary layer.'
# endif
# ifdef BULK_FLUX
     &  , wrthis(indxShflx_rlw), 'write shflx_rlw [W/m2]',
     &                                 'Long Wave heat flux.'
     &  , wrthis(indxShflx_lat), 'write shflx_lat [W/m2]',
     &                                 'Latent heat flux.'
     &  , wrthis(indxShflx_sen), 'write shflx_sen [W/m2]',
     &                                 'Sensible heat flux'
# endif
# ifdef MORPHODYN
     &  , wrthis(indxHm),   'write Hm [m]',
     &                                 'Time evolving bathymetry'
# endif
# if defined BIOLOGY && !defined PISCES
     &  , wrthis(indxHel),   'write Hel  ',
     &                       'Depth of the euphotic layer'
#  ifdef BIO_NChlPZD
     &  , wrthis(indxChC),   'write ChC  ',
     &                       'Chlorophyll to Carbon ratio'
#   ifdef OXYGEN
     &  , wrthis(indxU10),   'write u10 ',
     &                       'Wind speed at 10 m height'
     &  , wrthis(indxKvO2),  'write Kv_O2 ',
     &                       'Gas transfer coefficient for O2'
     &  , wrthis(indxO2sat), 'write O2sat ',
     &                       'Saturation concentration of O2'
#   endif
#  elif defined BIO_BioEBUS
     &  , wrthis(indxAOU), 'write AOU '
     &  ,                  'Apparent Oxygen Utilization.'
     &  , wrthis(indxWIND10), 'write wind10 ', 'wind speed at 10 m.'
#  endif
# endif /* BIOLOGY */
     &  , wrthis(indxBostr), 'write Bostr', 'Bottom Stress.'
     &  , wrthis(indxBustr), 'write Bustr', 'U-Bottom Stress.'
     &  , wrthis(indxBvstr), 'write Bvstr', 'V-Bottom Stress.'
     &  , wrthis(indxWstr),  'write Wstress', 'Wind Stress.'
     &  , wrthis(indxUWstr), 'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrthis(indxVWstr), 'write V-Wstress comp.', 'V-Wind Stress.'
# ifdef TEMPERATURE
     &  , wrthis(indxShflx), 'write Shflx [W/m2]',
     &                       'Surface net heat flux'
# endif
# ifdef SALINITY
     &  , wrthis(indxSwflx), 'write Swflx [cm/day]',
     &                       'Surface freshwater flux (E-P)'
# endif
# ifdef TEMPERATURE
     &  , wrthis(indxShflx_rsw),'write Shflx_rsw [W/m2]',
     &                          'Short-wave surface radiation'
# endif
# ifdef BHFLUX
     &  , wrthis(indxBhflx), 'write Bhflx [W/m2]',
     &                       'Bottom net heat flux'
# endif
# if defined BWFLUX && defined SALINITY
     &  , wrthis(indxBwflx), 'write Bwflx [cm/day]',
     &                       'Bottom freshwater flux (E-P)'
# endif

#endif /* SOLVE3D */

#ifdef AVERAGES
          MPI_master_only write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &  'Fields to be saved in averages file: (T/F)'
     &  , wrtavg(indxZ),  'write zeta ', 'free-surface.'
     &  , wrtavg(indxUb), 'write UBAR ', '2D U-momentum component.'
     &  , wrtavg(indxVb), 'write VBAR ', '2D V-momentum component.'
# ifdef SOLVE3D
     &  , wrtavg(indxU),  'write U    ', '3D U-momentum component.'
     &  , wrtavg(indxV),  'write V    ', '3D V-momentum component.'
#  ifdef TRACERS
          do itrc=1,NT
          MPI_master_only write(stdout,
     &                     '(6x,L1,2x,A,I2,A,2x,A,I2,A)')
     &                      wrtavg(indxV+itrc), 'write T(',
     &                      itrc,')', 'Tracer of index', itrc,'.'
          enddo
#  endif

          MPI_master_only write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrtavg(indxR),    'write RHO  ', 'Density anomaly'
     &  , wrtavg(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrtavg(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrtavg(indxAkv),  'write Akv  ', 'Vertical viscosity'
#  ifdef TEMPERATURE
     &  , wrtavg(indxAkt),  'write Akt  ',
     &                      'Vertical diffusivity for temperature.'
#  endif
#  ifdef SALINITY
     &  , wrtavg(indxAks),  'write Aks  ',
     &                         'Vertical diffusivity for salinity.'
#  endif
#  if defined ANA_VMIX || defined LMD_MIXING \
     || defined LMD_SKPP || defined LMD_BKPP \
     || defined GLS_MIXING
     &  , wrtavg(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
#  endif
#  ifdef VIS_COEF_3D
     &  , wrtavg(indxVisc),'write visc3d', 'Horizontal viscosity'
#  endif
#  ifdef DIF_COEF_3D
     &  , wrtavg(indxDiff),'write diff3d', 'Horizontal diffusivity'
#  endif
#  if defined LMD_SKPP || defined GLS_MIXING
     &  , wrtavg(indxHbl),  'write Hbl  ',
     &                          'Depth of model boundary layer'
#  endif
#  ifdef LMD_BKPP
     &  , wrtavg(indxHbbl),  'write Hbbl  ',
     &                    'Depth of the bottom planetary boundary layer'
#  endif
#  ifdef BULK_FLUX
     &  , wrtavg(indxShflx_rlw), 'write shflx_rlw [W/m2]',
     &                                 'Long Wave heat flux.'
     &  , wrtavg(indxShflx_lat), 'write shflx_lat[W/m2] ',
     &                                 'Latente heat flux.'
     &  , wrtavg(indxShflx_sen), 'write shflx_sen [W/m2]',
     &                                 'Sensible heat flux.'
#  endif
#  ifdef MORPHODYN
     &  , wrtavg(indxHm), 'write Hm [m]',
     &                                 'Time evolving bathymetry.'
#  endif
#  if defined BIOLOGY && !defined PISCES
     &  , wrtavg(indxHel),'write Hel  ',
     &                          'Depth of the euphotic layer'
#   ifdef BIO_NChlPZD
     &  , wrtavg(indxChC),'write ChC  ',
     &                          'Chlorophyll to Carbon ratio'
#    ifdef OXYGEN
     &  , wrtavg(indxU10),'write u10 ',
     &         'Wind speed at 10 m height'
     &  , wrtavg(indxKvO2),'write Kv_O2 ',
     &         'Gas transfer coefficient for O2'
     &  , wrtavg(indxO2sat),'write O2sat ',
     &         'Saturation concentration of O2'
#    endif

#   elif defined BIO_BioEBUS
     &  , wrtavg(indxAOU), 'write AOU '
     &  ,                      'Apparent Oxygen Utilization.'
     &  , wrtavg(indxWIND10), 'write wind10 ', 'wind speed at 10 m.'
#   endif
#  endif   /* BIOLOGY */
     &  , wrtavg(indxBostr),'write Bostr', 'Bottom Stress.'
     &  , wrtavg(indxBustr),'write Bustr', 'U-Bottom Stress.'
     &  , wrtavg(indxBvstr),'write Bvstr', 'V-Bottom Stress.'
     &  , wrtavg(indxWstr), 'write Wstr', 'Wind Stress.'
     &  , wrtavg(indxUWstr),'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrtavg(indxVWstr),'write V-Wstress comp.', 'V-Wind Stress.'
#  ifdef TEMPERATURE
     &  , wrtavg(indxShflx),'write Shflx [W/m2]',
     &                      'Surface net heat flux.'
#  endif
#  ifdef SALINITY
     &  , wrtavg(indxSwflx),'write Swflx [cm/day]',
     &                      'Surface freshwater flux (E-P)'
#  endif
#  ifdef TEMPERATURE
     &  , wrtavg(indxShflx_rsw),'write Shflx_rsw [W/m2]',
     &                      'Short-wave surface radiation.'
#  endif
#  ifdef BHFLUX
     &  , wrtavg(indxBhflx), 'write Bhflx [W/m2]',
     &       'Bottom net heat flux'
#  endif
#  if defined BWFLUX && defined SALINITY
     &  , wrtavg(indxBwflx), 'write Bwflx [cm/day]',
     &                       'Surface freshwater flux (E-P)'
#  endif
# endif /* SOLVE3D */
#endif /* AVERAGES */

#ifdef ANA_PSOURCE
# ifdef PSOURCE_NCFILE
          lstr = lenstr(qbarname)
          MPI_master_only write(stdout,'(2x,A,2x,A)')
     &             'Runoff Data File:', qbarname(1:lstr)

          MPI_master_only write(stdout,'(/6x,i6,2x,A,1x,A)')
     &                               Nsrc, 'Number of point sources'
          do is=1,Nsrc

          MPI_master_only write(stdout,'(3(/6x,i6,2x,A),(/6x,F6.0,2x,A))')
     &    Isrc(is), 'I point source indice'
     &  , Jsrc(is), 'J point source indice'
     &  , Dsrc(is), 'Orientation of point source flow'
     &  , qbardir(is), 'Direction of point source flow'

#  ifdef TRACERS
            do itrc=1,NT
              MPI_master_only write(stdout,
     &         '(6x,L1,2x,A,I2,A,2x,A,I2,A)')
     &         Lsrc(is,itrc), 'write Lsrc(',
     &         itrc,')', 'Tracer of index', itrc,'.'
            enddo
#  endif /* TRACERS */
#  ifndef PSOURCE_NCFILE_TS
#   ifdef TRACERS
            do itrc=1,NT
              MPI_master_only write(stdout,
     &         '(6x,1pe10.3,2x,A,I2,A,2x,A,I2,A)')
     &         Tsrc0(is,itrc), 'write Tsrc(',
     &         itrc,')', 'Tracer of index', itrc,'.'
            enddo
#   endif /* TRACERS */
#  endif
          enddo ! is
# else
!<- PSOURCE_NCFILE
          MPI_master_only write(stdout,'(/6x,i6,2x,A,1x,A)')
     &                               Nsrc, 'Number of point sources'
          do is=1,Nsrc
          MPI_master_only write(stdout,'(3(/6x,i6,2x,A))')
     &    Isrc(is), 'I point source indice'
     &  , Jsrc(is), 'J point source indice'
     &  , Dsrc(is), 'Direction of point source flow'
          MPI_master_only write(stdout,'(1pe10.3,2x,A)')
     &    Qbar(is), 'Total transport at point source'
#  ifdef MPI
          MPI_master_only write(*,*)'Isrc(is)=',Isrc(is)
          MPI_master_only write(*,*)'Jsrc(is)=',Jsrc(is)
#  endif
#  if defined SOLVE3D && defined TRACERS
          do itrc=1,NT
               MPI_master_only write(stdout,
     &                     '(6x,L1,2x,A,I2,A,2x,A,I2,A)')
     &                      Lsrc(is,itrc), 'write Lsrc(',
     &                      itrc,')', 'Tracer of index', itrc,'.'
          enddo
          do itrc=1,NT
               MPI_master_only write(stdout,
     &                     '(6x,1pe10.3,2x,A,I2,A,2x,A,I2,A)')
     &                      Tsrc0(is,itrc), 'write Tsrc(',
     &                      itrc,')', 'Tracer of index', itrc,'.'
          enddo
#  endif
          enddo
# endif   /*  define PSOURCE_NCFILE  */
#endif /* ANA_PSOURCE */




       
      endif ! ierr.eq.0  
      ! end of summary writing


      if (ierr.ne.0) then
        write(stdout,'(/1x,2A,I3,1x,A/)') 'READ_INP ERROR: ',
     & 'A total of', ierr, 'configuration errors discovered.'
        return
      endif
#ifdef MPI
!      call MPI_Barrier (MPI_COMM_WORLD, ierr)
#endif
      return
      end


