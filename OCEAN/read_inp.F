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


! auxiliary_history_fields and auxiliary_averages have been migrated to
! croco.nml namelist blocks. See croco_namelist.F90 / croco_namelist_init.F90.




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

      if (ierr.eq.0) then
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


