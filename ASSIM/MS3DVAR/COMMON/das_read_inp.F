#include "cppdefs.h"
                                        ! Read/report model input
      subroutine das_read_inp           ! parameters  from  startup 
                                        ! script file using keywords
      implicit none                     ! to recognize variables.
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "das_innov.h"
#include "das_covar.h"

      character fname*64, cfile*120
#ifdef MPI
      include 'mpif.h'
#endif
      integer ios, iargc, lstr, lenstr, input
      parameter (input=15)
!
#ifdef MPI
      if (mynode.eq.0 .and. iargc().eq.2) call getarg(2,fname) 
      call MPI_Bcast(fname,64,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#else
      if (iargc().eq.2) call getarg(2,fname)
#endif
! ==== == ======== ==== ======= ===== ======= == ======
!
      ios=0    ! <-- reset error counrter
      open(input, file=fname,status='old',
     &                  form='formatted', iostat=ios)
      if(ios.ne.0) then
        write(*, '(6x, A, /6x A)')
     &      '!!! ERROR: can not find the file namelist',
     &       fname
        stop 999
      else
        write(*, '(6x, A)') fname
      endif
!      read (input, *) time_level
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_tmi=cfile(1:lstr)

      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_mc=cfile(1:lstr)

      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_goes=cfile(1:lstr)

      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_swot=cfile(1:lstr)
 
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_js1=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_tp=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_sats=cfile(1:lstr)
!
! in situ
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_prf_TS_obs=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_sio_glider=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_whoi_glider=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_ptsur_ctd=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_martn_ctd=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_map_ctd=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_prof_ctd=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_flight_sst=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_hfradar_uv=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_hfradar_uv6=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_ship_sst=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_cal_auv=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_dor_auv=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_innov_sio=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_innov_whoi=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_innov_ana_sio=cfile(1:lstr)
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_innov_ana_whoi=cfile(1:lstr)
!
!ifdef DAS_BVAR_CORR
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_bvar=cfile(1:lstr)

      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_corr=cfile(1:lstr)

      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_corr_vert=cfile(1:lstr)
!
!endif
!
! distance to coast
!
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_dist_coast=cfile(1:lstr)
!
#if defined DAS_JASONSSH || defined DAS_TPSSH\
  || defined DAS_TPSSH || defined DAS_ERS2SSH\
  || defined DAS_GFOSSH || defined DAS_SWOTSSH
      read(input, '(A)') cfile
      lstr=lenstr(cfile)
      file_ssh_ref=cfile(1:lstr)
#endif
!
      close(input)

      write(*,*)
      write(*, '(6x, A,5(/6x,A))') 'OBSERVATION FILES:',
     &         file_tmi,file_mc,file_goes,file_swot,
     &         file_js1,file_tp,file_sats
      write(*,*)
      write(*, '(6x, A,18(/6x,A))') 'IN-SITU OBSERVATION FILES:',
     &         file_prf_TS_obs,
     &         file_sio_glider,file_whoi_glider,file_flight_sst,
     &         file_hfradar_uv,file_hfradar_uv6,file_ship_sst,
     &         file_ptsur_ctd,file_martn_ctd,file_map_ctd,
     &         file_prof_ctd,
     &         file_cal_auv,file_dor_auv,file_innov_sio,
     &         file_innov_whoi, file_innov_ana_sio,
     &         file_innov_ana_whoi

      write(*,*)
      write(*, '(6x, A, (/6x,A))') 'MODEL ERROR RMS FILES:',
     &         file_bvar
#ifdef DAS_BVAR_CORR
      write(*,*)
      write(*, '(6x, A,2(/6x,A))') 'CORRELATION FILES:',
     &         file_corr, file_corr_vert
      write(*,*)
#else
      write(*,'(6x, A)') 'NO CORRELATION IN DATA ASSIMILATION'
      write(*,*)
#endif
      write(*, '(6x, A, (/6x,A))') 'MODEL DIST-COAST FILES:',
     &         file_dist_coast
      write(*,*)

#if defined DAS_JASONSSH || defined DAS_TPSSH\
  || defined DAS_TPSSH || defined DAS_ERS2SSH\
  || defined DAS_GFOSSH || defined DAS_SWOTSSH
      write(*,'(6x, A, (/6x,A))') 'ZETA REFERENCE FILE:',
     &         file_ssh_ref
      write(*,*)
#endif
!
      return
      end

