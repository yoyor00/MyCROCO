#include "cppdefs.h"
                              ! Read initial conditions for the
      subroutine das_get_initial  ! primitive variables from NetCDF 
                              ! initialization file.
      implicit none
#include "netcdf.inc"
      real time_scale
      character vnts*4
      integer ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,
     &      start(2), count(2), ibuff(4),   das_nf_fread, checkdims
#ifdef SOLVE3D
     &                                      , itrc
#endif
      character units*64
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "ncscrum.h"
#include "das_ocean.h"
!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! used, provided that nrrec is within the available records; 
! (3) if the file has multiple records and nrrec<0, then THE LAST 
! available record is used.
!
      lstr=lenstr(ininame)
      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      if (ierr .ne. nf_noerr) then
        write(stdout,'(/1x,2A/15x,3A)') 'GET_INITIAL ERROR: Cannot ',
     &               'open netCDF file', '''', ininame(1:lstr) ,'''.'
        goto 99                                           !--> ERROR
      endif
!
! Read in evolving model variables:
! ---- -- -------- ----- ----------
      indx=1
!
! Free-surface.
!
      ierr=nf_inq_varid (ncid, 'zeta_da', varid)
      if (ierr .eq. nf_noerr) then
        ierr=das_nf_fread (zeta_s(START_2D_ARRAY), ncid, varid,
     &                                         indx, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'zeta_da',indx,ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'zeta_da', ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in the XI-direction.
!
      ierr=nf_inq_varid (ncid, 'u_da', varid)
      if (ierr .eq. nf_noerr) then
        ierr=das_nf_fread (u_s(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, u3dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'u_da', indx,
     &                                     ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'u_da', ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! 3D momentum component in the ETA-direction.
!
      ierr=nf_inq_varid (ncid, 'v_da', varid)
      if (ierr .eq. nf_noerr) then
        ierr=das_nf_fread (v_s(START_2D_ARRAY,1), ncid, varid,
     &                                        indx, v3dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'v_da', indx,
     &                                       ininame(1:lstr)
          goto 99                                        !--> ERROR
        endif
      else
        write(stdout,1) 'v_da', ininame(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Tracer variables.
!
      do itrc=1,NT
        if (itrc .eq. 1) then
           vnts='t_da'
        else if (itrc .eq. 2) then
           vnts='s_da'
        else
           vnts='BIO-VARIABLES'
        endif

        ierr=nf_inq_varid (ncid, vnts, varid)
        if (ierr .eq. nf_noerr) then
          ierr=das_nf_fread (t_s(START_2D_ARRAY,1,itrc), ncid,  varid,
     &                                               indx, r3dvar)
          if (ierr .ne. nf_noerr) then
            write(stdout,2) vnts, indx,
     &                                              ininame(1:lstr)
            goto 99                                       !--> ERROR
          endif
        else
          write(stdout,1) vnts, ininame(1:lstr)
          goto 99                                         !--> ERROR
        endif
      enddo
!
! mask
!
      ierr=nf_inq_varid (ncid, 'mask_da', varid)
      if (ierr .eq. nf_noerr) then           
        ierr=das_nf_fread (rmask_das(START_2D_ARRAY,1), ncid,  varid,
     &                                               indx, r3dvar)
        if (ierr .ne. nf_noerr) then                    
          write(stdout,2) 'mask_da', indx,
     &                                              ininame(1:lstr)
          goto 99                                       !--> ERROR
        endif                                           
      else
        write(stdout,1) 'mask_da', ininame(1:lstr)
        goto 99                                         !--> ERROR
      endif

  1   format(/1x,'GET_INITIAL - unable to find variable:',    1x,A,
     &                            /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'GET_INITIAL - error while reading variable:',1x, A,
     &    2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
!
!  Close input NetCDF file.
!
      ierr=nf_close(ncid)
      return
  99  may_day_flag=2
      return
      end
