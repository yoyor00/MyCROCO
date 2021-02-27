!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://roms.mpl.ird.fr
!======================================================================
!
#include "cppdefs.h"
                                      ! Read gridded bottom_drag parameters
      subroutine get_bottom_drag_2d   ! from bottom_drag_2d netCDF file. 
                                      ! Scale them. 
      ! author : Marie-Lou Bachely and Cristele Chevalier IRD
#if defined BOTTOM_DRAG_2D
      implicit none
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "bottom_drag_2d.h"
      real cff
      integer i,j, ierr, lstr,lvar, lenstr, nf_fread,
     &        status, varid, ncid, checkdims
#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif    
#include "netcdf.inc"
!
      lstr=lenstr(dragname)
      ierr=nf_open(dragname(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, dragname, lstr, varid)
        if (ierr. ne. nf_noerr) goto 99
      else
        write(stdout,'(/3(1x,A)/)') 'GET_BOTTOM_DRAG_2D ERROR:',
     &                       ' Cannot open ',
     &                       'input NetCDF file:', dragname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! bottom drag coefficient :
!
      ierr=nf_inq_varid (ncid, 'bottom_drag_coef', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread(RDRG2_2D (START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'bottom_drag_coef', dragname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'bottom_drag_coef', dragname(1:lstr)
        goto 99                                           !--> ERROR
      endif

      return

!
! Sort out error messages: The following portion of the code is 
!===== === ===== ========= not accessed unless something goes wrong.
!
      goto 99

  1   format(/1x, 'GET_BOTTOM_DRAG_2D - unable to find bottom_drag', 
     &                                                  'variable:',
     &           1x, A / 12x, 'in bottom_drag netCDF file:', 1x, A/)
  2   format(/1x, 'GET_BOTTOM_DRAG_2D - error while reading variable:',
     &           1x, A / 12x, 'in bottom_drag netCDF file:', 1x, A/)

 99   may_day_flag=2

#endif /* BOTTOM_DRAG_2D */
      return  
      end
