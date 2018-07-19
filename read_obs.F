c     this file is not differentiated

# include "cppdefs.h"

      subroutine read_obs()
      implicit none
# include "netcdf.inc"
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"

      character*(*) FILE_NAME
      parameter (FILE_NAME='obs.nc')

      integer ncid, varid, retval, nf_fread, k

      integer i,j

      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, 'zeta', varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      do k=ad_ast,ad_nt+1
         retval = nf_fread(ad_obs(START_2D_ARRAY,k-ad_ast+1),ncid,varid,
     &        k-ad_ast+1,r2dvar)
c         write (*,*) 'max_val ad_obs:',k,maxval(ad_obs)
      end do

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      write (*,*) 'max_val ad_obs:',maxval(ad_obs)

      
      end subroutine read_obs
