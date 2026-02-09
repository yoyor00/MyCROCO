#include "cppdefs.h"
#if defined DAS_GOES_SST

      subroutine das_get_goes

      implicit none

# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"

      integer status, nc_id, v_id, i, j
     
#include "netcdf.inc"
!
! READ goes sst
!
      flag_goes=.true.

      status = NF_OPEN(file_goes,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND OR OPENED',
     &          file_goes
        flag_goes=.false.
        goto 999
      endif
!
! sst
!
      status = NF_INQ_VARID(nc_id, 'sst',v_id)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring goes sst'
        flag_goes=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading goes sst'
        flag_goes=.false.
        goto 999
      endif
      do j=0,Mm+1
        do i=0,Lm+1
          sst_goes(i,j)=anc(i+1,j+1)
        enddo
      enddo
!
! mask
!
      do j=0,Mm+1
        do i=0,Lm+1
          if (sst_goes(i,j) .gt. 0.5 ) then
            goes_mask(i,j)=1.0
          else
            goes_mask(i,j)=0.0
          endif
        enddo
      enddo

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_goes) then
        write(*, '(6x, A,(/8x,A))') 'OK: GOES SST READ IN',
     &         file_goes
      endif
    
#else
      subroutine das_get_goes_empty
      return
      end
#endif /* DAS_GOES_SST */

