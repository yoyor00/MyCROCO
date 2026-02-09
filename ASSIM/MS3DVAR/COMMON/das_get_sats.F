#include "cppdefs.h"
#if defined DAS_SATSSSS

      subroutine das_get_sats

      implicit none

# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"

      integer status, nc_id, v_id, i, j
     
#include "netcdf.inc"
!
! READ sats sss
!
      flag_sats=.true.

      status = NF_OPEN(file_sats,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND',
     &          file_sats
        flag_sats=.false.
        goto 999
      endif
!
! sss
!
      status = NF_INQ_VARID(nc_id, 'sss',v_id)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring sats sss'
        flag_sats=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading sats sss'
        flag_sats=.false.
        goto 999
      endif
      do j=0,Mm+1
        do i=0,Lm+1
          sss_sats(i,j)=anc(i+1,j+1)
        enddo
      enddo
!
! mask
!
      do j=0,Mm+1
        do i=0,Lm+1
          if (sss_sats(i,j) .gt. -900) then
            sats_mask(i,j)=1.0
          else
            sats_mask(i,j)=0.0
          endif
        enddo
      enddo

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_sats) then
        write(*, '(6x, A,(/8x,A))') 'OK: SATS SSS READ IN',
     &         file_sats
      endif
    
      return
      end
#else
      subroutine das_get_sats_empty
      return
      end
#endif /* DAS_SATSSSS  */
