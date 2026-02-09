#include "cppdefs.h"
      subroutine das_get_sss

      implicit none

# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"

      integer status, nc_id, v_id, i, j
     
      include 'netcdf.inc'
!
! READ sss obs
!
      flag_sss=.true.

      status = NF_OPEN(file_sss,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: GOES NOT FOUND',
     &          file_sss
        flag_sss=.false.
        goto 999
      endif
!
! ssh
!
      status = NF_INQ_VARID(nc_id, 'ssh',v_id)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring sss'
        flag_sss=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading sss'
        flag_sss=.false.
        goto 999
      endif
      do j=0,Mm+1
        do i=0,Lm+1
          (i,j)=anc(i+1,j+1)
        enddo
      enddo
!
! mask
!
      do j=0,Mm+1
        do i=0,Lm+1
          if (ssh_sss(i,j) .gt. -900 
     &        .and. h(i,j) .gt. hcmin) then
            sss_mask(i,j)=1.0
          else
            sss_mask(i,j)=0.0
          endif
        enddo
      enddo

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_sss) then
        write(*, '(6x, A,(/8x,A))') 'OK: SSS READ IN',
     &         file_sss
      endif
    
      return
      end

