#include "cppdefs.h"
#if defined DAS_JASONSSH || defined DAS_TPSSH\
  || defined DAS_ERS2SSH\
  || defined DAS_GFOSSH || defined DAS_SWOTSSH

      subroutine das_get_ssh_ref
      implicit none
!
# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"
# include "das_ocean.h"

      integer status, nc_id, v_id, i, j

#include "netcdf.inc"
!
! READ ssh_ref ssh
!
      flag_ssh_ref=.true.

      status = NF_OPEN(file_ssh_ref,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: SSH REF NOT FOUND',
     &          file_ssh_ref

        do j=0,Mm+1
          do i=0,Lm+1
             ssh_ref(i,j)= 0.0
             ref_mask(i,j)= rmask_das(i,j,ndas)
          enddo
        enddo
!
!NOTE: ref_mask always applied to innovation, equivalent to zero bias at every gridpoint
!
        flag_ssh_ref=.false.
        goto 999
      endif
!
! ssh
!
      status = NF_INQ_VARID(nc_id, 'ssh',v_id)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring ssh ref'
        flag_swot=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading ssh ref'
        goto 999
      endif
      do j=0,Mm+1
        do i=0,Lm+1
          ssh_ref(i,j)=anc(i+1,j+1)
        enddo
      enddo
!
! mask
!
      do j=0,Mm+1
        do i=0,Lm+1
          if (ssh_ref(i,j) .gt. -900. ) then
           ref_mask(i,j)=rmask_das(i,j,ndas)
          else
            ref_mask(i,j)=0.0
          endif
        enddo
      enddo
      do j=0,Mm+1
        do i=0,Lm+1
           ssh_ref(i,j)= ssh_ref(i,j)*ref_mask(i,j)
        enddo
      enddo

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_ssh_ref) then
        write(*, '(6x, A,(/8x,A))') 'OK: SSH REF READ IN',
     &         file_ssh_ref
      endif


!
      return
      end
#else
      subroutine das_get_ssh_ref_empty
      return
      end
#endif
