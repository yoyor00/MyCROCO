#include "cppdefs.h"
#if defined DAS_SWOTSSH 

      subroutine das_get_swot

      implicit none

# include "param.h"
# include "das_param.h"
# include "grid.h"
# include "das_innov.h"
#if defined DAS_SWOTSSHMSK 
# include "das_buffer.h"
#endif


! BLXD ifdef cppkey DAS_SWOTD2C 
#if defined DAS_SWOTD2C 
! dist_coast at rho points read in das_get_d2c.F
! dist < threshold (var 'dist_rho' unit km) set the_mask to 0
! TODO introduce parameters to do dist to coast masking or not
!      real, parameter :: swot_d2c_thresh=4.
#endif
#if defined DAS_SWOTSSHMSK 
      logical :: flag_swot_mask
#endif

      integer status, nc_id, v_id, i, j
     
#include "netcdf.inc"
!
! READ swot ssh
!
      flag_swot=.true.

      status = NF_OPEN(file_swot,0,nc_id)
      if (status.ne.NF_NOERR) then
        write(*,'(6x,A,(/8x,A))') '??? WARNING: GOES NOT FOUND',
     &          file_swot
        flag_swot=.false.
        goto 999
      endif
!
! ssh
!
! BLXD rename netcdf var 'ssh' to 'adt'
!     status = NF_INQ_VARID(nc_id, 'ssh',v_id)
      status = NF_INQ_VARID(nc_id, 'adt',v_id)
      if (status.ne.NF_NOERR) then
          print*, '??? WARNING: error in inquiring adt'
          print*, '??? trying ssh'
          status = NF_INQ_VARID(nc_id,'ssh',v_id)
      endif
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in inquiring swot ssh'
        flag_swot=.false.
        goto 999
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
      if (status.ne.NF_NOERR) then
        print*, '??? WARNING: error in reading swot ssh'
        flag_swot=.false.
        goto 999
      endif
      do j=0,Mm+1
        do i=0,Lm+1
          ! BLXD shift indices OK 
          ssh_swot(i,j)=anc(i+1,j+1)
        enddo
      enddo

!
! BLXD nc mask ? TODO add a flag_swot_mask to use mask_src
!      if mask_src found then use it
!
#if defined DAS_SWOTSSHMSK 
      flag_swot_mask=.true.
      if(flag_swot_mask)then
       status = NF_INQ_VARID(nc_id, 'mask_src',v_id)
       if (status.ne.NF_NOERR) then
         print*, '??? WARNING: error in inquiring swot ssh mask '
         flag_swot_mask=.false.
       endif
      endif
      if(flag_swot_mask)then
          status = NF_GET_VAR_DOUBLE(nc_id, v_id, anc)
          if (status.ne.NF_NOERR) then
              print*, '??? WARNING: error in reading swot mask'
              flag_swot_mask=.false.
          endif
      endif
      if(flag_swot_mask)then
           do j=0,Mm+1
            do i=0,Lm+1
            ! BLXD shift indices :
            ! anc        array declared as 1:Lm+2,1:Mm+2
            ! das_buff2  array declared as 0:Lm+1,0:Mm+1
              das_buff2(i,j)=anc(i+1,j+1)
            enddo
           enddo
      endif
#endif

!
! mask
!
      do j=0,Mm+1
        do i=0,Lm+1
          if (ssh_swot(i,j) .gt. -900 
     &        .and. h(i,j) .gt. hcmin) then
            swot_mask(i,j)=1.0
          else
            swot_mask(i,j)=0.0
          endif
        enddo
      enddo

! BLXD
#if defined DAS_SWOTD2C 
! dist_coast at rho points read in das_get_d2c.F
! dist < threshold (var 'dist_rho' unit km) set the_mask to 0
! TODO introduce parameters to do dist to coast masking or not
!      swot_d2c_thresh
!     if(flag_d2c)then
!      do j=0,Mm+1
!        do i=0,Lm+1
!         if (dist_coast(i,j) .lt. swot_d2c_thresh) swot_mask(i,j)=0.0
!        enddo
!      enddo
!     endif
#endif

#if defined DAS_SWOTSSHMSK 
      if(flag_swot_mask)then
           do j=0,Mm+1
            do i=0,Lm+1
              swot_mask(i,j)=swot_mask(i,j)*das_buff2(i,j)
            enddo
           enddo
      endif
#endif

      status = NF_CLOSE(nc_id)

999   continue

      if (flag_swot) then
        write(*, '(6x, A,(/8x,A))') 'OK: SWOT SSH READ IN',
     &         file_swot
      endif
    
      return
      end
#else
      subroutine das_get_swot_empty
      return
      end
#endif /* DAS_SWOTSSH */
