#include "cppdefs.h"

      subroutine das_get_d2c
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
# include "das_buffer.h"
      integer i,j,ios, lenstr
      integer len_dist, status, nc_id, v_id
#include "netcdf.inc"
!
! READ distance to coast
!
      flag_d2c=.true.
      len_dist=lenstr(file_dist_coast)
      open(81, file=file_dist_coast(1:len_dist),form='unformatted'
     &               ,STATUS='old',                    IOSTAT=ios)
      close(81)
      if(ios.ne.0) then
        flag_d2c=.false.
        write(*,'(6x,A,3(/6x,A))')  ' ', '??? WARNING: can not find ',
     &               file_dist_coast, '??? no TMI SST used in DAS'
        return
      endif
      status = NF_OPEN(file_dist_coast(1:len_dist),0,nc_id)
      if (status.ne.NF_NOERR) then
        flag_d2c=.false.
        write(*,'(6x,A,3(/6x,A))')  ' ', '??? WARNING: OPEN ERROR ',
     &               file_dist_coast(1:len_dist), 
     &               '??? no TMI SST used in DAS'
        return
      endif
      status = NF_INQ_VARID(nc_id, 'dist_rho',v_id)
      if (status.ne.NF_NOERR) then
        flag_d2c=.false.
        write(*,'(6x,A,(/6x,A))')  
     &                  '??? WARNING: ERROR INQUIRE dist_rho',
     &                             '??? no TMI SST used in DAS'
        return
      endif
      status = NF_GET_VAR_DOUBLE(nc_id, v_id, das_buff2)
      if (status.ne.NF_NOERR) then
        flag_d2c=.false.
        write(*,'(6x,A,(/6x,A))') '??? WARNING: ERROR READ dist_coast',
     &                            '??? no TMI SST used in DAS'
        return
      endif
      status = NF_CLOSE(nc_id)
!
      do j=0,Mm+1
        do i=0,Lm+1
          dist_coast(i,j)=das_buff2(i,j)
        enddo
      enddo
!
! dist < threshold (CHECK UNITS dist_rhod) set the_mask to 0
!
!     do j=0,Mm+1
!       do i=0,Lm+1
!         if (dist_coast(i,j) .lt. 1.5) the_mask(i,j)=0.0
!       enddo
!     enddo
!
!      write(*,*) (dist_coast(i,20), i=0,Lm+1)
!
!
      return
      end
!#else
!      subroutine das_get_d2c_empty
!      return
!      end
!#endif


      
      

      

      


