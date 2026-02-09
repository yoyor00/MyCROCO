#include "cppdefs.h"
#ifdef DAS_TMISST

      subroutine das_get_tmi
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
# include "das_buffer.h"
      integer i,j,ios,ios1,ios2, lenstr,len_tmi, it
      integer len_dist, status, nc_id, v_id
      integer block_size
#include "netcdf.inc"
!
! WARNING: TMI have large warmer bias near coastal line
!          TMI SST is not used when it falls within 1.25 degrees
!          from coast.     
! NOTE
!     SST 00-06, 06-12, 12-18, 18-19 UTC
!     we pick 12-18
!
! READ distance to coast
!
! BLXD flag_d2c to .false. if no distance to coast data
!      Also abort TMI reading if no distance to coast

      len_dist=lenstr(file_dist_coast)
      open(81, file=file_dist_coast(1:len_dist),form='unformatted'
     &               ,STATUS='old',                    IOSTAT=ios)
      close(81)
      if(ios.ne.0) then
        write(*,'(6x,A,3(/6x,A))')  ' ', '??? WARNING: can not find ',
     &               file_dist_coast, '??? no TMI SST used in DAS'
        return
      endif
!     status = NF_OPEN(file_dist_coast(1:len_dist),0,nc_id)
!     if (status.ne.NF_NOERR) then
!       write(*,'(6x,A,3(/6x,A))')  ' ', '??? WARNING: OPEN ERROR ',
!    &               file_dist_coast(1:len_dist), 
!    &               '??? no TMI SST used in DAS'
!       return
!     endif
!     status = NF_INQ_VARID(nc_id, 'dist_rho',v_id)
!     if (status.ne.NF_NOERR) then
!       write(*,'(6x,A,(/6x,A))')  
!    &                  '??? WARNING: ERROR INQUIRE dist_rho',
!    &                             '??? no TMI SST used in DAS'
!       return
!     endif
!     status = NF_GET_VAR_DOUBLE(nc_id, v_id, das_buff2)
!     if (status.ne.NF_NOERR) then
!       write(*,'(6x,A,(/6x,A))') '??? WARNING: ERROR READ dist_coast',
!    &                            '??? no TMI SST used in DAS'
!       return
!     endif
!     status = NF_CLOSE(nc_id)
!
!      do j=0,Mm+1
!        do i=0,Lm+1
!          dist_coast(i,j)=das_buff2(i,j)
!        enddo
!      enddo
!
! READ TMI SST
!
      flag_tmi=.true.
      block_size=(Lm+2)*(Mm+2) *2   ! NOTE: 2 is due to double precision
                                    !    If not SGI and F77, it is 8
!      
      len_tmi=lenstr(file_tmi)
      open(81, file=file_tmi(1:len_tmi),form='unformatted',STATUS='old',
     &              IOSTAT=ios,ACCESS='direct',RECL=block_size)
      if(ios.ne.0) then
        flag_tmi=.false.
        write(*,'(6x,A,3(/6x,A))')  ' ', '??? WARNING: can not find ',
     &           file_tmi, '??? no TMI SST used in DAS'
        close(81)
        return
      endif
!
      it=time_level*2-1
      PRINT*, '        REC=',it
      read(81,REC=it,IOSTAT=ios1)
     &    ((sst_tmi(i,j),i=0,Lm+1),j=0,Mm+1)
      it=time_level*2
      read(81,REC=it,IOSTAT=ios2)
     &    ((tmi_mask(i,j),i=0,Lm+1),j=0,Mm+1)
      if(ios1.ne.0  .or. ios2.ne.0) then
        flag_tmi=.false.
        write(*,'(6x,A)')
     &        '??? WARNING: Error in reading TMI SST'
      else
        write(*,'(6x,A,(/8x,A))')  
     &         file_tmi, 'TMI SST READ IN'
      endif



!
! dist < 1.5 degrees, not use TMI SST
!
      if(flag_d2c)then

        do j=0,Mm+1
          do i=0,Lm+1
            if (dist_coast(i,j) .lt. 1.5) tmi_mask(i,j)=0.0
          enddo
        enddo

!      write(*,*) ' get tmi'
!      write(*,*) (tmi_mask(i,20), i=0,Lm+1)
!      write(*,*) ' get dist_coast'
!      write(*,*) (dist_coast(i,20), i=0,Lm+1)
!
      endif
! not use TMI SST around island Quadelupe
!
!     do j=1,11
!       do i=64,74
!         tmi_mask(i,j)=0.0
!       enddo
!     enddo
!     close(81)
!
      return
      end
#else
      subroutine das_get_tmi_empty
      return
      end
#endif


      
      

      

      


