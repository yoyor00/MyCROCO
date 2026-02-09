#include "cppdefs.h"
      subroutine das_get_auv
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i,j,k,kdas,ios,ios1,ios2,
     &        lenstr,len_file,it
      real hh
!NOTE
!     SSH 00-06, 06-12, 12-18, 18-19 UTC
!     we pick 12-18

!
! initilaizzation
!
      flag_cal=.true.
      flag_dor=.true.
      do k=1,max_prf_cal
        num_dor(k)=0
      enddo
      do k=1,max_prf_dor
        num_dor(k)=0
      enddo
      do k=1,max_prf_cal
        do i=1,max_cal
          lon_cal(i,k)=-9999.0
          lat_cal(i,k)=-9999.0
        enddo
      enddo
      do k=1,max_prf_dor
        do i=1,max_dor
          lon_dor(i,k)=-9999.0
          lat_dor(i,k)=-9999.0
        enddo
      enddo
!
! CalPoly auv
!
      len_file=lenstr(file_cal_auv)
      open(81, file=file_cal_auv(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_cal=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &         file_cal_auv
        goto 30
      endif
      do k=1,max_prf_cal
        do it=1,max_cal
          read(81,900,err=300)
     &         lon_cal(it,k),lat_cal(it,k),hh,
     &         t_cal_raw(it,k),s_cal_raw(it,k)
        enddo
300     continue
        num_cal(k)=it-1
        if ( it .ge. max_cal) then
          write(*,'(6x,A)')
     &      '??? WARNING: max_cal too small'
          goto 20
        endif
      enddo
20    continue
      write(*,'(6x,A,(/8x,A))')
     &         'OK: CALPOLY AUV READ IN',file_cal_auv
30    continue
      close(81)
!
! Dorado auv
!
      len_file=lenstr(file_dor_auv)
      open(81, file=file_dor_auv(1:len_file),form='formatted',
     &         STATUS='old',IOSTAT=ios)
      if(ios.ne.0) then
        flag_dor=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &         file_dor_auv
        goto 50
      endif
      do k=1,max_prf_dor
        do it=1,max_dor
          read(81,900,err=500)
     &         lon_dor(it,k),lat_dor(it,k),hh,
     &         t_dor_raw(it,k),s_dor_raw(it,k)
        enddo
500     continue
        num_dor(k)=it-1
        if ( it .ge. max_dor) then
          write(*,'(6x,A)')
     &      '??? WARNING: max_dor too small'
          goto 60
        endif
      enddo
60    continue
      write(*,'(6x,A,(/8x,A))')
     &         'OK: DORADO AUV READ IN',file_dor_auv
50    continue
      close(81)
!
900   format(5f10.4)
      return
      end
