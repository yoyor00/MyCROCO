#include "cppdefs.h"
#ifdef DAS_TPSSH

      subroutine das_get_tp
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i,j,ios,ios1,ios2, lenstr,len_tp, it
!NOTE
!
      flag_tp=.true.
      num_tp=0
!
      do i=1,max_tp
        lon_tp(i)=-9999.0
        lat_tp(i)=-9999.0
        ssh_tp_raw(i)=-9999.0
      enddo
!      
      len_tp=lenstr(file_tp)
      open(81, file=file_tp(1:len_tp),form='formatted',STATUS='old',
     &              IOSTAT=ios)
      if(ios.ne.0) then
        flag_tp=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &           file_tp
        goto 10
      endif
!
      do it=1,max_tp
        read(81,800,end=100)
     &   lon_tp(it),lat_tp(it),ssh_tp_raw(it)
      enddo
100   continue
      num_tp=it-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: TP SSH READ IN', file_tp
10    continue
      close(81)
!
800   format(3f10.4)
      return
      end
#else
      subroutine das_get_tp_empty
      return
      end
#endif


      
      

      

      


