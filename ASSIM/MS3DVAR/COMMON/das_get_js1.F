#include "cppdefs.h"
#ifdef DAS_JASONSSH

      subroutine das_get_js1
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_innov.h"
      integer i,j,ios,ios1,ios2, lenstr,len_js1, it
!NOTE
!
      flag_js1=.true.
      num_js1=0
!
      do i=1,max_js1
        lon_js1(i)=-9999.0
        lat_js1(i)=-9999.0
        ssh_js1_raw(i)=-9999.0
      enddo
!      
      len_js1=lenstr(file_js1)
      open(81, file=file_js1(1:len_js1),form='formatted',STATUS='old',
     &              IOSTAT=ios)
      if(ios.ne.0) then
        flag_js1=.false.
        write(*,'(6x,A,(/8x,A))') '??? WARNING: NOT FOUND ',
     &           file_js1
        goto 10
      endif
!
      do it=1,max_js1
        read(81,800,end=100)
     &   lon_js1(it),lat_js1(it),ssh_js1_raw(it)
      enddo
100   continue
      num_js1=it-1
      write(*,'(6x,A,(/8x,A))')
     &         'OK: JASON-1 SSH READ IN', file_js1
10    continue
      close(81)
!
800   format(3f10.4)
      return
      end
#else
      subroutine das_get_js1_empty
      return
      end
#endif


      
      

      

      


