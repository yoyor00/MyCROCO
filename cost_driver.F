c     direct call to cost function 

      subroutine cost_driver(iicroot)

      implicit none
# include "adinc.h"
# include "adparam.h"     
      integer iicroot

      double precision icost

      do k=1,ad_array_size/nnodes
         ad_g(k)=0.
         ad_x(k)=0.
      end do

      ad_x(1) = 0.1D-5
      
      call cost_fun(ad_x, icost)

      write(*,*) 'cost=', icost
      return
      end subroutine cost_driver
