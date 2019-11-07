c     this file is not differentiated

# include "cppdefs.h"      
      
      subroutine simul(indic,Sn,ad_x,cost,ad_g,izs,rzs,dzs)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer indic, Sn, i
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost, icost, costb
      double precision xbck(Sn), eps

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

         call save_croco_state()

         icost = 0
         call cost_fun(ad_x, icost)
         cost = icost
         call restore_croco_state()

         eps = 1.0D-10

         xbck(1) = ad_x(1) + eps
         icost = 0
         call cost_fun(xbck, icost)
         ad_g(1) = (icost-cost)/eps
         call restore_croco_state()

      end if
      return
      end

