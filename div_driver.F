c     this file is not differentiated

# include "cppdefs.h"      
      
      subroutine simul(indic,Sn,ad_x_f,f,ad_g_f,izs,rzs,dzs)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "adparam.h"

      integer indic, Sn, i
      double precision f
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost, costb
      double precision xbck(Sn), eps

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

         call save_croco_state()

         call cost_fun(ad_x_f, cost)
         f = cost
         call restore_croco_state()

         eps = 1.0D-10

         xbck(1) = ad_x_f(1) + eps
         call cost_fun(xbck, cost)
         ad_g_f(1) = (cost-f)/eps
         call restore_croco_state()

      end if
      return
      end

