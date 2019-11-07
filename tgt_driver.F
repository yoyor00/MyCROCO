c     this file is not differentiated

# include "cppdefs.h"

C     !! adjoint check only
C     !! x size = 1 ==> dfdx ok
      subroutine simul(indic,Sn,ad_x,cost,dfdx,izs,rzs,dzs)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer indic, Sn, i
      double precision dfdx(Sn)
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

         call save_croco_state()

c         call cost_fun(x, cost)
c         f = cost
c         call restore_croco_state()
c         call save_croco_state()
         dfdx(1) = 0
         ad_xd(1) = 1
         write (*,*) '--------------------------------------'
         write (*,*) 'start tangent'
         write (*,*) '--------------------------------------'
         call cost_fun_d(ad_x, ad_xd, cost, dfdx)

         write (*,*) '--------------------------------------'
         write (*,*) 'end tangent'
         write (*,*) '--------------------------------------'
         do i=1,Sn
            write(*,*) '1--->', i,ad_x(i),ad_xd(i),cost,dfdx(i)
         end do
         call restore_croco_state()

      end if
      return
      end

