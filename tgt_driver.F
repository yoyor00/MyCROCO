c     this file is not differentiated

# include "cppdefs.h"

C     !! adjoint check only
C     !! x size = 1 ==> dfdx ok
      subroutine simul(indic,Sn,ad_x,cost,grad,izs,rzs,dzs)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer indic, Sn, i, k
      double precision grad(Sn), idfdx
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost,icost

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

c         call cost_fun(x, cost)
c         f = cost
c         call restore_croco_state()
c     call save_croco_state()
         do i=1,Sn
            ad_xd(i) = 0
         end do
         icost = cost
         do i=1,Sn
            ad_xd(i) = 1
            cost = icost
            idfdx = grad(i)
            call save_croco_state()
            write (*,*) '--------------------------------------'
            write (*,*) 'start tangent:',i
            write (*,*) '--------------------------------------'
            call cost_fun_d(ad_x, ad_xd, cost, idfdx)

            write (*,*) '--------------------------------------'
            write (*,*) 'end tangent:',i
            write (*,*) '--------------------------------------'
            call restore_croco_state()
            ad_xd(i) = 0
            grad(i)  = idfdx
         end do
         write (*,*) '######################################'
         do i=1,Sn
            write(*,*) 'node:',mynode,'ad_x(',i,')',
     &           'grad(',i,')=',
     &           ad_x(i),grad(i)
         end do
         write (*,*) '######################################'

      end if
      return
      end

