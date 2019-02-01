c     this file is not differentiated

# include "cppdefs.h"

C     !! adjoint check only
C     !! x size = 1 ==> dfdx ok
      subroutine simul(indic,Sn,ad_x_f,f,dfdx,izs,rzs,dzs)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "adparam.h"

      integer indic, Sn, i
      double precision f, dfdx(Sn)
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
         call cost_fun_d(ad_x_f, ad_xd, cost, dfdx)

         write (*,*) '--------------------------------------'
         write (*,*) 'end tangent'
         write (*,*) '--------------------------------------'
         do i=1,Sn
            write(*,*) '1--->', i,ad_x_f(i),ad_xd(i),cost,dfdx(i)
         end do
         f = cost
         call restore_croco_state()

c         cost = f
c         dfdx(1) = 0
c         costd = 1
c         call cost_fun_d(x, dfdx, cost, costd)
c         do i=1,Sn
c            write(*,*) '2--->', i,ad_x(1),f,cost,dfdx(i),costd
c         end do
c         call restore_croco_state()

         
c         do i=1,Sn
c            write(*,*) '--->', i,f,cost,dfdx(i)
c         end do

      end if
      return
      end

