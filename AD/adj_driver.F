c     this file is not differentiated

# include "cppdefs.h"

      subroutine simul(indic,Sn,ad_x,cost,ad_g,izs,rzs,dzs)
      implicit none
# include "adinc.h"
# include "adparam.h"
      real t1, t2, xtime

      integer indic, Sn, k, i
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost, icost, costb
      double precision xbck(Sn), eps

      integer ierr

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

         call save_croco_state()

c        It is not trivial to get direct result from reverse mode with
c        tapenade cf Original statements and results disappearing from
c        the Reverse code:
c        https://www-sop.inria.fr/tropics/tapenade/faq.html#adjlive this
c        need -nooptim adjointliveness but it is not sufficient, and it
c        needs some modifications in generated code.  Here we call
c        directly the computation of cost_fun.
         write (*,*) '--------------------------------------'
         write (*,*) 'start cost'
         write (*,*) '--------------------------------------'
         t1 = xtime()
         call cost_fun(ad_x, icost)
         t2 = xtime()
         ad_dir_time=t2-t1
         cost = icost
         write (*,*) '--------------------------------------'
         write (*,*) 'end cost:',ad_dir_time
         write (*,*) '--------------------------------------'

c     call restore_croco_state()

c         eps = 1.0D-10

c         xbck(1) = ad_x(1) + eps
c         call cost_fun(xbck, cost)
c         write(*,*) '-||->', (cost-f)/eps
         call restore_croco_state()

         costb = 1
         icost = cost

         write (*,*) '--------------------------------------'
         write (*,*) 'start adjoint'
         write (*,*) '--------------------------------------'

         ad_g(:) = 0.
         t1 = xtime()
         call save_croco_state()
         call cost_fun_b(ad_x, ad_g, icost, costb)
         t2 = xtime()
         ad_adj_time=t2-t1
         write (*,*) '--------------------------------------'
         write (*,*) 'end adjoint:',ad_adj_time
         write (*,*) '--------------------------------------'

         call restore_croco_state()

         do i=1,Sn
            write(*,*) 'ad_x(',i,'),grad(',i,')=',
     &           ad_x(i),ad_g(i)
         end do

      end if
      return
      end
