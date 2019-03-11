c     this file is not differentiated

# include "cppdefs.h"      
      
      subroutine simul(indic,Sn,ad_x_f,f,ad_g_f,izs,rzs,dzs)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "adparam.h"
      real t1, t2, xtime
      
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
         call cost_fun(ad_x_f, cost)
         t2 = xtime()
         f = cost
         write (*,*) '--------------------------------------'
         write (*,*) 'end cost:',t2-t1
         write (*,*) '--------------------------------------'
c     call restore_croco_state()

c         eps = 1.0D-10

c         xbck(1) = ad_x(1) + eps
c         call cost_fun(xbck, cost)
c         write(*,*) '-||->', (cost-f)/eps
         call restore_croco_state()

         costb = 1
         cost = f

         write (*,*) '--------------------------------------'
         write (*,*) 'start adjoint'
         write (*,*) '--------------------------------------'
         t1 = xtime()
         call cost_fun_b(ad_x_f, ad_g_f, cost, costb)
         t2 = xtime()
         write (*,*) '--------------------------------------'
         write (*,*) 'end adjoint:',t2-t1
         write (*,*) '--------------------------------------'

         call restore_croco_state()

      end if
      return
      end

