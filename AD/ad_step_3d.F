      subroutine ad_step(nstep)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k,iifroot,nstep

      write(*,*) 'start ad_step'
C$AD BINOMIAL-CKP nstep+1 72 1
      do k=1,nstep
         do iifroot = 0,nfast+2
            call step()
         end do
c     write (*,*) 'k,iifroot,kstp,knew:',k,iifroot,kstp,knew
         ad_step_counter = ad_step_counter + 1
      end do
      write(*,*) 'end ad_step'

      end subroutine ad_step
