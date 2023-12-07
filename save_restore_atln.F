      subroutine save_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      CALL PUSHINTEGER4(ad_cost_counter)
c      CALL PUSHREAL8ARRAY(svstr, (lm+padd_x+4)*(mm+padd_e+4))
c      CALL PUSHREAL8ARRAY(sustr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL PUSHREAL8ARRAY(vbar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHREAL8ARRAY(ubar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHREAL8ARRAY(zeta, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL PUSHINTEGER4(may_day_flag)
      CALL PUSHBOOLEAN(synchro_flag)
      CALL PUSHREAL8(zobt)
      CALL PUSHINTEGER4(knew)
      CALL PUSHINTEGER4(iic)
      CALL PUSHREAL8ARRAY(zob, (lm+padd_x+4)*(mm+padd_e+4))

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
c      write(*,*) 'minzob',minval(zob_bck)
c      write(*,*) 'maxzob',maxval(zob_bck)
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "adinc.h"
# include "adparam.h"

      CALL POPREAL8ARRAY(zob, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPINTEGER4(iic)
      CALL POPINTEGER4(knew)
      CALL POPREAL8(zobt)
      CALL POPBOOLEAN(synchro_flag)
      CALL POPINTEGER4(may_day_flag)
      CALL POPREAL8ARRAY(zeta, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL POPREAL8ARRAY(ubar, (lm+padd_x+4)*(mm+padd_e+4)*4)
      CALL POPREAL8ARRAY(vbar, (lm+padd_x+4)*(mm+padd_e+4)*4)
c     CALL POPREAL8ARRAY(sustr, (lm+padd_x+4)*(mm+padd_e+4))
c     CALL POPREAL8ARRAY(svstr, (lm+padd_x+4)*(mm+padd_e+4))
      CALL POPINTEGER4(ad_cost_counter)

      write(*,*) 'restore+++'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '++++++++++'

      end subroutine restore_croco_state
