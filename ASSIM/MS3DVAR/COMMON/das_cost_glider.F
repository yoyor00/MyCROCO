#include "cppdefs.h"

      subroutine das_cost_glider
      implicit none
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
      integer num,k,kdas
      real cft,cfs
!
! sio and whoi gliders
!
!      sio_ot=1.0/(0.2*0.2)
!      sio_os=1.0/(0.5*0.5)
!
!      whoi_ot=1.0/(0.4*0.4)
!      whoi_os=1.0/(0.5*0.5)
!
!      ptsur_ot=1.0/(0.3*0.3)
!      ptsur_os=1.0/(0.3*0.3)
!
!      martn_ot=1.0/(0.3*0.3)
!      martn_os=1.0/(0.3*0.3)
!
! sio
!
!li      do num=1,prf_num_sio
!li        do k=1,max_prf
!li          kdas=ndas-k+1
!li          cft=t_s(Isio(num),Jsio(num),kdas,itemp)-t_sio(num,k)
!li          cfs=t_s(Isio(num),Jsio(num),kdas,isalt)-s_sio(num,k)
!li          costf=costf+0.5*(cft*cft*sio_ot + cfs*cfs*sio_os)
!li     &              *mask_sio(num,k)
!li        enddo
!li      enddo
!
! whoi
!
!li      do num=1,prf_num_whoi
!li        do k=1,max_prf
!li          kdas=ndas-k+1
!li          cft=t_s(Iwhoi(num),Jwhoi(num),kdas,itemp)-t_whoi(num,k)
!li          cfs=t_s(Iwhoi(num),Jwhoi(num),kdas,isalt)-s_whoi(num,k)
!li          costf=costf+0.5*(cft*cft*whoi_ot(k) + cfs*cfs*whoi_os(k))
!li     &              *mask_whoi(num,k)
!li        enddo
!li      enddo
!
! pt sur
!
!li      do num=1,prf_num_ptsur
!li        do k=1,max_prf
!li          kdas=ndas-k+1
!li          cft=t_s(Iptsur(num),Jptsur(num),kdas,itemp)
!li     &                                -t_ptsur(num,k)
!li          cfs=t_s(Iptsur(num),Jptsur(num),kdas,isalt)
!li     &                                -s_ptsur(num,k)
!li          costf=costf+0.5*(cft*cft*ptsur_ot 
!li     &                             + cfs*cfs*ptsur_os)
!li     &              *mask_ptsur(num,k)
!li        enddo
!li      enddo
!
! martin
!
      do num=1,prf_num_martn
        do k=1,max_prf
          kdas=ndas-k+1
          cft=t_s(Imartn(num),Jmartn(num),kdas,itemp)
     &                                -t_martn(num,k)
          cfs=t_s(Imartn(num),Jmartn(num),kdas,isalt)
     &                                -s_martn(num,k)
          costf=costf+0.5*(cft*cft*martn_ot
     &                             + cfs*cfs*martn_os)
     &              *mask_martn(num,k)
        enddo
      enddo
!
      return
      end
