#include "cppdefs.h"

      subroutine das_adj_cglider
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
!      ptsur_ot=1.0/(0.3*0.3)
!      ptsur_os=1.0/(0.3*0.3)
!
!      martn_ot=1.0/(0.3*0.3)
!      martn_os=1.0/(0.3*0.3)
!
! martn
!
      do num=1,prf_num_martn
        do k=1,max_prf
          kdas=ndas-k+1
          cft=t_s(Imartn(num),Jmartn(num),kdas,itemp)-t_martn(num,k)
          cfs=t_s(Imartn(num),Jmartn(num),kdas,isalt)-s_martn(num,k)
!
          t_s_adj(Imartn(num),Jmartn(num),kdas,itemp)
     &       =t_s_adj(Imartn(num),Jmartn(num),kdas,itemp)
     &               +cft*martn_ot*mask_martn(num,k)
          t_s_adj(Imartn(num),Jmartn(num),kdas,isalt)
     &       =t_s_adj(Imartn(num),Jmartn(num),kdas,isalt)
     &               +cfs*martn_os*mask_martn(num,k)
        enddo
      enddo
!
! ptsur
!
!LI      do num=1,prf_num_ptsur
!LI        do k=1,max_prf
!LI          kdas=ndas-k+1
!LI          cft=t_s(Iptsur(num),Jptsur(num),kdas,itemp)-t_ptsur(num,k)
!LI          cfs=t_s(Iptsur(num),Jptsur(num),kdas,isalt)-s_ptsur(num,k)
!LI!
!LI          t_s_adj(Iptsur(num),Jptsur(num),kdas,itemp)
!LI     &       =t_s_adj(Iptsur(num),Jptsur(num),kdas,itemp)
!LI     &               +cft*ptsur_ot*mask_ptsur(num,k)
!LI          t_s_adj(Iptsur(num),Jptsur(num),kdas,isalt)
!LI     &       =t_s_adj(Iptsur(num),Jptsur(num),kdas,isalt)
!LI     &               +cfs*ptsur_os*mask_ptsur(num,k)
!LI        enddo
!LI      enddo
!
! whoi
!
!LI      do num=1,prf_num_whoi
!LI        do k=1,max_prf
!LI          kdas=ndas-k+1
!LI          cft=t_s(Iwhoi(num),Jwhoi(num),kdas,itemp)-t_whoi(num,k)
!LI          cfs=t_s(Iwhoi(num),Jwhoi(num),kdas,isalt)-s_whoi(num,k)
!LI!
!LI          t_s_adj(Iwhoi(num),Jwhoi(num),kdas,itemp)
!LI    &       =t_s_adj(Iwhoi(num),Jwhoi(num),kdas,itemp)
!LI     &               +cft*whoi_ot(k)*mask_whoi(num,k)
!LI          t_s_adj(Iwhoi(num),Jwhoi(num),kdas,isalt)
!LI     &       =t_s_adj(Iwhoi(num),Jwhoi(num),kdas,isalt)
!LI     &               +cfs*whoi_os(k)*mask_whoi(num,k)
!LI        enddo
!LI      enddo
!
! sio
!
!LI      do num=1,prf_num_sio
!LI        do k=1,max_prf
!LI          kdas=ndas-k+1
!LI          cft=t_s(Isio(num),Jsio(num),kdas,itemp)-t_sio(num,k)
!LI          cfs=t_s(Isio(num),Jsio(num),kdas,isalt)-s_sio(num,k)
!LI!
!LI          t_s_adj(Isio(num),Jsio(num),kdas,itemp)
!LI     &       =t_s_adj(Isio(num),Jsio(num),kdas,itemp)
!LI     &               +cft*sio_ot*mask_sio(num,k)
!LI          t_s_adj(Isio(num),Jsio(num),kdas,isalt)
!LI     &       =t_s_adj(Isio(num),Jsio(num),kdas,isalt)
!LI     &               +cfs*sio_os*mask_sio(num,k)
!LI        enddo
!LI      enddo
!
      return
      end
