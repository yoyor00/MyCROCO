      subroutine cost_fun(x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"
      double precision delta_z, delta_r
      double precision cost, sqn_x, sqn_z, sqn_r
      double precision rdrg2_bck, h_bck(GLOBAL_2D_ARRAY)
      integer i
      integer j
      integer k
      integer Istr,Iend,Jstr,Jend
# include "compute_auxiliary_bounds.h"
      rdrg2_bck = rdrg2
      do j=JstrR,JendR
         do i=IstrR,IendR
            h_bck(i,j) = h(i,j)
         end do
      end do

      call set_state(x)
      cost = 0.
      sqn_x = 0.

      do k=1,120
         call step()
         sqn_z = 0.
         sqn_r = 0.
         do j=1,Mm
            do i=1,Lm
!               print*, '>>>',i,j,sim_iicroot+k, zeta(i,j,knew),
!     &              obsz(i,j,sim_iicroot+k)
               delta_z = (zeta(i,j,knew)-obsz(i,j,sim_iicroot+k))
               sqn_z = sqn_z + delta_z*delta_z
               sqn_r = sqn_r + (h(i,j)-h_bck(i,j))*(h(i,j)-h_bck(i,j))
            end do
         end do

         cost = cost + sqn_z
      end do
      sqn_x = sqn_x +
     &     (x(1) - rdrg2_bck)*(x(1) - rdrg2_bck)
      cost = cost + sqn_x
      return
      end subroutine cost_fun


      subroutine save_croco_state()
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"

      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real h_bck(GLOBAL_2D_ARRAY)
      real rdrg2_bck

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

      common /backup/ ubar_bck, vbar_bck, zeta_bck, h_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, rdrg2_bck

      ubar_bck(:,:,:) = ubar(:,:,:)
      vbar_bck(:,:,:) = vbar(:,:,:)
      zeta_bck(:,:,:) = zeta(:,:,:)
      h_bck(:,:) = h(:,:)
      rdrg2_bck = rdrg2

      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"

      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real h_bck(GLOBAL_2D_ARRAY)
      real rdrg2_bck

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

      common /backup/ ubar_bck, vbar_bck, zeta_bck, h_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, rdrg2_bck

      ubar(:,:,:) = ubar_bck(:,:,:)
      vbar(:,:,:) = vbar_bck(:,:,:)
      zeta(:,:,:) = zeta_bck(:,:,:)
      h(:,:) = h_bck(:,:)
      rdrg2 = rdrg2_bck

      kstp = kstp_bck
      krhs = krhs_bck
      knew = knew_bck
      iic = iic_bck

      write(*,*) 'restore+++'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '++++++++++'

      end subroutine restore_croco_state

      subroutine set_state(x)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"
      integer i
      integer j

      integer JS
      integer IS
      integer Istr,Iend,Jstr,Jend
# include "compute_auxiliary_bounds.h"

      JS = JendR - JstrR
      IS = IendR - IstrR

      write(*,*) 'set_state'
      write(*,*) 'krhs=',krhs

      do j=1,Mm
         do i=1,Lm
            h(i,j) = real(x(1+i+Lm*(j-1)))
         end do
      end do

      rdrg2 = real(x(1))

      write(777, 900) rdrg2
      write(778, 900) h

      write(*, *) 'rdrg2=',rdrg2
      write(*, *) 'h(4,4)=',h(4,4)

 900  format(1PE16.9)

      return
      end subroutine set_state
