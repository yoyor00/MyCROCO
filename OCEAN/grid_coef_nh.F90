#include "cppdefs.h"
#ifdef NBQ

      subroutine grid_coef_nh(                                                      &
#if defined NBQ_IJK
         Istr,Iend,Jstr,Jend  &
                ,Hzw_half_nbq_inv,Hzr_half_nbq_inv                                  &
		,Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v                             &
		,Hzu_half_qdmu, Hzv_half_qdmv                                       &
#endif
         )
      !           Hzu_half_qdmu,Hzv_half_qdmv                                        &

!**********************************************************************
!
!                 Pre-computations for NH / NBQ modes
!
!**********************************************************************

      use module_nh
      use module_nbq
# ifdef MPI
      use module_parallel_nbq, only : ierr,par,OUEST,EST
# endif
      implicit none

#if defined NBQ_IJK	  
	  integer :: Istr,Iend,Jstr,Jend
#endif
# include "param_F90.h"
# include "scalars_F90.h"
#include "private_scratch.h"
# include "grid.h"
# include "ocean3d.h"
# include "nbq.h"

#if defined NBQ_IJK
       real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)	   
       real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
       real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,0:N)
       real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,0:N)
#endif

      integer :: i,j,k,it
      double precision :: val1, val2

#if defined NBQ_IJK
#include "compute_auxiliary_bounds.h"
#endif
!
#ifndef NBQ_MASS
#  define Hzr_half_nbq Hz
#endif

!
!
!**********************************************************************
!    Initialisations and updates
!**********************************************************************
#if !defined NBQ_IJK
        do k=1,N
          do j=jstr_nh ,jend_nh
          do i=istru_nh,iend_nh+1
              gdepth_u(i,j,k) = zr_half_nbq(i,j,k)-zr_half_nbq(i-1,j,k)
              
              coefa_u(i,j,k)  = 0.25*pm_u(i,j)*                        &
                     (Hzr_half_nbq(i,j,k  )+Hzr_half_nbq(i-1,j,k ))/  &
                     (Hzw_half_nbq(i,j,k-1)+Hzw_half_nbq(i-1,j,k-1))
              coefb_u(i,j,k)  = 0.25*pm_u(i,j)*                        &
                     (Hzr_half_nbq(i,j,k  )+Hzr_half_nbq(i-1,j,k ))/   &
                     (Hzw_half_nbq(i,j,k  )+Hzw_half_nbq(i-1,j,k  ))
          enddo
          enddo

          do j=jstrv_nh,jend_nh+1
          do i=istr_nh ,iend_nh
              gdepth_v(i,j,k) = zr_half_nbq(i,j,k)-zr_half_nbq(i ,j-1,k)
 
              coefa_v(i,j,k)  = 0.25*pn_v(i,j)*                        &
                     (Hzr_half_nbq(i,j,k  )+Hzr_half_nbq(i,j-1,k ))/  &
                     (Hzw_half_nbq(i,j,k-1)+Hzw_half_nbq(i,j-1,k-1))
              coefb_v(i,j,k)  = 0.25*pn_v(i,j)*                        &
                     (Hzr_half_nbq(i,j,k  )+Hzr_half_nbq(i,j-1,k ))/  &
                     (Hzw_half_nbq(i,j,k  )+Hzw_half_nbq(i,j-1,k  ))
          enddo
          enddo
        enddo 
        
        do j = jstr_nh ,jend_nh
        do i = istru_nh,iend_nh+1
            gdepth_u(i,j,0)   = zw_half_nbq(i,j,0)-zw_half_nbq(i-1,j,0)
            gdepth_u(i,j,N+1) = zw_half_nbq(i,j,N)-zw_half_nbq(i-1,j,N)

            coefa_u(i,j,0)    = 0.5 * pm_u(i,j) * real (slip_nbq)  
            coefa_u(i,j,1)    = 0. 
            coefb_u(i,j,N)    = 0.     
            coefb_u(i,j,N+1)  = 0.5 * pm_u(i,j)
        enddo
        enddo

        do j = jstrv_nh,jend_nh+1
        do i = istr_nh ,iend_nh
            gdepth_v(i,j,0)   = zw_half_nbq(i,j,0)-zw_half_nbq(i,j-1,0)
            gdepth_v(i,j,N+1) = zw_half_nbq(i,j,N)-zw_half_nbq(i,j-1,N)
            
            coefa_v(i,j,0)    = 0.5 * pn_v(i,j) * real (slip_nbq) 
            coefa_v(i,j,1)    = 0.                    
            coefb_v(i,j,N)    = 0.        
            coefb_v(i,j,N+1)  = 0.5 * pn_v(i,j)
        enddo
        enddo
#endif

#if defined NBQ_IJK
        do k=1,N
        do j=jstr_nh-1,jend_nh+1
        do i=istr_nh-1,iend_nh+1
          Hzr_half_nbq_inv(i,j,k)=1.d0/max(1.e-30,Hzr(i,j,k))
# ifdef MASKING
          Hzr_half_nbq_inv(i,j,k)=Hzr_half_nbq_inv(i,j,k)*rmask(i,j)
# endif
        enddo
        enddo
        enddo

        do k=1,N
        do j=jstr_nh,jend_nh
        do i=istru_nh,iend_nh+1
          Hzu_half_qdmu(i,j,k)=0.5*(Hzr(i-1,j,k)+Hzr(i,j,k))*pm_u(i,j)   
#if defined MASKING
          Hzu_half_qdmu(i,j,k) = Hzu_half_qdmu(i,j,k) * umask(i,j)
#endif  
        enddo
        enddo
        enddo 
        do k=1,N
        do j=jstrv_nh,jend_nh+1
        do i=istr_nh,iend_nh
          Hzv_half_qdmv(i,j,k)=0.5*(Hzr(i,j-1,k)+Hzr(i,j,k))*pn_v(i,j)   
#if defined MASKING
          Hzv_half_qdmv(i,j,k) = Hzv_half_qdmv(i,j,k) * vmask(i,j)
#endif  
        enddo
        enddo
        enddo 

        do k=0,N
        do j=jstr_nh-1,jend_nh+1
        do i=istr_nh-1,iend_nh+1
          Hzw_half_nbq_inv(i,j,k)=1.d0/max(1.e-30,Hzw_half_nbq(i,j,k)) 
# ifdef MASKING
          Hzw_half_nbq_inv(i,j,k)=Hzw_half_nbq_inv(i,j,k)*rmask(i,j)
# endif
        enddo
        enddo
        enddo
		
        do k=0,N
        do j=jstr_nh,jend_nh
        do i=istru_nh,iend_nh+1
          Hzw_half_nbq_inv_u(i,j,k)=0.25d0*2.d0/max(1.e-30,Hzw_half_nbq(i,j,k)+Hzw_half_nbq(i-1,j,k)) 
# if defined MASKING
          Hzw_half_nbq_inv_u(i,j,k)=Hzw_half_nbq_inv_u(i,j,k) *umask(i,j)   
# endif      
        enddo
        enddo
        enddo
              
        do k=0,N
        do j=jstrv_nh,jend_nh+1
        do i=istr_nh,iend_nh
          Hzw_half_nbq_inv_v(i,j,k)=0.25d0*2.d0/max(1.e-30,Hzw_half_nbq(i,j,k)+Hzw_half_nbq(i,j-1,k)) 
# if defined MASKING	
          Hzw_half_nbq_inv_v(i,j,k)= Hzw_half_nbq_inv_v(i,j,k)*vmask(i,j)
# endif      	  
        enddo
        enddo
        enddo          

#endif

!**********************************************************************
!
! MPI NEEDS:
!
! gdepth_u, coefa_u and coefb_u: 
!                we need here ( istru_nh-1,[jstr_nh:jend_nh],[0:N+1] )
!                             ( iendu_nh+1,[jstr_nh:jend_nh],[0:N+1] )               
! gdepth_v, coefa_v and coefb_v: 
!                we need here ( [istr_nh:iend_nh],jstrv_nh-1,[0:N+1] )
!                             ( [istr_nh:iend_nh],jendv_nh+1,[0:N+1] ) 
!  
!**********************************************************************
        
     call grid_exchange()
        
      return
      end subroutine grid_coef_nh

#else
      subroutine grid_coef_nh_empty
      return
      end 
#endif


