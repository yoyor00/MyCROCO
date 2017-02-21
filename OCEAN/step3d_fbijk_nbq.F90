#include "cppdefs.h"
#if defined NBQ && defined NBQ_IJK
!
!======================================================================
!                      NBQ-Mode for NH-modeling
!                            Main Routine
!======================================================================
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief SNBQ driver : Non-hydrostatic algorithm with the 
!                       Non-boussinesq solver.
!
!> @details NBQ time step. See SNBQ web pages :
!  http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages
!    Algorithme_NBQ.htm      --> SNBQ algorithm:               
!    Algebrique_SNBQ.htm     --> SNBQ algebric representation
!    Couplage_Numerique.htm  --> Numerical coupling
!    
!    Couplage_Split_NBQ.htm  --> Coupling Splitting
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
      subroutine step3d_fbijk_nbq (Istr,Iend,Jstr,Jend, WORK, &
	    Hzw_half_nbq_inv,Hzr_half_nbq_inv,                    &
		Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v,               &
        zwrk1,zwrk2,zwrk3,zwrk4,zwrk5,FC,FX            &
		)

      use module_nh 
      use module_nbq

      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

      integer Istr,Iend,Jstr,Jend,i,j,k
      real WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real Hzw_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,0:N)	   
      real Hzr_half_nbq_inv(PRIVATE_2D_SCRATCH_ARRAY,N)
      real Hzw_half_nbq_inv_u(PRIVATE_2D_SCRATCH_ARRAY,0:N)
      real Hzw_half_nbq_inv_v(PRIVATE_2D_SCRATCH_ARRAY,0:N)	  
	  real zwrk1(PRIVATE_2D_SCRATCH_ARRAY,2)
	  real zwrk2(PRIVATE_2D_SCRATCH_ARRAY,2)
	  real zwrk3(PRIVATE_2D_SCRATCH_ARRAY,2)
	  real zwrk4(PRIVATE_2D_SCRATCH_ARRAY,2)	  
	  real zwrk5(PRIVATE_2D_SCRATCH_ARRAY)
	  real FC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
	  real FX(PRIVATE_2D_SCRATCH_ARRAY)
	   
# ifdef MPI
      include 'mpif.h'
# endif

      real :: dum_s
      double precision :: a_m,b_m
	  integer :: k1,k2
# ifndef MPI
      integer mynode
      mynode=0
# endif
# undef DEBUG

!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
        if (iif==1.and.iic==1) call initial_nh_tile (3,Istr,Iend,Jstr,Jend)
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equations:
!  ru+rubar (or lambda_ext+lambda_int)
!  dzdt*rhosurf
!-------------------------------------------------------------------
!
!     call ru_nbq(1)     ! MPI ! TBD
!
!------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!
# ifdef NBQ_IMP
      if (iif.eq.1.and.ifl_imp_nbq.eq.1) call implicitijk_nbq (1)
# endif
!
!*******************************************************************
!*******************************************************************
!              NBQ mode iteration (main loop)
!*******************************************************************
!*******************************************************************

      do iteration_nbq=1,iteration_nbq_max

!
!-------------------------------------------------------------------
!
!-------------------------------------------------------------------

        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh
              div_nbq(i,j,k)=-visc2_nbq*div_nbq(i,j,k)     &
                             +soundspeed2_nbq*rho_nbq(i,j,k)
            enddo
          enddo
        enddo

!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------
!
!  XI- and ETA-Directions:
!

#define ddiv_nbqdz_u zwrk1
#define ddiv_nbqdz_v zwrk2
#define ddiv_nbqdz zwrk5

        k2 = 1
		do k=0,N
		  k1=k2
		  k2=3-k1
		  if (k.eq.0 .or. k.eq.N) then ! Bottom/Top Boundary conditions
		    do j=JstrU_nh,JendU_nh
			  do i=IstrU_nh,IendU_nh
			    ddiv_nbqdz_u(i,j,k2)=0.
			  enddo
			enddo
			do j=JstrV_nh,JendV_nh
			  do i=IstrV_nh,IendV_nh
			    ddiv_nbqdz_v(i,j,k2)=0.
			  enddo
			enddo
		  else
		    do j=Jstr_nh,Jend_nh
			  do i=Istr_nh,Iend_nh
			    ddiv_nbqdz(i,j)=div_nbq(i  ,j,k+1) - div_nbq(i  ,j,k)
			  enddo
			enddo
			do j=JstrU_nh,JendU_nh
			  do i=IstrU_nh,IendU_nh
			    ddiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(ddiv_nbqdz(i,j)+ddiv_nbqdz(i-1,j))
			  enddo
			enddo
			do j=JstrV_nh,JendV_nh
			  do i=IstrV_nh,IendV_nh
			    ddiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(ddiv_nbqdz(i,j)+ddiv_nbqdz(i,j-1))
			  enddo
			enddo
		  endif
		  if (k.gt.0) then
		    do j=JstrU_nh,JendU_nh
			  do i=IstrU_nh,IendU_nh
			    dum_s=gdepth_u(i,j,k)*(ddiv_nbqdz_u(i,j,k2)+ddiv_nbqdz_u(i,j,k1)) &   ! dZdx * (d(delta p)dz)_u
				         -(div_nbq(i,j,k)-div_nbq(i-1,j,k))                           ! - d(delta p)dx
				dum_s=dum_s*(0.5*(Hzr_half_nbq(i-1,j,k)+Hzr_half_nbq(i,j,k))) * pm_u(i,j)

                rhssumu_nbq(i,j,k) = rhssumu_nbq(i,j,k) + dum_s
                qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k)                                             & 
                            + dtnbq * ( dum_s +  rho0*(ruint_nbq(i,j,k)+ruext_nbq(i,j,k)))
							
			  enddo
			enddo
		    do j=JstrV_nh,JendV_nh
			  do i=IstrV_nh,IendV_nh
			    dum_s=gdepth_v(i,j,k)*(ddiv_nbqdz_v(i,j,k2)+ddiv_nbqdz_v(i,j,k1)) &   ! dZdy * (d(delta p)dz)_v
				         -(div_nbq(i,j,k)-div_nbq(i,j-1,k))                           ! - d(delta p)dy
				dum_s=dum_s*(0.5*(Hzr_half_nbq(i,j-1,k)+Hzr_half_nbq(i,j,k))) * pn_v(i,j)

                rhssumv_nbq(i,j,k) = rhssumv_nbq(i,j,k) + dum_s
                qdmv_nbq(i,j,k) = qdmv_nbq(i,j,k)                                             & 
                            + dtnbq * ( dum_s +  rho0*(rvint_nbq(i,j,k)+rvext_nbq(i,j,k)))
							
			  enddo
			enddo
		  endif
		enddo

#undef ddiv_nbqdz_u
#undef ddiv_nbqdz_v
#undef ddiv_nbqdz

!
!  U-momentum open boundary conditions
!
# ifdef OBC_NBQ
!       call unbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)   ! TBD
!       call vnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)   ! TBD
# endif
          
!
!  Message passing: Send U (51) 
!
!       call parallele_nbq(51)   ! TBD
!
!  Message passing: Send V (52) 
!
!       call parallele_nbq(52)   ! TBD
!
!-------------------------------------------------------------------
!      Vertical Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
# ifndef NBQ_IMP
!
!  Z-Direction: Explicit
!
        do j=Jstr_nh,Jend_nh
          do k=1,N-1
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   div_nbq(i,j,k) - div_nbq(i,j,k+1)
               rhssumw_nbq(i,j,k) = rhssumw_nbq(i,j,k) + dum_s                              
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rho0 * rwint_nbq(i,j,k) )
             enddo             
           enddo
          k=N
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   div_nbq(i,j,k)
               rhssumw_nbq(i,j,k) = rhssumw_nbq(i,j,k) + dum_s                              
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rho0 * rwint_nbq(i,j,k) )
             enddo             		   
        enddo
        
# else

!       call parallele_nbq(151)  ! u only  ! TBD 
!       call parallele_nbq(152)  ! v only  ! TBD

        call implicitijk_nbq (2)    ! TBD

# endif
!
!      Vertical momentum open boundary conditions
!
# ifdef OBC_NBQ
!       call wnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!      Message passing 
!-------------------------------------------------------------------
!
!  Send
!       call parallele_nbq(53)    ! TBD
!
!  Receive
# ifndef NBQ_IMP
!       call parallele_nbq(151)   ! TBD 
!       call parallele_nbq(152)   ! TBD  
# endif 
!       call parallele_nbq(153)   ! TBD   
!
!-------------------------------------------------------------------
!      Message passing
!-------------------------------------------------------------------
!
!  Send
!       call parallele_nbq(7)     ! TBD

!  Receive
!       call parallele_nbq(17)    ! TBD

!
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
#  ifdef ACOUSTIC
!      call density_nbq(11)       ! TBD
#  endif
!
!-------------------------------------------------------------------
!      Mass equation: 
!-------------------------------------------------------------------
!

           do j=Jstr_nh,Jend_nh
             do i=Istr_nh,Iend_nh
			   WORK(i,j)=pm(i,j)*pn(i,j)
			 enddo
		   enddo
		   

#define dZdxq_u zwrk1
#define dZdxq_w zwrk2
#define dZdyq_v zwrk3
#define dZdyq_w zwrk4
#define FY zwrk5

        k2 = 1
		do k=0,N
		  k1=k2
		  k2=3-k1

		  if (k.lt.N) then
		     do j=Jstr_nh,Jend_nh
		       do i=Istr_nh,Iend_nh+1
		         dZdxq_u(i,j,k2)=gdepth_u(i,j,k+1)*qdmu_nbq(i,j,k+1)    ! (dZdx * (rho u))_u
		       enddo
		     enddo
		     do j=Jstr_nh,Jend_nh+1
		       do i=Istr_nh,Iend_nh
		         dZdyq_v(i,j,k2)=gdepth_v(i,j,k+1)*qdmv_nbq(i,j,k+1)    ! (dZdy * (rho v))_v
		       enddo
		     enddo			 
		  endif

		  if (k.eq.0 .or. k.eq.N) then	! Bottom/Top boundary conditions	  
		   do j=Jstr_nh,Jend_nh
		     do i=Istr_nh,Iend_nh+1
		       dZdxq_w(i,j,k2)=0.5*(zw_half_nbq(i,j,0)-zw_half_nbq(i-1,j,0))*qdmu_nbq(i,j,1)*Hzr_half_nbq_inv(i,j,1)
		     enddo
		   enddo
		   do j=Jstr_nh,Jend_nh+1
		     do i=Istr_nh,Iend_nh
		       dZdyq_w(i,j,k2)=0.5*(zw_half_nbq(i,j,0)-zw_half_nbq(i,j-1,0))*qdmv_nbq(i,j,1)*Hzr_half_nbq_inv(i,j,1)
		     enddo
		   enddo		   
		  else
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh+1
			     dZdxq_w(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) ! (dZdx * (rho u))_uw/Hzw_u
              enddo 
            enddo
			do j=Jstr_nh,Jend_nh+1
		      do i=Istr_nh,Iend_nh
			     dZdyq_w(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
              enddo 
            enddo
          endif

          if (k.gt.0) then
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh+1
			    FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
              enddo
            enddo
		    do j=Jstr_nh,Jend_nh+1
		      do i=Istr_nh,Iend_nh
			    FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
              enddo
            enddo			
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh+1
			    div_nbq(i,j,k)=FX(i,j)+FX(i+1,j)+FY(i,j)+FY(i,j+1)
              enddo
            enddo
          endif
		enddo		 

#undef dZdxq_u
#undef dZdxq_w
#undef zwrk5


#define FY zwrk5
			   
			   
           do j=Jstr_nh,Jend_nh
             do i=Istr_nh,Iend_nh
			   FC(i,0)=0.                 ! Bottom boundary condition
             enddo
		     do k=1,N
               do i=Istr_nh,Iend_nh
			     FC(i,k)=Hzw_half_nbq_inv(i,j,k) * qdmw_nbq(i,j,k)
               enddo	
               do i=Istr_nh,Iend_nh
			     div_nbq(i,j,k)=div_nbq(i,j,k)+FC(i,k)-FC(i,k-1)
               enddo
             enddo
           enddo 

           do k=1,N
             do j=Jstr_nh,Jend_nh		   
               do i=Istr_nh,Iend_nh+1
                 FX(i,j)=on_u(i,j)* qdmu_nbq(i,j,k)
#ifdef MASKING
                 FX(i,j) = FX(i,j) * umask(i,j)
#endif
               enddo
              enddo	
			  
             do j=Jstr_nh,Jend_nh+1		   
               do i=Istr_nh,Iend_nh
                 FY(i,j)=om_v(i,j)* qdmv_nbq(i,j,k)
#ifdef MASKING
                 FY(i,j) = FY(i,j) * vmask(i,j)
#endif
               enddo
              enddo

             do j=Jstr_nh,Jend_nh		   
               do i=Istr_nh,Iend_nh			   
                 div_nbq(i,j,k)=(div_nbq(i,j,k)                                    &
				                 +WORK(i,j)*(FX(i+1,j)-FX(i,j)+FY(i,j+1)-FY(i,j))  &
								 )*Hzr_half_nbq_inv(i,j,k)
                 rho_nbq(i,j,k) = rho_nbq(i,j,k)  - dtnbq * div_nbq(i,j,k)  
               enddo
             enddo
           enddo

!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
!        call rnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!*******************************************************************
!*******************************************************************
      enddo    ! NBQ loop

!*******************************************************************
!*******************************************************************
!
!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
      call ruijk_nbq(2)

#ifdef RVTK_DEBUG
!       call check_tab2d(rubar_nbq(:,:),'rubar_nbq step3d_nbq','u')
!       call check_tab2d(rvbar_nbq(:,:),'rvbar_nbq step3d_nbq','v')
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif      
      call densityijk_nbq(20)

   
      end subroutine step3d_fbijk_nbq

#else
      subroutine step3d_fbijk_nbq_empty
      end subroutine step3d_fbijk_nbq_empty
#endif
