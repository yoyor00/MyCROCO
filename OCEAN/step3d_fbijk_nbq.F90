#include "cppdefs.h"
# define NBQ_IMP
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
      subroutine step3d_fbijk_nbq (Istr,Iend,Jstr,Jend, WORK,     &
	    Hzw_half_nbq_inv,Hzr_half_nbq_inv,                    &
		Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v,            &
        zwrk1,zwrk2,zwrk3,zwrk4,zwrk5,FC,FX,                      &
        Hzu_half_qdmu,Hzv_half_qdmv                               &
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
      real Hzu_half_qdmu(PRIVATE_2D_SCRATCH_ARRAY,N)
      real Hzv_half_qdmv(PRIVATE_2D_SCRATCH_ARRAY,N)
      real DC(PRIVATE_1D_SCRATCH_ARRAY,N)
      real CF(PRIVATE_1D_SCRATCH_ARRAY,N)
	   
# ifdef MPI
      include 'mpif.h'
# endif

      real :: dum_s
      double precision :: a_m,b_m
	  integer :: k1, k2, kp1
      real :: cff1, cff2, cff3, cff

#include "compute_auxiliary_bounds.h"

# undef DEBUG
	   
!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
        if (iif==1.and.iic==1) call initial_nh_tile (3,Istr,Iend,Jstr,Jend)

        if (iic==1.and.(iif==1)) div_nbq=0.
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equations:
!  ru+rubar (or lambda_ext+lambda_int)
!  dzdt*rhosurf
!-------------------------------------------------------------------
!
!
!------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!
       do j=Jstr_nh,Jend_nh
         do i=Istr_nh,Iend_nh
			WORK(i,j)=pm(i,j)*pn(i,j)
		 enddo
       enddo

!*******************************************************************
!*******************************************************************
!              Initialize tendencies
!*******************************************************************
!*******************************************************************

         do k=1,N         
          do j=JstrU_nh,JendU_nh   
            do i=IstrU_nh,IendU_nh
              ru_nbq_ext (i,j,k) = qdmu_nbq(i,j,k)      
            enddo
          enddo
         enddo
         
         do k=1,N         
          do j=JstrV_nh,JendV_nh   
            do i=IstrV_nh,IendV_nh
              rv_nbq_ext (i,j,k) = qdmv_nbq(i,j,k)      
            enddo
          enddo
        enddo         
               
#ifdef M2FILTER_NONE
        if (LAST_2D_STEP) then
#endif
        do k=0,N 
          do j=Jstr_nh,Jend_nh             
            do i=Istr_nh,Iend_nh
              rw_nbq_ext (i,j,k) = qdmw_nbq(i,j,k)
            enddo
          enddo
        enddo
#ifdef M2FILTER_NONE        
        endif
#endif        
       
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

        do k=1,N
          do j=JstrV-1,Jend
            do i=IstrU-1,Iend
              div_nbq(i,j,k)=-visc2_nbq*div_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)      &
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
		  if (k.eq.0) then ! Bottom Boundary conditions
            do j=Jstr,Jend
              do i=IstrU,Iend
			    ddiv_nbqdz_u(i,j,k2)=0.
			  enddo
			enddo
            do j=JstrV,Jend
              do i=Istr,Iend
			    ddiv_nbqdz_v(i,j,k2)=0.
			  enddo
			enddo
		  else
            if (k.eq.N) then ! Top Boundary conditions
              do j=JstrV-1,Jend
                do i=IstrU-1,Iend
			      ddiv_nbqdz(i,j)= - div_nbq(i  ,j,k)
			    enddo
			  enddo
            else
              do j=JstrV-1,Jend
                do i=IstrU-1,Iend
                  ddiv_nbqdz(i,j)=div_nbq(i  ,j,k+1) - div_nbq(i  ,j,k)
                enddo
              enddo            
            endif
            do j=Jstr,Jend
              do i=IstrU,Iend
               ddiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(ddiv_nbqdz(i,j)+ddiv_nbqdz(i-1,j))              
              enddo
            enddo
            do j=JstrV,Jend
              do i=Istr,Iend
                ddiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(ddiv_nbqdz(i,j)+ddiv_nbqdz(i,j-1))
              enddo
            enddo            
            endif
            if (k.gt.0) then
            do j=Jstr,Jend
              do i=IstrU,Iend
			    dum_s=(zr_half_nbq(i,j,k)-zr_half_nbq(i-1,j,k)) &
                       *(ddiv_nbqdz_u(i,j,k2)+ddiv_nbqdz_u(i,j,k1)) &   ! dZdx * (d(delta p)dz)_u
				         -(div_nbq(i,j,k)-div_nbq(i-1,j,k))                           ! - d(delta p)dx

                dum_s=dum_s*Hzu_half_qdmu(i,j,k)
                qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k) + dtnbq * ( dum_s + ruint_nbq(i,j,k))               
              enddo
            enddo
            do j=JstrV,Jend
              do i=Istr,Iend
			    dum_s=(zr_half_nbq(i,j,k)-zr_half_nbq(i,j-1,k)) &
                       *(ddiv_nbqdz_v(i,j,k2)+ddiv_nbqdz_v(i,j,k1)) &   ! dZdy * (d(delta p)dz)_v
				         -(div_nbq(i,j,k)-div_nbq(i,j-1,k))                           ! - d(delta p)dy

                dum_s=dum_s*Hzv_half_qdmv(i,j,k)
                qdmv_nbq(i,j,k) = qdmv_nbq(i,j,k) + dtnbq * ( dum_s + rvint_nbq(i,j,k))					
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
       call unbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
       call vnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
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

!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
!

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,qdmv_nbq(START_2D_ARRAY,1))
#endif

#ifdef RVTK_DEBUG
       call check_tab3d(qdmu_nbq,'qdmu_nbq','u')
       call check_tab3d(qdmv_nbq,'qdmv_nbq','v')
#endif  
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
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rwint_nbq(i,j,k) )
#if defined MASKING
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#endif            
             enddo             
           enddo
          k=N
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   div_nbq(i,j,N)                              
               qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N)   &
                + dtnbq * ( dum_s + rwint_nbq(i,j,N) )
#if defined MASKING
                qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N) * rmask(i,j)
#endif               
             enddo             		   
          k=0
            do i=Istr_nh,Iend_nh                                                               
               dum_s =  -div_nbq(i,j,1)                              
               qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0)   &
                + dtnbq * ( dum_s + rwint_nbq(i,j,0) )
#if defined MASKING
                qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#endif               
             enddo  
        enddo

# endif
!
!      Vertical momentum open boundary conditions
!
# ifdef OBC_NBQ
        call wnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
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

#define dZdxq_u zwrk1
#define dZdxq_w zwrk2
#define dZdyq_v zwrk3
#define dZdyq_w zwrk4
#define FY zwrk5
 
! X -component 
        k2 = 1
		do k=0,N
		  k1=k2
		  k2=3-k1

		  if (k.lt.N) then
             kp1 = k + 1
		     do j=Jstr_nh,Jend_nh
		       do i=Istr_nh,Iend_nh+1
		         dZdxq_u(i,j,k2)=(zr_half_nbq(i,j,kp1)-zr_half_nbq(i-1,j,kp1)) &
                       *qdmu_nbq(i,j,kp1)    ! (dZdx * (rho u))_u
		       enddo
		     enddo	 
		  endif

		  if (k.eq.0) then	! Bottom boundary conditions	  
		   do j=Jstr_nh,Jend_nh
		     do i=Istr_nh,Iend_nh+1
                dZdxq_w(i,j,k2)=0.  
!                dZdxq_w(i,j,k2)= - (zw_half_nbq(i,j,0)-zw_half_nbq(i-1,j,0))*qdmu_nbq(i,j,1)  &
!                                 / (hzr_half_nbq(i,j,1)+hzr_half_nbq(i-1,j,1))            
		     enddo
		   enddo
          elseif (k==N) then ! Top boundary conditions
		   do j=Jstr_nh,Jend_nh
		     do i=Istr_nh,Iend_nh+1
		       dZdxq_w(i,j,k2)= (zw_half_nbq(i,j,N)-zw_half_nbq(i-1,j,N))   &
                                *qdmu_nbq(i,j,N)                            &                       
                              / (hzr_half_nbq(i,j,N)+hzr_half_nbq(i-1,j,N))
		     enddo
		   enddo     
          else
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh+1
			     dZdxq_w(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) ! (dZdx * (rho u))_uw/Hzw_u
              enddo 
            enddo
          endif

          if (k.gt.0) then
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh+1
			    FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
                
#if defined MASKING
                FX(i,j) = FX(i,j) * umask(i,j)
#endif                
              enddo
            enddo
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh
			    div_nbq(i,j,k)=FX(i,j)+FX(i+1,j)           
              enddo
            enddo
          endif
		enddo		 

! Y component        
        k2 = 1
		do k=0,N
		  k1=k2
		  k2=3-k1

		  if (k.lt.N) then
             kp1 = k + 1
		     do j=Jstr_nh,Jend_nh+1
		       do i=Istr_nh,Iend_nh
		         dZdyq_v(i,j,k2)=(zr_half_nbq(i,j,kp1)-zr_half_nbq(i,j-1,kp1)) &
                       *qdmv_nbq(i,j,kp1)    ! (dZdy * (rho v))_v
		       enddo
		     enddo			 
		  endif

		  if (k.eq.0) then	! Bottom boundary conditions
		   do j=Jstr_nh,Jend_nh+1
		     do i=Istr_nh,Iend_nh
               dZdyq_w(i,j,k2)=0.   
!               dZdyq_w(i,j,k2)= - (zw_half_nbq(i,j,0)-zw_half_nbq(i,j-1,0))*qdmv_nbq(i,j,1) &
!		                       / ( Hzr_half_nbq(i,j,1)+Hzr_half_nbq(i,j-1,1) )               
		     enddo
		   enddo
          elseif (k==N) then ! Top boundary conditions
		   do j=Jstr_nh,Jend_nh+1
		     do i=Istr_nh,Iend_nh
		       dZdyq_w(i,j,k2)= (zw_half_nbq(i,j,N)-zw_half_nbq(i,j-1,N))*qdmv_nbq(i,j,N) &
		                       / ( Hzr_half_nbq(i,j,N)+Hzr_half_nbq(i,j-1,N) )
		     enddo
		   enddo
          else
			do j=Jstr_nh,Jend_nh+1
		      do i=Istr_nh,Iend_nh
			     dZdyq_w(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
              enddo 
            enddo
          endif

          if (k.gt.0) then
		    do j=Jstr_nh,Jend_nh+1
		      do i=Istr_nh,Iend_nh
			    FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
#if defined MASKING
                FY(i,j) = FY(i,j) * vmask(i,j)
#endif                 
              enddo
            enddo			
		    do j=Jstr_nh,Jend_nh
		      do i=Istr_nh,Iend_nh
			    div_nbq(i,j,k)=div_nbq(i,j,k)+FY(i,j)+FY(i,j+1)              
              enddo
            enddo
          endif
		enddo		         

#undef dZdxq_u
#undef dZdxq_w
#undef FY

#define FY zwrk5

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
								 )
#ifdef MASKING
                 div_nbq(i,j,k) = div_nbq(i,j,k) * rmask(i,j)
#endif                              
               enddo
             enddo
           enddo

#ifdef NBQ_IMP
           
        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh
              FC(i,k)= soundspeed2_nbq*(rho_nbq(i,j,k)  - dtnbq * div_nbq(i,j,k) &
                     *Hzr_half_nbq_inv(i,j,k) )
            enddo
          enddo            
          do k=1,N-1
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   FC(i,k) - FC(i,k+1)            
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rwint_nbq(i,j,k) )
#if defined MASKING
                qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#endif               
             enddo             
           enddo
          k=N
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   FC(i,k)                              
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rwint_nbq(i,j,k) )
#if defined MASKING
                qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#endif              
             enddo             		   
        enddo
 
! Tridiagonal inversion 

        cff1=1./(soundspeed2_nbq*dtnbq**2)

        do j=Jstr_nh,Jend_nh
            do i=Istr_nh,Iend_nh
             cff=1.d0/(cff1+Hzw_half_nbq_inv(i,j,1)*(Hzr_half_nbq_inv(i,j,1)+Hzr_half_nbq_inv(i,j,2)))
             CF(i,1)=cff*(-Hzw_half_nbq_inv(i,j,2)*Hzr_half_nbq_inv(i,j,2))
             DC(i,1)=cff*qdmw_nbq(i,j,1)*cff1            
            enddo
            do k=2,N-1
              do i=Istr_nh,Iend_nh
                cff=1.d0/(cff1+                                                                   &
                    Hzw_half_nbq_inv(i,j,k)*(Hzr_half_nbq_inv(i,j,k)+Hzr_half_nbq_inv(i,j,k+1))   &
                    +Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)*CF(i,k-1))
                CF(i,k)=cff*(-Hzw_half_nbq_inv(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
                DC(i,k)=cff*(qdmw_nbq(i,j,k)*cff1+Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)*DC(i,k-1))             
              enddo            
            enddo
            k=N
            do i=Istr_nh,Iend_nh
              cff=1.d0/(cff1+Hzw_half_nbq_inv(i,j,k)*Hzr_half_nbq_inv(i,j,k) &
                        +Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)*CF(i,k-1))
              DC(i,k)=cff*(qdmw_nbq(i,j,k)*cff1+Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)*DC(i,k-1))             
            enddo 
            do i=Istr_nh,Iend_nh
              qdmw_nbq(i,j,k)=DC(i,k)           
            enddo
            do k=N-1,1,-1
              do i=Istr_nh,Iend_nh
                qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
              enddo            
            enddo                       
        enddo           
        
#endif
			   
           do j=Jstr_nh,Jend_nh
             do i=Istr_nh,Iend_nh
			   FC(i,0)=0.                 ! Bottom boundary condition
             enddo
		     do k=1,N
               do i=Istr_nh,Iend_nh
			     FC(i,k)=Hzw_half_nbq_inv(i,j,k) * qdmw_nbq(i,j,k)
               enddo
               do i=Istr_nh,Iend_nh           
			     div_nbq(i,j,k)=(div_nbq(i,j,k)+FC(i,k)-FC(i,k-1))             
               enddo
             enddo
           enddo
           
           do k=1,N
             do j=Jstr_nh,Jend_nh
               do i=Istr_nh,Iend_nh
                 rho_nbq(i,j,k) = rho_nbq(i,j,k)  - dtnbq * div_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k) 
               enddo
             enddo
           enddo
                 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,div_nbq(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,rho_nbq(START_2D_ARRAY,1))
# endif
           
#ifdef RVTK_DEBUG
       call check_tab3d(rho_nbq,'rho_nbq','r')
#endif           

!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
         call rnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
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
      call ruijk_nbq(2, Istr,Iend,Jstr,Jend,WORK)

#ifdef RVTK_DEBUG
       call check_tab2d(rubar_nbq,'rubar_nbq step3d_nbq','uint')
       call check_tab2d(rvbar_nbq,'rvbar_nbq step3d_nbq','vint')
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif      
      call densityijk_nbq(20)

   
      end subroutine step3d_fbijk_nbq

#else
      subroutine step3d_fbijk_nbq_empty
      end subroutine step3d_fbijk_nbq_empty
#endif
