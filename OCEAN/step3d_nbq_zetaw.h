!
! Resolve 3D NBQ equations
! This provides NBQ-filtered rhs terms for the barotropic equation.
! Then perform 2D filter of these terms to get forcing terms for 
! the total equation.
!
!*******************************************************************
!*******************************************************************
!  Indices & init
!*******************************************************************
!*******************************************************************
!

      if (istr.eq.1) then
        IstrR2=Istr-1
#ifdef MPI
        IstrU2=Istr  !+1
#else
        IstrU2=Istr +1
#endif
      else
        IstrR2=Istr
        IstrU2=Istr
      endif

      if (iend.eq.Lm) then
        IendR2=Iend+1
      else
        IendR2=Iend
      endif

      if (jstr.eq.1) then
        JstrR2=Jstr-1
#ifdef MPI
        JstrV2=Jstr  !+1
#else
        JstrV2=Jstr  +1
#endif
      else
        JstrR2=Jstr
        JstrV2=Jstr
      endif

      if (jend.eq.Mm) then
        JendR2=Jend+1
      else
        JendR2=Jend
      endif

#ifdef NBQ_DTDRHO2
       if (iic==1.and.iif==1) then
	  do j=JstrR2,JendR2
# ifdef NBQ_MASS
	    do k=1,N
	    do i=IstrR2,IendR2
	       zr_nbq(i,j,k,:)=z_r(i,j,k)
            enddo
            enddo
# endif
   	    do k=0,N
	    do i=IstrR2,IendR2
	       z_nbq (i,j,k,:)=z_w(i,j,k)
            enddo
            enddo
         enddo
       endif
#endif

# undef DEBUG
!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
!       if (iif==1.and.iic==1) call initial_nh_tile (3,Istr,Iend,Jstr,Jend)
!       if (iic==1.and.iif==1) thetadiv_nbq=0.
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
       do j=Jstr,Jend
         do i=Istr,Iend
	    WORK(i,j)=pm(i,j)*pn(i,j)
	 enddo
       enddo      
       do j=JstrU_nh,JendU_nh   
          do i=IstrU_nh,IendU_nh
               rubar_nbq(i,j) = 0 
          enddo
       enddo   
       do j=JstrV_nh,JendV_nh   
          do i=IstrV_nh,IendV_nh
               rvbar_nbq(i,j) = 0 
          enddo
       enddo
       do j=JstrU_nh,JendU_nh   
          do i=IstrU_nh,IendU_nh
               DU_nbq(i,j)=0.
          enddo
       enddo   
       do j=JstrV_nh,JendV_nh   
          do i=IstrV_nh,IendV_nh
               DV_nbq(i,j)=0.
          enddo
       enddo


#if defined NBQ_DTDRHO2 && defined NBQ_ZETAW && defined NBQ_MASS
       if (iic==1.and.iif==1) then
        do k=1,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
		 rho_bak(i,j,k)=rho(i,j,k)
	    enddo
	  enddo
	 enddo
       endif
#endif

	
!*******************************************************************
!*******************************************************************
!              Stores tendencies
!*******************************************************************
!*******************************************************************

        do k=0,N 
          do j=Jstr,Jend             
            do i=Istr,Iend
               rw_nbq_ext (i,j,k) = qdmw_nbq(i,j,k) 
            enddo
          enddo
        enddo
   
!*******************************************************************
!*******************************************************************
!              NBQ mode iteration (main loop)
!*******************************************************************
!*******************************************************************
!
!      do iteration_nbq=1,iteration_nbq_max
!
!-------------------------------------------------------------------
!       "Pressure - Viscosity" Variable (theta)
!               theta does not change
!-------------------------------------------------------------------
!       

        do k=1,N
          do j=JstrV2-1,Jend
            do i=IstrU2-1,Iend
              thetadiv_nbq(i,j,k)=(-visc2_nbq*thetadiv_nbq(i,j,k)
     &                                +soundspeed2_nbq*rho_nbq(i,j,k)) 
     &                               *Hzr_half_nbq_inv(i,j,k)  
            enddo
          enddo
        enddo

!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------

!---------------------------
!  XI- and ETA-Directions:
!---------------------------

!# define dthetadiv_nbqdz_u zwrk1
!# define dthetadiv_nbqdz_v zwrk2

#ifndef NBQ_NODS
!# define dthetadiv_nbqdz   zwrk5
#endif

        k2 = 1
        do k=0,N
          k1=k2
	  k2=3-k1

# ifdef NBQ_NODS
          if (NSTEP_DS) then
# endif
			
          if (k.eq.0) then ! Bottom Boundary conditions
	
	    do j=Jstr,Jend
              do i=IstrU2,Iend
	        dthetadiv_nbqdz_u(i,j,k2)=0. 
	      enddo
	    enddo

            do j=JstrV2,Jend
              do i=Istr,Iend
	        dthetadiv_nbqdz_v(i,j,k2)=0.
	      enddo
    	    enddo 

          else

            if (k.eq.N) then ! Top Boundary conditions

              do j=JstrV2-1,Jend
                do i=IstrU2-1,Iend
# ifndef NBQ_NODS
                  dthetadiv_nbqdz(i,j)    = - thetadiv_nbq(i  ,j,k)
# else
                  dthetadiv_nbqdz(i,j,k,1)= - thetadiv_nbq(i  ,j,k)
# endif
     	        enddo
	      enddo

            else

              do j=JstrV2-1,Jend
                do i=IstrU2-1,Iend
# ifndef NBQ_NODS
                  dthetadiv_nbqdz(i,j)    =thetadiv_nbq(i  ,j,k+1) 
     &             - thetadiv_nbq(i  ,j,k)
# else
                  dthetadiv_nbqdz(i,j,k,1)=thetadiv_nbq(i  ,j,k+1)
     &              - thetadiv_nbq(i  ,j,k)
# endif
                enddo
              enddo

            endif
  
            do j=Jstr,Jend
            do i=IstrU2,Iend
# ifndef NBQ_NODS
              dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &            dthetadiv_nbqdz(i,j)
     &           +dthetadiv_nbqdz(i-1,j)) 
# else
              dthetadiv_nbqdz_u(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dthetadiv_nbqdz(i,j,k,1)
     &          +dthetadiv_nbqdz(i-1,j,k,1))    
# endif          
            enddo
            enddo
            do j=JstrV2,Jend
              do i=Istr,Iend
# ifndef NBQ_NODS
                dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &         dthetadiv_nbqdz(i,j)
     &        +dthetadiv_nbqdz(i,j-1))
# else
                dthetadiv_nbqdz_v(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)*(
     &       dthetadiv_nbqdz(i,j,k,1)
     &      +dthetadiv_nbqdz(i,j-1,k,1))
# endif
              enddo
            enddo    
  
            endif    
# ifdef NBQ_NODS
          endif
# endif          

          if (k.gt.0) then

!...........U-momentum:
            do j=Jstr,Jend
              do i=IstrU2,Iend
                if (k.gt.1.and.k.lt.N) then 
# ifndef NBQ_NODS
	  	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                  *(dthetadiv_nbqdz_u(i,j,k2)+dthetadiv_nbqdz_u(i,j,k1))             ! dZdx * (d(delta p)dz)_u
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                         ! - d(delta p)dx  
# else
                  if (NSTEP_DS) then
	   	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &                    -z_r(i-1,j,k))       
     &                  *(dthetadiv_nbqdz_u(i,j,k2)+dthetadiv_nbqdz_u(i,j,k1))              ! dZdx * (d(delta p)dz)_u
                  endif 
                  dum_s=dthetadiv_nbqdz(i,j,k,1)                                        
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                         ! - d(delta p)dx  
# endif
                elseif (k.gt.1) then
# ifndef NBQ_NODS
	 	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                        *dthetadiv_nbqdz_u(i,j,k1)                                        ! dZdx * (d(delta p)dz)_u
     &	               -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                       ! - d(delta p)dx
     &                    +(z_w(i,j,N)-z_w(i-1,j,N))                      
     &                  *dthetadiv_nbqdz_u(i,j,k2)
# else
                  if (NSTEP_DS) then
	 	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &        -z_r(i-1,j,k))       
     &                   *dthetadiv_nbqdz_u(i,j,k1)                                       ! dZdx * (d(delta p)dz)_u
     &                  +(z_w(i,j,N)-z_w(i-1,j,N))                      
     &                  *dthetadiv_nbqdz_u(i,j,k2)                 
                  endif
	 	  dum_s=dthetadiv_nbqdz(i,j,k,1)                                          
     &	               -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                  
# endif
                else
# ifndef NBQ_NODS
	  	  dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                  *2.*dthetadiv_nbqdz_u(i,j,k2)                                      ! dZdx * (d(delta p)dz)_u
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                         ! - d(delta p)dx  
# else
                  if (NSTEP_DS) then
  	  	  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)
     &                     -z_r(i-1,j,k))       
     &                    *2.*dthetadiv_nbqdz_u(i,j,k2)                                    ! dZdx * (d(delta p)dz)_u
                  endif
	  	  dum_s=dthetadiv_nbqdz(i,j,k,1)                                          
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i-1,j,k))                         ! - d(delta p)dx  
# endif
                endif
                dum_s=dum_s   *Hzu_half_qdmu(i,j,k)
                qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k) + dtnbq * (
     &              dum_s + ru_int_nbq(i,j,k))  
                DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
                ru_nbq_ext (i,j,k) = dum_s / work(i,j) 
                rubar_nbq(i,j)=rubar_nbq(i,j)+ru_nbq_ext(i,j,k)
              enddo 
            enddo

!...........V-momentum:
            do j=JstrV2,Jend
              do i=Istr,Iend
                if (k.gt.1.and.k.lt.N) then 
# ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                     *(dthetadiv_nbqdz_v(i,j,k2)
     &                  +dthetadiv_nbqdz_v(i,j,k1))    ! dZdy * (d(delta p)dz)_v
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                ! - d(delta p)dy
# else
                  if (NSTEP_DS) then
	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &    -z_r(i,j-1,k)) 
     &                  *(dthetadiv_nbqdz_v(i,j,k2)+dthetadiv_nbqdz_v(i,j,k1))       ! dZdy * (d(delta p)dz)_v
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2)                         
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                ! - d(delta p)dy
# endif
                elseif (k.gt.1) then
# ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k))        
     &                  *dthetadiv_nbqdz_v(i,j,k1)                           ! dZdy * (d(delta p)dz)_v
     &	               -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                 ! - d(delta p)dy
     &                   +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                  *dthetadiv_nbqdz_v(i,j,k2)
# else
                  if (NSTEP_DS) then
     	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &                    -z_r(i,j-1,k)) 
     &                   *dthetadiv_nbqdz_v(i,j,k1)                               ! dZdy * (d(delta p)dz)_v
     &                  +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                  *dthetadiv_nbqdz_v(i,j,k2)
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2) 
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                ! - d(delta p)dy
# endif
                else
# ifndef NBQ_NODS
	          dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                  *2.*dthetadiv_nbqdz_v(i,j,k2)                             ! dZdy * (d(delta p)dz)_v
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                ! - d(delta p)dy
# else
                  if (NSTEP_DS) then
	          dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)
     &                      -z_r(i,j-1,k)) 
     &                  *2.*dthetadiv_nbqdz_v(i,j,k2)       ! dZdy * (d(delta p)dz)_v
                  endif
	          dum_s=dthetadiv_nbqdz(i,j,k,2) 
     &	              -(thetadiv_nbq(i,j,k)-thetadiv_nbq(i,j-1,k))                ! - d(delta p)dy
# endif
                endif
                
                dum_s=dum_s*Hzv_half_qdmv(i,j,k)
                qdmv_nbq(i,j,k) = qdmv_nbq(i,j,k) + dtnbq * (
     &                  dum_s + rv_int_nbq(i,j,k))    
                DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
                rv_nbq_ext (i,j,k) = dum_s / work(i,j)  
                rvbar_nbq(i,j)=rvbar_nbq(i,j)+rv_nbq_ext(i,j,k)	
              enddo
            enddo
          endif
        enddo        
        
!# undef dthetadiv_nbqdzu_
!# undef dthetadiv_nbqdz_v
!#ifndef NBQ_NODS
!# undef dthetadiv_nbqdz
!#endif

!---------------------------
!  U-momentum open boundary conditions
!---------------------------
       
# ifdef OBC_NBQ
       call unbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
       call vnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif


!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
!

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,qdmv_nbq(START_2D_ARRAY,1))
# endif
# ifdef RVTK_DEBUG
      call check_tab3d(qdmu_nbq,'qdmu_nbq','u')
      call check_tab3d(qdmv_nbq,'qdmv_nbq','v')
# endif  
!-------------------------------------------------------------------
!      Explicit Vertical Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
# ifndef NBQ_IMP
!---------------------------
!  Z-Direction: Explicit
!---------------------------
        do j=Jstr,Jend
          do k=1,N-1
            do i=Istr,Iend                                                               
               dum_s =   thetadiv_nbq(i,j,k) - thetadiv_nbq(i,j,k+1)   
               qdmw_nbq(i,j,k)   = qdmw_nbq(i,j,k)   
     &            + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
#ifdef NBQ_GRAV
     &            -0.25*(rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)
     &                 +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
     &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
     &                 *g*dtnbq
#endif
#ifdef MASKING
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
#endif            
            enddo             
          enddo
          k=N
          do i=Istr,Iend                                                               
               dum_s =   thetadiv_nbq(i,j,N)                              
               qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N)   
     &           + dtnbq * ( dum_s + rw_int_nbq(i,j,N) )
#ifdef NBQ_GRAV
     &            -rho_nbq(i,j,N)*0.5
     &                 *g*dtnbq
#endif
#ifdef MASKING
               qdmw_nbq(i,j,N) = qdmw_nbq(i,j,N) * rmask(i,j) 
#endif               
          enddo     
        		   
! Bottom boundary:        
!          do i=Istr,Iend                                                               
!                dum_s =  0.  !-thetadiv_nbq(i,j,1)                              
!                qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0)   &
!                 + dtnbq * ( dum_s + rw_int_nbq(i,j,0) )
!#ifdef MASKING
!                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
!#endif               
!                qdmw_nbq(i,j,0) = 0.
!           enddo  

         enddo

!---------------------------
! Vertical momentum open boundary conditions
!---------------------------
# ifdef OBC_NBQ
        call wnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,qdmw_nbq(START_2D_ARRAY,0))
# endif
# endif  /* NBQ_IMP */


!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
!#  ifdef ACOUSTIC
!      call density_nbq(11)       ! TBD
!#  endif
!
!-------------------------------------------------------------------
!      Mass equation (1): DX(p+1) ==> thetadiv_nbq
!-------------------------------------------------------------------
!

!# define dZdxq_u zwrk1
!# define dZdyq_v zwrk3
!#ifndef NBQ_NODS
!# define dZdxq_w zwrk2
!# define dZdyq_w zwrk4
!#endif
!#define FY zwrk5

!--------------------------- 
! X -component 
!---------------------------

        k2 = 1
        do k=0,N
          k1=k2
	  k2=3-k1

# ifdef NBQ_NODS
          if (NSTEP_DS) then
#endif

          if (k.lt.N) then
             kp1 = k + 1
	     do j=Jstr,Jend
             do i=Istr,Iend+1
	       dZdxq_u(i,j,k2)=(z_r(i,j,kp1)-z_r(i-1,j,kp1)) 
     &                  *qdmu_nbq(i,j,kp1)    ! (dZdx * (rho u))_u
             enddo
             enddo
          endif

	  if (k.eq.0) then	! Bottom boundary conditions

#  if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend
	    do i=Istr,Iend+1 
# ifndef NBQ_NODS
              dZdxq_w(i,j,k2)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                          *qdmu_nbq(i,j,1)  
     &                          / (Hzr(i,j,1)+Hzr(i-1,j,1))
# else
              dZdxq_w(i,j,k) = (z_w(i,j,0)-z_w(i-1,j,0))
     &                        *qdmu_nbq(i,j,1)  
     &                          / (Hzr(i,j,1)+Hzr(i-1,j,1))
# endif
	    enddo
	    enddo

 	    do j=Jstr,Jend
 	    do i=Istr,Iend    
# ifndef NBQ_NODS     
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k2) *pm_u(i,j) 
     &                +dZdxq_w(i+1,j,k2) *pm_u(i+1,j) ) 
     &                        * Hzr(i,j,1)      
# else
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k) *pm_u(i,j) 
     &                         +dZdxq_w(i+1,j,k) *pm_u(i+1,j) ) 
     &                        * Hzr(i,j,1)      
# endif

# if defined MASKING
              qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
# endif 
 	    enddo
 	    enddo 
#  else 
#   ifndef NBQ_NODS
 	    do j=Jstr,Jend
 	    do i=Istr,Iend +1
              dZdxq_w(i,j,k2)=0.  
              qdmw_nbq(i,j,0)=0.     
 	    enddo
 	    enddo 
!#   else	 
! 	    do j=Jstr,Jend
! 	    do i=Istr,Iend +1  
!              dZdxq_w(i,j,k)=0.     
 !             qdmw_nbq(i,j,0)=0.  
 !	    enddo
 !	    enddo 
#   endif	  
#  endif 

          elseif (k==N) then ! Top boundary conditions
           
            do j=Jstr,Jend
	    do i=Istr,Iend+1
# ifndef NBQ_NODS
	      dZdxq_w(i,j,k2)= (z_w(i,j,N)-z_w(i-1,j,N))   
     &                 *qdmu_nbq(i,j,N)                                                            
     &                 / (Hzr(i,j,N)+Hzr(i-1,j,N)) 
# else
	      dZdxq_w(i,j,k)= (z_w(i,j,N)-z_w(i-1,j,N))   
     &                 *qdmu_nbq(i,j,N)                                                            
     &                 / (Hzr(i,j,N)+Hzr(i-1,j,N))
# endif  
            enddo
            enddo  

!#  ifdef NBQ_SBBC
!            do j=Jstr,Jend
!	    do i=Istr,Iend
!# ifndef NBQ_NODS
!              qdmw_nbq(i,j,N+1)=qdmw_nbq(i,j,N+1)+0.5*(dZdxq_w(i,j,k2)+dZdxq_w(i+1,j,k2))
!# else
!              qdmw_nbq(i,j,N+1)=qdmw_nbq(i,j,N+1)+0.5*(dZdxq_w(i,j,k)+dZdxq_w(i+1,j,k))
!# endif
!#if defined MASKING
!              qdmw_nbq(i,j,N+1) = qdmw_nbq(i,j,N+1) * rmask(i,j)
!#endif 
!           enddo
!           enddo   
!#  endif
 
          else

            do j=Jstr,Jend
	    do i=Istr,Iend+1
# ifndef NBQ_NODS
	       dZdxq_w(i,j,k2)=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) 
# else
	       dZdxq_w(i,j,k )=Hzw_half_nbq_inv_u(i,j,k)*(
     &           dZdxq_u(i,j,k1)+dZdxq_u(i,j,k2)) 
# endif
            enddo 
            enddo

          endif

# ifdef NBQ_NODS

          else

	  if (k.eq.0) then	! Bottom boundary conditions

#  if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend
	    do i=Istr,Iend+1 
              dZdxq_w(i,j,k)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                           * qdmu_nbq(i,j,1)  
     &                           / (Hzr(i,j,1)+Hzr(i-1,j,1))
	    enddo
	    enddo

 	    do j=Jstr,Jend
 	    do i=Istr,Iend    
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i,j,k) *pm_u(i,j) +dZdxq_w(i+1,j,k) 
     &                         * pm_u(i+1,j) ) 
     &                         * Hzr(i,j,1)    

# if defined MASKING
              qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
# endif 
 	    enddo
 	    enddo 
!#  else 
! 	    do j=Jstr,Jend
! 	    do i=Istr,Iend    
!              dZdxq_w(i,j,k)=0.     
!              qdmw_nbq(i,j,0)=0.  
! 	    enddo
! 	    enddo 
#  endif 
          endif
          endif
# endif
          if (k.gt.0) then
           if (IstrU.le.Iend) then
            do j=Jstr,Jend
	    do i=Istr,Iend+1
# ifndef NBQ_NODS
	      FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
# else
	      FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k)-dZdxq_w(i,j,k-1))
# endif

#ifdef MASKING
              FX(i,j) = FX(i,j) * umask(i,j)
#endif                
            enddo
            enddo

            do j=Jstr,Jend
            do i=Istr,Iend
	      thetadiv_nbq(i,j,k)=FX(i,j)  +FX(i+1,j)           
            enddo
            enddo
           endif 
          endif
	enddo	
	 
!---------------------------
! Y component     
!---------------------------   

        k2 = 1
	do k=0,N
	  k1=k2
	  k2=3-k1

# ifdef NBQ_NODS
          if (NSTEP_DS) then
#endif

	  if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend+1
            do i=Istr,Iend
	      dZdyq_v(i,j,k2)=(z_r(i,j,kp1)-z_r(i,j-1,kp1)) 
     &                   *qdmv_nbq(i,j,kp1)    ! (dZdy * (rho v))_v
            enddo
	    enddo			 
          endif

          if (k.eq.0) then	! Bottom boundary conditions

#  if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend+1
            do i=Istr,Iend
# ifndef NBQ_NODS
               dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                     *qdmv_nbq(i,j,1) 
     &                     / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
# else
               dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                      *qdmv_nbq(i,j,1) 
     &                      / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
# endif 
	    enddo
	    enddo
 	    do j=Jstr,Jend
  	    do i=Istr,Iend   
# ifndef NBQ_NODS     
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k2)*pm_v(i,j)  
     &    +dZdyq_w(i,j+1,k2)*pm_v(i,j+1)  )    
     &                               * Hzr(i,j,1) 
# else
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k)*pm_v(i,j)  
     &       +dZdyq_w(i,j+1,k)*pm_v(i,j+1)  )    
     &                               * Hzr(i,j,1) 
# endif 
#if defined MASKING
                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#endif 
  	    enddo
            enddo
#  else  
# ifndef NBQ_NODS     
 	    do j=Jstr,Jend +1
 	    do i=Istr,Iend  
              dZdyq_w(i,j,k2)=0.  
              qdmw_nbq(i,j,0)=0.
 	    enddo
 	    enddo 
!#   else	 
! 	    do j=Jstr,Jend +1
! 	    do i=Istr,Iend    
!              dZdyq_w(i,j,k)=0.     
 !             qdmw_nbq(i,j,0)=0.  
 !	    enddo
 !	    enddo 
#  endif
#  endif

          elseif (k==N) then ! Top boundary conditions

            do j=Jstr,Jend+1
	    do i=Istr,Iend
# ifndef NBQ_NODS
              dZdyq_w(i,j,k2)= (z_w(i,j,N)-z_w(i,j-1,N))       
     &                        * qdmv_nbq(i,j,N)                                 
     &                        / ( Hzr(i,j,N)+Hzr(i,j-1,N) )
# else
              dZdyq_w(i,j,k)= (z_w(i,j,N)-z_w(i,j-1,N))       
     &                        * qdmv_nbq(i,j,N)                                 
     &                        / ( Hzr(i,j,N)+Hzr(i,j-1,N) )
# endif 
	    enddo
	    enddo

#  ifdef NBQ_SBBC
!           do j=Jstr,Jend
!	    do i=Istr,Iend
# ifndef NBQ_NODS
!              qdmw_nbq(i,j,N+1)=(qdmw_nbq(i,j,N+1)+0.5*(dZdyq_w(i,j,k2)+dZdyq_w(i,j+1,k2))) &
!                                * Hzw_half_nbq(i,j,N)
# else
!              qdmw_nbq(i,j,N+1)=(qdmw_nbq(i,j,N+1)+0.5*(dZdyq_w(i,j,k)+dZdyq_w(i,j+1,k)))&
!                                * Hzw_half_nbq(i,j,N)
# endif 
#   if defined MASKING
!              qdmw_nbq(i,j,N+1) = qdmw_nbq(i,j,N+1) * rmask(i,j)
#   endif 
!	    enddo
!	    enddo
#  endif
          else

      	    do j=Jstr,Jend+1
            do i=Istr,Iend
# ifndef NBQ_NODS
              dZdyq_w(i,j,k2)=Hzw_half_nbq_inv_v(i,j,k)
     &    *(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
# else
              dZdyq_w(i,j,k)=Hzw_half_nbq_inv_v(i,j,k)
     &     *(dZdyq_v(i,j,k1)+dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
# endif 
            enddo 
            enddo

          endif

# ifdef NBQ_NODS

          else

          if (k.eq.0) then	! Bottom boundary conditions

#  if defined NBQ_FREESLIP || defined NBQ_SBBC
	    do j=Jstr,Jend+1
            do i=Istr,Iend
               dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                   *qdmv_nbq(i,j,1) 
     &                   / ( Hzr(i,j,1)+Hzr(i,j-1,1) )
	    enddo
	    enddo
 	    do j=Jstr,Jend
  	    do i=Istr,Iend   
                 qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0) 	                             
     &                             +0.5*(dZdyq_w(i,j,k)*pm_v(i,j)  
     &                     +dZdyq_w(i,j+1,k)*pm_v(i,j+1)  )     
     &                               * Hzr(i,j,1) 
#   if defined MASKING
                 qdmw_nbq(i,j,0) = qdmw_nbq(i,j,0) * rmask(i,j)
#   endif 
  	    enddo
            enddo
!#  else  
! 	    do j=Jstr,Jend
! 	    do i=Istr,Iend    
!              dZdyq_w(i,j,k)=0. 
!              qdmw_nbq(i,j,0)=0.
! 	    enddo
! 	    enddo 
#  endif
          endif
          endif
# endif

          if (k.gt.0) then
           if (JstrV.le.Jend) then
            do j=Jstr,Jend+1
            do i=Istr,Iend 
# ifndef NBQ_NODS
	      FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
# else
	      FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k)-dZdyq_w(i,j,k-1))
# endif 
#ifdef MASKING
              FY(i,j) = FY(i,j) * vmask(i,j)
#endif                 
            enddo
            enddo
	    do j=Jstr,Jend
	    do i=Istr,Iend
   	      thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)+FY(i,j)+FY(i,j+1)              
            enddo
            enddo
           endif
          endif
	enddo		         

!#ifndef NBQ_NODS 
!#  undef dZdxq_u
!#  undef dZdxq_w
!#  undef dZdyq_v
!#  undef dZdyq_w
!#endif

!#undef FY
!#define FY zwrk5

        do k=1,N

         if (IstrU.le.Iend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend+1
            FX(i,j)=on_u(i,j)* qdmu_nbq(i,j,k)
#ifdef MASKING
            FX(i,j) = FX(i,j) * umask(i,j)
#endif
          enddo
          enddo	
         endif
			  
         if (JstrV.le.Jend) then
          do j=Jstr,Jend+1		   
          do i=Istr,Iend
            FY(i,j)=om_v(i,j)* qdmv_nbq(i,j,k)
#ifdef MASKING
            FY(i,j) = FY(i,j) * vmask(i,j)
#endif
          enddo
          enddo
         endif

!---------------------------
! Horizontal Divergence :
!     Dx(p+1)     
!---------------------------   
         if (IstrU.gt.Iend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FY(i,j+1)-FY(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         elseif (JstrV.gt.Jend) then
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FX(i+1,j)-FX(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         else
          do j=Jstr,Jend		   
          do i=Istr,Iend			   
            thetadiv_nbq(i,j,k)=(thetadiv_nbq(i,j,k)                         
     &		            +WORK(i,j)*(FX(i+1,j)-FX(i,j)+FY(i,j+1)-FY(i,j))  
     &		                       ) 
#ifdef MASKING
            thetadiv_nbq(i,j,k) = thetadiv_nbq(i,j,k) * rmask(i,j)
#endif                              
          enddo
          enddo
         endif  
        enddo

!---------------------------
! Time and Bp density variations
!---------------------------   
#if defined NBQ_DTDRHO2 && defined NBQ_ZETAW
        do j=Jstr,Jend
         do i=Istr,Iend
           FC(i,0)=0.              ! Bottom boundary condition
# ifdef NBQ_MASS
           DC(i,0)=rho(i,j,1)/rho0*0.
# endif
         enddo

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=   
     &          -(z_nbq(i,j,k,knew2)-z_nbq(i,j,k,kstp2))/dtfast 
     &          *0.5*( rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k  )  
     &                +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1) )
              DC(i,k)=0.5*(rho(i,j,k)+rho(i,j,k+1))/rho0
                

       thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)+(FC(i,k)-FC(i,k-1)) 
# ifdef NBQ_MASS
     &          -Hzr(i,j,k)*(rho(i,j,k)-rho_bak(i,j,k))/rho0/dt
     &          -(zr_nbq(i,j,k,knew2)-zr_nbq(i,j,k,kstp2))
     &                     /dtfast
     &                     *(DC(i,k)-DC(i,k-1))
# endif
            enddo
          enddo
          k=N
          do i=Istr,Iend
              FC(i,k)=    
     &          -(z_nbq(i,j,k,knew2)-z_nbq(i,j,k,kstp2))/dtfast 
     &               *(
     &                rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k)
     &              )
              DC(i,N)=rho(i,j,N) /rho0
          thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)+(FC(i,k)-FC(i,k-1)) 
# ifdef NBQ_MASS
     &          -Hzr(i,j,k)*(rho(i,j,k)-rho_bak(i,j,k))/rho0/dt
     &          -(zr_nbq(i,j,N,knew2)-zr_nbq(i,j,N,kstp2))
     &                     /dtfast
     &                     *(DC(i,N)-DC(i,N-1))
# endif
          enddo
        enddo

#endif

!
!-------------------------------------------------------------------
! Implicit Vertical Momentum equation: 
!-------------------------------------------------------------------
!
#ifdef NBQ_IMP
!  
        do j=Jstr,Jend
          do k=1,N
            do i=Istr,Iend
               FC(i,k)=  soundspeed2_nbq*(rho_nbq(i,j,k)- dtnbq 
     &                   * thetadiv_nbq(i,j,k) )  !XXX5
     &                   * Hzr_half_nbq_inv(i,j,k) 
            enddo
          enddo    
  
!.........Inner layers:
          do k=1,N-1
            do i=Istr,Iend                                                               
              dum_s =   FC(i,k) - FC(i,k+1)            
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &          + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
# if defined MASKING
              qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
# endif               
            enddo             
          enddo

!.........Surface BC:
          k=N
          do i=Istr,Iend                                                               
            dum_s =   FC(i,k)                              
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   
     &         + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
# if defined MASKING
             qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k) * rmask(i,j)
# endif              
          enddo   

!.........Bottom BC:       
!# ifdef NBQ_FREE_SLIP   		   
!          do i=Istr,Iend   
!             qdmw_nbq(i,j,0)=0.
!          enddo	   
!# else
!          do i=Istr,Iend   
!             qdmw_nbq(i,j,0)=0.
!          enddo
!# endif

        enddo


!--------------------------- 
! Gaussian Elimination:
!---------------------------

!.......Comptuts coef.
        cff1=1./(dtnbq*(soundspeed2_nbq*dtnbq+visc2_nbq)) 

        do j=Jstr,Jend

!..........Bottom BC:
           k=1
           do i=Istr,Iend
             cff=(cff1+Hzw_half_nbq_inv(i,j,1)*(Hzr_half_nbq_inv(i,j,1)
     &                 +Hzr_half_nbq_inv(i,j,2)))
             CF(i,1)=(-Hzw_half_nbq_inv(i,j,2)*Hzr_half_nbq_inv(i,j,2))/cff
             DC(i,1)=qdmw_nbq(i,j,1)*cff1/cff   
     &                +qdmw_nbq(i,j,0)/cff*Hzw_half_nbq_inv(i,j,0)
     &                  *Hzr_half_nbq_inv(i,j,1)
           enddo

!..........Inner layers:
           do k=2,N-1
             do i=Istr,Iend
               cff=(cff1+                                                                    
     &               Hzw_half_nbq_inv(i,j,k)*(Hzr_half_nbq_inv(i,j,k)
     &              +Hzr_half_nbq_inv(i,j,k+1))   
     &              +Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)
     &              *CF(i,k-1))
               CF(i,k)=(-Hzw_half_nbq_inv(i,j,k+1)
     &              *Hzr_half_nbq_inv(i,j,k+1))/cff
               DC(i,k)=(qdmw_nbq(i,j,k)*cff1+Hzw_half_nbq_inv(i,j,k-1)
     &              *Hzr_half_nbq_inv(i,j,k)*DC(i,k-1)) /cff           
             enddo            
           enddo

!..........Surface BC:
           k=N
           do i=Istr,Iend
             cff=(cff1+Hzw_half_nbq_inv(i,j,N)*Hzr_half_nbq_inv(i,j,N) 
     &                     +Hzw_half_nbq_inv(i,j,N-1)
     &                     *Hzr_half_nbq_inv(i,j,N)*CF(i,N-1))  
             CF(i,N)=0. 
             DC(i,k)=(qdmw_nbq(i,j,N)*cff1+Hzw_half_nbq_inv(i,j,N-1)
     &                  *Hzr_half_nbq_inv(i,j,N)*DC(i,N-1))/cff
           enddo 

!..........Solves tri-diag system:
           do i=Istr,Iend
             qdmw_nbq(i,j,N)=DC(i,k)   
!    &            -rho_nbq(i,j,N)*0.5
!    &                 *9.81*dtnbq        
           enddo
           do k=N-1,1,-1
             do i=Istr,Iend
               qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
!    &            -0.25*(rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)
!    &                 +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1))
!    &                *(Hzr(i,j,k)+Hzr(i,j,k+1))
!    &                 *9.81*dtnbq
             enddo            
           enddo                        
        enddo    
!---------------------------
! Vertical momentum open boundary conditions
!---------------------------
# ifdef OBC_NBQ
        call wnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
                
#endif /* NBQ_IMP */


#ifdef NBQ_ZETAW


!-------------------------------------------------------------------
!       Computes surface mean velocities (Zeta
!-------------------------------------------------------------------
       
        if (IstrU.le.Iend) then
         do j=Jstr,Jend
          do i=Istr,Iend+1     
               umean_nbq(i,j)=qdmu_nbq(i,j,N)  
#ifdef NBQ_MASS                        
     &           / ( rho_nbq(i,j,N)  *Hzr_half_nbq_inv(i,j,N)
     &              +rho(i,j,N)/rho0            
     &              +rho_nbq(i-1,j,N)*Hzr_half_nbq_inv(i-1,j,N)
     &              +rho(i-1,j,N)/rho0+2.)  
     &            / (Hzr(i,j,N)+Hzr(i-1,j,N)) * 4. 
#else
     &            / (Hzr(i,j,N)+Hzr(i-1,j,N)) * 2. 
#endif
#ifdef MASKING
     &            * umask(i,j) 
#endif
          enddo 
         enddo 
        endif

        if (JstrV.le.Jend) then
         do j=Jstr,Jend+1
          do i=Istr,Iend     
               vmean_nbq(i,j)=qdmv_nbq(i,j,N)              
#ifdef NBQ_MASS     
     &            / ( rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)
     &               +rho(i,j,N)/rho0   
     &               +rho_nbq(i,j-1,N)*Hzr_half_nbq_inv(i,j-1,N)
     &               +rho(i,j-1,N)/rho0   +2.)  
     &            / (Hzr(i,j,N)+Hzr(i,j-1,N)) * 4. 
#else
     &            / (Hzr(i,j,N)+Hzr(i,j-1,N)) * 2. 
#endif
#ifdef MASKING
     &            * vmask(i,j) 
#endif
          enddo
         enddo 
        endif

        do j=Jstr,Jend
          do i=Istr,Iend
               wmean_nbq(i,j,knew2)=qdmw_nbq(i,j,N)         
#ifdef NBQ_MASS     
     &     / (rho_nbq(i,j,N)*Hzr_half_nbq_inv(i,j,N)+1.+rho(i,j,N)/rho0)  
#endif
     &             * Hzw_half_nbq_inv(i,j,N)   
#ifdef MASKING
     &             * rmask(i,j) 
#endif
          enddo
        enddo 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
        if (IstrU.le.Iend) then
         call exchange_u2d_tile (Istr,Iend,Jstr,Jend,umean_nbq(START_2D_ARRAY))
        endif
        if (JstrV.le.Jend) then
         call exchange_v2d_tile (Istr,Iend,Jstr,Jend,vmean_nbq(START_2D_ARRAY))
        endif
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,wmean_nbq(START_2D_ARRAY,knew))
# endif

# ifdef NBQ_ZETA_OUT 
      i=4
      j=1 
      k=N
      write(300,*) qdmw_nbq(i,j,k)        
     &       / (rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)/rho0+1.+rho(i,j,k)) 
     &       /Hzw_half_nbq(i,j,k)
# endif

#endif
!
!-------------------------------------------------------------------
!      Mass equation (1)
!-------------------------------------------------------------------
!		
!
!.......Computes fluxes:  
!
        do j=Jstr,Jend
            
#ifdef NBQ_FREESLIP
         do i=Istr,Iend
           FC(i,0)=Hzw_half_nbq_inv(i,j,0) * qdmw_nbq(i,j,0)             ! Bottom boundary condition
         enddo
#else
         do i=Istr,Iend
           FC(i,0)=0.                                                    ! Bottom boundary condition
         enddo
#endif

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=Hzw_half_nbq_inv(i,j,k) * qdmw_nbq(i,j,k)   
	      thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                       +FC(i,k)-FC(i,k-1)  
            enddo
          enddo
            do i=Istr,Iend
              FC(i,N)=Hzw_half_nbq_inv(i,j,N) * qdmw_nbq(i,j,N)  
	      thetadiv_nbq(i,j,N)=thetadiv_nbq(i,j,N)
     &                       +FC(i,N)-FC(i,N-1)    
            enddo
        enddo

!.......Computes rho_nbq:
!                 
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!       call exchange_r3d_tile (Istr,Iend,Jstr,Jend,rho_nbq(START_2D_ARRAY,1))
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,thetadiv_nbq(START_2D_ARRAY,1))
# endif

         do k=1,N
         do j=Jstr-1,Jend+1
         do i=Istr-1,Iend+1
           rho_nbq(i,j,k) = rho_nbq(i,j,k)  - dtfast*thetadiv_nbq(i,j,k)    !*Hzr_half_nbq_inv(i,j,k) !XXX2
         enddo
         enddo
         enddo

#ifdef NBQ_DTDRHO2
	  do j=JstrR2,JendR2
# ifdef NBQ_MASS
	    do k=1,N
	    do i=IstrR2,IendR2
             zr_nbq(i,j,k,knew2)=z_r(i,j,k)
            enddo
            enddo
# endif
   	    do k=0,N
	    do i=IstrR2,IendR2
             z_nbq (i,j,k,knew2)=z_w(i,j,k)
            enddo
            enddo
         enddo
#endif

!
!-------------------------------------------------------------------
!      Grid !
!-------------------------------------------------------------------
!		
#  include "step2d_zetaw.h"   

!
!-------------------------------------------------------------------
!      
!-------------------------------------------------------------------
!

#ifdef NBQ_DTDRHO2B
         do j=Jstr-1,Jend
         do i=Istr-1,Iend
         dum_s=0.
         do k=1,N      ! A déplacer au pas de temps interne ?
           dum_s=dum_s
     &            +((Hz(i,j,k)-Hz_bak2(i,j,k))
     &            +( (qdmu_nbq(i+1,j,k)-qdmu_nbq(i,j,k))*pm(i,j)
     &              +(qdmv_nbq(i,j+1,k)-qdmv_nbq(i,j,k))*pn(i,j) )
     &              *dtnbq)
         enddo
         do k=1,N
           hz(i,j,k)=hz(i,j,k)
     &         -dum_s/(z_w(i,j,N)-z_w(i,j,0))
     &                *(z_w(i,j,k)-z_w(i,j,k-1))
         enddo
         enddo
         enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hz(START_2D_ARRAY,1))
#  endif 
#endif

#if defined NBQ_DTDRHO2 && defined NBQ_ZETAW && defined NBQ_MASS
          if (iif==nfast) then
            do k=1,N 
            do j=Jstr,Jend             
            do i=Istr,Iend
	       rho_bak(i,j,k)=rho(i,j,k)
            enddo
            enddo
            enddo
          endif
#endif
            
#ifdef RVTK_DEBUG
       call check_tab3d(rho_nbq,'rho_nbq','r')
#endif    
!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
!       call rnbqijk_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif

!
!*******************************************************************
!*******************************************************************
!      enddo    ! NBQ loop
!*******************************************************************
!*******************************************************************
!
   !     if (LAST_2D_STEP) then
!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
!          call ruijk_nbq(2, Istr,Iend,Jstr,Jend,WORK)


# ifdef M2FILTER_NONE
        if (LAST_2D_STEP) then
# endif 
          
        do k=0,N 
          do j=Jstr,Jend              
            do i=Istr,Iend
              rw_nbq_ext (i,j,k) = ((qdmw_nbq(i,j,k)-rw_nbq_ext(i,j,k))
     &              /dtnbq-ndtnbq*rw_int_nbq(i,j,k))/WORK(i,j)
            enddo
          enddo
        enddo
# ifdef M2FILTER_NONE        
        endif
# endif 

# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!      call exchange_r3d_tile (Istr,Iend,Jstr,Jend          ! TBD
!    &           ,  rw_nbq_ext(START_2D_ARRAY,0)) 
# endif
# ifdef RVTK_DEBUG
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext (ru_nbq)','v')
# endif    


#ifdef RVTK_DEBUG
          call check_tab2d(rubar_nbq,'rubar_nbq step3d_nbq','uint')
          call check_tab2d(rvbar_nbq,'rvbar_nbq step3d_nbq','vint')
!         call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif  
!    
#ifdef NBQ_MASS
         call densityijk_nbq(20)
#endif
!        endif
!
!  
!      end subroutine step3d_fbijk_nbq
!
!#else
!      subroutine step3d_fbijk_nbq_empty
!      end subroutine step3d_fbijk_nbq_empty
!#endif
! #endif
