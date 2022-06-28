
! FRANCIS LA
!
!
!=====================================================================
!      Fast-mode conservation of mass
!  
!  ... and integration of W-momentum with IMPLICIT scheme
!
!  From now on, thetadiv_nbq array is used for divergence (not theta):
!
!=====================================================================
!
      if (IstrU.gt.Iend) then
        do j=Jstr,Jend
          do i=Istr,Iend+1
            FX(i,j)=0.
          enddo
        enddo
      endif
      if (JstrV.gt.Jend) then
        do j=Jstr,Jend+1
          do i=Istr,Iend 
            FY(i,j)=0.
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!  X-component dZdx*qdmu for horizontal divergence
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif
          if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_u(i,j,k2)=(z_r(i,j,kp1)-z_r(i-1,j,kp1)) 
     &                           *qdmu_nbq(i,j,kp1)  ! (dZdx * (rho u))_u
              enddo
            enddo
          endif

          if (k.eq.0) then  ! Bottom boundary conditions

#  ifdef NBQ_FREESLIP
            do j=Jstr,Jend
              do i=Istr,Iend+1 
#   ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k ) = (z_w(i,j,0)-z_w(i-1,j,0))
     &                             *qdmu_nbq(i,j,1)  
     &                             /(Hzr(i,j,1)+Hzr(i-1,j,1))
#   else
                dZdxq_w(i,j,k2)= (z_w(i,j,0)-z_w(i-1,j,0))
     &                            *qdmu_nbq(i,j,1)  
     &                            /(Hzr(i,j,1)+Hzr(i-1,j,1))
#   endif
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend    
#   ifdef NBQ_GRID_SLOW
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k)*pm_u(i  ,j)
     &                              +dZdxq_w(i+1,j,k)*pm_u(i+1,j) )
     &                                                 * Hzr(i,j,1)
#   else
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k2)*pm_u(i  ,j)
     &                              +dZdxq_w(i+1,j,k2)*pm_u(i+1,j) )
     &                                                  * Hzr(i,j,1)
#   endif
#   ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#   endif 
              enddo
            enddo 

#  else /* NBQ_FREESLIP */

#   ifndef NBQ_GRID_SLOW
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_w(i,j,k2)=0.  
                qdmw_nbq(i,j,0)=0.     
              enddo
            enddo 
#   endif
 
#  endif /* NBQ_FREESLIP */

          elseif (k .eq. N) then ! Top boundary conditions
           
            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k )= (z_w(i,j,N)-z_w(i-1,j,N))
     &                            *qdmu_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i-1,j,N))
#  else
                dZdxq_w(i,j,k2)= (z_w(i,j,N)-z_w(i-1,j,N))
     &                            *qdmu_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i-1,j,N))
#  endif
              enddo
            enddo  
 
          else ! k<>0 & k<>N   ! Inner domain   

            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                dZdxq_w(i,j,k )=Hzw_nbq_inv_u(i,j,k)
     &                                 *(dZdxq_u(i,j,k1)+
     &                                   dZdxq_u(i,j,k2))
#  else
                dZdxq_w(i,j,k2)=Hzw_nbq_inv_u(i,j,k)
     &                                 *(dZdxq_u(i,j,k1)+
     &                                   dZdxq_u(i,j,k2))
#  endif
              enddo 
            enddo

          endif ! k<>0 , k<>N ,  Inner domain 

#  ifdef NBQ_GRID_SLOW

        else   ! NSTEP_DS: Update d./ds terms

          if (k.eq.0) then  ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend
              do i=Istr,Iend+1 
                dZdxq_w(i,j,k)=(z_w(i,j,0)-z_w(i-1,j,0))
     &                          *qdmu_nbq(i,j,1)  
     &                          /(Hzr(i,j,1)+Hzr(i-1,j,1))
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend    
                qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k)*pm_u(i  ,j) 
     &                              +dZdxq_w(i+1,j,k)*pm_u(i+1,j) ) 
     &                                                  *Hzr(i,j,1)    
#    ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
              enddo
            enddo 

#   endif /* NBQ_FREESLIP */

          endif ! k.eq.0
        endif ! NSTEP_DS: Update d./ds terms

#  else  /* NBQ_GRID_SLOW */

        if (k.eq.0) then ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
          do j=Jstr,Jend
            do i=Istr,Iend+1 
              dZdxq_w(i,j,k2)=(z_w(i,j,0)-z_w(i-1,j,0))
     &                         *qdmu_nbq(i,j,1)  
     &                         /(Hzr(i,j,1)+Hzr(i-1,j,1))
            enddo
          enddo

          do j=Jstr,Jend
            do i=Istr,Iend    
              qdmw_nbq(i,j,0)=0.5*(dZdxq_w(i  ,j,k2)*pm_u(i  ,j) 
     &                            +dZdxq_w(i+1,j,k2)*pm_u(i+1,j) ) 
     &                                                 *Hzr(i,j,1)    
#    ifdef MASKING
              qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
            enddo
          enddo 

#   endif /* NBQ_FREESLIP */ 

        endif ! k.eq.0

#  endif  /* NBQ_GRID_SLOW */

        if (k.gt.0) then
          if (IstrU.le.Iend) then
            do j=Jstr,Jend
              do i=Istr,Iend+1
#  ifdef NBQ_GRID_SLOW
                FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k )-dZdxq_w(i,j,k-1))
#  else
                FX(i,j)=-pm_u(i,j)*(dZdxq_w(i,j,k2)-dZdxq_w(i,j,k1))
#  endif
#  ifdef MASKING
                FX(i,j)=FX(i,j)*umask(i,j)
#  endif                
              enddo
            enddo

            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=FX(i,j)+FX(i+1,j)
              enddo
            enddo

          else ! IstrU.gt.Iend

            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=0.
              enddo
            enddo

          endif 
        endif
      enddo ! k=0,N
!
!-----------------------------------------------------------------------
!  Y-component dZdy*qdmv for horizontal divergence
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N  !<-- k loop
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif
          if (k.lt.N) then
            kp1 = k + 1
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_v(i,j,k2)=(z_r(i,j,kp1)-z_r(i,j-1,kp1)) 
     &                           *qdmv_nbq(i,j,kp1)  ! (dZdy * (rho v))_v
              enddo
            enddo
          endif

          if (k.eq.0) then  ! Bottom boundary conditions

#  ifdef NBQ_FREESLIP
            do j=Jstr,Jend+1
              do i=Istr,Iend
#   ifdef NBQ_GRID_SLOW
                 dZdyq_w(i,j,k )= (z_w(i,j,0)-z_w(i,j-1,0))
     &                             *qdmv_nbq(i,j,1) 
     &                             /(Hzr(i,j,1)+Hzr(i,j-1,1))
#   else
                 dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                             *qdmv_nbq(i,j,1) 
     &                             /(Hzr(i,j,1)+Hzr(i,j-1,1))
#   endif 
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend   
#   ifdef NBQ_GRID_SLOW
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                          +0.5*(dZdyq_w(i,j  ,k)*pn_v(i,j  )  
     &                               +dZdyq_w(i,j+1,k)*pn_v(i,j+1) )
     &                                                  * Hzr(i,j,1)
#   else
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                          +0.5*(dZdyq_w(i,j  ,k2)*pn_v(i,j  )  
     &                               +dZdyq_w(i,j+1,k2)*pn_v(i,j+1) )
     &                                                   * Hzr(i,j,1)
#   endif 
#   ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#   endif 
              enddo
            enddo

#  else  /* NBQ_FREESLIP */

#   ifndef NBQ_GRID_SLOW   
            do j=Jstr,Jend +1
              do i=Istr,Iend  
                dZdyq_w(i,j,k2)=0.  
                qdmw_nbq(i,j,0)=0.
              enddo
            enddo 
#   endif
#  endif /* NBQ_FREESLIP */

          elseif (k .eq. N) then ! Top boundary conditions

            do j=Jstr,Jend+1
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dZdyq_w(i,j,k )= (z_w(i,j,N)-z_w(i,j-1,N))
     &                           *qdmv_nbq(i,j,N)
     &                           /(Hzr(i,j,N)+Hzr(i,j-1,N))
#  else
                dZdyq_w(i,j,k2)= (z_w(i,j,N)-z_w(i,j-1,N))
     &                            *qdmv_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i,j-1,N))
#  endif 
              enddo
            enddo

          else

            do j=Jstr,Jend+1
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dZdyq_w(i,j,k )=Hzw_nbq_inv_v(i,j,k)
     &                           *(dZdyq_v(i,j,k1)+
     &                             dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#  else
                dZdyq_w(i,j,k2)=Hzw_nbq_inv_v(i,j,k)
     &                           *(dZdyq_v(i,j,k1)+
     &                             dZdyq_v(i,j,k2)) ! (dZdy * (rho v))_uw/Hzw_v
#  endif 
              enddo 
            enddo

          endif

#  ifdef NBQ_GRID_SLOW

        else  ! NSTEP_DS

          if (k.eq.0) then ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
            do j=Jstr,Jend+1
              do i=Istr,Iend
                dZdyq_w(i,j,k)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                           *qdmv_nbq(i,j,1) 
     &                           /(Hzr(i,j,1)+Hzr(i,j-1,1))
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend   
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                         +0.5*(dZdyq_w(i,j  ,k)*pn_v(i,j  )  
     &                              +dZdyq_w(i,j+1,k)*pn_v(i,j+1) )
     &                                                 * Hzr(i,j,1) 
#    ifdef MASKING
                qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
              enddo
            enddo
#   endif
         endif
        endif ! NSTEP_DS

#  else /* NBQ_GRID_SLOW */

        if (k.eq.0) then  ! Bottom boundary conditions

#   ifdef NBQ_FREESLIP
          do j=Jstr,Jend+1
            do i=Istr,Iend
              dZdyq_w(i,j,k2)= (z_w(i,j,0)-z_w(i,j-1,0))
     &                          *qdmv_nbq(i,j,1) 
     &                          /(Hzr(i,j,1)+Hzr(i,j-1,1))
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend   
               qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)
     &                         +0.5*(dZdyq_w(i,j  ,k2)*pn_v(i,j  )  
     &                              +dZdyq_w(i,j+1,k2)*pn_v(i,j+1) )
     &                                                  * Hzr(i,j,1) 
#    ifdef MASKING
               qdmw_nbq(i,j,0)=qdmw_nbq(i,j,0)*rmask(i,j)
#    endif 
            enddo
          enddo
#   endif  /* NBQ_FREESLIP */
        endif

#  endif /* NBQ_GRID_SLOW */

        if (k.gt.0) then
          if (JstrV.le.Jend) then
            do j=Jstr,Jend+1
              do i=Istr,Iend 
#  ifdef NBQ_GRID_SLOW
                FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k)-dZdyq_w(i,j,k-1))
#  else
                FY(i,j)=-pn_v(i,j)*(dZdyq_w(i,j,k2)-dZdyq_w(i,j,k1))
#  endif 
#  ifdef MASKING
                FY(i,j)=FY(i,j)*vmask(i,j)
#  endif                 
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                               +FY(i,j)+FY(i,j+1)
              enddo
            enddo
          endif
        endif

      enddo  !<-- k=0,N
!
!-----------------------------------------------------------------------
!  Compute total horizontal Divergence divH(qdmH(m+1))
!-----------------------------------------------------------------------
!
      do k=1,N !<-- k loop

        if (IstrU.le.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FX(i,j)=on_u(i,j)*qdmu_nbq(i,j,k)
            enddo
          enddo
        endif
        if (JstrV.le.Jend) then
          do j=Jstr,Jend+1  
            do i=Istr,Iend
              FY(i,j)=om_v(i,j)*qdmv_nbq(i,j,k)
            enddo
          enddo
        endif
        if (IstrU.gt.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend  
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                           +pm(i,j)*pn(i,j)*(FY(i,j+1)-FY(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        elseif (JstrV.gt.Jend) then
          do j=Jstr,Jend  
            do i=Istr,Iend   
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j)+
     &                                              FY(i,j+1)-FY(i,j))  
#  ifdef MASKING
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)*rmask(i,j)
#  endif                              
            enddo
          enddo
        endif  

      enddo ! <-- k=1,N
!
!-----------------------------------------------------------------------
!  thetadiv2_nbq: complet time-corrective term  (dh/dt included) 
!  thetadiv3_nbq: reduced time-corrective term  (no dh/dt)
!-----------------------------------------------------------------------
!
#  ifdef NBQ_GRID_SLOW
       if (FIRST_FAST_STEP) then
#  endif

        do j=JstrR,JendR
          do k=0,N
            do i=IstrR,IendR
              zw_nbq(i,j,k,knew)=z_w(i,j,k)
            enddo
          enddo
        enddo
        do j=Jstr,Jend !<-- j loop

          do i=Istr,Iend
            FC(i,0)=0.    ! Bottom boundary condition
            CF(i,0)=0.
          enddo

          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=   
     &          -(zw_nbq(i,j,k,knew)-zw_nbq(i,j,k,kstp))/dtgrid_nbq
     &           *0.5*( (1.+rho_grd(i,j,k  ))
     &                  +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  )  
     &                  +(1.+rho_grd(i,j,k+1))
     &                  +rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1))
              CF(i,k)=   
     &          -(zw_nbq(i,j,k,knew)-zw_nbq(i,j,k,kstp))/dtgrid_nbq 
     &          *0.5*(    rho_grd(i,j,k  ) 
     &                   +rho_nbq(i,j,k  )*Hzr_nbq_inv(i,j,k  ) 
     &                   +rho_grd(i,j,k+1)
     &                   +rho_nbq(i,j,k+1)*Hzr_nbq_inv(i,j,k+1) )

              thetadiv3_nbq(i,j,k)=CF(i,k)-CF(i,k-1) 
              thetadiv2_nbq(i,j,k)=FC(i,k)-FC(i,k-1) 
            enddo
          enddo
    
           do i=Istr,Iend
             FC(i,N)=   
     &        -(zw_nbq(i,j,N,knew)-zw_nbq(i,j,N,kstp))/dtgrid_nbq
     &         *( 1.+rho_grd(i,j,N)
     &          )

            CF(i,N)=  
     &        -(zw_nbq(i,j,N,knew)-zw_nbq(i,j,N,kstp))/dtgrid_nbq
     &         *( rho_grd(i,j,N)
     &          )
            thetadiv3_nbq(i,j,N)=CF(i,N)-CF(i,N-1) 
            thetadiv2_nbq(i,j,N)=FC(i,N)-FC(i,N-1) 
          enddo
        enddo !<-- j loop

#  ifdef NBQ_GRID_SLOW
       endif !<-- FIRST_FAST_STEP
#  endif
! !
! !====================================================================
! !         Exchange divergence (3)
! !====================================================================
! !	
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                thetadiv3_nbq(START_2D_ARRAY,-N_sl+1))

! FRANCIS LA
