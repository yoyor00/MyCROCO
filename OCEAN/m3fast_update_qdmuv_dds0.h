
! FRANCIS ICI

# ifdef NBQ
!
!-----------------------------------------------------------------------
!  Pressure-Viscosity forces in XI- and ETA-Directions
!-----------------------------------------------------------------------
!
      k2 = 1
      do k=0,N
        k1=k2
        k2=3-k1

#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif

          if (k.eq.0) then ! Bottom Boundary conditions

            do j=Jstr,Jend
              do i=Istr-1,Iend
                dthetadiv_nbqdz_u(i,j,k2)=0. 
              enddo
            enddo

            do j=Jstr-1,Jend
              do i=Istr,Iend
                dthetadiv_nbqdz_v(i,j,k2)=0.
              enddo
            enddo 

          else

            if (k.eq.N) then ! Top Boundary conditions

              do j=Jstr-1,Jend
                do i=Istr-1,Iend
#  ifdef NBQ_GRID_SLOW
                  dthetadiv_nbqdz(i,j,k,1)= - thetadiv_nbq(i,j,k)
#  else
                  dthetadiv_nbqdz(i,j)    = - thetadiv_nbq(i,j,k)
#  endif
                enddo
              enddo

            else

              do j=Jstr-1,Jend
                do i=Istr-1,Iend
#  ifdef NBQ_GRID_SLOW
                  dthetadiv_nbqdz(i,j,k,1)=thetadiv_nbq(i,j,k+1)
     &                                   - thetadiv_nbq(i,j,k)
#  else
                  dthetadiv_nbqdz(i,j)    =thetadiv_nbq(i,j,k+1) 
     &                                   - thetadiv_nbq(i  ,j,k)
#  endif
                enddo
              enddo

            endif
  
            do j=Jstr,Jend
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                dthetadiv_nbqdz_u(i,j,k2)=Hzw_nbq_inv_u(i,j,k)*(
     &                                    dthetadiv_nbqdz(i,j,k,1)
     &                                   +dthetadiv_nbqdz(i-1,j,k,1))
#  else
                dthetadiv_nbqdz_u(i,j,k2)=Hzw_nbq_inv_u(i,j,k)*(
     &                                    dthetadiv_nbqdz(i,j)
     &                                   +dthetadiv_nbqdz(i-1,j))
#  endif          
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
#  ifdef NBQ_GRID_SLOW
                 dthetadiv_nbqdz_v(i,j,k2)=Hzw_nbq_inv_v(i,j,k)*(
     &                                     dthetadiv_nbqdz(i,j,k,1)
     &                                    +dthetadiv_nbqdz(i,j-1,k,1))
#  else
                 dthetadiv_nbqdz_v(i,j,k2)=Hzw_nbq_inv_v(i,j,k)*(
     &                                     dthetadiv_nbqdz(i,j)
     &                                    +dthetadiv_nbqdz(i,j-1))
#  endif
              enddo
            enddo    
  
          endif    

#  ifdef NBQ_GRID_SLOW
        endif ! NSTEP_DS
#  endif

        if (k.gt.0) then
!
!-----------------------------------------------------------------------
!  Fast-mode U-momentum: qdmu_nbq
!-----------------------------------------------------------------------
!
          do j=Jstr,Jend
            do i=IstrU,Iend
              if (k.gt.1.and.k.lt.N) then 
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i  ,j,k)
     &                                     -z_r(i-1,j,k))  
     &              *(dthetadiv_nbqdz_u(i,j,k2)+
     &                dthetadiv_nbqdz_u(i,j,k1)) ! dZdx * (d(delta p)dz)_u
                endif 
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                   
     &              *(dthetadiv_nbqdz_u(i,j,k2)+
     &                dthetadiv_nbqdz_u(i,j,k1)) ! dZdx * (d(delta p)dz)_u
#  endif
              elseif (k.gt.1) then
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i ,j,k)
     &                                     -z_r(i-1,j,k))       
     &                 *dthetadiv_nbqdz_u(i,j,k1) ! dZdx * (d(delta p)dz)_u
     &                 +(z_w(i,j,N)-z_w(i-1,j,N))
     &                 *dthetadiv_nbqdz_u(i,j,k2)
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                *dthetadiv_nbqdz_u(i,j,k1) ! dZdx * (d(delta p)dz)_u
     &                +(z_w(i,j,N)-z_w(i-1,j,N))   
     &                *dthetadiv_nbqdz_u(i,j,k2)
#  endif
              else
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i ,j,k)
     &                                     -z_r(i-1,j,k))       
     &             *2.*dthetadiv_nbqdz_u(i,j,k2) ! dZdx * (d(delta p)dz)_u
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#  else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &              *2.*dthetadiv_nbqdz_u(i,j,k2) ! dZdx * (d(delta p)dz)_u
#  endif
              endif
#  ifdef MASKING
              if (umask(i-1,j)*umask(i+1,j) .eq. 0.) then
                dum_s=dum_s-(thetadiv_nbq(i  ,j,k)-
     &                       thetadiv_nbq(i-1,j,k))  ! - d(delta p)dx
              else
#  endif
                dum_s=dum_s
     &                -(gammau  *thetadiv_nbq(i  ,j,k)+
     &                  gammau_2*thetadiv_nbq(i+1,j,k)-
     &                  gammau  *thetadiv_nbq(i-1,j,k)-
     &                  gammau_2*thetadiv_nbq(i-2,j,k)) ! - d(delta p)dx
#  ifdef MASKING
              endif
#  endif
              dum_s=dum_s*Hzu_qdmu(i,j,k)
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcoru(i,j,k)
#  endif
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.5*omega*U0*cos(omega*time)
     &                       *(Hz(i-1,j,k)+Hz(i,j,k))
#  endif
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bustr(i,j)
#  endif
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtnbq * (
     &                        dum_s + ru_int_nbq(i,j,k) )  
#  ifdef MASKING
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*umask(i,j)
#  endif
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
              ru_nbq(i,j,k)=dum_s/work(i,j) 

#  if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*(1.-NBQnudgcof(i,j))
     &                        +u(i,j,k,nrhs)*Hzu_qdmu(i,j,k)
     &                                           *NBQnudgcof(i,j)
#  endif
            enddo 
          enddo
!
!-----------------------------------------------------------------------
!  Fast-mode V-momentum: qdmv_nbq
!-----------------------------------------------------------------------
!
          do j=JstrV,Jend
            do i=Istr,Iend
              if (k.gt.1.and.k.lt.N) then 
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k2)+
     &                  dthetadiv_nbqdz_v(i,j,k1)) ! dZdy * (d(delta p)dz)_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k2)+
     &                  dthetadiv_nbqdz_v(i,j,k1)) ! dZdy * (d(delta p)dz)_v
#  endif
              elseif (k.gt.1) then
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *dthetadiv_nbqdz_v(i,j,k1) ! dZdy * (d(delta p)dz)_v
     &                +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                *dthetadiv_nbqdz_v(i,j,k2)
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))            
     &                *dthetadiv_nbqdz_v(i,j,k1) ! dZdy * (d(delta p)dz)_v
     &                +(z_w(i,j,N)-z_w(i,j-1,N))      
     &                *dthetadiv_nbqdz_v(i,j,k2)
#  endif
              else
#  ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)
     &                                     -z_r(i,j-1,k)) 
     &                *2.*dthetadiv_nbqdz_v(i,j,k2) ! dZdy * (d(delta p)dz)_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2) 
#  else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *2.*dthetadiv_nbqdz_v(i,j,k2) ! dZdy * (d(delta p)dz)_v
#  endif
              endif
                
#  ifdef MASKING
              if (vmask(i,j-1)*vmask(i,j+1) .eq. 0.) then
              dum_s=dum_s
     &             -(thetadiv_nbq(i,j  ,k)-
     &               thetadiv_nbq(i,j-1,k)) ! - d(delta p)dy
              else
#  endif
              dum_s=dum_s
     &             -(gammau  *thetadiv_nbq(i,j  ,k)+
     &               gammau_2*thetadiv_nbq(i,j+1,k)-
     &               gammau  *thetadiv_nbq(i,j-1,k)-
     &               gammau_2*thetadiv_nbq(i,j-2,k)) ! - d(delta p)dy
#  ifdef MASKING
              endif
#  endif
              dum_s=dum_s*Hzv_qdmv(i,j,k)
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcorv(i,j,k)
#  endif
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.25*(f(i,j)+f(i,j-1))*
     &                       U0*sin(omega*time)*
     &                       (Hz(i,j,k)+Hz(i,j-1,k))
#  endif
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bvstr(i,j)
#  endif
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtnbq * (
     &                        dum_s + rv_int_nbq(i,j,k) )
#  ifdef MASKING
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*vmask(i,j)
#  endif
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
              rv_nbq(i,j,k)=dum_s/work(i,j)  

#  if defined NBQ_NUDGING && defined NBQCLIMATOLOGY
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*(1.-NBQnudgcof(i,j))
     &                        +v(i,j,k,nrhs)*Hzv_qdmv(i,j,k)
     &                                           *NBQnudgcof(i,j)
#  endif
            enddo
          enddo
        endif  !<-- k>0
      enddo   !<-- k=0,N   

# else /* NBQ */
!
!********************************************************************
!
!  Advance fast 3D u,v for the hydrostatic case *****
!
!********************************************************************
!
#  ifdef BSTRESS_FAST
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar_nbq(i,j)=rubar_nbq(i,j)-bustr(i,j)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend      
          rvbar_nbq(i,j)=rvbar_nbq(i,j)-bvstr(i,j)
        enddo
      enddo
#  endif
      if ((mod(iif-1,inc_faststep) .eq. inc_faststep-1) .OR.
     &    (LAST_FAST_STEP)) then
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              dum_s=0.
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bustr(i,j)
#  endif
              qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)+dtfast * (
     &                      nb_faststep*(dum_s + ru_int_nbq(i,j,k))
     &                     +rubar_sum(i,j)*(Hz(i-1,j,k)+Hz(i,j,k)))    
#  ifdef MASKING
     &                                        *umask(i,j)
#  endif
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
              if (LAST_FAST_STEP) ru_nbq(i,j,k)=dum_s/work(i,j) 
            enddo 
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              dum_s=0.
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bvstr(i,j)
#  endif
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)+dtfast * (
     &                       nb_faststep*(dum_s + rv_int_nbq(i,j,k))
     &                      +rvbar_sum(i,j)*(Hz(i,j-1,k)+Hz(i,j,k)))
#  ifdef MASKING
     &                                        *vmask(i,j)
#  endif
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
              if (LAST_FAST_STEP) rv_nbq(i,j,k)=dum_s/work(i,j)
            enddo
          enddo
        enddo
      endif

# endif /* NBQ */ 
!
!-----------------------------------------------------------------------
! Apply point sources for river runoff simulations
!-----------------------------------------------------------------------
!
# ifdef PSOURCE
      do is=1,Nsrc 
#  ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#  else
        i=Isrc(is)
        j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            DU_nbq(i,j)=0.
#  ifdef WET_DRY
            umask_wet(i,j)=1
#  endif
            do k=1,N
              qdmu_nbq(i,j,k)=Qsrc(is,k)*pn_u(i,j)
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
            enddo
          else
            DV_nbq(i,j)=0.
#  ifdef WET_DRY
            vmask_wet(i,j)=1
#  endif
            do k=1,N
              qdmv_nbq(i,j,k)=Qsrc(is,k)*pm_v(i,j)
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
            enddo
          endif
        endif
      enddo
# endif
!
!-----------------------------------------------------------------------
!  U & V momentum wet mask
!-----------------------------------------------------------------------
!
# ifdef WET_DRY
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DU_nbq(i,j))*umask_wet(i,j)
          umask_wet(i,j)=0.5*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DU_nbq(i,j)=DU_nbq(i,j)*umask_wet(i,j)
!#  ifdef MRL_WCI
!          ust2d(i,j)=ust2d(i,j)*umask_wet(i,j)
!#  endif
          do k=1,N
            qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*umask_wet(i,j)
          enddo
        enddo 
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.)
          cff2_WD=0.5+SIGN(0.5,DV_nbq(i,j))*vmask_wet(i,j)
          vmask_wet(i,j)=0.5*vmask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.-cff1_WD)
          DV_nbq(i,j)=DV_nbq(i,j)*vmask_wet(i,j)
!#  ifdef MRL_WCI
!          vst2d(i,j)=vst2d(i,j)*vmask_wet(i,j)
!#  endif
          do k=1,N
            qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*vmask_wet(i,j)
          enddo
        enddo 
      enddo
# endif

# ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(qdmu_nbq,'qdmu_nbqint','uint')
      call check_tab3d(qdmv_nbq,'qdmv_nbqint','vint')
C$OMP END MASTER
# endif
!
!-----------------------------------------------------------------------
!  U & V momentum open boundary conditions
!-----------------------------------------------------------------------
!
      do j=Jstr,Jend
        do i=IstrU,Iend
          ubar(i,j,knew)=urhs(i,j)
        enddo 
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          vbar(i,j,knew)=vrhs(i,j)
        enddo 
      enddo
# ifdef NBQ
      M2bc_nbq_flag=.true. ! apply boundary wet/dry conditions
                           ! and compute DU_nbq
      call u2dbc_tile   (Istr,Iend,Jstr,Jend, work)
      call v2dbc_tile   (Istr,Iend,Jstr,Jend, work)
# endif
      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
!
!-----------------------------------------------------------------------
! Exchange periodic boundaries and computational margins
!-----------------------------------------------------------------------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  
#  ifndef NBQ
      if ((mod(iif-1,inc_faststep) .eq. inc_faststep-1) .OR.
     &     (LAST_FAST_STEP)) then
#  endif
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmv_nbq(START_2D_ARRAY,1))
#  ifndef NBQ
      endif
#  endif
!      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
!     &                        DU_nbq(START_2D_ARRAY))
!      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
!     &                        DV_nbq(START_2D_ARRAY))
# endif
# ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
      call check_tab3d(qdmu_nbq,'qdmu_nbq','u')
      call check_tab3d(qdmv_nbq,'qdmv_nbq','v')
C$OMP END MASTER
# endif  
!   FRANCIS ICI
