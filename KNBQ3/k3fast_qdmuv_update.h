! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_qdmuv_update.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! !********************************
! ! U and V fast-momentum
! !********************************
! !********************************
! !
! !
! !********************************
! ! Define local variables:
! !********************************
! !
#define dthetadiv_nbqdz_u FC3D
#define dthetadiv_nbqdz_v DC3D
#ifdef NBQ_GRID_SLOW
# define dthetadiv_nbqdz_ij dthetadiv_nbqdz(i,j,k,1)
# define dthetadiv_nbqdz_im1j dthetadiv_nbqdz(i-1,j,k,1)
# define dthetadiv_nbqdz_ijm1 dthetadiv_nbqdz(i,j-1,k,1)
#else
# define dthetadiv_nbqdz_ij dthetadiv_nbqdz(i,j)
# define dthetadiv_nbqdz_im1j dthetadiv_nbqdz(i-1,j)
# define dthetadiv_nbqdz_ijm1 dthetadiv_nbqdz(i,j-1)
#endif 
! !
# ifdef K3FAST_UV
!$acc kernels if(compute_on_device) default(present)

       do k=-N_sl,N    ! Loop is on w-levels
#  ifdef K3FAST_RHO
! !
! !********************************
! !  Pressure-Viscosity forces 
! !   in XI- and ETA-Directions
! !                     \Delta k
! !********************************
! !
#  ifdef NBQ_GRID_SLOW
        if (NSTEP_DS) then
#  endif
! !
! !.......k=-N_sl..................
! ! 
          if (k.eq.-N_sl) then  
            do j=Jstr,Jend
              do i=Istr-1,Iend
                dthetadiv_nbqdz_u(i,j,k)=0. 
              enddo
            enddo
            do j=Jstr-1,Jend
              do i=Istr,Iend
                dthetadiv_nbqdz_v(i,j,k)=0.
              enddo
            enddo 
          else
            if (k.eq.N) then
! !
! !......k=N.......................
! ! Mathematica
              do j=Jstr-1,Jend
                do i=Istr-1,Iend
                  dthetadiv_nbqdz_ij=-thetadiv_nbq(i,j,k)
                enddo
              enddo
            else            
! !
! !......0<k<N.....................
! ! Mathematica
              do j=Jstr-1,Jend
                do i=Istr-1,Iend
                  dthetadiv_nbqdz_ij= thetadiv_nbq(i,j,k+1)
     &                               -thetadiv_nbq(i,j,k)
                enddo
              enddo
            endif
! !
! !********************************
! !  Pressure-Viscosity forces 
! !     in XI- and ETA-Directions
! !                     \Sigma k
! !********************************
! !
            do j=Jstr,Jend
              do i=Istr,Iend
                dthetadiv_nbqdz_u(i,j,k)=Hzw_nbq_inv_u(i,j,k)*
     &                                   ( dthetadiv_nbqdz_ij
     &                                    +dthetadiv_nbqdz_im1j )
#   ifdef MASKING
     &              *umask(i,j)
#   endif
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                 dthetadiv_nbqdz_v(i,j,k)=Hzw_nbq_inv_v(i,j,k)*
     &                                    ( dthetadiv_nbqdz_ij
     &                                     +dthetadiv_nbqdz_ijm1 )
#   ifdef MASKING
     &              *vmask(i,j)
#   endif
              enddo
            enddo    
          endif   
#  ifdef NBQ_GRID_SLOW
        endif ! NSTEP_DS
#  endif
        
#  endif /* K3FAST_RHO */
      enddo   !<-- k=0,   
!$acc end kernels
!$acc kernels if(compute_on_device) default(present)

! !
      do k=-N_sl+1,N    ! Loop is on w-levels
! !
! !********************************
! !  Fast-mode U-momentum: qdmu_nbq
! !********************************
! !
          do j=Jstr,Jend
            do i=IstrU,Iend
! !
! !................................
! !  Comp. pressure gradient & 2nd visc.
! !         U-component: \Sigma  k
! !................................
! !
#  ifdef K3FAST_RHO
! !
! !....1<k<N (-Nsl+1<k<N)..........
! !
#   ifndef K3FAST_SEDLAYERS
              if (k.gt.1      .and.k.lt.N) then 
#   else
              if (k.gt.-N_sl+1.and.k.lt.N) then 
#   endif
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)-z_r(i-1,j,k))  
     &              *(dthetadiv_nbqdz_u(i,j,k)+
     &                dthetadiv_nbqdz_u(i,j,k-1)) 
                endif 
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#   else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                   
     &              *(dthetadiv_nbqdz_u(i,j,k)+
     &                dthetadiv_nbqdz_u(i,j,k-1)) 
#   endif
! !
! !.......k=N......................
! ! Mathematica
              elseif (k.eq.N) then      
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=
     &                  (z_r(i ,j,k)-z_r(i-1,j,k))       
     &                 *dthetadiv_nbqdz_u(i,j,k-1) 
     &                 +(z_w(i,j,N)-z_w(i-1,j,N))
     &                 *dthetadiv_nbqdz_u(i,j,k)
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#   else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                *dthetadiv_nbqdz_u(i,j,k-1) 
     &                +(z_w(i,j,N)-z_w(i-1,j,N))   
     &                *dthetadiv_nbqdz_u(i,j,k) ! BC on dthetadiv_nbqdz_u
#   endif
! !
! !......k=-Nsl+1..................
! !
              else    
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,1)=(z_r(i,j,k)-z_r(i-1,j,k))       
     &                                   *2.*dthetadiv_nbqdz_u(i,j,k) 
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,1)
#   else
                dum_s=(z_r(i,j,k)-z_r(i-1,j,k))                      
     &                *dthetadiv_nbqdz_u(i,j,k)
#   endif
              endif                     ! -------- elseif 1<k<N
! !
! !................................
! ! Horiz. PG
! !................................
! !
#   ifdef K3FAST_PG2
#    ifdef EW_PERIODIC 
              dum_s=dum_s
     &                -( gammau  *thetadiv_nbq(i  ,j,k)
     &                  +gammau_2*thetadiv_nbq(i+1,j,k)
     &                  -gammau  *thetadiv_nbq(i-1,j,k)
     &                  -gammau_2*thetadiv_nbq(i-2,j,k)) ! - d(delta p)dx
#    else
#     if defined MPI
         if ((WEST_INTER.or.i.ne.IstrU).and.(EAST_INTER.or.i.ne.Iend)
#      ifdef MASKING
     &   .and.(umask(i-1,j)*umask(i+1,j) .ne. 0.)
#      endif
     &      ) then
#     else
         if (i.ne.IstrU.and.i.ne.Iend
#      ifdef MASKING
     &   .and.(umask(i-1,j)*umask(i+1,j) .ne. 0.)
#      endif
     &      ) then
#     endif
              dum_s=dum_s
     &                -( gammau  *thetadiv_nbq(i  ,j,k)
     &                  +gammau_2*thetadiv_nbq(i+1,j,k)
     &                  -gammau  *thetadiv_nbq(i-1,j,k)
     &                  -gammau_2*thetadiv_nbq(i-2,j,k)) ! - d(delta p)dx
         else
              dum_s=dum_s
     &                -( thetadiv_nbq(i  ,j,k)
     &                  -thetadiv_nbq(i-1,j,k))
	     endif
#    endif /* EW_PERIODIC	*/     
#   else /* !K3FAST_PG2 */
              dum_s=dum_s
     &                -( thetadiv_nbq(i  ,j,k)
     &                  -thetadiv_nbq(i-1,j,k))
#   endif /* K3FAST_PG2 */

              dum_s=dum_s*0.5*(Hzr(i-1,j,k)+Hzr(i,j,k))*pm_u(i,j)
#   ifdef MASKING
     &                   *umask(i,j)
#   endif  
#  else  /* ! K3FAST_RHO */
              dum_s=0.
#  endif /* K3FAST_RHO */
! !
! !................................
! !  Non-Trad. Coriolis
! !................................
! !
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcoru(i,j,k)
#  endif
! !
! !................................
! !  Body-tide
! !................................
! !
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.5*omega*U0*cos(omega*time)
     &                       *(Hzr(i-1,j,k)+Hzr(i,j,k))
#  endif
! !
! !................................
! !  Bottom-friction
! !................................
! !
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bustr(i,j)
#  endif
! !
! !................................
! !  Nudging
! !................................
! !
#  ifdef NBQ_NUDGING 
              cff=NBQnudgcof(i,j)/dtfast
#  elif defined KNHINT_CORR
              cff=0
#  endif
#  ifdef KNHINT_CORR
               cff=cff+alphaw_nbq/dtfast
     &             *exp(-(z_r(i,j,k)            -z_r(i,j,N))**2
     &                  /(z_r(i,j,N-alphaNw_nbq)-z_r(i,j,N))**2)
#  endif
#  if (defined NBQ_NUDGING && defined NBQCLIMATOLOGY) || defined KNHINT_CORR
              dum_s=dum_s-qdmu_nbq(i,j,k)*cff
     &                   +u(i,j,k,nrhs)
     &                    *0.5*(Hzr(i-1,j,k)+Hzr(i,j,k))
     &                    *cff
#   ifdef MASKING
     &                    *umask(i,j)
#   endif  
#  endif
! !
! !................................
! !  Compute qdmu_nbq
! !................................
! !
#  ifndef K3FAST_SEDLAYERS 
               qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtfast * (
     &                         dum_s
#   ifdef K3FAST_C3D_UVSF         
     &                        +ru_int_nbq(i,j,k) 
#   endif    
     &                        )
#   ifdef MASKING
     &                          *umask(i,j)
#   endif
#  else   /* K3FAST_SEDLAYERS */          
              if (k.gt.0) then   
               qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtfast * (
     &                         dum_s
#   ifdef K3FAST_C3D_UVSF         
     &                        +ru_int_nbq(i,j,k) 
#   endif    
     &                        )
#   ifdef MASKING
     &                          *umask(i,j)
#   endif
              else 
               qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtfast * (
     &                         dum_s
     &                        )
#   ifdef MASKING
     &                          *umask(i,j)
#   endif
              endif
#  endif  /* K3FAST_SEDLAYERS */
#  ifdef K3FAST_SEDLAYERS 
              if (k.gt.0) then
#  endif
#  ifdef K3FAST_ZETAW                  
              DU_nbq(i,j)=DU_nbq(i,j)+qdmu_nbq(i,j,k)
#  endif
              
              if (LAST_FAST_STEP) ru_nbq(i,j,k)=dum_s/work(i,j)               
#  ifdef K3FAST_SEDLAYERS 
              endif
#  endif

            enddo 
          enddo
! !
! !********************************
! !  Fast-mode V-momentum: qdmv_nbq
! !********************************
! !
          do j=JstrV,Jend
            do i=Istr,Iend
! !
! !................................
! !  Comp. pressure gradient & 2nd visc.
! !         v-component: \Sigma  k
! !................................
! !
#  ifdef K3FAST_RHO
! !
! !....1<k<N (-Nsl+1<k<N)..........
! !
#   ifndef K3FAST_SEDLAYERS
              if (k.gt.1      .and.k.lt.N) then 
#   else
              if (k.gt.-N_sl+1.and.k.lt.N) then 
#   endif
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)-z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k)+
     &                  dthetadiv_nbqdz_v(i,j,k-1)) ! dZdy * (d(delta p)dz)_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#   else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *(dthetadiv_nbqdz_v(i,j,k)+
     &                  dthetadiv_nbqdz_v(i,j,k-1)) ! dZdy * (d(delta p)dz)_v
#   endif
! !
! !.....k=N........................
! ! Mathematica
              elseif (k.eq.N) then      
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j  ,k)-z_r(i,j-1,k)) 
     &                *dthetadiv_nbqdz_v(i,j,k-1) 
     &                +(z_w(i,j,N)-z_w(i,j-1,N))        
     &                *dthetadiv_nbqdz_v(i,j,k)  ! BC on dthetadiv_nbqdz_v
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2)
#   else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))            
     &                *dthetadiv_nbqdz_v(i,j,k-1) 
     &                +(z_w(i,j,N)-z_w(i,j-1,N))      
     &                *dthetadiv_nbqdz_v(i,j,k)
#   endif
! !
! !......k=-Nsl+1..................
! !
              else      
#   ifdef NBQ_GRID_SLOW
                if (NSTEP_DS) then
                  dthetadiv_nbqdz(i,j,k,2)=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                                  *2.*dthetadiv_nbqdz_v(i,j,k) 
                endif
                dum_s=dthetadiv_nbqdz(i,j,k,2) 
#   else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k)) 
     &                *dthetadiv_nbqdz_v(i,j,k) 
#   endif
              endif                     ! -------- elseif 1<k<N
! !
! !................................
! ! Horiz. PG
! !................................
! !
#   ifdef K3FAST_PG2
#    ifdef NS_PERIODIC
              dum_s=dum_s
     &             -(gammau  *thetadiv_nbq(i,j  ,k)+
     &               gammau_2*thetadiv_nbq(i,j+1,k)-
     &               gammau  *thetadiv_nbq(i,j-1,k)-
     &               gammau_2*thetadiv_nbq(i,j-2,k)) ! - d(delta p)dy
#    else
#     if defined MPI
         if ((SOUTH_INTER.or.j.ne.JstrV).and.(NORTH_INTER.or.j.ne.Jend)
#      ifdef MASKING
     &      .and.(vmask(i,j-1)*vmask(i,j+1).ne.0.) 
#      endif              
     &      ) then
#     else
         if (j.ne.JstrV.and.j.ne.Jend
#      ifdef MASKING
     &      .and.(vmask(i,j-1)*vmask(i,j+1).ne.0.) 
#      endif              
     &      ) then
#     endif
              dum_s=dum_s
     &             -(gammau  *thetadiv_nbq(i,j  ,k)+
     &               gammau_2*thetadiv_nbq(i,j+1,k)-
     &               gammau  *thetadiv_nbq(i,j-1,k)-
     &               gammau_2*thetadiv_nbq(i,j-2,k)) ! - d(delta p)dy
         else
              dum_s=dum_s
     &             -( thetadiv_nbq(i,j  ,k)
     &               -thetadiv_nbq(i,j-1,k)) 
        endif
#    endif /* NS_PERIODIC	*/     
#   else /* !K3FAST_PG2 */
              dum_s=dum_s
     &                -( thetadiv_nbq(i,j  ,k)
     &                  -thetadiv_nbq(i,j-1,k))
#   endif /* K3FAST_PG2 */
     
              dum_s=dum_s*0.5*(Hzr(i,j-1,k)+Hzr(i,j,k))*pn_v(i,j)
#   ifdef MASKING
     &                   *vmask(i,j)
#   endif  
#  else   /* ! K3FAST_RHO */
              dum_s = 0.
#  endif  /* K3FAST_RHO */
! !
! !................................
! !  Non-Trad. Coriolis
! !................................
! !
#  ifdef UV_COR_NT
              dum_s=dum_s+ntcorv(i,j,k)
#  endif
! !
! !................................
! !  Body-tide
! !................................
! !
#  if defined INTERNAL || defined BODYTIDE
              dum_s=dum_s+0.25*(f(i,j)+f(i,j-1))*
     &                       U0*sin(omega*time)*
     &                       (Hzr(i,j,k)+Hzr(i,j-1,k))
#  endif
! !
! !................................
! !  Bottom-friction
! !................................
! !
#  ifdef BSTRESS_FAST
              if (k.eq.1) dum_s=dum_s-bvstr(i,j)
#  endif
! !
! !................................
! !  Nudging
! !................................
! !   
#  ifdef NBQ_NUDGING 
              cff=NBQnudgcof(i,j)/dtfast
#  elif defined KNHINT_CORR
              cff=0
#  endif
#  ifdef KNHINT_CORR
               cff=cff+alphaw_nbq/dtfast 
     &             *exp(-(z_r(i,j,k)            -z_r(i,j,N))**2
     &                  /(z_r(i,j,N-alphaNw_nbq)-z_r(i,j,N))**2)
#  endif
#  if (defined NBQ_NUDGING && defined NBQCLIMATOLOGY) || defined KNHINT_CORR
              dum_s=dum_s-qdmv_nbq(i,j,k)*cff
     &                   +v(i,j,k,nrhs)*cff
     &                    *0.5*(Hzr(i,j-1,k)+Hzr(i,j,k))
#   ifdef MASKING
     &                    *vmask(i,j)
#   endif   
#  endif   
! !
! !................................
! !  Compute qdmv_nbq
! !................................
! !
#  ifndef K3FAST_SEDLAYERS 
               qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtfast * (
     &                         dum_s
#   ifdef K3FAST_C3D_UVSF         
     &                        +rv_int_nbq(i,j,k) 
#   endif    
     &                        )
#   ifdef MASKING
     &                          *vmask(i,j)
#   endif
#  else   /* ! K3FAST_SEDLAYERS */          
              if (k.gt.0) then   
               qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtfast * (
     &                         dum_s
#   ifdef K3FAST_C3D_UVSF         
     &                        +rv_int_nbq(i,j,k) 
#   endif    
     &                        )
#   ifdef MASKING
     &                          *vmask(i,j)
#   endif
              else 
               qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtfast * (
     &                         dum_s
     &                        )
#   ifdef MASKING
     &                          *vmask(i,j)
#   endif
              endif
#  endif  /* K3FAST_SEDLAYERS */
     
#  ifdef K3FAST_SEDLAYERS 
              if (k.gt.0) then
#  endif
#  ifdef K3FAST_ZETAW  
              DV_nbq(i,j)=DV_nbq(i,j)+qdmv_nbq(i,j,k)
#  endif              
              if (LAST_FAST_STEP) rv_nbq(i,j,k)=dum_s/work(i,j) 
#  ifdef K3FAST_SEDLAYERS 
              endif
#  endif
            enddo
          enddo
!        endif  !<-- k>0
      enddo   !<-- k=0,N   
  
!$acc end kernels
! !
      if (LAST_FAST_STEP) then
!$acc update host( ru_nbq,rv_nbq )   ! iif=last
      endif
# endif   /*   K3FAST_UV  */
! !
! !********************************
! ! Apply point sources for river 
! !    runoff simulations
! !********************************
! !
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
! !********************************
! !  U & V momentum wet mask
! !********************************
! !
# if defined WET_DRY && defined K3FAST_ZETAW
!$acc kernels if(compute_on_device) default(present)
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
!$acc end kernels      
# endif

# ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
# ifdef K3FAST_SEDLAYER      
      call check_tab3d_sedlay(qdmu_nbq,'qdmu_nbqint_a','uint',N_sl+1,N,
     &  ondevice=.TRUE.)
      call check_tab3d_sedlay(qdmv_nbq,'qdmv_nbqint_a','vint',N_sl+1,N,
     &  ondevice=.TRUE.)
# else      
       call check_tab3d(qdmu_nbq,'qdmu_nbqint_a','uint',
     &  ondevice=.TRUE.)
       call check_tab3d(qdmv_nbq,'qdmv_nbqint_a','vint',
     &  ondevice=.TRUE.)
# endif       
c C$OMP END MASTER
# endif
# ifdef K3FAST_UV
! !
! !********************************
! !  U & V momentum open 
! !   boundary conditions
! !********************************
! !
      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
! !
! !********************************
! ! Exchange periodic boundaries 
! !   and computational margins
! !********************************
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI  

!      if ((mod(iif-1,inc_faststep) .eq. inc_faststep-1) .OR.
!     &     (LAST_FAST_STEP)) then
#  ifndef K3FAST_SEDLAYERS
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmu_nbq(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmv_nbq(START_2D_ARRAY,1))
#  else
      call exchange_u3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        qdmu_nbq(START_2D_ARRAY,-N_sl+1))
      call exchange_v3d_sedlay_tile (Istr,Iend,Jstr,Jend,
     &                        qdmv_nbq(START_2D_ARRAY,-N_sl+1))

#  endif
# endif
! !
! !********************************
! ! Debug
! !********************************
! !
# ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
# ifdef K3FAST_SEDLAYER      
      call check_tab3d_sedlay(qdmu_nbq,'qdmu_nbq','uint',-N_sl+1,N,
     &  ondevice=.TRUE.)
      call check_tab3d_sedlay(qdmv_nbq,'qdmv_nbq','vint',-N_sl+1,N,
     &  ondevice=.TRUE.)
# else      
      call check_tab3d(qdmu_nbq,'qdmu_nbq','uint',
     &  ondevice=.TRUE.)
      call check_tab3d(qdmv_nbq,'qdmv_nbq','vint',
     &  ondevice=.TRUE.)
# endif
C$OMP END MASTER
# endif
# endif /* K3FAST_UV */
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_qdmuv_update.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
