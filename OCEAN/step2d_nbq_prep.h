!
!-----------------------------------------------------------------------
! Pre-computes Grid coef.
!-----------------------------------------------------------------------
!
#ifdef NBQ_ZETAW
!      if (iic.eq.1.and.iif==1) flag_grid=1
!
!      if (flag_grid.eq.1) then
!       flag_grid=0
!        call grid_coef_nh(
!     &   Istr,Iend,Jstr,Jend,
!     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
!     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
!     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
!     &   )
!     endif
#endif

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      if (FIRST_2D_STEP ) then	
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

#ifndef NBQ_ZETAW
        call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,1),
     &   work3d_nbq(PRIVATE_2D_SCRATCH_ARRAY,1,2)      
!     &   Hzu_half_qdmu, Hzv_half_qdmv                               
     &   )
#endif
!
!-----------------------------------------------------------------------
! DTDRHO...
!-----------------------------------------------------------------------
!
!# if defined NBQ_DTDRHO && !defined NBQ_ZETAW
!        do k=1,N
!          do j=jstrq_nh-1,jendq_nh+1
!            do i=istrq_nh-1,iendq_nh+1
!              hrho_nbq(i,j,k,nnew)=Hzr_half_nbq(i,j,k)*rho(i,j,k)/rho0
!            enddo  
!          enddo  
!        enddo
!        if (iic.le.2) hrho_nbq(:,:,:,nstp)=hrho_nbq(:,:,:,nnew)
!        do k=0,N
!          do j=jstrq_nh-1,jendq_nh+1
!            do i=istrq_nh-1,iendq_nh+1
!              z_nbq(i,j,k,nnew)=zw_half_nbq(i,j,k)
!            enddo
!          enddo  
!        enddo
!        if (iic.eq.1) z_nbq(:,:,:,nstp)=z_nbq(:,:,:,nnew)
!# endif

# ifndef NBQ_IJK
        call mat_mom_nh
        call mat_cont_nh
# endif
!
!-----------------------------------------------------------------------
! Initializes rho_bar
!-----------------------------------------------------------------------    
!
# ifndef NBQ_ZETAW
#  ifdef NBQ_MASS
         do j=jstrq_nh-1,jendq_nh+1
           do i=istrq_nh-1,iendq_nh+1
             rhobar_nbq_int(i,j)= 0.
           enddo
         enddo
		  
         do k=1,N
           do j=jstrq_nh-1,jendq_nh+1
             do i=istrq_nh-1,iendq_nh+1
               rhobar_nbq_int(i,j)= rhobar_nbq_int(i,j)
     &              +(1./rho0)*rho(i,j,k)*Hzr_half_nbq(i,j,k)
             enddo  
           enddo  
         enddo
#  endif
# endif /* NBQ_ZETAW */

!
!-----------------------------------------------------------------------
! Initializes Ru-avg2 variables
!-----------------------------------------------------------------------
!     
# if defined M2FILTER_NONE 
          do k=1,N
            do j=Jstr,Jend
              do i=IstrU,Iend
#  if defined NBQ_IJK
                ru_nbq_avg2(i,j,k)=qdmu_nbq(i,j,k)
#  else
                l_nbq=ijk2lmom_nh(i,j,k,1)
                ru_nbq_avg2(i,j,k)=qdm_nbq_a(l_nbq)
#  endif
              enddo
            enddo 
          enddo
          
     
          do k=1,N
            do j=JstrV,Jend
              do i=Istr,Iend
#  if defined NBQ_IJK
                rv_nbq_avg2(i,j,k)=qdmv_nbq(i,j,k)
#  else
                l_nbq=ijk2lmom_nh(i,j,k,2)
                rv_nbq_avg2(i,j,k)=qdm_nbq_a(l_nbq)
#  endif
              enddo
            enddo 
          enddo
          do k=0,N
            do j=Jstr,Jend
              do i=Istr,Iend
#  if defined NBQ_IJK
                rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
#  else
                l_nbq=ijk2lmom_nh(i,j,k,3)
                rw_nbq_avg2(i,j,k)=qdm_nbq_a(l_nbq)
#  endif
              enddo
            enddo 
          enddo
# endif  

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      endif    ! FIRST_2D_STEP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!# ifdef NBQ
!#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
!      call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,
!     &                                 ru_int_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,
!     &                                 rv_int_nbq(START_2D_ARRAY,1))
!      call exchange_u3d_tile (Istr_nh,Iend_nh,Jstr_nh,Jend_nh,
!     &                                 rw_int_nbq(START_2D_ARRAY,0))
!#  endif     
!# endif     

