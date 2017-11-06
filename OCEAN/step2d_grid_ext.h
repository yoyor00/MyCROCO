!***********************************************************************
! s-grid for NH
!***********************************************************************
!***********************************************************************
! s-grid evolves with external mode (NH)
!    step 1: extrapolation at n+1/2
!***********************************************************************

!#define zwrk UFx
!#define rzeta  UFe
!#define rzeta2  VFe
!#define rzetaSA VFx
!
!-----------------------------------------------------------------------
! Computes grid at m
!-----------------------------------------------------------------------
!

!==================================
# ifdef NBQ_MASS
!==================================
!      do k=1,N
!        do j=JstrV-1,Jend
!          do i=IstrU-1,Iend
!            Hzr(i,j,k)=Hz(i,j,k)
!     &                         /rho_nbq_ext(i,j,k)
!          enddo
!        enddo
!      enddo
!
!  Set lateral boundary conditions for Hz
!
      call hnbq_bc_tile (Istr,Iend,Jstr,Jend)

#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        Hzr(START_2D_ARRAY,1))
#  endif

!==================================
# endif /* NBQ_MASS */
!==================================


!  Set extended range
!
#  ifdef EW_PERIODIC
      imin=Istr-2
      imax=Iend+2
#  else
      if (WESTERN_EDGE) then
        imin=Istr-1
      else
        imin=Istr-2
      endif
      if (EASTERN_EDGE) then
        imax=Iend+1
      else
        imax=Iend+2
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=Jstr-2
      jmax=Jend+2
#  else
      if (SOUTHERN_EDGE) then
        jmin=Jstr-1
      else
        jmin=Jstr-2
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend+1
      else
        jmax=Jend+2
      endif

#  endif    

!
!-----------------------------------------------------------------------
!  Compute other vertical grid variables at m
!-----------------------------------------------------------------------
!	  
      do k=1,N-1
        do j=jmin,jmax
          do i=imin,imax
            Hzw_half_nbq(i,j,k)=z_r(i,j,k+1)-z_r(i,j,k)
          enddo
        enddo
      enddo
      
      do j=jmin,jmax
        do i=imin,imax
          Hzw_half_nbq(i,j,0)=z_r(i,j,1)-z_w(i,j,0)
          Hzw_half_nbq(i,j,N)=z_w(i,j,N)-z_r(i,j,N)
        enddo
      enddo 

!# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
!     &                      z_w(START_2D_ARRAY,0))
!      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
!     &                      z_r(START_2D_ARRAY,1))
!# endif

        call grid_coef_nh(
     &   Istr,Iend,Jstr,Jend,
     &   Hzw_half_nbq_inv,Hzr_half_nbq_inv,
     &   Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &   Hzu_half_qdmu, Hzv_half_qdmv                                     
     &   )

