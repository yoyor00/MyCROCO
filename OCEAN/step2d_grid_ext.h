!***********************************************************************
! s-grid for NH
!***********************************************************************
!***********************************************************************
! s-grid evolves with external mode (NH)
!    step 1: extrapolation at n+1/2
!***********************************************************************

#define zwrk UFx
#define rzeta  UFe
#define rzeta2  VFe
#define rzetaSA VFx
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
c LAURENT: should be put in an OBC tile -with zero gradient) routine for Hzr (or Hz)
c since it is used at several places in the code

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=Jstr,Jend
            Hzr(0,j,k)=Hzr(1,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=Jstr,Jend
            Hzr(LOCALLM+1,j,k)=Hzr(LOCALLM,j,k)
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=Istr,Iend
            Hzr(i,0,k)=Hzr(i,1,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=Istr,Iend
            Hzr(i,LOCALMM+1,k)=Hzr(i,LOCALMM,k)
          enddo
        enddo
      endif
#  endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          Hzr(0,0,k)=Hzr(1,1,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          Hzr(0,LOCALMM+1,k)=Hzr(1,LOCALMM,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          Hzr(LOCALLM+1,0,k)=Hzr(LOCALLM,1,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          Hzr(LOCALLM+1,LOCALMM+1,k)=
     &                     Hzr(LOCALLM,LOCALMM,k)
        enddo
      endif
#  endif
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
#if TOTO
c LAURENT: next lines have probably to be removed
c  step2d_grid_ext.h is always called after a set_depth
c where z_r, z_w are properly computed
      do j=jmin,jmax
        do i=imin,imax
           z_w(i,j,0)=-h(i,j)
           z_r(i,j,1)=z_w(i,j,0)
     &                        +0.5*Hzr(i,j,1)
           z_w(i,j,1)   =z_w(i,j,0)+Hzr(i,j,1)
        enddo
      enddo
      do k=2,N
        do j=jmin,jmax
          do i=imin,imax
            z_w(i,j,k)=z_w(i,j,k-1)+Hzr(i,j,k)
            z_r(i,j,k)=z_r(i,j,k-1)
     &                 +0.5*(Hzr(i,j,k)+Hzr(i,j,k-1))
          enddo
        enddo
      enddo
#endif
	  
	  
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

