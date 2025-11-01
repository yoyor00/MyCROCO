!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
#include "cppdefs.h"

#if defined SOLVE3D && defined BACKSCATTER
!
!=======================================================================

      SUBROUTINE backscatter (tile)

!=======================================================================
!
      IMPLICIT NONE
      INTEGER         :: tile, trd
      INTEGER         :: omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "ocean3d.h"
# include "compute_tile_bounds.h"

      trd=omp_get_thread_num()

      call backscatter_tile ( Istr,Iend,Jstr,Jend,
     &                        A3d(1, 1,trd), A3d(1, 2,trd),
     &                        A2d(1, 1,trd), A2d(1, 2,trd),
     &                        A2d(1, 3,trd), A2d(1, 4,trd))

      RETURN
      END

      SUBROUTINE backscatter_tile ( Istr,Iend,Jstr,Jend,
     &                                 uke_new, rhs_uke,
     &                                    FX,FE,work,FC)
!
!======================================================================
!                 ***  SUBROUTINE  backscatter  ***
!
!  This routine computes antiviscosities for backscattering small-scale
!  numerical dissipation of energy into large-scale resolved motions. 
!  First, it solves the closure equation for the unresolved (subgrid) 
!  kinetic energy, uke:
!
!        	d(uke)/dt = - Cdiss*Ediss - Eback
!
!    with Ediss numerical dissipation of resolved energy (uke source)
!         Eback energy backscattered to resolved motions (uke sink)
!           
!    Eback = Nub * u.L(u)  
!       with the antiviscosity Nub (called visc3d_BS here) :
!
!       	Nub = - Cback * dx * sqrt(uke)
!
!       L(u) is a harmonic diffusion operator used to set the 
!       antiviscous term added to the momentum equation (in uv3dmix_BS):
!
!          	B(u) = Nub * L(u) 
!
! *** References :
!
! Jansen, M.F., Held, I.M., 2014. Parameterizing subgrid-scale eddy 
! effects using energetically consistent backscatter. Ocean Model, 
! 80, 36–48
!
! Juricke S., Danilov S., Kutsenkob A., Oliver M., 2019: Ocean kinetic 
! energy backscatter parametrizations on unstructured grids. Ocean Model,
! 138, 51-67
!
!======================================================================
!
      IMPLICIT NONE
# include "param.h"
# include "mixing.h"
# include "scalars.h"
# include "ocean3d.h"
# include "forces.h"
# include "grid.h"
# include "penalization.h"

      ! Local integers
      INTEGER         ::  Istr,Iend,Jstr,Jend
      INTEGER         ::  i,j,k, ismoo
      INTEGER         ::  imin,imax,jmin,jmax
      ! Local arrays
      REAL            ::  uke_new  (PRIVATE_2D_SCRATCH_ARRAY,1:N)
      REAL            ::  rhs_uke  (PRIVATE_2D_SCRATCH_ARRAY,1:N)
      REAL            ::  FX       (PRIVATE_2D_SCRATCH_ARRAY)
      REAL            ::  FE       (PRIVATE_2D_SCRATCH_ARRAY)
      REAL            ::  work     (PRIVATE_2D_SCRATCH_ARRAY)
      REAL            ::  FC       (PRIVATE_1D_SCRATCH_ARRAY,0:N)
      ! Local scalars
      REAL            ::  Cdiss,Ro,Drate
      REAL            ::  dudx,dvdy,dudy,dvdx, dzdx,dzdy,slope
      REAL            ::  cff,cff1,cff2

      ! Parameter values
      REAL,    PARAMETER ::  Cvisc   = 0.1   ! Antiviscosity coef
      REAL,    PARAMETER ::  Rdiss   = 1.0    ! Filter scale (vs Ro)
                                              ! for KE --> UKE
      INTEGER, PARAMETER ::  Nsmoo   = 1
      REAL,    PARAMETER ::  smoo_a  = 1./12. ! 9-point isotropic 
      REAL,    PARAMETER ::  smoo_b  = 3./16. ! laplacian filter

      REAL,    PARAMETER ::  eps     = 1.0E-12

# define UKE_HADV_UP1     /* CPP keys for */
# define UKE_VADV_UP1     /* advection    */
# define UKE_SMOOTH       /* & smoothing  */
# undef  RUKE_SMOOTH

# include "compute_auxiliary_bounds.h"

# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
!--------------------------------------------------------------------
! Compute extended bounds (needed for UKE smoothing)
!--------------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+1
      imin = Istr-1
      imax = Iend+1
# else
      if (WESTERN_EDGE) then
        imin=Istr
      else
        imin=Istr-1
      endif
      if (EASTERN_EDGE) then
        imax=Iend
      else
        imax=Iend+1
      endif
#  define I_EXT_RANGE imin,imax
# endif
# ifdef NS_PERIODIC
#  define J_EXT_RANGE Jstr-1,Jend+1
      jmin = Jstr-1
      jmax = Jend+1
# else
      if (SOUTHERN_EDGE) then
        jmin=Jstr
      else
        jmin=Jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=Jend
      else
        jmax=Jend+1
      endif
#  define J_EXT_RANGE jmin,jmax
# endif
!
!--------------------------------------------------------------------
! Horizontal advection (first-order upwind)
!--------------------------------------------------------------------
!
!$acc kernels if(compute_on_device) default(present)
      
# ifdef UKE_HADV_UP1
!$acc loop private(FX,FE)
      DO k=1,N
        DO j=jmin,jmax
          DO i=imin,imax+1
            cff=Huon(i,j,k) SWITCH umask(i,j)
            FX(i,j)=uke(i-1,j,k)*max(cff,0.)
     &             +uke(i  ,j,k)*min(cff,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax+1
          DO i=imin,imax
            cff=Hvom(i,j,k) SWITCH vmask(i,j)
            FE(i,j)=uke(i,j-1,k)*max(cff,0.)
     &             +uke(i,j  ,k)*min(cff,0.)
          ENDDO
        ENDDO
        DO j=jmin,jmax
          DO i=imin,imax
            cff=pm(i,j)*pn(i,j)/Hz(i,j,k)
            uke_new(i,j,k)=uke(i,j,k) -dt*cff
     &                       *(  FX(i+1,j)-FX(i,j)
     &                         + FE(i,j+1)-FE(i,j) )
          ENDDO
        ENDDO
      ENDDO
# else
      DO k=1,N
        DO j=jmin,jmax
          DO i=imin,imax
            uke_new(i,j,k)=uke(i,j,k)
          ENDDO
        ENDDO
      ENDDO
# endif /* UKE_HADV_UP1 */
!
!--------------------------------------------------------------------
! Vertical advection (first-order upwind)
!--------------------------------------------------------------------
!
# ifdef UKE_VADV_UP1
      DO j=jmin,jmax
        DO k=2,N
          DO i=imin,imax
            cff=We(i,j,k-1)
            FC(i,k)=uke(i,j,k-1)*max(cff,0.)
     &             +uke(i,j,k  )*min(cff,0.)
          ENDDO
        ENDDO
        DO i=imin,imax
           FC(i,1)=0.
           FC(i,N)=0.
        ENDDO
        DO k=1,N-1
          DO i=imin,imax
            cff=pm(i,j)*pn(i,j)/Hz(i,j,k)
            uke_new(i,j,k)=uke_new(i,j,k)
     &                   - dt*cff*(FC(i,k+1)-FC(i,k))
          ENDDO
        ENDDO
      ENDDO
# endif /* UKE_VADV_UP1 */
!
!--------------------------------------------------------------------
! Compute UKE RHS (source and sink)
!
! Source and sink terms (KE_Hdiss, KE_Hback) are calculated
! in rhs3d and uv3dmix_BS, respectively
!--------------------------------------------------------------------
!
      DO k=1,N
        DO j=jmin,jmax
          DO i=imin,imax
            dudx=(u(i+1,j,k,nstp)-u(i,j,k,nstp))*pm(i,j)
            dvdy=(v(i,j+1,k,nstp)-v(i,j,k,nstp))*pn(i,j)
            dudy=(u(i,j+1,k,nstp)+u(i+1,j+1,k,nstp)
     &           -u(i,j-1,k,nstp)-u(i+1,j-1,k,nstp))*0.25*pn(i,j)
            dvdx=(v(i+1,j,k,nstp)+v(i+1,j+1,k,nstp)
     &           -v(i-1,j,k,nstp)-v(i-1,j+1,k,nstp))*0.25*pm(i,j)
            Drate=sqrt((dudx-dvdy)**2+(dudy+dvdx)**2)

            Ro=Drate/max(eps,f(i,j))   ! Local Rossby nb
            Cdiss=1./(1.+Ro/Rdiss)     ! Frac. of KE diss.
                                       ! going to UKE
            rhs_uke(i,j,k)= Cdiss*
     &         ( 0.5*(KE_Hdiss_U(i,j,k)+KE_Hdiss_U(i+1,j,k))
     &         + 0.5*(KE_Hdiss_V(i,j,k)+KE_Hdiss_V(i,j+1,k)) ) ! source
     &       + ( 0.5*(KE_Hback_U(i,j,k)+KE_Hback_U(i+1,j,k))
     &         + 0.5*(KE_Hback_V(i,j,k)+KE_Hback_V(i,j+1,k)) ) ! sink
# ifdef PERMEABILITY
            rhs_uke(i,j,k)=rhs_uke(i,j,k)*(1.-pena_r(i,j,k))
# endif
          ENDDO
        ENDDO
      ENDDO
!
! Smooth UKE RHS
!
# ifdef RUKE_SMOOTH
#  define WORK rhs_uke
#  define FE1 work
      DO ismoo=1,Nsmoo
#  include "lap_filter.h"
      ENDDO
#  undef FE1
#  undef WORK
# endif
# undef HTRUE
!
!--------------------------------------------------------------------
! Compute time integration
!--------------------------------------------------------------------
!
      DO k=1,N
        DO j=jmin,jmax
          DO i=imin,imax
            cff=pm(i,j)*pn(i,j)
            uke_new(i,j,k)=( Hz_bak(i,j,k)*uke_new(i,j,k)
     &                            - dt*cff*rhs_uke(i,j,k) )
     &                                          / Hz(i,j,k)
# ifdef PERMEABILITY
     &                                  *(1.-pena_r(i,j,k))
# endif
          ENDDO
        ENDDO
      ENDDO
!
!--------------------------------------------------------------------
! Smooth UKE
!--------------------------------------------------------------------
!
# ifdef UKE_SMOOTH
#  define WORK uke_new
#  define FE1 work
      DO ismoo=1,Nsmoo
#  include "lap_filter.h"
      ENDDO
#  undef FE1
#  undef WORK
# endif /* UKE_SMOOTH */

# undef I_EXT_RANGE
# undef J_EXT_RANGE
!
!--------------------------------------------------------------------
! Store new UKE
!--------------------------------------------------------------------
!
      DO k=1,N
        DO j=Jstr,Jend
          DO i=Istr,Iend
            uke(i,j,k)=uke_new(i,j,k)
#  ifdef PERMEABILITY
     &                      *(1.-pena_r(i,j,k))
#  endif
          ENDDO
        ENDDO
      ENDDO
!
!--------------------------------------------------------------------
! Compute visc3d_BS 
!--------------------------------------------------------------------
!
!$acc kernels if(compute_on_device) default(present)

      DO k=1,N
        DO j=Jstr,Jend
          DO i=Istr,Iend
             cff1=sqrt(MAX(2.*uke(i,j,k),0.))
             cff2=sqrt(om_r(i,j)*on_r(i,j))
             visc3d_BS_r(i,j,k) = - Cvisc*cff1*cff2
     &                                  SWITCH rmask(i,j)
          ENDDO
        ENDDO
      ENDDO
!
!------------------------------------------------
! Apply boundary conditions
!------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,N
            uke(Istr-1,j,k)=uke(Istr,j,k)
            visc3d_BS_r(Istr-1,j,k)=visc3d_BS_r(Istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr,Jend
          do k=1,N
            uke(Iend+1,j,k)=uke(Iend,j,k)
            visc3d_BS_r(Iend+1,j,k)=visc3d_BS_r(Iend,j,k)
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,N
            uke(i,Jstr-1,k)=uke(i,Jstr,k)
            visc3d_BS_r(i,Jstr-1,k)=visc3d_BS_r(i,Jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr,Iend
          do k=1,N
            uke(i,Jend+1,k)=uke(i,Jend,k)
            visc3d_BS_r(i,Jend+1,k)=visc3d_BS_r(i,Jend,k)
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          uke(Istr-1,Jstr-1,k)=uke(Istr,Jstr,k)
          visc3d_BS_r(Istr-1,Jstr-1,k)=visc3d_BS_r(Istr,Jstr,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          uke(Istr-1,Jend+1,k)=uke(Istr,Jend,k)
          visc3d_BS_r(Istr-1,Jend+1,k)=visc3d_BS_r(Istr,Jend,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          uke(Iend+1,Jstr-1,k)=uke(Iend,Jstr,k)
          visc3d_BS_r(Iend+1,Jstr-1,k)=visc3d_BS_r(Iend,Jstr,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          uke(Iend+1,Jend+1,k)=uke(Iend,Jend,k)
          visc3d_BS_r(Iend+1,Jend+1,k)=visc3d_BS_r(Iend,Jend,k)
        enddo
      endif
#  endif
# endif
!
! Viscosity at PSI points
!
      do k=1,N
        do j=Jstr,JendR
          do i=Istr,IendR
            visc3d_BS_p(i,j,k)=0.25*
     &            ( visc3d_BS_r(i,j  ,k)+visc3d_BS_r(i-1,j  ,k)
     &             +visc3d_BS_r(i,j-1,k)+visc3d_BS_r(i-1,j-1,k) )
          enddo
        enddo
      enddo
!
!$acc end kernels
!
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
!
# if defined BD_EXCHANGE
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend, uke)
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend, visc3d_BS_r)
      call exchange_p3d_tile (Istr,Iend,Jstr,Jend, visc3d_BS_p)
# endif
!======================================================================
      return
      end
!======================================================================
#else
      subroutine backscatter_empty (tile)
      end
#endif  /* SOLVE3D & BACKSCATTER */

