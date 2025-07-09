! $Id: gls_mixing.F 1524 2014-04-14 17:00:06Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
#include "cppdefs.h"

#if defined SOLVE3D && defined TKE3D_MIXING
!
# define TKE_HADV_UP1
# define TKE_VADV_UP1
# define TKE_SFC_NEUMANN
# define TKE_BOT_NEUMANN
!
!=======================================================================
      SUBROUTINE tke3d_mixing (tile)
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

      call tke_mixing_tile ( Istr, Iend, Jstr, Jend,
     &                    A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd),
     &                    A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd))

      RETURN
      END

!
!===============================================================================
      SUBROUTINE tke_mixing_tile ( Istr, Iend, Jstr, Jend,
     &              ustar_sfc_sq,ustar_bot_sq, DC, FC, CF, RH)
!===============================================================================
!
!
!===============================================================================
!                 ***  SUBROUTINE  tke3d_mixing  ***
!* Compute TKE : the corresponding prognostic equation is
!  d TKE/dt = - d/dx(u TKE) - d/dy(v TKE) + d/dz(Kz d/dz(TKE)) + Sprod3d - Bprod - espilon
!       with Sprod3d a 3D TKE production term (computed in uv3dmix)
!            Bprod a TKE destruction/production by stratification
!       end epsilon a dissipation paramterized as epsilon = ceps TKE sqrt(TKE) / Leps
!       Horizontal advection is achieved using a first-order upwind scheme
!* The boundary conditions for the vertical diffusion of TKE are homogeneous Neumann
!    at top and bottom
!
!* Future developments : at this stage the following terms are missing and/or could be added
!    - Vertical advection by the mean vertical velocity
!    - Explicit horizontal diffusion of TKE (achieved implicitly by the 1st-order upwind for the moment)
!    - Dirichlet boundary conditions
!    - Non-homogeneous Neumann BC to account for wave effects (e.g. Craig and Banner 1994)
!    - More advanced calculation of Turbulent Prandtl number ...
!
!* Note that to simplify the computation of the shear production term and
!         the use of a FV approach to TKE advection, TKE is cell-centered
!
!  References :
!     - Deardorff, J.W., 1980. Stratocumulus-capped mixing layers derived from a
!           three-dimensional model. Bound.-Lay. Meteorol. 18, 495–527.
!
!     - Sullivan, P., McWilliams, J., Moeng, C., 1994. A subgrid-scale model for
!           large-eddy simulation of planetary boundary-layer flows 71, 247–276.
!
!     - Zhang, X., Bao, J.W., Chen, B., Grell, E.D., 2018. A three-dimensional
!            scale-adaptive turbulent kinetic energy scheme in the wrf-arw model.
!            Mon. Weather Rev. 146(7), 2023 – 2045
!===============================================================================
!
      IMPLICIT NONE
# include "param.h"
# include "mixing.h"
# include "scalars.h"

      ! Local integers
      INTEGER         ::   Istr, Iend, Jstr, Jend
      INTEGER         ::   i,       j,    k, kref
      ! Local arrays
      REAL            ::  diss        (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  ustar_sfc_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  ustar_bot_sq(PRIVATE_2D_SCRATCH_ARRAY      )
      REAL            ::  DC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL            ::  FC          (PRIVATE_1D_SCRATCH_ARRAY,0:N  )
      REAL            ::  CF          (PRIVATE_1D_SCRATCH_ARRAY,1:N-1)
      REAL            ::  RH          (PRIVATE_1D_SCRATCH_ARRAY,1:N  )
      ! Local scalars
      REAL            ::  Ake, bvf_r
      REAL            ::  cff , cff1 , cff2
      REAL            ::  Bprod, Sprod
      REAL            ::  cm0inv2,invPrt
      REAL            ::  flux_top, flux_bot
      REAL            ::  trb_sfc, trb_bot, ld80, rn2
      REAL            ::  HUon_w, HVom_w, delta
      REAL            ::  su_r,sv_r,dissl,wstar2
      ! Parameter values
      REAL, PARAMETER ::  tke_min =  1.0E-08 ! min value for TKE
      REAL, PARAMETER ::  rsmall  =  1.0E-16 ! minimum stratification
      REAL, PARAMETER ::  nuws    =  1.0E-06 ! background value for diffusivities
      REAL, PARAMETER ::  nuwm    =  1.0E-06 ! background value for viscosities
      REAL, PARAMETER ::  oneThird =  1./3.
!
!--------------------------------------------------
# include "ocean3d.h"
# include "forces.h"
# include "grid.h"
# include "averages.h"
# include "compute_auxiliary_bounds.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif
!
!-------------------------------------------------------------------------------
! Initialization of the mixing length
!-------------------------------------------------------------------------------
!
      DO k = 1,N
        DO j = jstr,jend
          DO i = istr,iend
            ! harmonic avergaing to estimate a cell-centered value of bvf
            cff = bvf(i,j,k)*bvf(i,j,k-1)
            IF(cff > 0.D0) THEN
              bvf_r=cff/(bvf(i,j,k)+bvf(i,j,k-1))
            ELSE
              bvf_r=0.D0
            ENDIF
            rn2             = MAX( bvf_r, rsmall )
            ! Deardorff 1980 mixing length
            ld80            = SQRT( 2.*tke(i,j,k,nstp)/rn2 )
            ! Grid box size (should be recomputed as Hz changes over time)
            delta           = ( om_r(i,j)*on_r(i,j)*Hz(i,j,k) )**oneThird
            ! Length-scale as in Sullivan et al. (1994) or Zhang et al. (2018)
            Lscale(i,j,k)   = MIN( delta, ld80 )
          ENDDO
        ENDDO
      ENDDO
!
!-------------------------------------------------------------------------------
! Horizontal advection of TKE (1st-order upwind and no explicit diffusion)
!-------------------------------------------------------------------------------
!
#  define FX ustar_sfc_sq
#  define FE ustar_bot_sq
      DO k=1,N
        !==
        DO j=jstr,jend
          DO i=istr,iend+1
#ifdef TKE_HADV_UP1
            IF( Huon(i,j,k) > 0. ) THEN
              FX(i,j)=  tke(i-1,j,k,nstp)*Huon(i,j,k) SWITCH umask(i,j)
            ELSE
              FX(i,j)=  tke(i  ,j,k,nstp)*Huon(i,j,k) SWITCH umask(i,j)
            ENDIF
#else
            FX(i,j) = 0.
#endif
          ENDDO
        ENDDO
        !==
        DO j=jstr,jend+1
          DO i=istr,iend
#ifdef TKE_HADV_UP1
            IF( Hvom(i,j,k) > 0. ) THEN
              FE(i,j)=  tke(i,j-1,k,nstp)*Hvom(i,j,k) SWITCH vmask(i,j)
            ELSE
              FE(i,j)=  tke(i,j  ,k,nstp)*Hvom(i,j,k) SWITCH vmask(i,j)
            ENDIF
#else
            FE(i,j)= 0.
#endif
          ENDDO
        ENDDO
        !==
        DO j=jstr,jend
          DO i=istr,iend
            tke(i,j,k,nnew) = Hz_bak(i,j,k)*tke(i,j,k,nstp) - dt*pm(i,j)*pn(i,j)
     &                    *(  FX(i+1,j)-FX(i,j) + FE(i,j+1)-FE(i,j) )           ! tke contains Hz x TKE at this point
          ENDDO
        ENDDO
        !==
      ENDDO
#undef FX
#undef FE
!
!--------------------------------------------------------------------
!  Compute ustar squared at the surface and at the bottom
!      useful for Dirichlet boundary conditions under the assumption
!             of logarithmic boundary layers (yet to be implemented)
!--------------------------------------------------------------------
!
#if defined TKE_SFC_NEUMANN || defined TKE_BOT_NEUMANN

#else
      DO j=jstr,jend
        DO i=istr,iend
          su_r=0.5*(sustr(i,j)+sustr(i+1,j))
          sv_r=0.5*(svstr(i,j)+svstr(i,j+1))
          ustar_sfc_sq( i, j ) = SQRT(su_r**2+sv_r**2)
          ustar_bot_sq( i, j ) =
     &          SQRT( (0.5*(bustr(i,j)+bustr(i+1,j)))**2
     &               +(0.5*(bvstr(i,j)+bvstr(i,j+1)))**2 )
         ENDDO
      ENDDO
#endif
!
!--------------------------------------------------
! Tridiagonal inversion
!--------------------------------------------------
!
      DO j=jstr,jend         !<-- j-outer loop
        ! Off-diagonal terms for the tridiagonal problem
        cff=-0.5*dt*(tke_ce/tke_cm)
        DO k=1,N-1
          DO i=istr,iend
            FC(i,k)= cff*(Akv(i,j,k+1)+Akv(i,j,k)) /
     &                         (z_r(i,j,k+1)-z_r(i,j,k))
          ENDDO
        ENDDO
        ! Top and bottom boundary conditions
        DO i=istr,iend
# ifdef TKE_SFC_NEUMANN
          FC(i,N)=0.
# endif
# ifdef TKE_BOT_NEUMANN
          FC(i,0)=0.
# endif
        ENDDO
        ! Production/Dissipation terms and diagonal term
        DO k=1,N
          DO i=istr,iend
            ! Shear and buoyancy production
            Sprod =  Sprod3d(i,j,k)
# if defined TEMPERATURE || defined SALINITY
            ! Harmonic averaging to determine bvf
            cff = bvf(i,j,k)*bvf(i,j,k-1)
            IF(cff > 0.D0) THEN
              bvf_r=cff/(bvf(i,j,k)+bvf(i,j,k-1))
            ELSE
              bvf_r=0.D0
            ENDIF
            Bprod = -Akt(i,j,k,itemp)*bvf_r
# else
            Bprod = 0.
# endif
            RH(i,k) = tke(i,j,k,nnew) + dt*Hz(i,j,k)*(Bprod+Sprod)
            dissl   = tke_ceps*SQRT(tke(i,j,k,nstp))/Lscale(i,j,k)      ! <-- parameterization of the epsilon term
            DC(i,k) = Hz(i,j,k)*(1.+dt*dissl) - FC(i,k) - FC(i,k-1)     !<-- diagonal term
          ENDDO
        ENDDO
        ! tridiagonal resolution
        DO i=istr,iend
          cff       =  1./DC(i,1)
          CF(i,1)   = cff*FC(i,1)  ! q(1) = c(1)/b(1)
          RH(i,1)   = cff*RH(i,1)  ! f(1) = f(1)/b(1)
        ENDDO

        DO k=2,N
          DO i=istr,iend
            cff     = 1./(DC(i,k)-CF(i,k-1)*FC(i,k-1))     ! p = 1/(b(k)+a(k)*q(k-1))
            CF(i,k) = cff*FC(i,k)                          ! q(k) = c(k)*p
            RH(i,k) = cff*( RH(i,k)-FC(i,k-1)*RH(i,k-1) )  ! f(k) = ( f(k)-a(k)*f(k-1) )*p
          ENDDO
        ENDDO

        DO i=istr,iend
          tke(i,j,N,nnew) = MAX( RH(i,N), tke_min )
        ENDDO

        DO k=N-1,1,-1
          DO i=istr,iend
            RH(i,k)         = RH(i,k)-CF(i,k)*RH(i,k+1)
            tke(i,j,k,nnew) = MAX( RH(i,k), tke_min )
          ENDDO
        ENDDO
      !--------------------------------------------------
      ENDDO     ! j loop
      !--------------------------------------------------

!
!--------------------------------------------------
! Finalize by calculating Akv/Akt
!--------------------------------------------------
!
      !--------------------------------------------------
      DO j=jstr,jend  !<-- j-outer loop
      !--------------------------------------------------
        DO k=1,N ! Since TKE is cell-centered, Akv and Akt are also cell-centered
          DO i=istr,iend
            ! Compute inverse of Turbulent Prandtl number (as in Sullivan et al. (1994))
            delta      = ( om_r(i,j)*on_r(i,j)*Hz(i,j,k) )**oneThird
            invPrt     = 1.+2.*(Lscale(i,j,k)/delta)
            !
            cff        = tke_cm * SQRT(tke(i,j,k,nnew)) * Lscale(i,j,k)
            Akv(i,j,k) = MAX( cff , nuwm) SWITCH rmask(i,j)
# ifdef TEMPERATURE
            Akt(i,j,k,itemp)= MAX( invPrt * cff, nuws ) SWITCH rmask(i,j)
# endif
# ifdef SALINITY
            Akt(i,j,k,isalt)= MAX( invPrt * cff, nuws ) SWITCH rmask(i,j)
# endif
#ifdef MASKING
            Lscale(i,j,k) = Lscale(i,j,k) * rmask(i,j)
#endif
          ENDDO
        ENDDO
      !--------------------------------------------------
      ENDDO  !<-- end j-outer loop
      !--------------------------------------------------
!
!------------------------------------------------
!  Apply lateral boundary conditions
!------------------------------------------------
!
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          do k=1,N
            tke(istr-1,j,k,nnew)=tke(istr,j,k,nnew)
            Akv(istr-1,j,k     )=Akv(istr,j,k    )
#  ifdef TEMPERATURE
            Akt(istr-1,j,k,itemp)=Akt(istr,j,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(istr-1,j,k,isalt)=Akt(istr,j,k,isalt)
#  endif
            Lscale(istr-1,j,k)=Lscale(istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          do k=1,N
            tke(iend+1,j,k,nnew)=tke(iend,j,k,nnew)
            Akv(iend+1,j,k      )=Akv(iend,j,k      )
#  ifdef TEMPERATURE
            Akt(iend+1,j,k,itemp)=Akt(iend,j,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(iend+1,j,k,isalt)=Akt(iend,j,k,isalt)
#  endif
            Lscale(iend+1,j,k)=Lscale(iend,j,k)
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          do k=1,N
            tke(i,jstr-1,k,nnew)=tke(i,jstr,k,nnew)
            Akv(i,jstr-1,k     )=Akv(i,jstr,k     )
#  ifdef TEMPERATURE
            Akt(i,jstr-1,k,itemp)=Akt(i,jstr,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(i,jstr-1,k,isalt)=Akt(i,jstr,k,isalt)
#  endif
            Lscale(i,jstr-1,k)=Lscale(i,jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          do k=1,N
            tke(i,jend+1,k,nnew)=tke(i,jend,k,nnew)
            Akv(i,jend+1,k)=Akv(i,jend,k)
#  ifdef TEMPERATURE
            Akt(i,jend+1,k,itemp)=Akt(i,jend,k,itemp)
#  endif
#  ifdef SALINITY
            Akt(i,jend+1,k,isalt)=Akt(i,jend,k,isalt)
#  endif
            Lscale(i,jend+1,k)=Lscale(i,jend,k)
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          tke(istr-1,jstr-1,k,nnew)=tke(istr,jstr,k,nnew)
          Akv(istr-1,jstr-1,k      )=Akv(istr,jstr,k      )
#  ifdef TEMPERATURE
          Akt(istr-1,jstr-1,k,itemp)=Akt(istr,jstr,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(istr-1,jstr-1,k,isalt)=Akt(istr,jstr,k,isalt)
#  endif
          Lscale(istr-1,jstr-1,k)=Lscale(istr,jstr,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          tke(istr-1,jend+1,k,nnew)=tke(istr,jend,k,nnew)
          Akv(istr-1,jend+1,k      )=Akv(istr,jend,k      )
#  ifdef TEMPERATURE
          Akt(istr-1,jend+1,k,itemp)=Akt(istr,jend,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(istr-1,jend+1,k,isalt)=Akt(istr,jend,k,isalt)
#  endif
          Lscale(istr-1,jend+1,k)=Lscale(istr,jend,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          tke(iend+1,jstr-1,k,nnew)=tke(iend,jstr,k,nnew)
          Akv(iend+1,jstr-1,k      )=Akv(iend,jstr,k      )
#  ifdef TEMPERATURE
          Akt(iend+1,jstr-1,k,itemp)=Akt(iend,jstr,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(iend+1,jstr-1,k,isalt)=Akt(iend,jstr,k,isalt)
#  endif
          Lscale(iend+1,jstr-1,k)=Lscale(iend,jstr,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          tke(iend+1,jend+1,k,nnew)=tke(iend,jend,k,nnew)
          Akv(iend+1,jend+1,k)=Akv(iend,jend,k)
#  ifdef TEMPERATURE
          Akt(iend+1,jend+1,k,itemp)=Akt(iend,jend,k,itemp)
#  endif
#  ifdef SALINITY
          Akt(iend+1,jend+1,k,isalt)=Akt(iend,jend,k,isalt)
#  endif
          Lscale(iend+1,jend+1,k)=Lscale(iend,jend,k)
        enddo
      endif
#  endif
# endif
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- -----------
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        tke(START_2D_ARRAY,1,nnew))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Akv)
!      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Shear3d)
#  ifdef TEMPERATURE
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,itemp))
#  endif
#  ifdef SALINITY
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Akt(START_2D_ARRAY,0,isalt))
#  endif
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        Lscale(START_2D_ARRAY,1))
# endif
!======================================================================
      return
      end
!======================================================================
#endif  /* SOLVE3D & TKE3D_MIXING */
