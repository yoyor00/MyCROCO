#include "cppdefs.h"

#if defined DIAGNOSTICS_DIAPYCNAL
      subroutine compute_buoyancy_gradient (tile)
!
      implicit none
      integer tile, trd, omp_get_thread_num
#  include "param.h"
#  include "private_scratch.h"
#  include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call compute_buoyancy_gradient_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine compute_buoyancy_gradient_tile (Istr,Iend,Jstr,Jend)
!
!======================================================================
! Compute adiabatic buoyancy gradients
! used to diagnose diapycnal fluxes and diapycnal velocity 
!
!====================================================================== 
!
      implicit none
#  include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j,k, kbl_tmp,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc
# endif
      real cff,cff1,cff2,dpth
#  include "grid.h"
#  include "ocean3d.h"
#  include "coupling.h"
#  include "scalars.h"
#  include "diags_pv.h"
!
#  include "compute_auxiliary_bounds.h"
!#  include "compute_extended_bounds.h"

!
#  include "compute_buoyancy_gradient.h"

      
      return
      end
#else
      subroutine compute_buoyancy_gradient_empty
      end
#endif /* DIAGNOSTICS_DIAPYCNAL*/


!======================================================================

#if defined DIAGNOSTICS_PV

      subroutine set_diags_pv (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_pv_tile (istr,iend,jstr,jend,tile)
      return
      end

      subroutine set_diags_pv_tile (istr,iend,jstr,jend,tile)
!
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
      real eps
      parameter (eps=1.D-20)
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "diagnostics.h"

# if defined DIAGNOSTICS_DIAPYCNAL
#  include "mixing.h"
# endif

# ifdef DIAGNOSTICS_EK
#  include "diags_ek.h"
# endif

# include "diags_pv.h"

      real alfa(PRIVATE_2D_SCRATCH_ARRAY),
     &     beta(PRIVATE_2D_SCRATCH_ARRAY)

# if defined DIAGNOSTICS_DIAPYCNAL
      real epsil,dpth,cfr,
     &     dR(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     FC(PRIVATE_2D_SCRATCH_ARRAY),
     &     rx(PRIVATE_2D_SCRATCH_ARRAY),
     &     dZx(PRIVATE_2D_SCRATCH_ARRAY),
     &     gradbuoy(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     TF_xHmix_rot(PRIVATE_2D_SCRATCH_ARRAY,N,2),
     &     TF_yHmix_rot(PRIVATE_2D_SCRATCH_ARRAY,N,2),
     &     TF_zHmix_rot(PRIVATE_2D_SCRATCH_ARRAY,N,2),
     &     TF_zVmix_rot(PRIVATE_2D_SCRATCH_ARRAY,N,2)
      parameter (epsil=0.) !1.E-10)

# endif
# include "compute_auxiliary_bounds.h"

# if ! defined DIAGNOSTICS_DIAPYCNAL
        do itrc=1,NTA
          do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
                Trhs(i,j,k,itrc) = TForc(i,j,k,itrc)
     &                           + THmix(i,j,k,itrc)
     &                           + TVmix(i,j,k,itrc)
              enddo
            enddo
          enddo
        enddo

      do k=1,N
          call alfabeta_k_tile(Istr,Iend,Jstr,Jend,k,alfa,beta)
          do j=JstrR,JendR
              do i=IstrR,IendR
                  Trhs(i,j,k,1) = -alfa(i,j)*Trhs(i,j,k,1)
#  ifdef SALINITY
                  Trhs(i,j,k,2) = beta(i,j)*Trhs(i,j,k,2)
#  endif
              enddo
          enddo
      enddo

# else





! compute horizontal dbdx, dbdy, dbdz (on rho-grid)
! # include "compute_buoyancy_gradient.h"

! product of fluxes and buoyancy gradients (on sigma-levels)
# include "compute_diapycnal_fluxes.h"

! compute horizontal dbdx, dbdy, dbdz (on rho-grid)
# include "compute_buoyancy_gradient_horizontal.h"

      do k=1,N

          call alfabeta_k_tile(Istr,Iend,Jstr,Jend,k,alfa,beta)

          do j=Jstr,Jend
              do i=Istr,Iend

                  gradbuoy(i,j,k) = sqrt( dbdx(i,j,k)**2
     &                                  + dbdy(i,j,k)**2
     &                                  + dbdz(i,j,k)**2 )




                  ! r.h.s. of buoyancy equation
                  Trhs(i,j,k,1)=  (  alfa(i,j)*TVmix(i,j,k,1)
     &                             + alfa(i,j)*TForc(i,j,k,1)
     &                             + alfa(i,j)*THmix(i,j,k,1)
# ifdef SALINITY
     &                             - beta(i,j)*TVmix(i,j,k,2)
     &                             - beta(i,j)*TForc(i,j,k,2)
     &                             - beta(i,j)*THmix(i,j,k,2)
# endif
     &                             ) * g

                  ! divide by the density gradient amplitude to get diapycnal velocity
                  ! the density gradient is computed in prsgrd.F
                  Trhs(i,j,k,2)=  Trhs(i,j,k,1)
     &                         / max(eps, gradbuoy(i,j,k) )


                  ! rate of buoyancy equation
                  Trhs(i,j,k,3)=  (  alfa(i,j)*Trate(i,j,k,1)
# ifdef SALINITY
     &                             - beta(i,j)*Trate(i,j,k,2)
# endif
     &                             ) * g

                  ! rate of buoyancy equation
                  Trhs(i,j,k,4)=  (  alfa(i,j)*TXadv(i,j,k,1)
     &                             + alfa(i,j)*TYadv(i,j,k,1)
     &                             + alfa(i,j)*TVadv(i,j,k,1)
# ifdef SALINITY
     &                             - beta(i,j)*TXadv(i,j,k,2)
     &                             - beta(i,j)*TYadv(i,j,k,2)
     &                             - beta(i,j)*TVadv(i,j,k,2)
# endif
     &                             ) * g


                  bflux(i,j,k) = - g * (
     &                            alfa(i,j)* TF_xHmix_rot(i,j,k,1)
     &                          + alfa(i,j)* TF_yHmix_rot(i,j,k,1)
     &                          + alfa(i,j)* TF_zHmix_rot(i,j,k,1)
     &                          + alfa(i,j)* TF_zVmix_rot(i,j,k,1)
# ifdef SALINITY
     &                          - beta(i,j)* TF_xHmix_rot(i,j,k,2) 
     &                          - beta(i,j)* TF_yHmix_rot(i,j,k,2) 
     &                          - beta(i,j)* TF_zHmix_rot(i,j,k,2)
     &                          - beta(i,j)* TF_zVmix_rot(i,j,k,2)
# endif
     &                              )


                  diffusivity(i,j,k)  = bflux(i,j,k)
     &                          / max(eps, gradbuoy(i,j,k)**2 )

                  
                  ! debug
                  !Trhs(i,j,k,1) =  TF_xHmix(i,j,k,1)
                  !Trhs(i,j,k,2) =  TF_yHmix(i,j,k,1)
                  !Trhs(i,j,k,3) =  TF_zHmix(i,j,k,1)
                  !Trhs(i,j,k,4) =  TF_zVmix(i,j,k,1)
                  !Trhs(i,j,k,5) =  TF_xHmix(i,j,k,2)
                  !Trhs(i,j,k,6) =  TF_yHmix(i,j,k,2)
                  !Trhs(i,j,k,7) =  TF_zHmix(i,j,k,2)
                  !Trhs(i,j,k,8) =  TF_zVmix(i,j,k,2)

                  !bflux(i,j,k)     =  gradbuoy(i,j,k)

                  !diffusivity(i,j,k)    =  TF_zVmix(i,j,k,1)



              enddo
          enddo
      enddo
# endif /* DIAGNOSTICS_DIAPYCNAL */



      return
      end
#else /* DIAGNOSTICS_PV*/
      subroutine set_diags_pv_empty
      end
#endif /* DIAGNOSTICS_PV */










