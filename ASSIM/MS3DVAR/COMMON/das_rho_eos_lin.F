#include "cppdefs.h"

      subroutine das_lin_rho_eos (tile)
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "das_param.h"
# include "das_private_scratch.h"
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call das_lin_rho_eos_tile (Istr,Iend,Jstr,Jend, 
     &                          A2d(1,1,trd), A2d(1,2,trd),
     &                          A2d(1,3,trd), A2d(1,4,trd),
     &                  A2d(1,5,trd), A2d(1,6,trd),
     &                  A2d(1,7,trd), A2d(1,8,trd))
      return
      end

      subroutine das_lin_rho_eos_tile (Istr,Iend,Jstr,Jend, 
     &                                   rho1,K0,K1,K2,
     &                   rho1_adj,K0_adj,K1_adj,K2_adj)
!
! Computes density increment via Equation Of State (EOS) for seawater.
! If so prescribed, linearized non-linear EOS of Jackett and McDougall (1995)
! is used.
!
! Tt potential temperature [deg Celsius].
! Ts salinity [PSU].
! Tz pressure/depth, [depth in meters and negative].
!
! K0, K1 and K2 are the pressure polynomial coefficients for secant
! bulk modulus, so that
!
!               bulk = K0 - K1 * z + K2 * z**2 ;
!
! while rho1 is sea-water density [kg/m^3] at standard pressure
! of 1 Atm, so that the density anomaly at in-sity pressure is
!
!               rho = rho1 / (1 + z / bulk) - 1000
! 
! If so prescribed, it also computes the Brunt-Vaisala frequency
! [1/s^2] at horizontal RHO-points and vertical W-points,
!
!                   bvf = - g/rho0 d(rho)/d(z).
!
!  In computation of bvf the density anomaly difference is computed
!  by adiabatically lowering/rising the water parcel from RHO point
!  above/below to the W-point depth at "z_w".
!
!  Reference:
!
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of
!  Hydrostatic Profiles to Achieve Static Stability, Journ of Atmos.
!  and Oceanic Techn., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
!    and McDougall (1992), unpublished manuscript, CSIRO, Australia.
!    It computes in-situ density anomaly as a function of potential
!    temperature (Celsius) relative to the surface, salinity (PSU),
!    and depth (meters).  It assumes  no  pressure  variation along 
!    geopotential  surfaces,  that  is,  depth  and  pressure  are  
!    interchangeable. >>                                             
!                                          John Wilkin, 29 July 92   
!
      implicit none
# include "param.h"
# include "das_param.h"
      integer Istr,Iend,Jstr,Jend, i,j,k 
      real rho1(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K0(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K1(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K2(PRIVATE_1D_SCRATCH_ARRAY,NDAS)
      real rho1_adj(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K0_adj(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K1_adj(PRIVATE_1D_SCRATCH_ARRAY,NDAS),
     &       K2_adj(PRIVATE_1D_SCRATCH_ARRAY,NDAS)
# include "grid.h"
# include "das_ocean.h"
# include "das_ocean9.h"
# include "scalars.h"
!
# ifdef NONLIN_EOS
      real A00, A01, A02, A03, A04, A10, A11, A12, A13, Tt,
     &     AS0, AS1, AS2, B00, B01, B02, B03, B10, B11,   Ts,
     &     B12, BS1, E00, E01, E02, E10, E11, E12, sqrtTs, cff
      real Tt_adj, Ts_adj, sqrtTs_adj
      parameter(A00=+19092.56 ,  A01=+209.8925   , A02=-3.041638 ,
     &          A03=-1.852732e-3, A04=-1.361629e-5, A10=104.4077   ,
     &          A11=-6.500517   , A12=+0.1553190  , A13=2.326469e-4 ,
     &          AS0=-5.587545   , AS1=+0.7390729  , AS2=-1.909078e-2,
     &          B00=+4.721788e-1, B01=+1.028859e-2, B02=-2.512549e-4,
     &          B03=-5.939910e-7, B10=-1.571896e-2, B11=-2.598241e-4,
     &          B12=+7.267926e-6, BS1=+2.042967e-3,
     &          E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,
     &          E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
      real QR, Q01, Q02, Q03, Q04, Q05, Q10, Q11,
     &         Q12, Q13, Q14, QS0, QS1, QS2, Q20
      parameter(QR=+999.842594 , Q01=+6.793952e-2, Q02=-9.095290e-3,
     &          Q03=+1.001685e-4, Q04=-1.120083e-6, Q05=+6.536332e-9,
     &          Q10=+0.824493   , Q11=-4.08990e-3 , Q12=+7.64380e-5 ,
     &          Q13=-8.24670e-7 , Q14=+5.38750e-9 , QS0=-5.72466e-3 ,
     &          QS1=+1.02270e-4 , QS2=-1.65460e-6 , Q20=+4.8314e-4 )
      real cccc,cccc_adj
      real cccc1,cccc1_adj
# else
      real cff
# endif
!
# include "compute_extended_bounds.h"
!
        do j=JstrR,JendR
# ifdef NONLIN_EOS
          do k=1,NDAS                        
            do i=IstrR,IendR                
              Tt_adj=t_s(i,j,k,itemp) 
                                                           !error standar deviation. This 
                                                           !is introduced because of preconditioning 
              Tt=t_das9(i,j,k,itemp)   
#  ifdef SALINITY
              Ts_adj=t_s(i,j,k,isalt)
              Ts=t_das9(i,j,k,isalt)
              if (Ts.lt.0.0000000001) then
                sqrtTs_adj=0.0
              else
                sqrtTs_adj=0.5/sqrt(Ts)*Ts_adj
              endif
              sqrtTs=sqrt(Ts)
#  else
              Ts_adj=0.0
              Ts=0.0
              sqrtTs_adj=0.0
              sqrtTs=0.0
#  endif
              K0_adj(i,k)=Tt_adj*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
     &                   +Tt*Tt_adj*(A02+2.0*Tt*A03
     &                       +3.0*Tt*Tt*A04)
     &                   +Ts_adj*(A10+Tt*(A11+Tt*(A12+Tt*A13))
     &                       +sqrtTs*(AS0+Tt*(AS1+Tt*AS2)))

     &                   +Ts*( Tt_adj*(A11+Tt*(2.0*A12
     &                       +3.0*Tt*A13))
     &                   +sqrtTs_adj*(AS0+Tt*(AS1+Tt*AS2))
     &                   +sqrtTs*Tt_adj*(AS1+2.0*Tt*AS2) )

              K0(i,k)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
     &                   +Ts*(A10+Tt*(A11+Tt*(A12+Tt*A13))
     &                       +sqrtTs*(AS0+Tt*(AS1+Tt*AS2)))
              K1_adj(i,k)=Tt_adj*(B01+Tt*(B02+Tt*B03))
     &                   +Tt*Tt_adj*(B02+2.0*Tt*B03)
     &                   +Ts_adj*(B10+Tt*(B11+Tt*B12)+sqrtTs*BS1)
     &                   +Ts*(Tt_adj*(B11+2.0*Tt*B12)
     &                           +sqrtTs_adj*BS1)
              K1(i,k)=B00+Tt*(B01+Tt*(B02+Tt*B03))
     &                   +Ts*(B10+Tt*(B11+Tt*B12)+sqrtTs*BS1)
              K2_adj(i,k)=Tt_adj*(E01+2.0*Tt*E02)
     &                   +Ts_adj*(E10+Tt*(E11+Tt*E12))
     &                   +Ts*Tt_adj*(E11+2.0*Tt*E12)
              K2(i,k)=E00+Tt*(E01+Tt*E02)
     &                   +Ts*(E10+Tt*(E11+Tt*E12))
              rho1_adj(i,k)=Tt_adj*(Q01+Tt*(Q02+Tt
     &                   *(Q03+Tt*(Q04+Tt*Q05))))
     &                  +Tt*(Q02+Tt*(2.0*Q03
     &                    +Tt*(3.0*Q04+4.0*Tt*Q05)))*Tt_adj
     &                  +Ts_adj*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))
     &                           +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20)
     &                  +Ts*( (Q11+Tt*(2.0*Q12
     &                    +Tt*(3.0*Q13+4.0*Tt*Q14)))*Tt_adj
     &                  +sqrtTs_adj*(QS0+Tt*(QS1+Tt*QS2))
     &                  +sqrtTs*Tt_adj*(QS1+2.0*Tt*QS2)+Ts_adj*Q20 )
              rho1(i,k)=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
     &                    +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))
     &                          +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20)
            enddo
                       /* R8000: 35 clock cycles, 53%(64%) of peak */
            do i=IstrR,IendR
              cccc_adj=-0.1*z_das(k)/
     &             ( K0(i,k)-z_das(k)*(K1(i,k)-z_das(k)*K2(i,k)) )/
     &             ( K0(i,k)-z_das(k)*(K1(i,k)-z_das(k)*K2(i,k)) )
     &           * ( K0_adj(i,k)-z_das(k)*(K1_adj(i,k)-z_das(k)
     &           * K2_adj(i,k)) )
              cccc=1.+0.1*z_das(k)/( K0(i,k)
     &                  -z_das(k)*(K1(i,k)-z_das(k)*K2(i,k)) )

              rho_s(i,j,k)=( rho1_adj(i,k)
     &                      -rho1(i,k)*cccc_adj/cccc
     &                       )/cccc

            enddo   
          enddo
# else
          do k=1,NDAS                                       ! LINEAR
            do i=IstrR,IendR                             ! EQUATION
              rho_s(i,j,k)=Tcoef*t_s(i,j,k,itemp)    ! OF STATE
#  ifdef SALINITY
     &                  +Scoef*t_s(i,j,k,isalt)
#  endif
            enddo
          enddo
# endif /* NONLIN_EOS */
        enddo
      return
      end

