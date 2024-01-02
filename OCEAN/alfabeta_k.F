! alfabeta_k.F -- compute 2D fields alpha, beta (lin EOS coeff.) at 1 sigma level
! adapted from alfabeta.F ; NJAL 2017-07-12
! only used if defined DIAGNOSTICS_DIAPYCNAL
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined DIAGNOSTICS_DIAPYCNAL
      subroutine alfabeta_k_tile (Istr,Iend,Jstr,Jend,k, alpha,beta)
!
!--------------------------------------------------------------------
!  This routine computes thermal expansion and saline contraction
!  coefficients as a function of potential temperature, salinity,
!  and pressure from a polynomial expression (Jackett & McDougall,
!  1992). The polynomial expression was found from fitting to 248
!  values  in the  oceanographic  ranges of  salinity,  potential
!  temperature,  and pressure.  It  assumes no pressure variation
!  along geopotential surfaces,  that is,  depth and pressure are
!  interchangeable. The coefficients are evaluated at the surface.
!
!  On Output:
!
!     alpha    Thermal expansion coefficient [kg/m^3/Celsius].
!     beta     Saline contraction coefficient [kg/m^3/PSU].
!
!  Adapted from original "rati" and "beta" routines.
!
!  Copyright (c) 1996 Rutgers University
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
      integer Istr,Iend,Jstr,Jend, i,j, k
      integer imin,imax,jmin,jmax
      real    alpha(PRIVATE_2D_SCRATCH_ARRAY),
     &        beta(PRIVATE_2D_SCRATCH_ARRAY)
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# ifdef NONLIN_EOS
#  define DUKO_2001
      real      Q00, Q01, Q02, Q03, Q04, Q05, U00, U01, U02, U03,
     &          U04, V00, V01, V02, W00
      parameter(Q00=+999.842594 , Q01=+6.793952E-2, Q02=-9.095290E-3,
     &          Q03=+1.001685E-4, Q04=-1.120083E-6, Q05=+6.536332E-9,
     &          U00=+0.824493   , U01=-4.08990E-3 , U02=+7.64380E-5 ,
     &          U03=-8.24670E-7 , U04=+5.38750E-9 , V00=-5.72466E-3 ,
     &          V01=+1.02270E-4 , V02=-1.65460E-6 , W00=+4.8314E-4  )
      real      cfq, cff_rho1, 
     &          drho1dT, drho1dS, dK0dT, dK0dS, dqp1dT, dqp1dS

      real K0,K1,K2,
     &     r00,r01,r02,r03,r04,r05,    K00,K01,K02,K03,K04,  dr00,
     &     r10,r11,r12,r13,r14,        K10,K11,K12,K13,      Ts, Tt,
     &     rS0,rS1,rS2,                KS0,KS1,KS2,          sqrtTs,
     &     r20,                                              dpth,
     &             B00,B01,B02,B03,    E00,E01,E02,
     &             B10,B11,B12,        E10,E11,E12,          cff,
     &                 BS1,                            cff1,cff2

      parameter(r00=999.842594,   r01=6.793952E-2,  r02=-9.095290E-3,
     &                            r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,

     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2,


     &          B00=0.4721788,    B01=0.01028859,   B02=-2.512549e-4,
     &                                              B03=-5.939910e-7,
     &          B10=-0.01571896,  B11=-2.598241e-4, B12=7.267926e-6,
     &                            BS1=2.042967e-3,

     &          E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,
     &          E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
#  ifdef DUKO_2001
      real K0_Duk
#  endif
#  endif
!
# include "compute_auxiliary_bounds.h"
!
!---------------------------------------------------------------
!  Compute thermal expansion and saline contraction coefficients
!  at surface
!
!  Ts salinity [PSU].
!  Tt potential temperature [deg Celsius].
!  den(Ts,Tt,0) surface density [kg/m^3]
!  rho1(Ts,Tt,0)=den(Ts,Tt,0)-1000. , computed from Jackett &
!                                               McDougall, 1992)
!  alpha(Ts,Tt,0)=-d(rho1(Ts,Tt,0))/d(Tt) / den(Ts,Tt,0)
!  beta(Ts,Tt,0) = d(rho1(Ts,Tt,0))/d(Ts) / den(Ts,Tt,0)
!---------------------------------------------------------------
!
# ifdef EW_PERIODIC
#  define I_EXT_RANGE Istr-1,Iend+1
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

#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)
      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif

# ifdef NONLIN_EOS
!
!  compute rho as a perturbation to rho0 (at the surface)
!
      dr00=r00-rho0
!
      cff=1./rho0
      cfq=0.1D0/(K00+K0_Duk)
#  endif

      do j=JstrR,JendR                                 ! NONLINEAR
          do i=IstrR,IendR                             !  EQUATION
# ifdef NONLIN_EOS
#  ifdef TEMPERATURE
            Tt=t(i,j,k,nstp,itemp)                     !  OF STATE
#  else
            Tt=25
#  endif
#  ifdef SALINITY
            Ts=max(t(i,j,k,nstp,isalt), 0.)
            sqrtTs=sqrt(Ts)
#  else
            Ts=33.64
            sqrtTs=5.8
#  endif
            cff_rho1 =( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                            r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                    rS1+Tt*rS2 ))+Ts*r20 ))
#  ifdef MASKING
     &                                                *rmask(i,j)
#  endif
            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))



            drho1dT = r01+Tt*( 2.*r02+Tt*( 3.*r03+Tt*(
     &                                  4.*r04 +Tt*5.*r05 )))
     &                         +Ts*( r11+Tt*( 2.*r12+Tt*(
     &                                   3.*r13 +Tt*4.*r14 ))
     &                              +sqrtTs*(rS1+Tt*2.*rS2))

            drho1dS =  r10+Tt*(r11+Tt*(r12+Tt*(r13+Tt*r14)))
     &                         +1.5*(rS0+Tt*(rS1+Tt*rS2))*sqrtTs
     &                                              +2.*r20*Ts 
            
            dK0dT = K01+Tt*( 2.*K02+Tt*( 3.*K03+ 4.*Tt*K04))
     &                         +Ts*((K11+Tt*( 2.*K12+ 3.*Tt*K13))
     &                              +sqrtTs*(KS1+Tt*2.*KS2))

            dK0dS =  K10+Tt*(K11+Tt*(K12+Tt*K13))
     &                         +1.5*(KS0+Tt*(KS1+Tt*KS2))*sqrtTs
                                                        
            dqp1dT = cfq*(drho1dT*(K0_Duk*K00+K0*(
     &                                     K0_Duk-K00-K0))
     &                    -dK0dT*((rho0+cff_rho1)*(K00+K0_Duk)))
     &                  /((K00+K0)*(K00+K0))    
     
            dqp1dS = cfq*(drho1dS*(K0_Duk*K00+K0*(
     &                                     K0_Duk-K00-K0))
     &                    -dK0dS*((rho0+cff_rho1)*(K00+K0_Duk)))
     &                  /((K00+K0)*(K00+K0))

            dpth=z_w(i,j,N)-z_r(i,j,k)
            alpha(i,j) = -cff*(drho1dT+dqp1dT*dpth*(1.-qp2*dpth))
            beta(i,j) = cff*(drho1dS+dqp1dS*dpth*(1.-qp2*dpth))


# else
!
!   Linear Equation of state thermal expansion and saline
!   contraction coefficients:
!
#  ifdef TEMPERATURE
          alpha(i,j)=abs(Tcoef)/rho0
#  else
          alpha(i,j)=0.
#  endif
#  ifdef SALINITY
          beta(i,j)=abs(Scoef)/rho0
#  else
          beta(i,j)=0.
#  endif
# endif /* NONLIN_EOS */
        enddo
      enddo

c*      do i=Istr,Iend,4
c*        write(6,15) i,alpha(i,10)*1.e7,beta(i,10)*1.e7
c*      enddo
c*  15  format(1x,'i = ',i4,' alpha = ',f6.1,1x,'  beta = ',f6.1)

#else
      subroutine alfabeta_k_empty
#endif /* DIAGNOSTICS_DIAPYCNAL */
      return
      end
