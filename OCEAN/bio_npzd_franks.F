!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined BIOLOGY && defined BIO_NPZD_FRANKS

      subroutine biology_tile (Istr,Iend,Jstr,Jend)
!
! Compute biological forcing functions
!
! In this particular implementation there is 4 compartments:
! NO3, PHYTOplankton, ZOOplanknton, DETritus.
!
      implicit none
      integer Istr,Iend,Jstr,Jend
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "diagnostics.h"
# include "scalars.h"
# include "forces.h"
# include "mixing.h"

      real K_ext, K_NO3, K_phy, Vm_NO3, PhyMR, ZooGR, ZooMR,
     &     ZooMD, ZooGA, ZooEC, DetRR, wDet
      integer ITERMAX
      integer nsink
# ifdef DIAGNOSTICS_BIO

      real trend_no3,trend_phy,trend_zoo,trend_det,somme
      real bilan_no3,bilan_phy,bilan_zoo,bilan_det, sinking_loss,
     &     trend_total
      integer l, iflux
      real ThisVSinkFlux(N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(N, NumFluxTerms)
# endif
      parameter (
     &  ITERMAX = 3,      ! number of small implicit time step
     &  nsink   = NumVSinkTerms, ! Vertical sink terms

     &  K_ext = 0.06,     ! light extinction coefficient (1/m), default = 0.06
     &  K_NO3 = 0.1,      ! half-saturation concentration for phytoplankton 
                          ! nitrate uptake (mmol N/m3), default = 0.1
     &  K_phy = 0.4,      ! phytoplankton saturation coefficient (mmol N/m3), 
                          ! default = 0.4
     &  Vm_NO3 = 2.0,     ! nitrate uptake rate (1/d), default = 2.0
     &  PhyMR = 0.05,     ! phytoplankton senescence/mortality rate (1/d), 
                          ! default = 0.05
     &  ZooGR = 0.6,      ! zooplankton maximum growth rate (1/d), default = 0.6
     &  ZooMR = 0.2,      ! zooplankton mortality rate (1/d), default = 0.2
     &  ZooMD = 0.05,     ! zooplankton death bits rate (1/d), default = 0.05
     &  ZooGA = 0.3,      ! zooplankton grazing inefficiency (1/d), default = 0.3
     &  ZooEC = 0.15,     ! zooplankton excreted fraction (1/d), default = 0.15
     &  DetRR = 0.1,      ! detritus remineralization rate (1/d), default = 0.1
     &  wDet = 5.0  )     ! detrital sinking rate (m/d), default = 5.0

      integer i,j,k, ITER, iB
      real    NO3(N), Phyt(N), Zoo(N), Det(N),
     &        aJ(N),FC(0:N),
     &        PAR, PARsup, attn, Vp, Epp, cu, aL,aR, dtdays, dtsec, L_NO3,
     &        E_NO3,cff,cff1,cff2,cff3,cff6,
     &        SB(N,nsink),dSB(0:N,nsink),wSB(nsink)

# if defined DIAGNOSTICS_BIO
     &      , LastVSinkFlux,ColumnMassOld(NumVSinkTerms)
     &      , ColumnMassNew(NumVSinkTerms)
# endif /* DIAGNOSTICS_BIO */
!
# include "compute_auxiliary_bounds.h"
!
      dtdays=dt/(24.*3600.*float(ITERMAX))  ! time step as fraction of day.
# if defined DIAGNOSTICS_BIO
      dtsec = dt / float(ITERMAX)           ! time step in seconds
# endif /* DIAGNOSTICS_BIO */
!
!
! Since the following solver is iterative to achieve implicit
! discretization of the biological interaction, two time slices are
! required, BIO where BIO is understood as vector of
! biological state variables: BIO=[NO3,Phyt,Zoo,Det]. Assume
! that the iterations converge, the newly obtained state variables
! satisfy equations
!
!           BIO = BIO + dtdays * rhs(BIO)
!
! where rhs(BIO) is the vector of biological r.h.s. computed at
! the new time step. During the iterative procedure a series of
! fractional time steps is performed in a chained mode (splitting
! by different biological conversion processes) in sequence NO3 --
! Phyt -- Zoo -- Det, that is the main food chain. In all
! stages the concentration of the component being consumed is
! treated in fully implicit manner, so that the algorithm guarantees
! non-negative values, no matter how strong is the concentration of
! active consuming component (Phyto or Zoo).
!
! The overall algorithm, as well as any stage of it is formulated
! in conservative form (except explicit sinking) in sense that the
! sum of concentration of all four components is conserved.
!

/*
#  ifdef EW_PERIODIC
#   define I_RANGE Istr,Iend
#  else
#   define I_RANGE IstrR,IendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE Jstr,Jend
#  else
#   define J_RANGE JstrR,JendR
#  endif
*/
# define I_RANGE Istr,Iend
# define J_RANGE Jstr,Jend



      do j=J_RANGE
        do i=I_RANGE
# ifdef DIAGNOSTICS_BIO
! Reset the biogeochemical fluxes. This is necessary because the
! biological routine uses multiple. time steps for each physical time
! step.
          do k=1,N
            do l=1,NumFluxTerms
              bioFlux(i,j,k,l) = 0.0
            enddo
          end do
          do k=0,N
            do l=1,NumVSinkTerms
              bioVSink(i,j,k,l) = 0.0
            enddo
          enddo
# endif /* DIAGNOSTICS_BIO */
!
! Extract biological variables from tracer arrays; place them into
! scratch variables; restrict their values to be positive definite.
!

      do k=1,N
        NO3(k) =max(t(i,j,k,nnew,iNO3_)  ,0.)   ! Nitrate
        Phyt(k)=max(t(i,j,k,nnew,iPhy1)  ,0.)   ! Phytoplankton
        Zoo(k) =max(t(i,j,k,nnew,iZoo1)  ,0.)   ! Zooplankton
        Det(k) =max(t(i,j,k,nnew,iDet1)  ,0.)   ! Detritus
      enddo

      DO ITER=1,ITERMAX      !--> Start internal iterations to achieve
                             !    nonlinear backward-implicit solution.
        cff1=dtdays*Vm_NO3
        do k=N,1,-1     !<-- irreversible
! (1) NO3 uptake by Phyto
!
          cff=Phyt(k)*cff1*EXP(K_ext*z_r(i,j,k))/(K_NO3+NO3(k))
          NO3(k)=No3(k)/(1.+cff)
          Phyt(k)=Phyt(k)+NO3(k)*cff
# ifdef DIAGNOSTICS_BIO
          ThisFlux(k, NFlux_NewProd) = cff*NO3(k)
# endif /* DIAGNOSTICS_BIO */
        enddo
!
! (1) Phytoplankton grazing by Zooplankton to Zoo and Detr
! (2) Phytoplankton mortality to Detr (mu_P_D)
!
        cff1=dtdays*ZooGR
        cff2=dtdays*PhyMR
        cff3=K_phy*K_phy
        do k=1,N
          cff=Zoo(k)*Phyt(k)*cff1/(cff3+Phyt(k)*Phyt(k))
          Phyt(k)=Phyt(k)/(1.+cff+cff2)
          Zoo(k)=Zoo(k)+Phyt(k)*cff*(1.-ZooGA)
          Det(k)=Det(k)+Phyt(k)*(cff2+cff*(ZooGA-ZooEC))
          NO3(k)=NO3(k)+Phyt(k)*cff*ZooEC
          
# ifdef DIAGNOSTICS_BIO
          ThisFlux(k, NFlux_Grazing)=Phyt(k)*cff*(1.-ZooGA)
          
          ThisFlux(k, NFlux_SlopFeed) = Phyt(k)*cff*ZooGA
          ThisFlux(k, NFlux_Pmort) = Phyt(k) * cff2
# endif /* DIAGNOSTICS_BIO */
        enddo
!
! (1) Zoo excretion to NO3  (rate mu_Z_A)
! (2) Zoo mortality to Det (rate mu_Z_D)
!
        cff1=1./(1.+dtdays*(ZooMR+ZooMD))
        cff2=dtdays*ZooMR
        cff3=dtdays*ZooMD
        do k=1,N
          Zoo(k)=Zoo(k)*cff1
          NO3(k)=NO3(k)+Zoo(k)*cff2
          Det(k)=Det(k)+Zoo(k)*cff3

# ifdef DIAGNOSTICS_BIO
          ThisFlux(k, NFlux_Zmetab)=Zoo(k)*cff2
          ThisFlux(k, NFlux_Zmort)=Zoo(k)*cff3
# endif /* DIAGNOSTICS_BIO */
        enddo
!
! (1) Det remineralization to N03
!
        cff1=dtdays*DetRR
        cff2=1./(1.+cff1)
        do k=1,N
          Det(k)=Det(k)*cff2
          NO3(k)=NO3(k)+Det(k)*cff1

# ifdef DIAGNOSTICS_BIO
          ThisFlux(k, NFlux_ReminD)=Det(k)*cff1
# endif /* DIAGNOSTICS_BIO */
        enddo

!
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
!
        do k=1,N
          SB(k,1)=Det(k)
        enddo
        wSB(1)=wDet

        do iB=1,nsink

! Part (i): Construct parabolic splines: compute vertical derivatives
! of the fields SB. The derivatives are located at W-points;
! Neumann boundary conditions are assumed on top and bottom.
!
          dSB(0,iB)=0.
          FC(0)=0.
          cff6=6.
          do k=1,N-1
            cff=1./(2.*Hz(i,j,k+1)+Hz(i,j,k)*(2.-FC(k-1)))
            FC(k)=cff*Hz(i,j,k+1)
            dSB(k,iB)=cff*(cff6*(SB(k+1,iB)-SB(k,iB))
     &                             -Hz(i,j,k)*dSB(k-1,iB))
          enddo
          dSB(N,iB)=0.
          do k=N-1,1,-1     !<-- irreversible
            dSB(k,iB)=dSB(k,iB)-FC(k)*dSB(k+1,iB)
          enddo
!
! Part (ii): Convert dSB [which are now vertical derivatives
! of fields SB at the grid box interfaces] into field values
! at these interfaces, assuming parabolic profiles within each grid
! box. Restrict these values to lie between bounds determined from
! box-averaged values of grid boxes adjscent from above and below.
! (This restriction is part of PPM-like monotonization procedure.)
!
          cff=1./3.
          dSB(0,iB)=SB(1,iB) !-cff*Hz(1)*(dSB(0,iB)+0.5*dSB(1,iB))
          dSB(N,iB)=SB(N,iB) !+cff*Hz(N)*(dSB(N,iB)+0.5*dSB(N-1,iB))
          do k=2,N          !<-- irreversible
            dSB(k-1,iB)=SB(k,iB)
     &                 -cff*Hz(i,j,k)*(0.5*dSB(k,iB)+dSB(k-1,iB))
            dSB(k-1,iB)=max(dSB(k-1,iB),min(SB(k-1,iB),SB(k,iB)))
            dSB(k-1,iB)=min(dSB(k-1,iB),max(SB(k-1,iB),SB(k,iB)))
          enddo
!
! Part (iii): Convert dSB into flux-integrated values,
! complete PPM flux limiting. This procedure starts from assigning
! Left and Right (aR,aL) values of the interpolating parabolae, then
! monotonicity conditions are checked and aL,aR are modified to fit.
! Overall, from this moment and further on it follows Colella--
! --Woodward, 1984 bombmaking code almost exactly.
!
          do k=1,N           !<-- irreversible
            FC(k)=dtdays/Hz(i,j,k)
            aR=dSB(k,iB)
            aL=dSB(k-1,iB)
            cff1=(aR-aL)*6.*(SB(k,iB)-.5*(aR+aL))
            cff2=(aR-aL)**2
            if ((aR-SB(k,iB))*(SB(k,iB)-aL).lt.0.) then
              aL=SB(k,iB)
              aR=SB(k,iB)
            elseif (cff1.gt.cff2) then
              aL=3.*SB(k,iB)-2.*aR
            elseif (cff1.lt.-cff2) then
              aR=3.*SB(k,iB)-2.*aL
            endif
            cu=wSB(iB)*FC(k)
            dSB(k-1,iB)=SB(k,iB)-(1.-cu)*(.5*(aR-aL)-(.5*(aR+aL)
     &                                   -SB(k,iB) )*(1.-2.*cu))
          enddo
          dSB(N,iB)=0.   ! Set no-flux boundary conditions at top.
!
! Apply fluxes:
!
          do k=1,N
            SB(k,iB)=SB(k,iB)+wSB(iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
          enddo
        enddo  ! <-- iB

# ifdef DIAGNOSTICS_BIO
        do iflux = 1, NumVSinkTerms
          ColumnMassOld(iflux) = 0.0
          ColumnMassNew(iflux) = 0.0
        end do
# endif /* DIAGNOSTICS_BIO */

        do k=1,N
          ! detritus
# ifdef DIAGNOSTICS_BIO
          ColumnMassOld(1)=ColumnMassOld(1)
     &                    +Det(k)
          ThisVSinkFlux(k, NFlux_VSinkD1)=Det(k)-SB(k,1)
# endif /* DIAGNOSTICS_BIO */
          Det(k) = SB(k,1)
# ifdef DIAGNOSTICS_BIO
          ColumnMassNew(1)=ColumnMassNew(1)
     &                    +Det(k)
# endif /* DIAGNOSTICS_BIO */
        enddo
!
# ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
! for computational efficiency, divide now by dtsec to get the correct units
        do iflux = 1, NumFluxTerms
          do k = 1, N
            bioFlux(i,j,k,iflux) = ( bioFlux(i,j,k,iflux) +
     &                               ThisFlux(k, iflux) / dt )
#  ifdef MASKING
     &             * rmask(i,j)
#  endif /* MASKING */
          end do
        end do
        do iflux = 1, NumVSinkTerms
! Compute the vertical sinking flux into the sediment by comparing
! previous and current mass in this (i,j) column
! The flux is positive if upward, so usually it will be
! negative, i.e. into the sediment.
          LastVSinkFlux = ( ColumnMassNew(iflux) -
     &                      ColumnMassOld(iflux) )
          bioVSink(i,j,0,iflux) = (  bioVSink(i,j,0,iflux) +
     &                               LastVSinkFlux / dt )
#  ifdef MASKING
     &          * rmask(i,j)
#  endif /* MASKING */
          do k = 1, N
            LastVSinkFlux = LastVSinkFlux +
     &                      ThisVSinkFlux(k,iflux)
            bioVSink(i,j,k,iflux) =(  bioVSink(i,j,k,iflux) +
     &                                LastVSinkFlux / dt )
#  ifdef MASKING
     &          * rmask(i,j)
#  endif /* MASKING */
          end do
        end do
# endif /* DIAGNOSTICS_BIO */
!
      ENDDO  ! <-- ITER
!
! Write back
!
!          print*,'N=',N
# undef DEBUG_BIO
# if defined DIAGNOSTICS_BIO && defined DEBUG_BIO
!
      k=N
      if ((i.eq.13).and.(j.eq.15)) then
        bilan_no3 = bioFlux(i,j,k,NFlux_Zmetab)
     &           + bioFlux(i,j,k,NFlux_ReminD)
     &           - bioFlux(i,j,k,NFlux_NewProd)
!
        bilan_phy = bioFlux(i,j,k,NFlux_NewProd)
     &           - bioFlux(i,j,k,NFlux_Pmort)
     &           - bioFlux(i,j,k,NFlux_Grazing)
     &           - bioFlux(i,j,k,NFlux_SlopFeed)
     &           - ( bioVSink(i,j,k,NFlux_VSinkP1)
     &               - bioVSink(i,j,k-1,NFlux_VSinkP1) )
!
        bilan_zoo = bioFlux(i,j,k,NFlux_Grazing)
     &            - bioFlux(i,j,k,NFlux_Zmetab)
     &            - bioFlux(i,j,k,NFlux_Zmort)
!
        bilan_det = bioFlux(i,j,k,NFlux_SlopFeed)
     &            + bioFlux(i,j,k,NFlux_Zmort)
     &            + bioFlux(i,j,k,NFlux_Pmort)
     &            - bioFlux(i,j,k,NFlux_ReminD)
     &            - ( bioVSink(i,j,k,NFlux_VSinkD1)
     &                - bioVSink(i,j,k-1,NFlux_VSinkD1) )
!
        somme = bilan_no3+bilan_phy+bilan_zoo+bilan_det
        trend_no3 = ( (min(t(i,j,k,nnew,iNO3_),0.) +NO3(k))
     &               - t(i,j,k,nnew,iNO3_) ) / dt
        trend_phy = ( (min(t(i,j,k,nnew,iPhy1),0.) +Phyt(k))
     &               - t(i,j,k,nnew,iPhy1) )/ dt
        trend_zoo = ( (min(t(i,j,k,nnew,iZoo1),0.) +Zoo(k))
     &               - t(i,j,k,nnew,iZoo1) )/ dt
        trend_det = ( (min(t(i,j,k,nnew,iDet1),0.) +Det(k))
     &       - t(i,j,k,nnew,iDet1) )/ dt
!
        trend_total =  trend_no3 + trend_phy + trend_zoo + trend_det
!
        sinking_loss = - ( bioVSink(i,j,k,NFlux_VSinkD1)
     &                     - bioVSink(i,j,k-1,NFlux_VSinkD1) )
     &                 - ( bioVSink(i,j,k,NFlux_VSinkP1)
     &                     - bioVSink(i,j,k-1,NFlux_VSinkP1) )
!
        print*, '=================='
        print*, 'i=',i,' j=',j,' k=',k
        print*, 'somme SMS(of N)=         ',somme
        print*, 'trend_total of N=        ',trend_total
        print*, 'Sinking_loss of N=       ',sinking_loss
        print*, 'These three values have to be the same!'
        print*, 'Error for N = trend_total-somme=',trend_total-somme
        print*, "-----------------"
        print*, 'bilan_no3 - trend_no3=', bilan_no3-trend_no3
        print*, 'bilan_phy-trend_phy=  ', bilan_phy-trend_phy
        print*, 'bilan_zoo-trend_zoo=  ', bilan_zoo-trend_zoo
        print*, 'bilan_det-trend_det=  ', bilan_det-trend_det
        print*, '-----------------'
        print*, 'bioFlux(i,j,k,NFlux_NewProd)= ',
     &                               bioFlux(i,j,k,NFlux_NewProd)
        print*, 'bioFlux(i,j,k,NFlux_Grazing)= ',
     &                               bioFlux(i,j,k,NFlux_Grazing)
        print*, 'bioFlux(i,j,k,NFlux_SlopFeed)=',
     &                               bioFlux(i,j,k,NFlux_SlopFeed)
        print*, 'bioFlux(i,j,k,NFlux_Pmort)=   ',
     &                               bioFlux(i,j,k,NFlux_Pmort)
        print*, 'bioFlux(i,j,k,NFlux_Zmetab)=  ',
     &                               bioFlux(i,j,k,NFlux_Zmetab)
        print*, 'bioFlux(i,j,k,NFlux_Zmort)=   ',
     &                               bioFlux(i,j,k,NFlux_Zmort)
        print*, 'bioFlux(i,j,k,NFlux_ReminD)=  ',
     &                               bioFlux(i,j,k,NFlux_ReminD)
        print*, 'bioVSink(i,j,k,NFlux_VSinkP1)=',
     &                                bioVSink(i,j,k,NFlux_VSinkP1)
        print*, 'bioVSink(i,j,k,NFlux_VSinkD1)=',
     &                                bioVSink(i,j,k,NFlux_VSinkD1)
        print*, '=================='
      endif
# endif

          do k=1,N
            t(i,j,k,nnew,iNO3_)=min(t(i,j,k,nnew,iNO3_),0.) +NO3(k)
            t(i,j,k,nnew,iPhy1)=min(t(i,j,k,nnew,iPhy1),0.) +Phyt(k)
            t(i,j,k,nnew,iZoo1)=min(t(i,j,k,nnew,iZoo1),0.) +Zoo(k)
            t(i,j,k,nnew,iDet1)=min(t(i,j,k,nnew,iDet1),0.) +Det(k)
          enddo
        enddo
      enddo


#else
      subroutine biology_empty ()
#endif
      return
      end
