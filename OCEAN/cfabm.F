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
#if defined BIOLOGY && defined CFABM

      subroutine biology_tile (Istr,Iend,Jstr,Jend,trd)
!
! Compute biological forcing functions
!
! In this particular implementation there is 4 compartments:
! NO3, PHYTOplankton, ZOOplanknton, DETritus.
!
        use fabm
        use fabm_config

        implicit none
        integer Istr,Iend,Jstr,Jend,trd
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "diagnostics.h"
# include "scalars.h"
# include "forces.h"
# include "mixing.h"

# ifdef DIAGNOSTICS_BIO
      integer l, iflux
      real ThisVSinkFlux(Istr:Iend, 0:N, NumVSinkTerms),  ! [mmol m-2 s-1], upward flux is positive
     &     ThisFlux(Istr:Iend, N, NumFluxTerms)
     &     ThisGasExcFlux(Istr:Iend, NumGasExcTerms)
# endif

      logical, parameter :: repair = .TRUE.
      logical :: valid_int, valid_sf, valid_bt
      integer i, j, k, iB
      integer NT, itrc     !!! FABM related edit
      integer nx, ny    !!! FABM related edit
      integer :: i1, j1 !These indices range strictly over (1,nx), (1,ny), while (i,j) range over (Istr,Iend), (Jstr,Jend)

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

!!! FABM related edit start
      nx = Iend-Istr+1
      ny = Jend-Jstr+1
      NT = ntrc_bio
      real    state1(Istr:Iend,Jstr:Jend,N,NBT), 
     &        state_sf(Istr:Iend,Jstr:Jend,3,NSAT), 
     &        state_bt(Istr:Iend,Jstr:Jend,3,NBAT), 
     &        state_old(Istr:Iend,N,NBT),                   !Old state (accounting for clipping)
     &        dBdt(Istr:Iend,N,NBT),                        !FABM SMS term (units [C/s])
     &        sms_sf(Istr:Iend,NBT),                        !FABM surface-attached SMS terms (units [Matter/m2/s])
     &        sms_bt(Istr:Iend,NBT),                        !FABM bottom-attached SMS terms (units [Matter/m2/s])
     &        flux_sf(Istr:Iend,NBT),                       !FABM surface fluxes (units [Matter/m2/s])
     &        flux_bt(Istr:Iend,NBT),                       !FABM bottom fluxes (units [Matter/m2/s])
     &        FC(Istr:Iend,0:N),
     &        Hz_inv(Istr:Iend,N),
     &        Hz_inv2(Istr:Iend,N),
     &        Hz_inv3(Istr:Iend,N),
     &        SB(Istr:Iend,N,NT),
     &        wSB(Istr:Iend,N,NT)

#ifdef FABM_NONNEG_S
      real    Snn(Istr:Iend,Jstr:Jend,N)                          !Non-negative salinity for input to FABM
#endif

!!! FABM related edit end
      real    cu, aL,aR
     &        E_NO3,cff,cff1,cff2,cff6

# if defined DIAGNOSTICS_BIO
     &      , LastVSinkFlux,ColumnMassOld(Istr:Iend,NumVSinkTerms)
     &      , ColumnMassNew(Istr:Iend,NumVSinkTerms)
# endif /* DIAGNOSTICS_BIO */
!
# include "compute_auxiliary_bounds.h"
!-----------------------------------------------------------------------
! Relink environmental data as required (due to variable nnew).
!-----------------------------------------------------------------------
!Link environmental data to FABM (for variables that must be relinked for each time step since nnew varies)
!Note: We cannot use variable_needs_values to parse these calls, since these will return
!      .false. after the first linking calls in roms_fabm.F.
      IF (FMODEL%needsT(trd)) THEN
        CALL FMODEL%f(trd)%model%link_interior_data(               &
      &                    FMODEL%id_temp(trd),                     &
      &                    t(Istr:Iend,Jstr:Jend,1:N,nnew,1))
      END IF
      IF (FMODEL%needsS(trd)) THEN
      #ifdef FABM_NONNEG_S
        DO k=1,N
          DO j=J_RANGE
            DO i=I_RANGE
              Snn(i,j,k) = MAX(0., t(i,j,k,nnew,2))
            END DO
          END DO
        END DO
        CALL FMODEL%f(trd)%model%link_interior_data(               &
      &                    FMODEL%id_salt(trd),                     &
      &                    Snn(Istr:Iend,Jstr:Jend,1:N))
      #else      
        CALL FMODEL%f(trd)%model%link_interior_data(               &
      &                    FMODEL%id_salt(trd),                     &
      &                    t(Istr:Iend,Jstr:Jend,1:N,nnew,2))
      #endif
      END IF
!-----------------------------------------------------------------------
! Relink surface/bottom state data as required (due to variable nnew).
!-----------------------------------------------------------------------
      DO itrc=1,NSAT
        CALL FMODEL%f(trd)%model%link_surface_state_data(          &
     &           itrc,state_sf(Istr:Iend,Jstr:Jend,nnew,itrc))
      END DO

      DO itrc=1,NBAT
        CALL FMODEL%f(trd)%model%link_bottom_state_data(           &
     &           itrc,state_bt(Istr:Iend,Jstr:Jend,nnew,itrc))
      END DO
!-----------------------------------------------------------------------------------
! Copy interior data from "t" to "state", capping if required using FABM check_state
!-----------------------------------------------------------------------------------
!  
!  Note: The combined FABM sms + ROMS sinking increment will be calculated as (state - state_old).
!        state_old is assigned within the J loop to save memory.
!  Note: If capping is applied here using FABM (cpp FABM_CHECK_STATE) then this capping
!        will be applied to BOTH "state_old" and "state".  This is consistent with the
!        capping applied to Bio and Bio_old in larger ROMS models (e.g. fennel.h, ecosim.h). 
!        However, in smaller ROMS biology modules e.g. npzd_Franks.h, nonnegativity is imposed 
!        more loosely as a correction that conserves total mass, and is not applied to Bio_old.
!        To reproduce the latter behaviour, FABM_CHECK_STATE should be deactivated (or parameter
!        "repair" set to "false") and the mass-conserving correction is done within the FABM module 
!        (e.g. niva_roms_npzd_Franks.F90).
      DO itrc=1,NBT
        ibio=itrc_bio + itrc - 1    ! because itrc starts from 1, we need to subtract 1.
        DO k=1,N
          DO j=J_RANGE
            DO i=I_RANGE
              state1(i,j,k,itrc)=t(i,j,k,nnew,ibio) !Note: FABM is permanently linked to state1 (in roms_fabm.F)
            END DO
          END DO
        END DO
      END DO
#ifdef FABM_CHECK_STATE
      DO j1=1,ny
        DO k=1,N
          CALL FMODEL%f(trd)%model%check_interior_state(1,nx,j1, &
     &        k,repair,valid_int)   !This will cap state1(Istr:Iend,Jstr+j1,k,1:NBT)
        END DO
      END DO
      IF (NSAT.gt.0) THEN
        DO j1=1,ny
          CALL FMODEL%f(trd)%model%check_surface_state(1,nx,j1,  &
     &        repair,valid_sf)      !This will cap state_sf(Istr:Iend,Jstr+j1,nnew,1:NSAT)
        END DO
      END IF
      IF (NBAT.gt.0) THEN
        DO j1=1,ny
          CALL FMODEL%f(trd)%model%check_bottom_state(1,nx,j1,   &
     &        repair,valid_bt)      !This will cap state_bt(Istr:Iend,Jstr+j1,nnew,1:NBAT)
        END DO
      END IF
#endif
!-----------------------------------------------------------------------
!  Prepare all fields FABM needs to compute SMS terms (e.g., light)
!-----------------------------------------------------------------------
      CALL FMODEL%f(trd)%model%prepare_inputs()
!
!
!      dtdays=dt/(24.*3600.*float(ITERMAX))  ! time step as fraction of day.
# if defined DIAGNOSTICS_BIO
!      dtsec = dtFABM           ! time step in seconds
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
      do j=J_RANGE
        j1 = j-Jstr+1 !FABM chunks are always indexed 1:chunksize
                      !Therefore, all calls to FABM APIs should use index j1 instead of j
# ifdef DIAGNOSTICS_BIO
! Reset the biogeochemical fluxes. This is necessary because the
! biological routine uses multiple time steps for each physical time step.
        do i=I_RANGE
          do k=1,N
            do l=1,NumFluxTerms
              ThisFlux(i,k,l) = 0.0
              bioFlux(i,j,k,l) = 0.0
            enddo
          enddo
          do k=0,N
            do l=1,NumVSinkTerms
              bioVSink(i,j,k,l) = 0.0
              ThisVSinkFlux(i,k,l) = 0.0
            enddo
          enddo
          do l=1,NumGasExcTerms
            ThisGasExcFlux(i,l) = 0.0
            GasExcFlux(i,j,l) = 0.0
          enddo
        enddo
# endif /* DIAGNOSTICS_BIO */
!-----------------------------------------------------------------------
!  Compute inverse thickness to avoid repeated divisions.
!-----------------------------------------------------------------------
        DO k=1,N
          DO i=I_RANGE
            Hz_inv(i,k)=1./Hz(i,j,k)
          END DO
        END DO
!-----------------------------------------------------------------------
!  Store the old interior state.
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO k=1,N
            DO i=I_RANGE
              state_old(i,k,itrc)=state1(i,j,k,itrc)
            END DO
          END DO
        END DO
!-----------------------------------------------------------------------
!  Generate SMS terms for pelagic variables
!-----------------------------------------------------------------------
        DO k=1,N
          dBdt(Istr:Iend,k,1:NBT)=0.
          CALL FMODEL%f(trd)%model%get_interior_sources(1,nx,j1,k, &
     &                                        dBdt(Istr:Iend,k,1:NBT))
!Note: the biological tracer indices idbio range over (NAT+NPT+NCS+NNS)+1:NT = NT-NBT+1:NT, see rfabm_mod.h
        END DO
!-----------------------------------------------------------------------
!  Generate surface fluxes and SMS for surface-attached variables
!-----------------------------------------------------------------------
        flux_sf = 0.
        sms_sf = 0.
        !Note: surface fluxes may still be needed even if NSAT = 0
        CALL FMODEL%f(trd)%model%get_surface_sources(1,nx,j1,      &
     &         flux_sf(Istr:Iend,1:NBT),sms_sf(Istr:Iend,1:NSAT))
!-----------------------------------------------------------------------
!  Generate bottom fluxes and SMS terms for bottom-attached variables
!-----------------------------------------------------------------------
        flux_bt = 0.
        sms_bt = 0.
        !Note: bottom fluxes may still be needed even if NBAT = 0
        CALL FMODEL%f(trd)%model%get_bottom_sources(1,nx,j1,       &
     &        flux_bt(Istr:Iend,1:NBT),sms_bt(Istr:Iend,1:NBAT))
!-----------------------------------------------------------------------
!  Add contributions from surface and bottom fluxes
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO i=I_RANGE
            dBdt(i,N,itrc) = dBdt(i,N,itrc) +                   &
      &              flux_sf(i,itrc)*Hz_inv(i,N)
            dBdt(i,1,itrc) = dBdt(i,1,itrc) +                           &
      &              flux_bt(i,itrc)*Hz_inv(i,1)
          !E.g surface/bottom attached variables may have [mass/m2] while pelagic variables may have [mass/m3]
          END DO
        END DO
!-----------------------------------------------------------------------
!  Update tracers with rates of change from FABM (Euler step)
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO k=1,N
            DO i=I_RANGE
              state1(i,j,k,itrc) = state1(i,j,k,itrc) +                 &
     &                      dBdt(i,k,itrc)*dtFABM
# ifdef DIAGNOSTICS_BIO
              ThisFlux(i,k,itrc) = dBdt(i,k,itrc)*dtFABM
# endif /* DIAGNOSTICS_BIO */
            END DO
          END DO
        END DO
!-----------------------------------------------------------------------
!  Update surface states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NSAT
          DO i=I_RANGE
            state_sf(i,j,nnew,itrc) = MAX(state_sf(i,j,nnew,itrc) +     &
     &                         sms_sf(i,itrc)*dtFABM, 0.)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
# ifdef DIAGNOSTICS_BIO
            ThisGasExcFlux(itrc) = sms_sf(i,itrc)*dtFABM
# endif /* DIAGNOSTICS_BIO */
          END DO
        END DO
!-----------------------------------------------------------------------
!  Update bottom states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NBAT
          DO i=I_RANGE
            state_bt(i,j,nnew,itrc) = MAX(state_bt(i,j,nnew,itrc) +     &
     &                         sms_bt(i,itrc)*dtFABM, 0.)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
# ifdef DIAGNOSTICS_BIO
            ThisFlux(i,0,itrc) = sms_bt(i,itrc)*dtFABM
# endif /* DIAGNOSTICS_BIO */
          END DO
        END DO
!-----------------------------------------------------------------------
!  Get vertical sinking velocities from FABM
!-----------------------------------------------------------------------
        DO k=1,N
          CALL FMODEL%f(trd)%model%get_vertical_movement(          &
     &           1,nx,j1,k,wSB(Istr:Iend,k,NT-NBT+1:NT))
        END DO
!        wSB(Istr:Iend,1:N,NT-NBT+1:NT) = -1*                  &
!     &        wSB(Istr:Iend,1:N,NT-NBT+1:NT)
        !FABM outputs vertical velocities in [m/s] positive upward
        !commented out ROMS code expected wSB in [m/s] positive downward
        !CROCO code below expects wSB in [m/s] positive upward as FABM
!
!-----------------------------------------------------------------------
! Vertical sinking: Vertical advection algorithm based on monotonic,
! continuous conservative parabolic splines.
!-----------------------------------------------------------------------
        do k=1,N
          do i=I_RANGE
            do itrc=1,NBT 
              iB=NT-NBT+itrc
              SB(i,k,iB)=state1(i,j,k,itrc)
            enddo
          enddo
        enddo
!
!!! FABM edit start
!        do iB=1,nsink       !!!FABM: This is not convenient for arbitrary fabm model
                            !!!FABM: Instead we loop over all biol tracers and use IF statement
        do itrc=1,NBT
          iB= NT-NBT+itrc
          IF (MAXVAL(ABS(wSB(Istr:Iend,1:N,iB))).gt.1.E(-12)) THEN
            do i=I_RANGE
!!!Note: 1.0E-12 m/s => 0.3 mm/decade (safely negligible)
!!! FABM edit end

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
                dSB(k,iB)=cff*(cff6*(SB(i,k+1,iB)-SB(i,k,iB))
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
              dSB(0,iB)=SB(i,1,iB) !-cff*Hz(1)*(dSB(0,iB)+0.5*dSB(1,iB))
              dSB(N,iB)=SB(i,N,iB) !+cff*Hz(N)*(dSB(N,iB)+0.5*dSB(N-1,iB))
              do k=2,N          !<-- irreversible
                dSB(k-1,iB)=SB(i,k,iB)
     &                 -cff*Hz(i,j,k)*(0.5*dSB(k,iB)+dSB(k-1,iB))
                dSB(k-1,iB)=max(dSB(k-1,iB),min(SB(i,k-1,iB),SB(i,k,iB)))
                dSB(k-1,iB)=min(dSB(k-1,iB),max(SB(i,k-1,iB),SB(i,k,iB)))
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
                FC(k)=dtFABM/Hz(i,j,k)
                aR=dSB(k,iB)
                aL=dSB(k-1,iB)
                cff1=(aR-aL)*6.*(SB(i,k,iB)-.5*(aR+aL))
                cff2=(aR-aL)**2
                if ((aR-SB(i,k,iB))*(SB(i,k,iB)-aL).lt.0.) then
                  aL=SB(i,k,iB)
                  aR=SB(i,k,iB)
                elseif (cff1.gt.cff2) then
                  aL=3.*SB(i,k,iB)-2.*aR
                elseif (cff1.lt.-cff2) then
                  aR=3.*SB(i,k,iB)-2.*aL
                endif
                cu=wSB(i,k,iB)*FC(k)
                dSB(k-1,iB)=SB(i,k,iB)-(1.-cu)*(.5*(aR-aL)-(.5*(aR+aL)
     &                                   -SB(i,k,iB) )*(1.-2.*cu))
              enddo
              dSB(N,iB)=0.   ! Set no-flux boundary conditions at top.
!
! Apply fluxes:
!
              do k=1,N
                SB(i,k,iB)=SB(i,k,iB)+wSB(i,k,iB)*FC(k)*(dSB(k,iB)-dSB(k-1,iB))
              enddo
            enddo  !  <-- I_RANGE
        enddo  ! <-- iB
!
        do i=I_RANGE
# ifdef DIAGNOSTICS_BIO
          do iflux = 1, NBT
            iB= NT-NBT+iflux
            ColumnMassOld(i,iflux) = 0.0
            ColumnMassNew(i,iflux) = 0.0
# endif /* DIAGNOSTICS_BIO */
          ! pelagic biotracer i
            do k=1,N
# ifdef DIAGNOSTICS_BIO
              ColumnMassOld(i,iflux)=ColumnMassOld(i,iflux)
     &                    +state1(i,j,k,iflux)
              ThisVSinkFlux(i, k, iflux)=state1(i,j,k,iflux)-SB(i,k,iB)
# endif /* DIAGNOSTICS_BIO */
              state1(i,j,k,iflux)=SB(i,k,iB)
# ifdef DIAGNOSTICS_BIO
              ColumnMassNew(i,iflux)=ColumnMassNew(i,iflux)
     &                    +state1(i,j,k,iflux)
# endif /* DIAGNOSTICS_BIO */
!
# ifdef DIAGNOSTICS_BIO
! Transfer fluxes to global arrays at the end of each biological time step
! for computational efficiency, divide now by dtsec to get the correct units
            ! TODO
#  ifdef MASKING
            ! TODO
#  endif /* MASKING */
            end do
! Compute the vertical sinking flux into the sediment by comparing
! previous and current mass in this (i,j) column
! The flux is positive if upward, so usually it will be
! negative, i.e. into the sediment.
            LastVSinkFlux = ( ColumnMassNew(i,iflux) -
     &                      ColumnMassOld(i,iflux) )
            bioVSink(i,j,0,iflux) = (  bioVSink(i,j,0,iflux) +
     &                               LastVSinkFlux / dtFABM )
#  ifdef MASKING
     &          * rmask(i,j)
#  endif /* MASKING */
            do k = 1, N
              LastVSinkFlux = LastVSinkFlux +
     &                      ThisVSinkFlux(i,k,iflux)
              bioVSink(i,j,k,iflux) =(  bioVSink(i,j,k,iflux) +
     &                                LastVSinkFlux / dtFABM )
#  ifdef MASKING
     &          * rmask(i,j)
#  endif /* MASKING */
            end do
# endif /* DIAGNOSTICS_BIO */
          end do    ! <== NBT pelagic biotracers
        enddo     ! <== I_RANGE
!
! Write back
!
        do itrc=1,NBT
          ibio=itrc_bio + itrc - 1    ! because itrc starts from 1, we need to subtract 1.
          do i=I_RANGE
            do k=1,N
              t(i,j,k,nnew,ibio)=min(t(i,j,k,nnew,ibio),0.) +state1(i,j,k,itrc)
            enddo
        enddo
      enddo
#else
      subroutine biology_empty ()
#endif
      return
      end
