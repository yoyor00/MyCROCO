/*
!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================

! This routine couple fabm sources and sinks to the global arrays

*/
#include "cppdefs.h"
#if defined CFABM

      subroutine biology_tile (Istr, Iend, Jstr, Jend)
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

! -- hes 2024
! For FABM coupling
# include "fabm"
# include "fabm_config"
! --

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

# define I_RANGE Istr,Iend
# define J_RANGE Jstr,Jend

! local variables
      real(dp) :: yday ! Modified for updated ARANGO code (30/07/2019)
      real(r8) :: cff, cffL, cffR, cu, dltL, dltR

      ! == SINTEF 2024
      ! These initially used wider arrays (Istr-4:Iend+3), but I could not find a CROCO equivalent, so I'm trying this
      ! ==
      real(r8), dimension(Istr:Iend, N(ng), NT(ng)) :: w      ! ROMS vertical sinking velocity [m/s] positive downward
      integer, dimension(Istr:Iend, N(ng)) :: ksource
      real(r8), dimension(Istr:Iend, N(ng), NBT) :: state_old ! Old state (accounting for clipping)
      real(r8), dimension(Istr:Iend, N(ng), NBT) :: dBdt      ! FABM term (units [C/s])
      real(r8), dimension(Istr:Iend, NBT) :: flux_sf          ! FABM surface fluxes (units [Matter/m2/s])
      real(r8), dimension(Istr:Iend, NSAT) :: sms_sf          ! FABM surface-attached terms (units [Matter/m2/s])
      real(r8), dimension(Istr:Iend, NBT) :: flux_bt          ! FABM bottom fluxes (units [Matter/m2/s])
      real(r8), dimension(Istr:Iend, NBAT) :: sms_bt          ! FABM bottom-attached terms (units [Matter/m2/s])
      real(r8), dimension(Istr:Iend, 0:N(ng)) :: FC
      real(r8), dimension(Istr:Iend, N(ng)) :: Hz_inv
      real(r8), dimension(Istr:Iend, N(ng)) :: Hz_inv2
      real(r8), dimension(Istr:Iend, N(ng)) :: Hz_inv3
      real(r8), dimension(Istr:Iend, N(ng)) :: WL
      real(r8), dimension(Istr:Iend, N(ng)) :: WR
      real(r8), dimension(Istr:Iend, N(ng)) :: bL
      real(r8), dimension(Istr:Iend, N(ng)) :: bR
      real(r8), dimension(Istr:Iend, N(ng)) :: qc

! Chunksize variables for calls to FABM APIs
      nx = Iend-Istr+1
      ny = Jend-Jstr+1

!-----------------------------------------------------------------------
!  If wetting and drying, re-point FABM to full mask on first time step
!-----------------------------------------------------------------------
# ifdef WET_DRY
      IF (iic(ng).eq.1) THEN
        CALL FMODEL(ng)%f(tile)%model%set_mask(                         &
            & rmask_full(Istr:Iend,Jstr:Jend))
      END IF
# endif

! Relink environmental data as required (due to variable nstp).
!-----------------------------------------------------------------------
!Link environmental data to FABM (for variables that must be relinked for each time step since nstp varies)
!Note: We cannot use variable_needs_values to parse these calls, since these will return
!      .false. after the first linking calls in roms_fabm.F.
      IF (FMODEL(ng)%needsT(tile)) THEN
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
      &                    FMODEL(ng)%id_temp(tile),                     &
      &                    t(Istr:Iend, Jstr:Jend, 1:N(ng), nstp, 1)) ! Temperature is indezed 1st
      END IF

      IF (FMODEL(ng)%needsS(tile)) THEN
#ifdef FABM_NONNEG_S
        DO k=1,N(ng)
          DO j=J_range
            DO i=I_range
              Snn(i,j,k) = MAX(0.0_r8, t(i, j, k, nstp, 2)) ! Salinity is indexed 2nd
            END DO
          END DO
        END DO
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
      &                    FMODEL(ng)%id_salt(tile),                     &
      &                    Snn(Istr:Iend,Jstr:Jend,1:N(ng)))
#else
        CALL FMODEL(ng)%f(tile)%model%link_interior_data(               &
      &                    FMODEL(ng)%id_salt(tile),                     &
      &                    t(Istr:Iend,Jstr:Jend,1:N(ng),nstp,2))
#endif
      END IF

#ifdef FABM_TSS_ONLINE
      !-----------------------------------------------------------------------
      ! If needed, update the total suspended sediment.
      ! (Done here to ensure use of latest values (operator splitting))
      !-----------------------------------------------------------------------
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tss(i,j,k) = 0.0_r8
          END DO
        END DO
      END DO

      ! 
      DO itrc=1,NST
        ivar=idsed(itrc)
        DO k=1,N(ng)
          DO j=J_range
            DO i=I_range
#  ifdef MASKING
              tss(i,j,k) = tss(i,j,k)+rmask_full(i,j)*t(i,j,k,nstp,ivar) 
#  else
              tss(i,j,k) = tss(i,j,k)+t(i,j,k,nstp,ivar)
#  endif
            END DO
          END DO
        END DO
      END DO
#endif

!-----------------------------------------------------------------------
! Relink surface/bottom state data as required (due to variable nstp).
!-----------------------------------------------------------------------
      DO itrc=1,NSAT
        CALL FMODEL(ng)%f(tile)%model%link_surface_state_data(          &
      &           itrc,state_sf(Istr:Iend,Jstr:Jend,nstp,itrc))
      END DO

      DO itrc=1,NBAT
        CALL FMODEL(ng)%f(tile)%model%link_bottom_state_data(           &
      &           itrc,state_bt(Istr:Iend,Jstr:Jend,nstp,itrc))
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
        ibio=idbio(itrc)
        DO k=1,N(ng)
          DO j=J_range
            DO i=I_range
              state1(i,j,k,itrc)=t(i,j,k,nstp,ibio) !Note: FABM is permanently linked to state1 (in roms_fabm.F)
            END DO
          END DO
        END DO
      END DO

#ifdef FABM_CHECK_STATE
      DO j1=1,ny
        DO k=1,N(ng) ! vertical loop
          CALL FMODEL(ng)%f(tile)%model%check_interior_state(1,nx,j1, &
      &        k,repair,valid_int)   !This will cap state1(Istr:Iend,Jstr+j1,k,1:NBT)
        END DO
      END DO
      IF (NSAT.gt.0) THEN
        DO j1=1,ny
          CALL FMODEL(ng)%f(tile)%model%check_surface_state(1,nx,j1,  &
      &        repair,valid_sf)      !This will cap state_sf(Istr:Iend,Jstr+j1,nstp,1:NSAT)
        END DO
      END IF
      IF (NBAT.gt.0) THEN
        DO j1=1,ny
          CALL FMODEL(ng)%f(tile)%model%check_bottom_state(1,nx,j1,   &
      &        repair,valid_bt)      !This will cap state_bt(Istr:Iend,Jstr+j1,nstp,1:NBAT)
        END DO
      END IF
#endif

!-----------------------------------------------------------------------
!  Set date variables
!-----------------------------------------------------------------------
CALL caldate (tdays(ng), yd_dp=yday) ! Modified for updated ARANGO code (30/07/2019)
      ydayc(ng) = yday-1.0_r8
#ifdef DEBUGFABM
      IF (iic(ng).le.icheckmax(ng)) THEN
        write(*,*) "tdays = ", tdays(ng)
        write(*,*) "yday, ydayc = ", yday, ydayc(ng)
      END IF
      !     dBdt1max = 1e3 !Maximum absolute rate of change from FABM model (def = 1e3)
                    !Set to zero to examine initial FABM output
      !     icheckmax = 0  !Maximum time step counter value (iig) for which MIN/MAX dBdt is checked
                    !Set to zero to switch off all checks
                    !Set to large value to check all time steps
#endif

!-----------------------------------------------------------------------
!  Prepare all fields FABM needs to compute SMS terms (e.g., light)
!-----------------------------------------------------------------------
      CALL FMODEL(ng)%f(tile)%model%prepare_inputs()

!-----------------------------------------------------------------------
!  Enter J loop to compute SMS terms
!-----------------------------------------------------------------------
      J_LOOP : DO j=J_range
          j1 = j-Jstr+1 !FABM chunks are always indexed 1:chunksize
                      !Therefore, all calls to FABM APIs should use index j1 instead of j

!-----------------------------------------------------------------------
!  Compute inverse thickness to avoid repeated divisions.
!-----------------------------------------------------------------------
          DO k=1,N(ng)
            DO i=I_range
              Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
            END DO
          END DO
          DO k=1,N(ng)-1
            DO i=I_range
              Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
            END DO
          END DO
          DO k=2,N(ng)-1
            DO i=I_range
              Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
            END DO
          END DO

!-----------------------------------------------------------------------
!  Store the old interior state.
!-----------------------------------------------------------------------
          DO itrc=1,NBT
            DO k=1,N(ng)
              DO i=I_range
                state_old(i,k,itrc)=state1(i,j,k,itrc)
              END DO
            END DO
          END DO

!-----------------------------------------------------------------------
!  Output grid parameters for this tile if req'd
!-----------------------------------------------------------------------
#ifdef DEBUGFABM
          IF (iic(ng).le.icheckmax(ng)) THEN
            write(*,*) "Istr,Iend,Jstr,Jend = ",Istr,Iend,Jstr,Jend
            write(*,*) "iic,icheckmax,j,j1 = ",iic(ng),icheckmax(ng),j,j1
# ifdef MASKING
            write(*,*) "MIN,MAX(rmask_full(Istr:Iend,j)) = ",             &
      &  MINVAL(rmask_full(Istr:Iend,j)), MAXVAL(rmask_full(Istr:Iend,j))
# endif
          END IF
# endif

!-----------------------------------------------------------------------
!  Generate SMS terms for pelagic variables
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          dBdt(Istr:Iend,k,1:NBT)=0.0_r8
          CALL FMODEL(ng)%f(tile)%model%get_interior_sources(1,nx,j1,k, &
        &                                        dBdt(Istr:Iend,k,1:NBT))
  !Note: the biological tracer indices idbio range over (NAT+NPT+NCS+NNS)+1:NT = NT-NBT+1:NT, see rfabm_mod.h
        END DO

#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do"
          write(*,*) "MIN,MAX(dBdt(Istr:Iend,1:N(ng),1:NBT)) = ", &
        &      MINVAL(dBdt(Istr:Iend,1:N(ng),1:NBT)), &
        &      MAXVAL(dBdt(Istr:Iend,1:N(ng),1:NBT))
        END IF
#endif

!-----------------------------------------------------------------------
!  Generate surface fluxes and SMS for surface-attached variables
!-----------------------------------------------------------------------
        flux_sf = 0.0_r8
        sms_sf = 0.0_r8

        !Note: surface fluxes may still be needed even if NSAT = 0
        CALL FMODEL(ng)%f(tile)%model%get_surface_sources(1,nx,j1,      &
        &         flux_sf(Istr:Iend,1:NBT),sms_sf(Istr:Iend,1:NSAT))

#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_surface"
          write(*,*) "MIN,MAX(flux_sf(Istr:Iend,:)) = ",                &
        &      MINVAL(flux_sf(Istr:Iend,:)), MAXVAL(flux_sf(Istr:Iend,:))
          IF (NSAT.gt.0) THEN
            write(*,*) "MIN,MAX(sms_sf(Istr:Iend,:)) = ",               &
        &        MINVAL(sms_sf(Istr:Iend,:)), MAXVAL(sms_sf(Istr:Iend,:))
          END IF
        END IF
#endif

!-----------------------------------------------------------------------
!  Generate bottom fluxes and SMS terms for bottom-attached variables
!-----------------------------------------------------------------------
        flux_bt = 0.0_r8
        sms_bt = 0.0_r8

        !Note: bottom fluxes may still be needed even if NBAT = 0
        CALL FMODEL(ng)%f(tile)%model%get_bottom_sources(1,nx,j1,       &
        &        flux_bt(Istr:Iend,1:NBT),sms_bt(Istr:Iend,1:NBAT))

#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm_do_bottom"
          write(*,*) "MIN,MAX(flux_bt(Istr:Iend,:)) = ",                &
        &      MINVAL(flux_bt(Istr:Iend,:)), MAXVAL(flux_bt(Istr:Iend,:))
          IF (NBAT.gt.0) THEN
            write(*,*) "MIN,MAX(sms_bt(Istr:Iend,:)) = ",               &
        &      MINVAL(sms_bt(Istr:Iend,:)), MAXVAL(sms_bt(Istr:Iend,:))
          END IF
        END IF
#endif

!-----------------------------------------------------------------------
!  Add contributions from surface and bottom fluxes
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO i=I_range
            dBdt(i,N(ng),itrc) = dBdt(i,N(ng),itrc) +                   &
        &              flux_sf(i,itrc)*Hz_inv(i,N(ng))
            dBdt(i,1,itrc) = dBdt(i,1,itrc) +                           &
        &              flux_bt(i,itrc)*Hz_inv(i,1)
          !E.g surface/bottom attached variables may have [mass/m2] while pelagic variables may have [mass/m3]
          END DO
        END DO


!-----------------------------------------------------------------------
!  Check the model state and FABM computations in detail if req'd
!-----------------------------------------------------------------------
!#ifdef DEBUGFABM
        !IF (iic(ng).le.icheckmax(ng)) THEN
          ! This is just an example of some old debugging code, to be adapted as required
          !dBdt1 = MAXVAL(ABS(dBdt(Istr:Iend,:,1:NBT)))

          !IF (dBdt1.ge.dBdt1max(ng)) THEN
          !  write(*,*) "dBdt1 = MAX(ABS(dBdt)) = ", dBdt1
          !END IF
        !END IF
!#endif

!-----------------------------------------------------------------------
!  Update tracers with rates of change from FABM (Euler step)
!-----------------------------------------------------------------------
        DO itrc=1,NBT
          DO k=1,N(ng)
            DO i=I_range
              state1(i,j,k,itrc) = state1(i,j,k,itrc) + dBdt(i,k,itrc)*dtFABM(ng)
            END DO
          END DO
        END DO

#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Updated water column state with FABM dBdt"
          write(*,*) "MIN,MAX(state1(Istr:Iend,j,1:N(ng),1:NBT))=",     &
            MINVAL(state1(Istr:Iend,j,1:N(ng),1:NBT)),                  &
        &      MAXVAL(state1(Istr:Iend,j,1:N(ng),1:NBT))
        END IF
#endif

!-----------------------------------------------------------------------
!  Update surface states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NSAT
          DO i=I_range
            state_sf(i,j,nnew,itrc) = MAX(state_sf(i,j,nstp,itrc) +     &
        &                         sms_sf(i,itrc)*dtFABM(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO

#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          IF (NSAT.gt.0) THEN
            write(*,*) "Updated surface states with FABM sms_sf"
            write(*,*) "MIN,MAX(state_sf(Istr:Iend,j,nnew,1:NSAT))=",   &
              MINVAL(state_sf(Istr:Iend,j,nnew,1:NSAT)),                &
        &        MAXVAL(state_sf(Istr:Iend,j,nnew,1:NSAT))
          END IF
        END IF
#endif

!-----------------------------------------------------------------------
!  Update bottom states with rate of change from FABM
!-----------------------------------------------------------------------
        DO itrc=1,NBAT
          DO i=I_range
            state_bt(i,j,nnew,itrc) = MAX(state_bt(i,j,nstp,itrc) +  &
        &                                 sms_bt(i,itrc)*dtFABM(ng), 0.0_r8)
          !Note: we DO impose a zero lower bound on non-tracer state variables (cf. "t" below)
          END DO
        END DO
#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          IF (NBAT.gt.0) THEN
            write(*,*) "Updated bottom states with FABM sms_bt"
            write(*,*) "MIN,MAX(state_bf(Istr:Iend,j,nnew,1:NBAT))=",   &
              MINVAL(state_bt(Istr:Iend,j,nnew,1:NBAT)),                &
        &        MAXVAL(state_bt(Istr:Iend,j,nnew,1:NBAT))
          END IF
        END IF
#endif

!-----------------------------------------------------------------------
!  Get vertical sinking velocities from FABM
!-----------------------------------------------------------------------
        DO k=1,N(ng)
          CALL FMODEL(ng)%f(tile)%model%get_vertical_movement(          &
        &           1,nx,j1,k,w(Istr:Iend,k,NT(ng)-NBT+1:NT(ng)))
        END DO

        w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)) = -1*                  &
        &        w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))
        !FABM outputs vertical velocities in [m/s] positive upward
        !ROMS code below expects w in [m/s] positive downward
#ifdef DEBUGFABM
        IF (iic(ng).le.icheckmax(ng)) THEN
          write(*,*) "Done fabm get_vertical_movement"
          write(*,*) "MIN,MAX(w(Istr:Iend,1:N(ng),NT-NBT+1:NT)) = ",    &
            MINVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng))),           &
        &      MAXVAL(w(Istr:Iend,1:N(ng),NT(ng)-NBT+1:NT(ng)))
        END IF
#endif


!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
! Modified by Paul Wallhead in ROMS, 08/03/2017
!  SINK_LOOP: DO isink=1,Nsink is not convenient for arbitrary fabm model
!  Instead we loop over all biol tracers and use IF statement

        SINK_LOOP: DO itrc=1,NBT
        ibio=idbio(itrc)
        IF (MAXVAL(ABS(w(Istr:Iend,1:N(ng),ibio))).gt.1.0E-12_r8) THEN
    !!!Note: 1.0E-12 m/s => 0.3 mm/decade (safely negligible)
    !!!End modification PWA 08/03/2017
    !
    !  Copy concentration of biological particulates into scratch array
    !  "qc" (q-central, restrict it to be positive) which is hereafter
    !  interpreted as a set of grid-box averaged values for biogeochemical
    !  constituent concentration.
    !
            DO k=1,N(ng)
              DO i=I_range
                qc(i,k)=state1(i,j,k,itrc) !PWA: Note this uses the bgc-updated state (hence operator splitting), consistent with ROMS biology modules
              END DO
            END DO
            DO k=N(ng)-1,1,-1
              DO i=I_range
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=I_range
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
    !
    !  Apply PPM monotonicity constraint to prevent oscillations within the
    !  grid box.
    !
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
    !
    !  Compute right and left side values (bR,bL) of parabolic segments
    !  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
    !
    !  NOTE: Although each parabolic segment is monotonic within its grid
    !        box, monotonicity of the whole profile is not guaranteed,
    !        because bL(k+1)-bR(k) may still have different sign than
    !        qc(i,k+1)-qc(i,k).  This possibility is excluded,
    !        after bL and bR are reconciled using WENO procedure.
    !
                  cff=(dltR-dltL)*Hz_inv3(i,k)
                  dltR=dltR-cff*Hz(i,j,k+1)
                  dltL=dltL+cff*Hz(i,j,k-1)
                  bR(i,k)=qc(i,k)+dltR
                  bL(i,k)=qc(i,k)-dltL
                  WR(i,k)=(2.0_r8*dltR-dltL)**2
                  WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO

            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=I_range
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO

            DO i=I_range
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))

#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif

#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)

#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)

#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO ! Istr,Iend
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=I_range
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO

!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            DO k=1,N(ng)
              DO i=I_range
!!!Modified PWA 21/06/2022
#ifdef CAP_WSINK_CFL
                    cff = MIN(maxCFLwsink(ng)*Hz(i,j,k), dtFABM(ng)*w(i,k,ibio)) !Cap vertical velocity using maximum CFL number
!Note: Imposing the limit here rather than in a separate loop after
!      retrieval from FABM gives a small performance advantage (~5% speedup in tests).
#else
                    cff=dtFABM(ng)*w(i,k,ibio) !Modified 12/11/2020 to allow negative w => buoyant particles
#endif
    !!!End modified PWA 21/06/2022
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
      !
      !  Finalize computation of flux: add fractional part.
      !
            DO k=1,N(ng)
              DO i=I_range
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
            &                    Hz(i,j,ks)*cu*                                &
            &                    (bL(i,ks)+                                    &
            &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
            &                         (1.5_r8-cu)*                             &
            &                         (bR(i,ks)+bL(i,ks)-                      &
            &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=I_range
                state1(i,j,k,itrc)=qc(i,k)+                             &
            &           (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
          END IF
        END DO SINK_LOOP

!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "state_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "state"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of nutrients
!  when advection causes tracer concentration to go negative.
!-----------------------------------------------------------------------

        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff = state1(i,j,k,itrc) - state_old(i,k,itrc)
#ifdef MASKING
              cff = cff*rmask_full(i,j)
!It seems safer to ensure that no cells are not updated when dry.
!(rmask_full = rmask_wet*rmask, where rmask is the imposed time-independent mask (0=dry, 1=wet))
#endif
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
            END DO
          END DO
        END DO
      END DO J_LOOP

#else
      subroutine biology_empty ()
#endif
      return
      end
