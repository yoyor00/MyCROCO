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
      module croco_fabm

#if defined BIOLOGY && defined CFABM
      implicit none

      PRIVATE
      PUBLIC  :: croco_fabm_init

      CONTAINS

      SUBROUTINE croco_fabm_init (tile)     !Istr,Iend,Jstr,Jend,tile)
        ! == Added by SINTEF 2024
        USE fabm
        USE fabm_config
        USE module_param

# include "grid.h"
# include "ocean3d.h"
# include "ocean2d.h"
# include "private_scratch.h"
# ifdef MPI
      include 'mpif.h'
# endif
# include "diagnostics.h"
# include "scalars.h"
# include "forces.h"
# include "mixing.h"

      integer :: tile, nx, ny, i, j, k
      integer :: itrc, ibio, NDbio2dmax, NDbio3dmax 
      integer trd, omp_get_thread_num
#include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      FMODEL%f(trd)%model => fabm_create_model()
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_create_model"
# endif
      nx = Iend-Istr+1
      ny = Jend-Jstr+1

      ! FABM time step
      dtFABM = DBLE(nFABM)*dt

      CALL FMODEL%f(trd)%model%set_domain(nx,ny,N)
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm set_domain"
# endif

# ifdef MASKING
        CALL FMODEL%f(trd)%model%set_mask(                         &
     &                  rmask(Istr:Iend,Jstr:Jend))
        ! Note: Only the permanent land mask rmask should be passed to FABM at this point
#  ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_set_mask"
#  endif
# endif

        !! Link the FABM state variables
        ! First link the interior state variables to FABM (using temporary variable "state1")
        ! This is a permanent link.
        DO itrc=1,NBT
          CALL FMODEL%f(trd)%model%link_interior_state_data(       &
     &     itrc,state1(Istr:Iend,Jstr:Jend,1:N,itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_interior_state_data"
# endif

        ! Next link the surface-attached state variables (using "state_sf" for nstp=1, cf. ana_initial.h)
        ! This link will have to be adjusted every time step in rfabm.h as nstp varies
        DO itrc=1,NSAT
          CALL FMODEL%f(trd)%model%link_surface_state_data(        &
     &     itrc,state_sf(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_surface_state_data"
# endif

        ! Next link the bottom-attached (benthic) state variables (using "state_bt" for nstp=1, cf. ana_initial.h)
        ! This link will have to be adjusted every time step in rfabm.h as nstp varies
        DO itrc=1,NBAT
          CALL FMODEL%f(trd)%model%link_bottom_state_data(         &
     &     itrc,state_bt(Istr:Iend,Jstr:Jend,1,itrc))
        END DO
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm link_bottom_state_data"
# endif

        ! Set the pointers for cell thickness and bgc forcings
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%cell_thickness)) THEN
          CALL FMODEL%f(trd)%model%link_interior_data(             & !cell thickness (m)
     &     fabm_standard_variables%cell_thickness,                      &
     &     Hz(Istr:Iend,Jstr:Jend,1:N))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%depth)) THEN
          CALL FMODEL%f(trd)%model%link_interior_data(             & !depth of layer midpoints (m)
     &     fabm_standard_variables%depth,                               &
     &     z_r(Istr:Iend,Jstr:Jend,1:N))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%temperature)) THEN
          FMODEL%needsT(trd) = .true.                                  !This will be used to parse subsequent linking calls
          CALL FMODEL%f(trd)%model%link_interior_data(             & !temperature [degC]
     &     fabm_standard_variables%temperature,                         &
     &     t(Istr:Iend,Jstr:Jend,1:N,1,1))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%practical_salinity)) THEN
          FMODEL%needsS(trd) = .true.                                  !This will be used to parse subsequent linking calls
          CALL FMODEL%f(trd)%model%link_interior_data(             & !salinity [psu]
     &     fabm_standard_variables%practical_salinity,                  &
     &     t(Istr:Iend,Jstr:Jend,1:N,1,2))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%density)) THEN
          CALL FMODEL%f(trd)%model%link_interior_data(             & !density [kg/m3]
     &     fabm_standard_variables%density,                             &
     &     rho(Istr:Iend,Jstr:Jend,1:N))
        END IF
!        IF (FMODEL%f(trd)%model%variable_needs_values(             &
!     &   fabm_standard_variables%pressure)) THEN
!          CALL FMODEL%f(trd)%model%link_interior_data(             & !pressure due to water [dbar]
!     &     fabm_standard_variables%pressure,                            &
!     &     OCEAN(ng)%pres(Istr:Iend,Jstr:Jend,1:N))
!        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &fabm_standard_variables%surface_downwelling_shortwave_flux)) THEN
          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !surface downwelling shortwave flux [W/m2]
     &     fabm_standard_variables%surface_downwelling_shortwave_flux,  &
     &     radsw(Istr:Iend,Jstr:Jend))
        END IF
!        IF (FMODEL%f(trd)%model%variable_needs_values(             &
!     &   fabm_standard_variables%wind_speed)) THEN
!          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !wind speed [m/s]
!     &     fabm_standard_variables%wind_speed,                          &
!     &     wspd(Istr:Iend,Jstr:Jend))
!        END IF
!        IF (FMODEL%f(trd)%model%variable_needs_values(             &
!     &   fabm_standard_variables%bottom_stress)) THEN
!          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !bottom stress [Pa]
!     &     fabm_standard_variables%bottom_stress,                       &
!     &     FORCES(ng)%bstrPa(Istr:Iend,Jstr:Jend))
!        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%bottom_depth_below_geoid)) THEN
          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !bottom depth below geoid [m]
     &     fabm_standard_variables%bottom_depth_below_geoid,            &
     &     h(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%longitude)) THEN
          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !longitude [degrees]
     &     fabm_standard_variables%longitude,                           &
     &     lonr(Istr:Iend,Jstr:Jend))
        END IF
        IF (FMODEL%f(trd)%model%variable_needs_values(             &
     &   fabm_standard_variables%latitude)) THEN
          CALL FMODEL%f(trd)%model%link_horizontal_data(           & !latitude [degrees]
     &     fabm_standard_variables%latitude,                            &
     &     latr(Istr:Iend,Jstr:Jend))
        END IF
!        IF (FMODEL(ng)%f(trd)%model%variable_needs_values(             &
!     &  fabm_standard_variables%number_of_days_since_start_of_the_year)) THEN
!          CALL FMODEL(ng)%f(trd)%model%link_scalar(                    & !year day (0-365) [days]
!     &fabm_standard_variables%number_of_days_since_start_of_the_year,   &
!     &     ydayc(ng))
!        END IF

!!! Additional forcing variables that require cpp activation
! # ifdef FABM_ASILT
!         IF (FMODEL(ng)%f(trd)%model%variable_needs_values(             &
!      &  fabm_standard_variables%absorption_of_silt)) THEN
!           CALL FMODEL(ng)%f(trd)%model%link_interior_data(             & !absorption due to silt [m-1]
!      &     fabm_standard_variables%absorption_of_silt,                  &
!      &     FORCES(ng)%Asilt(Istr:Iend,Jstr:Jend,1:N(ng)))
!         END IF
! # endif
!! Note: absorption of silt (needed by light_iop) is a custom variable for ERSEM, not a standard variable
!! It is actually the PAR attenuation due to silt
!! Constant values can be implemented more simply by adding to the fabm.yaml:
! silt_abs:
!   model: bulk_constant
!   parameters:
!     standard_name: absorption_of_silt
!     value: 0.1
! # if defined FABM_TSS || defined FABM_TSS_ONLINE
!         CALL FMODEL(ng)%f(trd)%model%link_interior_data(               & !total suspended sediment [kg/m3]
!      &   'tss',FORCES(ng)%tss(Istr:Iend,Jstr:Jend,1:N(ng)))
! # endif

!# ifdef FABM_PCO2ATM
!        IF (FMODEL(ng)%f(trd)%model%variable_needs_values(             &
!     &  fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air)) THEN
!          CALL FMODEL(ng)%f(trd)%model%link_horizontal_data(           & !partial pressure of atmospheric CO2 [uatm]
!     &  fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, & !note: the name in FABM is not quite correct (only equivalent if total atmospheric pressure = 1 atm)
!     &     FORCES(ng)%pCO2atm(Istr:Iend,Jstr:Jend))
!        END IF
!# endif
! # ifdef FABM_ADYTRACER
!         IF (FMODEL(ng)%f(trd)%model%variable_needs_values(             &
!      &  fabm_standard_variables%gelbstoff_absorption_satellite)) THEN
!           CALL FMODEL(ng)%f(trd)%model%link_horizontal_data(           & !absorption due to dissolved and detrital matter (adg) [m-1]
!      &  fabm_standard_variables%gelbstoff_absorption_satellite,         &
!      &     FORCES(ng)%ADY_0(Istr:Iend,Jstr:Jend))
!         END IF
! # endif
! # ifdef FABM_N3ATMDEPO
!         CALL FMODEL(ng)%f(trd)%model%link_horizontal_data(             & !atmospheric deposition of oxidized nitrogen [mmolN/m2/s]
!      &   'n3_deposition/flux',FORCES(ng)%N3atmd(Istr:Iend,Jstr:Jend))
! # endif
! # ifdef FABM_N4ATMDEPO
!         CALL FMODEL(ng)%f(trd)%model%link_horizontal_data(             & !atmospheric deposition of reduced nitrogen [mmolN/m2/s]
!      &   'n4_deposition/flux',FORCES(ng)%N4atmd(Istr:Iend,Jstr:Jend))
! # endif
! # ifdef FABM_ISOURCES
!         IF (FMODELS(ng)%fmodel(trd)%variable_needs_values(             &
!      &  standard_variables%absorption_of_silt)) THEN
!           CALL FMODELS(ng)%fmodel(trd)%link_interior_data(             & !internal sources via FABM (work in progress, not yet usable)
!      &     standard_variables%absorption_of_silt,                       &
!      &     FORCES(ng)%Asilt(Istr:Iend,Jstr:Jend,1:N(ng)))
!         END IF
! # endif
! # ifdef FABM_AICE
!         CALL FMODEL(ng)%f(trd)%model%link_horizontal_data(             & !fraction of cell covered by ice from internal sea-ice module (Hedstrom code only) []
!      &   'aice',ICE(ng)%ai(Istr:Iend,Jstr:Jend,1))
! # endif
! # if defined CALC_TDR && defined CALC_TSR
! #  ifdef FABM_TDR
!         CALL FMODEL(ng)%f(trd)%model%link_interior_data(               & !turbulent energy dissipation rate [m2/s3]
!      &   'tdr',MIXING(ng)%tdr(Istr:Iend,Jstr:Jend,1:N(ng)))
! #  endif
! #  ifdef FABM_MU
!         CALL FMODEL(ng)%f(trd)%model%link_interior_data(               & !dynamic viscosity [kg/m/s]
!      &   'mu',MIXING(ng)%mu(Istr:Iend,Jstr:Jend,1:N(ng)))
! #  endif
! #  ifdef FABM_NU
!         CALL FMODEL(ng)%f(trd)%model%link_interior_data(               & !kinematic viscosity [m2/s]
!      &   'nu',MIXING(ng)%nu(Istr:Iend,Jstr:Jend,1:N(ng)))
! #  endif
! #  ifdef FABM_TSR
!         CALL FMODEL(ng)%f(trd)%model%link_interior_data(               & !turbulent shear rate [s-1]
!      &   'tsr',MIXING(ng)%tsr(Istr:Iend,Jstr:Jend,1:N(ng)))
! #  endif
! # endif
# ifdef FABM_DT
        CALL FMODEL%f(trd)%model%link_scalar('dtFABM', dtFABM)   !FABM time step [s]
# endif

# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm_link_data calls"
# endif

        ! Set variable ids for pointers that will change depending on time step (via nstp)
        !Note: We cannot use variable_needs_values to parse these calls, since these will return
        !      .false. after the first linking calls above.
        IF (FMODEL%needsT(trd)) THEN
          FMODEL%id_temp(trd) =                                    &
     &      FMODEL%f(trd)%model%get_interior_variable_id(          &
     &      fabm_standard_variables%temperature)
        END IF
        IF (FMODEL%needsS(trd)) THEN
          FMODEL%id_salt(trd) =                                    &
     &      FMODEL%f(trd)%model%get_interior_variable_id(          &
     &       fabm_standard_variables%practical_salinity)
        END IF
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done setting runtime variable ids"
# endif

        ! Check this fabm model is ready
!        CALL fabm_check_ready(FMODEL(ng)%f(trd)%model)
        CALL FMODEL%f(trd)%model%start()
# ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done fabm start"
# endif

# ifdef FABM_INITIAL
        ! Call fabm_initialize_state for each layer and j coordinate
        ! NOTE: This will OVERRIDE any ROMS tracer initialization via .nc input
        !       --- call to ROMS_initialise occurs before call to roms_fabm_init in ocean.h
        DO k=1,N
          DO j=1,ny
            CALL FMODEL%f(trd)%model%initialize_interior_state(    &
     &                                 1,nx,j,k)
          END DO
        END DO
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N
            DO j=Jstr,Jend
              DO i=Istr,Iend
                t(i,j,k,1,ibio) = state1(i,j,k,itrc)
              END DO
            END DO
          END DO
        END DO
#  ifdef DEBUGFABM
        IF (Master) write(stdout,*) "Done initialize_interior_state"
#  endif
# endif

# if defined FABM_INITIAL || defined FABM_INITIAL_SB
        IF (NSAT.gt.0) THEN
          DO j=1,ny
            CALL FMODEL%f(trd)%model%initialize_surface_state(     &
     &                                        1,nx,j)
          END DO
#  ifdef DEBUGFABM
          IF (Master) write(stdout,*) "Done initialize_surface_state"
#  endif
        END IF
        IF (NBAT.gt.0) THEN
          DO j=1,ny
            CALL FMODEL%f(trd)%model%initialize_bottom_state(      &
     &                                        1,nx,j)
          END DO
#  ifdef DEBUGFABM
          IF (Master) write(stdout,*) "Done initialize_bottom_state"
#  endif
        END IF
# endif

      END DO

      RETURN
      END SUBROUTINE croco_fabm_init
      ! ==
#endif /* BIOLOGY & CFABM */

      end module croco_fabm
      