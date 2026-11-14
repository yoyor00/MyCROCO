!
! OLS (overlapping Schwarz latency hiding) for distributed memory communication.
!
! Defining `MPI_OLS_2D` results in
! - exchanging more halo layers than needed for a single time step,
! - running multiple barotropic time steps without any data exchange,
! - computing also in halo layers to compensate for the missing data.
!

#ifdef MPI_OLS_2D

/*
  OLS overlaps the halo communication of distributed memory runs. Without MPI
  there is no communication to overlap at all, so the method degenerates to the
  plain scheme. Disabling it here keeps serial builds working (the
  `exchange_*_ols_tile` routines only exist in `MPI` builds). This header is
  reached via `cppdefs.h` -> `set_global_definitions.h`, i.e. before any
  `#ifdef MPI_OLS_2D` in the Fortran sources.
*/
# ifndef MPI
#  undef MPI_OLS_2D
# endif

# ifndef MPI_OLS_2D_ADD_LAYERS
#  error "MPI_OLS_2D_ADD_LAYERS not defined, but MPI_OLS_2D activated"
# endif

# ifndef MPI_OLS_2D_COMM_N_TIMES
#  define MPI_OLS_2D_COMM_N_TIMES (MPI_OLS_2D_ADD_LAYERS)
#  error "MPI_OLS_2D_COMM_N_TIMES not defined, but MPI_OLS_2D activated"
# endif

#else

# define MPI_OLS_2D_ADD_LAYERS 0

#endif

