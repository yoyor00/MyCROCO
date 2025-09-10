!
! Control of latency hiding for distributed memory communication.
!
! Defining `MPI_OVERLAPPING_SCHWARZ_2D` will result in:
!
! - exchanging more layers than necessary for a single time step.
!
! - running multiple time steps without data exchange.
!
! - doing computations also in halo layers to compensate for not communicated data.
!

#ifdef MPI_OVERLAPPING_SCHWARZ_2D

# ifndef MPI_OVERLAPPING_SCHWARZ_2D_NUM_LAYERS
#  error "MPI_OVERLAPPING_SCHWARZ_2D_NUM_LAYERS not defined, but MPI_OVERLAPPING_SCHWARZ_2D activated"
# endif

# ifndef MPI_OVERLAPPING_SCHWARZ_2D_COMM_N_TIMES
#  define MPI_OVERLAPPING_SCHWARZ_2D_COMM_N_TIMES (MPI_OVERLAPPING_SCHWARZ_2D_NUM_LAYERS)
#  error "MPI_OVERLAPPING_SCHWARZ_2D_COMM_N_TIMES not defined, but MPI_OVERLAPPING_SCHWARZ_2D activated"
# endif

#else

# define MPI_OVERLAPPING_SCHWARZ_2D_NUM_LAYERS 0

#endif

