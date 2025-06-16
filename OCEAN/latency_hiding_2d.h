!
! Control of latency hiding for distributed memory communication.
!
! Defining `MPI_LAT_HID_2D` will result in:
!
! - exchanging more layers than necessary for a single time step.
!
! - running multiple time steps without data exchange.
!
! - doing computations also in halo layers to compensate for not communicated data.
!

#ifdef MPI_LAT_HID_2D

# ifndef MPI_LAT_HID_2D_ADD_LAYERS
#  error "MPI_LAT_HID_2D_ADD_LAYERS not defined, but MPI_LAT_HID_2D activated"
# endif

# ifndef MPI_LAT_HID_2D_COMM_N_TIMES
#  define MPI_LAT_HID_2D_COMM_N_TIMES (MPI_LAT_HID_2D_ADD_LAYERS)
#  error "MPI_LAT_HID_2D_COMM_N_TIMES not defined, but MPI_LAT_HID_2D activated"
# endif

#endif
