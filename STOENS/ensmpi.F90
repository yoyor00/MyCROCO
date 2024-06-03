#include "cppdefs.h"
MODULE ensmpi
   !!======================================================================
   !!                       ***  MODULE  ensmpi  ***
   !! Ensemble simulation : run several members in parallel using MPI
   !!=====================================================================
   !! History :  2024-06 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   ens_comm_set : set up communicators for ensemble run
   !!----------------------------------------------------------------------
   ! include parameters from CROCO
   USE scalars

   IMPLICIT NONE
   PRIVATE

   ! include definition of local CROCO communicator
   ! (previosuly defined for OA_COUPLING)
#  include "mpi_cpl.h"

   ! Global communicator with all ensemble members
   INTEGER, PUBLIC :: mpi_comm_all

   ! Parameter of ensemble simulation
   INTEGER, PUBLIC :: nn_ens_size = 1         ! ensemble size
   INTEGER, PUBLIC :: nn_ens_start = 1        ! index of the first ensemble member
   LOGICAL, PUBLIC :: ln_ens_rst_in = .FALSE. ! use ensemble (T) or single (F) input restart file

   ! Arrays to store communicators for ensemble members
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE         :: ncomm_member   ! communicator for every ensemble member
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE :: ncomm_area     ! communicator for every subdomain

   ! Public routines
   PUBLIC :: ens_comm_set

CONTAINS

   SUBROUTINE ens_comm_set ( )
      !!----------------------------------------------------------------------
      !!                  ***  routine ens_comm_set  ***
      !!
      !! ** Purpose :   Set up communicators for ensemble run 
      !!
      !! ** Method  :   The global communicator is divided in nmember x mpi_com_croco
      !!----------------------------------------------------------------------
      INTEGER                            :: jpnij          ! number of MPI subdomains in each ensemble
      INTEGER                            :: jmember, jjproc
      INTEGER                            :: imember, impprank, ierr
      INTEGER, DIMENSION(:), ALLOCATABLE :: irank_member   ! list of processors for current member (dimension jpnij)
      INTEGER, DIMENSION(:), ALLOCATABLE :: igrp_member    ! group ID for every ensemble member     
      INTEGER, DIMENSION(:), ALLOCATABLE :: irank_area     ! list of processors for current subdomain (dimension nn_ens_size)
      INTEGER, DIMENSION(:), ALLOCATABLE :: igrp_area      ! group ID for every subdomain
      INTEGER ::                            ngrp_world     ! group ID for the world processors
      ! Get number of subdomains (for each member)
      jpnij = NNODES

      ! Read parameters for ensemble simulation
      CALL ens_param_read

      ! Create global group
      CALL mpi_comm_group( mpi_comm_all, ngrp_world, ierr )

      ! Arrays to store local communicators for ensemble members
      ALLOCATE( ncomm_member(nn_ens_size)  )     ! communicator for members
      ALLOCATE( irank_member(jpnij), igrp_member(nn_ens_size)  )   ! temporary array of domain rank for a specific member

      ! Create one communicator for each ensemble member
      DO jmember = 1, nn_ens_size
          irank_member(:) = (/ (jjproc, jjproc=0,jpnij-1) /)
          irank_member(:) = irank_member(:) + ( jmember - 1 ) * jpnij
          ! Create the group for current member from the world group
          CALL mpi_group_incl( ngrp_world, jpnij, irank_member, igrp_member(jmember), ierr )
          ! Create the communicator for current member
          CALL mpi_comm_create( mpi_comm_all, igrp_member(jmember), ncomm_member(jmember), ierr )
      ENDDO

      ! Get rank of processor in global communicator
      ! and identify to which member communicator it belongs
      CALL mpi_comm_rank( mpi_comm_all, impprank, ierr )
      imember = 1 + impprank / jpnij

      ! Set CROCO communicator to communicator for current ensemble member
      ocean_grid_comm = ncomm_member(imember)

      ! Return index of current ensemble member
      ! This must be used in CROCO to produce a different behaviour in different members, e.g.:
      ! - a different name for parameter files
      ! - a different name for input files
      ! - a different name for output files (needed)
      ! - a different name for XIOS context (needed if XIOS is used)
      ! - a different seed for the random number generator (needed if STOCHASTIC is used)
      kmember = nn_ens_start + imember - 1

      ! Deallocate arrays
      DEALLOCATE( irank_member , igrp_member )

   END SUBROUTINE ens_comm_set


   SUBROUTINE ens_param_read
      !!----------------------------------------------------------------------
      !!                  ***  routine ens_param_read  ***
      !!
      !! ** Purpose :   Read parameters for ensemble simulation
      !!
      !!----------------------------------------------------------------------

      ! Get ensemble size, index of first member, and restart option

   END SUBROUTINE ens_param_read

END MODULE ensmpi

