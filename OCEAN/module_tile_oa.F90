#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS

      module module_tile_oa

!======================================================================
!
!> @brief Online Analysis (OA=OnlineAnalysis)
!! - Derived data type to handle the vectorization of
!!   the spatial structure of the OA state vector
!!   in the case of dual OpenMP-MPI parallelization of the 
!!   horizontal domain.
!
!! @details required when applying croco tiles in the case of tile-threads
!!  with or without dual OpenMP-MPI parallelization of the horizontal domain.
!!  Fields of the derived type inherits of pieces of code 
!!  from older OA version (2006).
!!  More history/info in source module_interface_oa.F90
!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Croco Tile-thread compliant version (2021) 
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo
!
!  REFERENCE:
!  
!======================================================================

      implicit none

      integer :: ntiles

      type tile_space_str

        integer                            :: imin,imax,jmin,jmax
        integer                            :: kmin,kmax  ! vertical dim added even though not tile dependent
        integer, dimension(:), allocatable :: begvs_oa   ! nzv_oa 
        integer, dimension(:), allocatable :: begvs3d_oa, kmin3d_oa ! nzvs_oa + 1
        integer, dimension(:), allocatable :: l2i_oa, l2j_oa ! nzvs_oa
        integer, dimension(:,:,:), allocatable :: ij2l_oa    ! imin:imax,jmin:jmax, nzv_oa

        ! BLD TODO remove from module_oa_variables since it's the 2D and 3D dimension
        ! of the 2D and 3D spatial structures - dependent from tiles 
        integer :: nzvs3d_oa, nzvs_oa

      end type tile_space_str

      type(tile_space_str), allocatable, dimension(:), target :: st
      type(tile_space_str), pointer :: pst => null()

      private

      public :: st, pst, ntiles, allocate_tile_space_str, deallocate_tile_space_str, &
                allocate_begvs_oa, allocate_begvs3d_oa, deallocate_tile_varspace_oa


      CONTAINS 

      subroutine allocate_tile_space_str( tile_size )

        implicit none

        integer, intent(in) :: tile_size

        if ( .not. allocated(st) ) then
            allocate(st(1:tile_size))
            ntiles = tile_size
        end if

        return
      end subroutine allocate_tile_space_str
                                                                                        

      subroutine allocate_begvs_oa( imin, imax, jmin, jmax, kmin, kmax, tile )

        use module_oa_variables, only : nzv_oa

        implicit none
        
        integer, intent(in) :: tile
        integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax

        st(tile)%imin = imin
        st(tile)%imax = imax
        st(tile)%jmin = jmin
        st(tile)%jmax = jmax
        st(tile)%kmin = kmin
        st(tile)%kmax = kmax

        if ( .not. allocated( st(tile)%begvs_oa ) ) then
            allocate( st(tile)%begvs_oa( nzv_oa + 1 ) )
        end if 

        if ( .not. allocated( st(tile)%ij2l_oa ) ) then
            allocate( st(tile)%ij2l_oa( imin:imax, jmin:jmax, nzv_oa ) )
        end if 

        return 
      end subroutine allocate_begvs_oa

      subroutine allocate_begvs3d_oa(tile)
        ! BLXD not true anymore tile dependent size
        !use module_oa_variables, only : nzvs_oa
        implicit none
        integer, intent(in) :: tile
        integer :: nzvs_oa

        nzvs_oa = st(tile)%nzvs_oa

        if ( .not. allocated( st(tile)%begvs3d_oa ) ) then
            allocate( st(tile)%begvs3d_oa( nzvs_oa + 1 ) )
        end if 

        ! BLXD check why kmin3d_oa as size nzvs_oa PLUS ONE
        if ( .not. allocated( st(tile)%kmin3d_oa ) ) then
            allocate( st(tile)%kmin3d_oa( nzvs_oa + 1 ) )
        end if 

        if ( .not. allocated( st(tile)%l2i_oa ) ) then
            allocate( st(tile)%l2i_oa( nzvs_oa ) )
        end if 

        if ( .not. allocated( st(tile)%l2j_oa ) ) then
            allocate( st(tile)%l2j_oa( nzvs_oa ) )
        end if 

        return
      end subroutine allocate_begvs3d_oa

      subroutine deallocate_tile_varspace_oa(tile)
        implicit none
        integer, intent(in) :: tile
 
        if ( allocated( st(tile)%begvs_oa ) ) then
            deallocate( st(tile)%begvs_oa )
        end if 
        if ( allocated( st(tile)%begvs3d_oa ) ) then
            deallocate( st(tile)%begvs3d_oa )
        end if 
        if ( allocated( st(tile)%kmin3d_oa ) ) then
            deallocate( st(tile)%kmin3d_oa )
        end if 
        if ( allocated( st(tile)%ij2l_oa ) ) then
            deallocate( st(tile)%ij2l_oa )
        end if 
        if ( allocated( st(tile)%l2i_oa ) ) then
            deallocate( st(tile)%l2i_oa )
        end if 
        if ( allocated( st(tile)%l2j_oa ) ) then
            deallocate( st(tile)%l2j_oa )
        end if 
        return
      end subroutine deallocate_tile_varspace_oa

      subroutine deallocate_tile_space_str( )
      implicit none
        ! BLXD removed since in online_spectral_diags after
        ! if ( allocated( st(tile)%begvs_oa ) ) then
        !    stop
        ! else
        if ( allocated( st ) ) then
            deallocate(st)
        end if
        return
      end subroutine deallocate_tile_space_str


end module module_tile_oa

#else /* ONLINE_ANALYSIS */
      module module_tile_oa_empty
      end module
#endif /* ONLINE_ANALYSIS */

