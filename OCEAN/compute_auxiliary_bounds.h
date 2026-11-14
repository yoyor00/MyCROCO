!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
/* Auxiliary module "compute_auxiliary_bounds.h":
---------------------------------------------------------------
  Compute derived bounds for the loop indices over a subdomain
 "tile". The extended bounds [labelled by suffix R] are designed
 to cover also the outer ghost points, if the subdomain "tile" is
 adjacent to the PHYSICAL boundary.

 NOTE: IstrR,IendR,JstrR,JendR computed by this module DO NOT COVER
 ghost points associated with the internal computational boundaries
 of MPI subdomains.

   This module also computes loop-bounds for U- and V-type variables
 which belong to the interior of the computational domain. These
 are labelled by suffixes U,V and they step one grid point inward
 from the side of the subdomain adjacent to the physical boundary.
   Conversely, for an internal subdomain [which does not have
 segments of the physical boundary] all variables with suffixes
 R,U,V are set to the same values are the corresponding non-suffixed
 variables.
   Because this module also contains type declarations for these
 bounds, it must be included just after the last type declaration
 inside a subroutine, but before the first executable statement.
*/


      integer IstrR,IendR,JstrR,JendR

#ifdef EW_PERIODIC
! `IstrO` is the OLS-extended tile (and plain `Istr` without OLS)
# define IstrU IstrO
#else
      integer IstrU
#endif
#ifdef NS_PERIODIC
# define JstrV JstrO
#else
      integer JstrV
#endif


#ifdef MPI_OLS_2D
# include "overlapping_schwarz_2d.h"
      integer :: IstrO, IendO, JstrO, JendO
      integer :: ols_barotropic_timestep

      ! Phase within a block of `MPI_OLS_2D_COMM_N_TIMES` barotropic steps.
      ! With `SOLVE3D`, `step2d` runs `nfast` times per baroclinic step while
      ! `iic` advances only once, hence use `iif` (1..nfast).
# ifdef SOLVE3D
      if (iif.ge.1) then
        ols_barotropic_timestep = MOD(iif-1,
     &              MPI_OLS_2D_COMM_N_TIMES)
      else
        ols_barotropic_timestep = 0
      endif
# else
      ols_barotropic_timestep = MOD(iic-1,
     &              MPI_OLS_2D_COMM_N_TIMES)
# endif

      ! Reduce the extended halo layers within the block by
      ! `MPI_OLS_2D_TS_HALO_REDUCTION` per barotropic step.
# ifdef MPI_OLS_2D_TS_HALO_REDUCTION
#  define MPI_OLS_2D_ADD_AUX_LAYERS (MPI_OLS_2D_ADD_LAYERS-MPI_OLS_2D_TS_HALO_REDUCTION*ols_barotropic_timestep)
# else
#  define MPI_OLS_2D_ADD_AUX_LAYERS MPI_OLS_2D_ADD_LAYERS
# endif


# ifdef MPI_OLS_2D_EXTEND_IJ
      ! Bounds of the OLS-extended tile. `Istr`..`Jend` keep their original
      ! meaning; the extended bounds are used at the exchange call sites.
#  ifdef EW_PERIODIC
      IstrO = Istr - MPI_OLS_2D_ADD_AUX_LAYERS
      IendO = Iend + MPI_OLS_2D_ADD_AUX_LAYERS
#  else
      IstrO = Istr
      IendO = Iend
      if (.not. (WESTERN_EDGE)) then
        IstrO = IstrO - MPI_OLS_2D_ADD_AUX_LAYERS
      endif
      if (.not. (EASTERN_EDGE)) then
        IendO = IendO + MPI_OLS_2D_ADD_AUX_LAYERS
      endif
#  endif

#  ifdef NS_PERIODIC
      JstrO = Jstr - MPI_OLS_2D_ADD_AUX_LAYERS
      JendO = Jend + MPI_OLS_2D_ADD_AUX_LAYERS
#  else
      JstrO = Jstr
      JendO = Jend
      if (.not. (SOUTHERN_EDGE)) then
        JstrO = JstrO - MPI_OLS_2D_ADD_AUX_LAYERS
      endif
      if (.not. (NORTHERN_EDGE)) then
        JendO = JendO + MPI_OLS_2D_ADD_AUX_LAYERS
      endif
#  endif
# else
      IstrO = Istr
      IendO = Iend
      JstrO = Jstr
      JendO = Jend
# endif

#else

#  define IstrO Istr
#  define IendO Iend
#  define JstrO Jstr
#  define JendO Jend

#endif


      ! With OLS the derived bounds (`IstrR`, `IstrU`, ...) must be based on
      ! the extended tile, exactly as if `Istr`..`Jend` had been widened.
#ifdef MPI_OLS_2D
# define Istr_ IstrO
# define Iend_ IendO
# define Jstr_ JstrO
# define Jend_ JendO
#else
# define Istr_ Istr
# define Iend_ Iend
# define Jstr_ Jstr
# define Jend_ Jend
#endif

      if (WESTERN_EDGE) then
#ifdef EW_PERIODIC
        IstrR=Istr_-2
#else
        IstrR=Istr_-1
        IstrU=Istr_+1
#endif
      else
        IstrR=Istr_
#ifndef EW_PERIODIC
        IstrU=Istr_
#endif
      endif

      if (EASTERN_EDGE) then
#ifdef EW_PERIODIC
        IendR=Iend_+2
#else
        IendR=Iend_+1
#endif
      else
        IendR=Iend_
      endif

      if (SOUTHERN_EDGE) then
#ifdef NS_PERIODIC
        JstrR=Jstr_-2
#else
        JstrR=Jstr_-1
        JstrV=Jstr_+1
#endif
      else
        JstrR=Jstr_
#ifndef NS_PERIODIC
        JstrV=Jstr_
#endif
      endif

      if (NORTHERN_EDGE) then
#ifdef NS_PERIODIC
        JendR=Jend_+2
#else
        JendR=Jend_+1
#endif
      else
        JendR=Jend_
      endif

#undef Istr_
#undef Iend_
#undef Jstr_
#undef Jend_
