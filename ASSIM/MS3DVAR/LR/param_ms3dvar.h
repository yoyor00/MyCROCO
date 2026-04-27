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
! param_ms3dvar.h - MS3DVAR Grid Parameters
! ==========================================
!
! This file contains ONLY the grid parameters specific to MS3DVAR.
! It does NOT include the full CROCO param.h (1155 lines) which has
! configurations for all test cases and physics modules.
!
! MS3DVAR-specific grid parameters:
! ----------------------------------
! MS3DVAR works with a HIGH-RESOLUTION geopotential grid (z-levels)
! and optionally a LOW-RESOLUTION grid for multi-scale analysis.
!
! Key Parameters:
! - LLmH, MMmH, NnH : High-resolution grid dimensions (geopotential levels)
! - nratio, nhalf   : Downsampling ratio (high-res to low-res)
! - LLm_lr, MMm_lr  : Low-resolution derived dimensions
! - NSUB_XL, NSUB_EL: Tiling for OpenMP parallelization ; L = Low reso. grid LR
!
! CPP Keys:
! - DAS          : Enable data assimilation (defined in cppdefs_ms3dvar.h)
! - DAS_IN_MBFGS : Special flag used only in das_lbfgs.F minimizer
!
! IMPORTANT:
! - DAS is incompatible with CROCO's AUTOTILING + OpenMP combination
! - DAS is incompatible with MPI in its current implementation
!
! USAGE:
! ------
! This file should be included by variant-specific param.h files
! (LR/param.h, MR/param.h, etc.) which can then override or extend
! these settings as needed.
!
!======================================================================

!----------------------------------------------------------------------
! MS3DVAR High-Resolution Grid Configuration
!----------------------------------------------------------------------
! This is the geopotential grid used for data assimilation.
! It typically matches your high-resolution ocean model grid.
!
! LLmH, MMmH : Horizontal dimensions (interior points)
! NnH        : Number of vertical geopotential levels
!----------------------------------------------------------------------

      integer  LLmH,  MMmH, NnH, NSUB_XL, NSUB_EL, NPPL

! CONFIGURATION: Edit this section for your domain
! -------------------------------------------------
! Example 1: WMED (West Mediterranean) - DEFAULT
!    parameter (LLmH=629,  MMmH=537, NnH=50)

! Example 2: High-resolution variant (commented out)
!     parameter (LLmH=384,  MMmH=390, NnH=66)

! Example 3: Low-resolution variant (commented out)
      parameter (LLmH=209,  MMmH=179, NnH=50)

! Example 4: Your custom configuration
!     parameter (LLmH=xxx,  MMmH=yyy, NnH=zz)

!----------------------------------------------------------------------
! MS3DVAR Low-Resolution Grid Configuration
!----------------------------------------------------------------------
! For multi-scale analysis, MS3DVAR can use a downsampled low-res grid.
!
! nratio : Downsampling ratio (1 = no downsampling, 2 = half resolution, etc.)
! nhalf  : Helper parameter for grid alignment
! LLm_lr, MMm_lr : Computed low-resolution dimensions
!
! Example: If LLmH=629, nratio=2 -> LLm_lr ≈ 314 (half resolution)
!----------------------------------------------------------------------

#define NRATIO_DEFINED
      integer  nratio, nhalf,  LLm_lr, Lm_lr,  MMm_lr, Mm_lr

! CONFIGURATION: Set downsampling ratio
! --------------------------------------
      parameter (nratio=1, nhalf=nratio/2+1)

! Compute low-resolution grid dimensions
      parameter (LLm_lr=(LLmH+2-nhalf)/nratio-1,
     &           MMm_lr=(MMmH+2-nhalf)/nratio-1)
      parameter (Lm_lr=LLm_lr, Mm_lr=MMm_lr)

!----------------------------------------------------------------------
! OpenMP Tiling Configuration
!----------------------------------------------------------------------
! For OpenMP parallelization, the domain is split into tiles (LR grid).
!
! NSUB_XL : Number of tiles in X direction
! NSUB_EL : Number of tiles in Y (Eta) direction
! NPPL    : Total number of tiles (= NSUB_XL * NSUB_EL)
!
! Adjust based on your number of OpenMP threads and domain size.
!----------------------------------------------------------------------

! CONFIGURATION: Set tiling
! --------------------------
!      parameter (NSUB_XL=6, NSUB_EL=6, NPPL=NSUB_XL*NSUB_EL)
!      parameter (NSUB_XL=4, NSUB_EL=4, NPPL=16)  ! Alternative: 4x4 = 16 tiles
      parameter (NSUB_XL=5, NSUB_EL=4, NPPL=20)

!----------------------------------------------------------------------
! Special Case: Minimizer Grid Configuration
!----------------------------------------------------------------------
! The L-BFGS minimizer (das_lbfgs.F) temporarily defines DAS_IN_MBFGS
! and needs to work with the high-resolution grid directly.
!----------------------------------------------------------------------

#if defined DAS_IN_MBFGS
! When DAS_IN_MBFGS is defined (only in das_lbfgs.F):
! Use high-resolution grid for minimization
      integer NSUB_X, NSUB_E, NPP
      integer  LLm, Lm,  MMm, Mm
      parameter (LLm=LLm_lr,  MMm=MMm_lr)
      parameter (NSUB_X=NSUB_XL, NSUB_E=NSUB_EL, NPP=NPPL,
     &           Lm=LLm, Mm=MMm)

#else
!----------------------------------------------------------------------
! Standard CROCO Grid Parameters
!----------------------------------------------------------------------
! These parameters are required by CROCO's standard source files.
! They are derived from the configuration above.
!----------------------------------------------------------------------

      integer  LLm, Lm, MMm, Mm, N, LLm0, MMm0

!----------------------------------------------------------------------
! Grid Dimensions Based on Configuration
!----------------------------------------------------------------------
! Set LLm0, MMm0, N based on your configuration name defined in
! cppdefs_ms3dvar.h (e.g., WMED, BENGUELA, etc.)
!----------------------------------------------------------------------

#if defined REGIONAL

# if defined WMED
! West Mediterranean configuration
! Use low-resolution grid derived from high-res
      parameter (LLm0=LLm_lr,  MMm0=MMm_lr,  N=NnH)

# elif defined BENGUELA
! Benguela upwelling system
      parameter (LLm0=41,   MMm0=42,   N=32)

# elif defined YOUR_CONFIG
! Your custom regional configuration
!     parameter (LLm0=xxx,  MMm0=yyy,  N=zz)

# else
! Default/generic regional configuration
! YOU MUST SET THESE VALUES FOR YOUR DOMAIN
      parameter (LLm0=LLm_lr,  MMm0=MMm_lr,  N=NnH)
!     parameter (LLm0=xx,   MMm0=xx,   N=xx)   ! Placeholder

# endif  /* Regional configuration selection */

#else
! For non-REGIONAL configurations (test cases, etc.)
! Set explicit dimensions here if needed
      parameter (LLm0=xx,   MMm0=xx,   N=xx)

#endif  /* REGIONAL */

!----------------------------------------------------------------------
! Derived CROCO Parameters
!----------------------------------------------------------------------
! These are standard CROCO parameters computed from LLm0, MMm0, N
!----------------------------------------------------------------------

      parameter (LLm=LLm0,  MMm=MMm0)
      parameter (Lm=LLm, Mm=MMm)

! Standard output unit
      integer stdout
#ifdef LOGFILE
      common /stdout/stdout
#else
      parameter (stdout=6)
#endif

! Padding parameters for staggered grid
      integer padd_X, padd_E
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)

! Vertical grid parameters
      integer Np, NpHz, N_sl
      parameter (N_sl=0)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))

! Horizontal grid parameters
      integer NSUB_X, NSUB_E, NPP
      integer size_XI, size_ETA
      integer sse, ssz, se, sz, N2d, N3d, N3dHz
      integer NSA
! NSA: number of arrays in the workspace (28 without NBQ, 35 with NBQ)
# ifdef NBQ
      parameter (NSA=35)
# else
      parameter (NSA=28)
# endif
! NWEIGHT: maximum number of weights for barotropic mode
      integer NWEIGHT
      parameter (NWEIGHT=1000)

! Point sources
#if defined PSOURCE || defined PSOURCE_MASS || defined PSOURCE_NCFILE
      integer Msrc
# ifdef RIVER
      parameter (Msrc=2)
# elif defined SEAGRASS
      parameter (Msrc=1)
# elif defined VILAINE
      parameter (Msrc=2)
# elif defined ESTUARY
      parameter (Msrc=1)
# else
      parameter (Msrc=30)
# endif
#endif

! Tiling configuration (default = low-res tiling)
      parameter (NSUB_X=NSUB_XL, NSUB_E=NSUB_EL, NPP=NPPL)

! Tile size computation
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)

! 2D/3D array size computation
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)

#endif    /* DAS_IN_MBFGS */

!----------------------------------------------------------------------
! End of param_ms3dvar.h grid section.
! Include CROCO's own param.h for tracer/physics parameters.
! The MS3DVAR guard tells param.h to skip its grid section
! (lines 1-425) since we already defined our own grid above.
! This avoids duplicating CROCO's param.h physics section and ensures
! MS3DVAR automatically tracks CROCO code evolution.
!----------------------------------------------------------------------

#define MS3DVAR
#include "croco_ocean_param.h"
