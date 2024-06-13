MODULE stoexternal
   !!======================================================================
   !!                       ***  MODULE  stoexternal  ***
   !! Purpose        : external resources provide by the model (user supplied)
   !!=====================================================================
   !!   lbc_lnk      : generic interface for lbc_lnk_3d and lbc_lnk_2d
   !!   ctl_nam      : generate error message if failed to read namelist
   !!----------------------------------------------------------------------

   ! include parameters from CROCO
#include "cppdefs.h"
   USE scalars
   IMPLICIT NONE
   PRIVATE
   ! include definition of grid and mask from CROCO
#include "grid.h"

   ! import MPI include file                                         |
   include 'mpif.h'

   ! Type of variables
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                           !: working precision
   INTEGER, PUBLIC, PARAMETER ::   i4 = SELECTED_INT_KIND( 9)        !: single precision (integer 4)
   INTEGER, PUBLIC, PARAMETER ::   i8 = SELECTED_INT_KIND(14)        !: double precision (integer 8)
   INTEGER, PUBLIC, PARAMETER ::   lc = 256                          !: Length of Character strings

   ! Problem dimension
   INTEGER, PUBLIC ::   jpi        !: size dim 1 of MPI tiles, INCLUDING ghost cells
   INTEGER, PUBLIC ::   jpj        !: size dim 2 of MPI tiles, INCLUDING ghost cells
   INTEGER, PUBLIC ::   jpk        !: size dimension 3
   INTEGER, PUBLIC ::   narea      !: index of local domain
   INTEGER, PUBLIC ::   mppsize    !: number of processes
   INTEGER, PUBLIC ::   jpiglo     !: size of global domain (dim 1)
   INTEGER, PUBLIC ::   jpjglo     !: size of global domain (dim 2)
   ! Starting indices of MPI tiles, EXCLUDING ghost cells
   INTEGER, PUBLIC ::   Istr2, Iend2, Jstr2, Jend2

   ! Description of the grid
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:), POINTER :: glamt, gphit   ! longitude and latitude
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:), ALLOCATABLE :: glamtglo   ! global longitude
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:), ALLOCATABLE :: gphitglo   ! global latitude
   INTEGER, PUBLIC, SAVE, DIMENSION(:), ALLOCATABLE    :: mig        ! index of grid point in global grid
   INTEGER, PUBLIC, SAVE, DIMENSION(:), ALLOCATABLE    :: mjg        ! index of grid point in global grid

   ! Description of the mask
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: rmask_sto  ! land/ocean mask at T-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: umask_sto  ! land/ocean mask at U-points
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), POINTER :: vmask_sto  ! land/ocean mask at V-points

   ! I/O parameters
   INTEGER, PUBLIC ::   numout      =    6      !: logical unit for output print; set to stdout; do not change
   LOGICAL, PUBLIC ::   lwm         = .TRUE.    !: true on the 1st processor only (always)
   LOGICAL, PUBLIC ::   lwp         = .TRUE.    !: true on the 1st processor only .OR. ln_ctl
   INTEGER, PUBLIC ::   numnam_ref  =   -1      !: logical unit for reference namelist
   INTEGER, PUBLIC ::   numnam_cfg  =   -1      !: logical unit for configuration specific namelist
   INTEGER, PUBLIC ::   numond      =   -1      !: logical unit for Output Namelist Dynamics

   ! Ensemble parameters
   CHARACTER(len=3), PUBLIC :: cn_mem='001'            !: charcater string with ensemble member index
   LOGICAL, PUBLIC          :: ln_ensemble = .FALSE.   !: control of ensemble simulations
   LOGICAL, PUBLIC          :: ln_ens_rst_in = .FALSE. !: use ensemble (T) or single (F) input restart file
   INTEGER, PUBLIC          :: nn_ens_size = 1         !: ensemble size
   INTEGER, PUBLIC          :: nn_ens_start = 1        !: index of the first ensemble member
   INTEGER, PUBLIC          :: nmember = 1             !: index of current ensemble member

   ! Public routines
   INTERFACE lbc_lnk
      MODULE PROCEDURE lbc_lnk_2d, lbc_lnk_3d
   END INTERFACE

   PUBLIC ctl_nam, ocean_2_stogen, broadcast_array, lbc_lnk

CONTAINS

   SUBROUTINE lbc_lnk_2d( pt2d, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   cd_type : type of grid (T, U, V, F)
      !!                psign not used in CROCO (no grid folding)
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign

#    if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      IF (cd_type=='T') THEN
        call exchange_r2d_tile (Istr2,Iend2,Jstr2,Jend2, pt2d)
      ELSEIF (cd_type=='U') THEN
        call exchange_u2d_tile (Istr2,Iend2,Jstr2,Jend2, pt2d)
      ELSEIF (cd_type=='V') THEN
        call exchange_v2d_tile (Istr2,Iend2,Jstr2,Jend2, pt2d)
      ELSE
        call exchange_r2d_tile (Istr2,Iend2,Jstr2,Jend2, pt2d)
      ENDIF
#     endif

   END SUBROUTINE lbc_lnk_2d


   SUBROUTINE lbc_lnk_3d( pt3d, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 3D array (non mpp case)
      !!
      !! ** Method  :   cd_type : type of grid (T, U, V, F)
      !!                psign not used in CROCO (no grid folding)
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pt3d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign

#    if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      IF (cd_type=='T') THEN
        call exchange_r3d_tile (Istr2,Iend2,Jstr2,Jend2, pt3d)
      ELSEIF (cd_type=='U') THEN
        call exchange_u3d_tile (Istr2,Iend2,Jstr2,Jend2, pt3d)
      ELSEIF (cd_type=='V') THEN
        call exchange_v3d_tile (Istr2,Iend2,Jstr2,Jend2, pt3d)
      ELSE
        call exchange_r3d_tile (Istr2,Iend2,Jstr2,Jend2, pt3d)
      ENDIF
#     endif

   END SUBROUTINE lbc_lnk_3d


   SUBROUTINE ctl_nam ( kios, cdnam, ldwp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Information when error while reading a namelist
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(inout) ::   kios      ! IO status after reading the namelist
      CHARACTER(len=*) , INTENT(in   ) ::   cdnam     ! group name of namelist for which error occurs
      CHARACTER(len=5)                 ::   clios     ! string to convert iostat in character for print
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      !!----------------------------------------------------------------------
      WRITE (clios, '(I5.0)') kios
      IF( kios < 0 ) THEN
         print *, 'W A R N I N G:  end of record or file while reading namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios)
      ENDIF

      IF( kios > 0 ) THEN
         print *, 'E R R O R :   misspelled variable in namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) 
      ENDIF
      kios = 0
      RETURN
   END SUBROUTINE ctl_nam


   SUBROUTINE ocean_2_stogen (tile)
      implicit none
      integer tile
#ifdef  ALLOW_SINGLE_BLOCK_MODE
C$    integer  trd, omp_get_thread_num
#endif
# include "compute_tile_bounds.h"
      call ocean_2_stogen_tile (Istr,Iend,Jstr,Jend)
      return
      end


   SUBROUTINE ocean_2_stogen_tile (Istr,Iend,Jstr,Jend)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ocean_2_stogen  ***
      !!
      !! ** Purpose :   provide CROCO parameters to STOGEN
      !!----------------------------------------------------------------------
      INTEGER :: ji1, jj1, jk1
      integer :: Istr,Iend,Jstr,Jend

      ! Open namelist files
      numnam_ref = 10 ; numnam_cfg = 11 ; lwm = .FALSE.
      OPEN(UNIT=numnam_ref,FILE='namelist_sto_ref',STATUS='OLD',FORM='FORMATTED',ACCESS='SEQUENTIAL')
      OPEN(UNIT=numnam_cfg,FILE='namelist_sto_cfg',STATUS='OLD',FORM='FORMATTED',ACCESS='SEQUENTIAL')

      ! Define grid size (for local subdomain) -- 
      ! QJ: should correcpond to GLOBAL_2D_ARRAY defined in set_global_definitions.h
      ! follow the rule: #if undef THREE_GHOST_POINTS & defined MPI
      ! need some updates for other options ...
      jpi = Lm+4+padd_X  !size_XI
      jpj = Mm+4+padd_E  !size_ETA
      jpk = N

      ! Define global grid size (with all subdomains)
      jpiglo = LLm
      jpjglo = MMm

      ! Define number of subdomain and index of local subdomain
      mppsize = NNODES
      narea = mynode + 1

      ! Define starting and ending indices of MPI tiles, excuding ghost cells
      Istr2 = Istr+2
      Iend2 = Iend+2
      Jstr2 = Jstr+2
      Jend2 = Jend+2

      ! Associate the grid pointers with the CROCO common block arrays
      glamt(1:jpi,1:jpj) => lonr
      gphit(1:jpi,1:jpj) => latr

      ! Associate the mask pointers with the CROCO common block arrays
      ALLOCATE(rmask_sto(1:jpi,1:jpj,1:jpk))
      ALLOCATE(umask_sto(1:jpi,1:jpj,1:jpk))
      ALLOCATE(vmask_sto(1:jpi,1:jpj,1:jpk))
      DO jk1 = 1, jpk
         rmask_sto(1:jpi,1:jpj,jk1:jk1) => rmask
         umask_sto(1:jpi,1:jpj,jk1:jk1) => umask
         vmask_sto(1:jpi,1:jpj,jk1:jk1) => vmask
      ENDDO

      ! Define index of grid points (of local subdomain) in global grid (all subdomains)
      ALLOCATE(mig(1:jpi))
      ALLOCATE(mjg(1:jpj))
      DO ji1 = 1, jpi
        mig(ji1) = ji1 + ii * Lm
      ENDDO
      DO jj = 1, jpj
        mjg(jj1) = jj1 + jj * Mm
      ENDDO

      ! Warning: glamtglo and gphitglo are left unallocated.
      ! Options using them (in stokernel) will fail
      ! (with error message: 'Incorrect grid in stokernel').

   END SUBROUTINE ocean_2_stogen_tile


   SUBROUTINE broadcast_array( ptab )
      REAL(wp), DIMENSION(:), INTENT(in) :: ptab   ! array to broadcast

      INTEGER :: ierr

#ifdef MPI
      CALL MPI_bcast(ptab, size(ptab,1), MPI_DOUBLE_PRECISION, &
     &                     0, MPI_COMM_WORLD, ierr)
#endif
 
   END SUBROUTINE broadcast_array

END MODULE stoexternal
