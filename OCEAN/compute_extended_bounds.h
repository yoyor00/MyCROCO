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
/* Auxiliary module "compute_extended_bounds.h":
------------------------------------------------------
 Bounds designed to cover interior points of an array and
 ghost points associated with PHYSICAL side boundaries AND
 also ghost points of internal computational boundaries of
 MPI-subdomains.

 In the case of shared memory code (non-MPI version)
 IstrR, IendR,JstrR,JendR are equivalent to that computed
 by "compute_auxiliary_bounds.h".
*/
      integer IstrR,IendR,JstrR,JendR
#ifdef MPI_OLS_2D
      integer :: IstrO, IendO, JstrO, JendO
#endif

      if (Istr.eq.1) then
#ifdef EW_PERIODIC
        IstrR=Istr-2
#else
# ifdef MPI
        if (WEST_INTER) then
          IstrR=Istr-2
        else
          IstrR=Istr-1
        endif
# else
        IstrR=Istr-1
# endif
#endif
      else
        IstrR=Istr
      endif

#ifdef MPI
      if (Iend.eq.Lmmpi) then
#else
      if (Iend.eq.Lm) then
#endif
#ifdef EW_PERIODIC
        IendR=Iend+2
#else
# ifdef MPI
        if (EAST_INTER) then
          IendR=Iend+2
        else
          IendR=Iend+1
        endif
# else
        IendR=Iend+1
# endif
#endif
      else
        IendR=Iend
      endif

      if (Jstr.eq.1) then
#ifdef NS_PERIODIC
        JstrR=Jstr-2
#else
# ifdef MPI
        if (SOUTH_INTER) then
          JstrR=Jstr-2
        else
          JstrR=Jstr-1
        endif
# else
        JstrR=Jstr-1
# endif
#endif
      else
        JstrR=Jstr
      endif

#ifdef MPI
      if (Jend.eq.Mmmpi) then
#else
      if (Jend.eq.Mm) then
#endif
#ifdef NS_PERIODIC
        JendR=Jend+2
#else
# ifdef MPI
        if (NORTH_INTER) then
          JendR=Jend+2
        else
          JendR=Jend+1
        endif
# else
        JendR=Jend+1
# endif
#endif
      else
        JendR=Jend
      endif


#ifdef MPI_OLS_2D
      ! Bounds of the OLS-extended tile; `IstrR`..`JendR` are extended as well,
      ! while `Istr`..`Jend` keep their original meaning.
# ifdef EW_PERIODIC
      IstrO = Istr - MPI_OLS_2D_ADD_LAYERS
      IendO = Iend + MPI_OLS_2D_ADD_LAYERS
      IstrR = IstrR - MPI_OLS_2D_ADD_LAYERS
      IendR = IendR + MPI_OLS_2D_ADD_LAYERS
# else
      IstrO = Istr
      IendO = Iend
      if (.not. (WESTERN_EDGE)) then
        IstrO = IstrO - MPI_OLS_2D_ADD_LAYERS
        IstrR = IstrR - MPI_OLS_2D_ADD_LAYERS
      endif
      if (.not. (EASTERN_EDGE)) then
        IendO = IendO + MPI_OLS_2D_ADD_LAYERS
        IendR = IendR + MPI_OLS_2D_ADD_LAYERS
      endif
# endif

# ifdef NS_PERIODIC
      JstrO = Jstr - MPI_OLS_2D_ADD_LAYERS
      JendO = Jend + MPI_OLS_2D_ADD_LAYERS
      JstrR = JstrR - MPI_OLS_2D_ADD_LAYERS
      JendR = JendR + MPI_OLS_2D_ADD_LAYERS
# else
      JstrO = Jstr
      JendO = Jend
      if (.not. (SOUTHERN_EDGE)) then
        JstrO = JstrO - MPI_OLS_2D_ADD_LAYERS
        JstrR = JstrR - MPI_OLS_2D_ADD_LAYERS
      endif
      if (.not. (NORTHERN_EDGE)) then
        JendO = JendO + MPI_OLS_2D_ADD_LAYERS
        JendR = JendR + MPI_OLS_2D_ADD_LAYERS
      endif
# endif
#else
#  define IstrO Istr
#  define IendO Iend
#  define JstrO Jstr
#  define JendO Jend
#endif
