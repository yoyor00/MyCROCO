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

      integer :: IstrR,IendR,JstrR,JendR

#ifdef MPI_OVERLAPPING_SCHWARZ_2D
      integer :: Istr_orig, Iend_orig, Jstr_orig, Jend_orig
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


#ifdef MPI_OVERLAPPING_SCHWARZ_2D
      Istr_orig = Istr
      Iend_orig = Iend
      Jstr_orig = Jstr
      Jend_orig = Jend

# ifdef EW_PERIODIC

      IstrR = IstrR - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      IendR = IendR + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS

      Istr = Istr - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      Iend = Iend + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS

# else

      if (.not. (WESTERN_EDGE)) then
        IstrR = IstrR - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
        Istr = Istr - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      endif

      if (.not. (EASTERN_EDGE)) then
        IendR = IendR + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
        Iend = Iend + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      endif

# endif


# ifdef NS_PERIODIC

      JstrR = JstrR - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      JendR = JendR + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS

      Jstr = Jstr - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      Jend = Jend + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS

# else

      if (.not. (SOUTHERN_EDGE)) then
        JstrR = JstrR - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
        Jstr = Jstr - MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      endif

      if (.not. (NORTHERN_EDGE)) then
        JendR = JendR + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
        Jend = Jend + MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS
      endif

# endif
#else
#  define Istr_orig Istr
#  define Iend_orig Iend
#  define Jstr_orig Jstr
#  define Jend_orig Jend
#endif


#define IJ_ORIG Istr_orig,Iend_orig,Jstr_orig,Jend_orig
