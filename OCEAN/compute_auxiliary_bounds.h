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
      integer IstrUm1,Iendp1
      integer JstrVm1,Jendp1
#ifdef EW_PERIODIC
# define IstrU Istr
#else
      integer IstrU
#endif
#ifdef NS_PERIODIC
# define JstrV Jstr
#else
      integer JstrV
#endif

      if (WESTERN_EDGE) then
# ifdef EW_PERIODIC
        IstrR=Istr-2
# else
        IstrR=Istr-1
        IstrU=Istr+1
# endif
        IstrUm1 = max(IstrU-1,2)
      else
!       IstrR=Istr
#if defined MPI && defined K3FAST
        if (istr.eq.1 .and. WEST_INTER) then
          IstrR = Istr - 2
        else
          IstrR = Istr
        endif
#else
        IstrR=Istr
#endif
# ifndef EW_PERIODIC
        IstrU=Istr
# endif
        IstrUm1 = IstrU-1
      endif

      if (EASTERN_EDGE) then
# ifdef EW_PERIODIC
        IendR=Iend+2
# else
        IendR=Iend+1
# endif
        Iendp1=min(Iend+1,Lmmpi)
      else
!        IendR=Iend
#if defined MPI && defined K3FAST
        if (iend.eq.Lmmpi .and. EAST_INTER) then
          IendR = Iend + 2
        else
          IendR = Iend
        endif
#else
        IendR=Iend
#endif
        Iendp1=Iend+1
      endif

      if (SOUTHERN_EDGE) then
# ifdef NS_PERIODIC
        JstrR=Jstr-2
# else
        JstrR=Jstr-1
        JstrV=Jstr+1
# endif
        JstrVm1 = max(JstrV-1,2)
      else
!       JstrR=Jstr
#if defined MPI && defined K3FAST
        if (jstr.eq.1 .and. SOUTH_INTER) then
          JstrR = Jstr - 2
        else
          JstrR = Jstr
        endif
#else
        JstrR=Jstr
#endif
# ifndef NS_PERIODIC
        JstrV=Jstr
# endif
        JstrVm1 = JstrV-1
      endif

      if (NORTHERN_EDGE) then
# ifdef NS_PERIODIC
        JendR=Jend+2
# else
        JendR=Jend+1
# endif
        Jendp1=min(Jend+1,Mmmpi)
      else
!       JendR=Jend
#if defined MPI && defined K3FAST
        if (jend.eq.Mmmpi .and. NORTH_INTER) then
          JendR = Jend + 2
        else
          JendR = Jend
        endif
#else
        JendR=Jend
#endif
        Jendp1=Jend+1
      endif
