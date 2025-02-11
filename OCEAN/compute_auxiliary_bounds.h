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
#ifdef MPI
      integer IstrUm1, IendUp1
      integer JstrVm1, JendVp1
      integer Istrm1, Iendp1
      integer Jstrm1, Jendp1
#endif

      if (WESTERN_EDGE) then
#ifdef EW_PERIODIC
        IstrR=Istr-2
#else
        IstrR=Istr-1
        IstrU=Istr+1
#endif
      else
        IstrR=Istr
#ifndef EW_PERIODIC
        IstrU=Istr
#endif
      endif

      if (EASTERN_EDGE) then
#ifdef EW_PERIODIC
        IendR=Iend+2
#else
        IendR=Iend+1
#endif
      else
        IendR=Iend
      endif

      if (SOUTHERN_EDGE) then
#ifdef NS_PERIODIC
        JstrR=Jstr-2
#else
        JstrR=Jstr-1
        JstrV=Jstr+1
#endif
      else
        JstrR=Jstr
#ifndef NS_PERIODIC
        JstrV=Jstr
#endif
      endif

      if (NORTHERN_EDGE) then
#ifdef NS_PERIODIC
        JendR=Jend+2
#else
        JendR=Jend+1
#endif
      else
        JendR=Jend
      endif

#ifdef MPI
      if (WEST_INTER) then
        IstrUm1 = IstrU-1
	Istrm1  = Istr-1
      else
        IstrUm1 = max(IstrU-1,2)
	Istrm1  = max(Istr-1,1)
      endif
      if (EAST_INTER) then
        IendUp1 = Iend+1
	Iendp1  = Iend+1
      else
        IendUp1 = min(Iend+1,Lmmpi)
	Iendp1  = min(Iend+1,Lmmpi)
      endif

      if (SOUTH_INTER) then
        JstrVm1 = JstrV-1
	Jstrm1  = Jstr-1
      else
        JstrVm1 = max(JstrV-1,2)
	Jstrm1  = max(Jstr-1,1)
      endif
      if (NORTH_INTER) then
        JendVp1 = Jend+1
	Jendp1  = Jend+1
      else
        JendVp1 = min(Jend+1,Mmmpi)
	Jendp1  = min(Jend+1,Mmmpi)
      endif
#endif
