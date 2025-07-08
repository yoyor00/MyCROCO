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
#include "cppdefs.h"
#ifdef AGRIF
AGRIFISDEFINED
#endif
#ifdef OPENMP
OPENMPISDEFINED
#endif
#ifdef MPI
MPIISDEFINED
#endif
#if defined OA_COUPLING || defined OW_COUPLING
OACPLISDEFINED
#endif
#ifdef XIOS
XIOSISDEFINED
#endif
#ifdef KNBQ
KNBQISDEFINED
#endif
#ifdef KNHINT
KNHINTISDEFINED
#endif
#ifdef KHCOMP
KHCOMPISDEFINED
#endif
#ifdef KH3D
KH3DISDEFINED
#endif
#ifdef OPENACC
OPENACCISDEFINED
#endif
#ifdef KNBQ3 
KNBQ3ISDEFINED
#endif
#ifdef key_pisces_npzd
P2ZISDEFINED
#endif
#ifdef key_pisces_quota
P5ZISDEFINED
#endif
#ifdef KXX 
KXXISDEFINED
#endif
