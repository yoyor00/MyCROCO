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
/* This is include file "ocean2d.h".
--------------------------------------------------------------------
 zeta,rheta     Free surface elevation [m] and its time tendency;
 ubar,rubar     Vertically integrated  2D velocity components in
 vbar,rvbar     XI- and ETA-directions and their time tendencies;
*/
      real zeta(GLOBAL_2D_ARRAY,4)
      real ubar(GLOBAL_2D_ARRAY,4)
      real vbar(GLOBAL_2D_ARRAY,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar

#if !defined SOLVE3D && defined M2_HADV_UP3
      real urhs(GLOBAL_2D_ARRAY)
      real vrhs(GLOBAL_2D_ARRAY)
      real Duon(GLOBAL_2D_ARRAY)
      real DVom(GLOBAL_2D_ARRAY)
#endif

# ifdef STO_XI2D
      real sto_xi2d(GLOBAL_2D_ARRAY)
      common /stogen2d/ sto_xi2d
# endif

