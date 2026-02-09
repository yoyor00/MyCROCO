/* This is include file "das_covar.h". 
  --------------------------------------------
*/
! z-coord
!
      real costf, cost0, costf_adj
      real xmin(NDIM)
      real xgrd(NDIM)
      real smin(NXM)
      real ymin(NXM)
      real diag(NDIM)
      real wmin(NWORK)
      common /das_miniz/xmin,xgrd,smin,ymin,diag,wmin
      common /das_cost/costf, cost0, costf_adj
      
