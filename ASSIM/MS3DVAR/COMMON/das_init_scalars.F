#include "cppdefs.h"
      subroutine das_init_scalars
!
! prepare all the variables that are initialized
! by one single processor
!

      implicit none
#include "param.h"
#include "das_param.h"
#include "das_ocean.h"
#include "das_innov.h"
#include "das_lbgfs.h"
      real LCT,LCS, alpha, exphT,exphS,alphaT, alphaS
      real tclm1d(ndas), sclm1d(ndas)

      integer k
                         ! One-dimesional profiles of T/S are need
                         ! for extrapolating the background to every
                         ! model grid. Note: 1) data is given from surface
                         ! to bottom, but the DA grid from bottom to surface;
                         ! 2) the DA result is little influenced by them.
      data tclm1d/
     &  15.49977, 15.51811, 15.56084, 15.61235, 15.65109, 15.77505, 
     &  15.81649, 15.83611, 15.8651, 15.92468, 15.98445, 15.97575, 
     &  15.89624, 
     &  15.77042, 15.61898, 15.4514, 15.26752, 15.09449, 14.9323, 
     &  14.78689, 
     &  14.66382, 14.55353, 14.45278, 14.3566, 14.26343, 14.17232, 
     &  14.08381, 
     &  13.99825, 13.91781, 13.84177, 13.77105, 13.70481, 13.64334, 
     &  13.5863, 
     &  13.53486, 13.48823, 13.44806, 13.41255, 13.3812, 13.35467, 
     &  13.33107, 
     &  13.31229 /
      data sclm1d/ 
     &   38.30985, 38.33209, 38.36414, 38.38966, 38.40676, 38.42036, 
     &   38.42034, 38.42803, 38.44218, 38.46166, 38.48612, 38.51652,
     &   38.56033, 
     &   38.60695, 38.65437, 38.70253, 38.74897, 38.78949, 38.82145,
     &   38.84357, 
     &   38.85528, 38.85835, 38.85538, 38.84804, 38.83733, 38.82408,
     &   38.80991, 
     &   38.79455, 38.77901, 38.76347, 38.74842, 38.73407, 38.7205, 
     &   38.70749, 
     &   38.69595, 38.68532, 38.67633, 38.6681, 38.66053, 38.65361, 
     &   38.64771, 
     &   38.644 /

!      data tclm1d/30.125, 30.097, 29.970, 29.744, 29.428, 29.004,
!     &            28.456, 27.886, 27.393, 26.551, 26.037, 25.706,
!     &            25.590, 25.523, 25.410, 24.432, 23.079, 21.792,
!     &            19.600, 17.615, 15.506, 11.818, 7.813, 6.255/
!      data sclm1d/36.546, 36.505, 36.472, 36.455, 36.443, 36.432,
!     &            36.419, 36.409, 36.406, 36.425, 36.440, 36.456,
!     &            36.446, 36.467, 36.505, 36.592, 36.633, 36.574,
!     &            36.243, 36.009, 35.773, 35.287, 35.024, 34.896/
!

      do k=1,ndas
        z_das(k)=-z_dep(ndas-k+1)
        tc1d_das(k) = tclm1d(ndas-k+1)-0.5   ! help stabalize the stratification
        sc1d_das(k) = sclm1d(ndas-k+1)-0.3
      enddo

      costf=0.
      cost0=0.

      LCT=30.0
      LCS=35.0
      alpha=0.5
      do k=1,max_prf
        exphT=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCT*LCT)))
        exphS=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCS*LCS)))
        sio_ot(k) = sio_ot_raw/(exphT*exphT)
        sio_os(k) = sio_os_raw/(exphS*exphS)
      enddo
!whoi
      LCT=30.0
      LCS=35.0
      alpha=0.25
      do k=1,max_prf
        exphT=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCT*LCT)))
        exphS=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCS*LCS)))
        whoi_ot(k) = whoi_ot_raw/(exphT*exphT)
        whoi_os(k) = whoi_os_raw/(exphS*exphS)
      enddo
! ptsur
      LCT=30.0
      LCS=35.0
      alpha=0.25
      do k=1,max_prf
        exphT=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCT*LCT)))
        exphS=(1.0+ alpha*exp(-z_dep(k)*z_dep(k)/(LCS*LCS)))
        ptsur_ot(k) = ptsur_ot_raw/(exphT*exphT)
        ptsur_os(k) = ptsur_os_raw/(exphS*exphS)
      enddo
!
      LCT=50.0
      LCS=60.0
      alphaT=0.2
      alphaS=0.1
      do k=1,max_prf
        exphT=(1.0+ alphaT*exp(-z_dep(k)*z_dep(k)/(LCT*LCT)))
        exphS=(1.0+ alphaS*exp(-z_dep(k)*z_dep(k)/(LCS*LCS)))
        map_ot(k) = map_ot_raw/(exphT*exphT)
        map_os(k) = map_os_raw/(exphS*exphS)
      enddo


      return
      end
