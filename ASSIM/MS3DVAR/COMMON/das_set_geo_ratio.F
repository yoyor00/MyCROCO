#include "cppdefs.h"
      subroutine das_set_geo_ratio  
      implicit none
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
      integer i,j,isp,ibnd
      real Lg2, ssc, crr
!
! compute the entire model domain
! but the lower triangular matrices
!
# ifdef EW_PERIODIC
#  define IR_ENTIRE 1,LLm
#  define IU_ENTIRE 1,LLm
# else
#  define IR_ENTIRE 0,LLm+1
#  define IU_ENTIRE 1,LLm+1
# endif

# ifdef NS_PERIODIC
#  define JR_ENTIRE 1,MMm
#  define JV_ENTIRE 1,MMm
# else
#  define JR_ENTIRE 0,MMm+1
#  define JV_ENTIRE 1,MMm+1
# endif

      isp=4
      Lg2=80.0*80.0
      crr=5.0       !< 3.5, agoe?
      ssc=sin(3.1415926*crr/180.0)

      do j=JR_ENTIRE
        do i=IR_ENTIRE

          ibnd=LLm*MMm
          ibnd=min(ibnd,i)
          ibnd=min(ibnd,LLm+1-i)
          ibnd=min(ibnd,j)
          ibnd=min(ibnd,MMm+1-j)
          

          georatio(i,j)=1.0
          
          if ( abs(latr(i,j)) .lt. crr ) then
             georatio(i,j)=1.0 *
     &                 (abs(sin(3.1415926*latr(i,j)/180.)))**1.35 /ssc
          endif
             
!          georatio(i,j)=0.618*( 1.0-exp(-h(i,j)*h(i,j)/Lg2) )
!    &     *(1.0-exp(-float(ibnd)*float(ibnd)/(float(isp)*float(isp))))

        enddo
      enddo
        
!      
# undef IR_ENTIRE
# undef IU_ENTIRE
# undef JR_ENTIRE
# undef JV_ENTIRE

      return
      end
