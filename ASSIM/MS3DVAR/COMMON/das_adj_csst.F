#include "cppdefs.h"
#if defined DAS_MCSST || defined DAS_TMISST\
  || defined DAS_GOES_SST

      subroutine das_adj_csst (tile)
      implicit none
      integer tile
#include "param.h"
#include "das_param.h"
#include "compute_tile_bounds.h"
      call das_adj_csst_tile (Istr,Iend,Jstr,Jend) 
      return
      end

      subroutine das_adj_csst_tile (Istr,Iend,Jstr,Jend) 
      implicit none
      integer Istr,Iend,Jstr,Jend,    i,j,k
      real cft, cfm, cfg
#if defined DAS_MCSST && defined DAS_FDN_MCSST
      integer ifndlev      
#endif
#include "param.h"
# include "das_param.h"
#include "grid.h"
#include "scalars.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_lbgfs.h"
# include "das_innov.h"
!
# include "compute_auxiliary_bounds.h"

# ifdef EW_PERIODIC
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE  Istr,IendR
# endif

# ifdef NS_PERIODIC
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE  Jstr,JendR
# endif
!
      do j=JR_RANGE
         do i=IR_RANGE
# ifdef DAS_TMISST
           cft=t_s(i,j,NDAS,itemp)-sst_tmi(i,j)
           t_s_adj(i,j,NDAS,itemp)=t_s_adj(i,j,NDAS,itemp)
     &                     +cft*tmi_mask(i,j)*tmi_oin(i,j)
# endif
# ifdef DAS_MCSST
#  if !defined DAS_FDN_MCSST
           cfm=t_s(i,j,NDAS,itemp)-sst_mc(i,j)
           t_s_adj(i,j,NDAS,itemp)=t_s_adj(i,j,NDAS,itemp)
     &                       +cfm*mc_mask(i,j)*mc_oin(i,j)
#  endif 
# endif
# ifdef DAS_GOES_SST
           cfg=t_s(i,j,NDAS,itemp)-sst_goes(i,j)
           t_s_adj(i,j,NDAS,itemp)=t_s_adj(i,j,NDAS,itemp)
     &                       +cfg*goes_mask(i,j)*goes_oin(i,j)
# endif
         enddo
       enddo
!
#if defined DAS_MCSST && defined DAS_FDN_MCSST
      do ifndlev=1,nfndlev
       do j=JR_RANGE
        do i=IR_RANGE
           cfm=t_s(i,j,NDAS-(nfndlev-ifndlev),itemp)
     &         -sst_mc(i,j,ifndlev)
           t_s_adj(i,j,NDAS-(nfndlev-ifndlev),itemp)=
     &          t_s_adj(i,j,NDAS-(nfndlev-ifndlev),itemp)
     &         +cfm*mc_mask(i,j,ifndlev)*mc_oin(i,j,ifndlev)
        enddo
       enddo
      enddo
#endif
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
!
      return
      end
#else
      subroutine das_adj_csst_empty
      return
      end
#endif /* DAS_TMISST DAS_MCSST */

