#include "cppdefs.h"
#if defined DAS_TMISST || defined DAS_MCSST\
  || defined DAS_GOES_SST

      subroutine das_innov_sst (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call das_innov_sst_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_innov_sst_tile (Istr,Iend,Jstr,Jend)
!
!--------------------------------------------------------------------
!  Set initial conditions for momentum and tracer variables using
!  analytical expressions.
!--------------------------------------------------------------------
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k, itrc
      real otmi, omc, ogoes,cft, crt, crt_h, crt_l
      real crt_hg, crt_lg
#if defined DAS_MCSST && defined DAS_FDN_MCSST
      integer ifndlev      
#endif
# include "param.h"
# include "das_param.h"
# include "scalars.h"
# include "grid.h"
# include "das_ocean.h"
# include "das_covar.h"
# include "das_innov.h"
!
!
# include "compute_auxiliary_bounds.h"
!
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
!  forecast error variance
!------------------------------------------------------------------
!
! Douglas A. May,* Michelle M. Parmeter, Daniel S. Olszewski, and Bruce 
! D. McKenzie, 1998: Operational Processing of Satellite
! Sea Surface Temperature Retrievals at the Naval Oceanographic Office.
! Bulletin of the American Meteorological Society, Vol. 79, No. 3, 397-406.
!   ::SST accuracy maintains a root-mean-square difference (rmsd) error of 
!     less than 0.7
!
      otmi=1.0*1.0
      omc=0.8*0.8
      ogoes=1.2*1.2
      crt_l=-2.5    ! mcsst, cloud contamination
      crt_h=2.80
      crt_lg=-2.0    ! mcsst, cloud contamination
      crt_hg=2.00
      crt=2.80
!
# ifdef DAS_TMISST
      do j=JR_RANGE
        do i=IR_RANGE
          cft=sst_tmi(i,j)-t_das(i,j,NDAS,itemp)
          if (abs(cft) .gt. crt ) tmi_mask(i,j)=0.0
          sst_tmi(i,j)=cft * tmi_mask(i,j)
          tmi_oin(i,j)=1.0/otmi
        enddo
      enddo
# endif
# ifdef DAS_MCSST
# if !defined DAS_FDN_MCSST
      do j=JR_RANGE
        do i=IR_RANGE
          cft=sst_mc(i,j)-t_das(i,j,NDAS,itemp)
          if ( (cft .gt. crt_h) .or. (cft .lt. crt_l)) 
     &                                mc_mask(i,j)=0.0
          sst_mc(i,j)=cft * mc_mask(i,j)
          mc_oin(i,j)=1.0/omc
        enddo
      enddo
# else
      do ifndlev=1,nfndlev
        do j=JR_RANGE
          do i=IR_RANGE
            cft=sst_mc(i,j,ifndlev)-t_das(i,j,NDAS-(nfndlev-ifndlev),itemp)
            if ( (cft .gt. crt_h) .or. (cft .lt. crt_l)) 
     &                                mc_mask(i,j,ifndlev)=0.0
            sst_mc(i,j,ifndlev)=cft * mc_mask(i,j,ifndlev)
            mc_oin(i,j,ifndlev)=fndlev(ifndlev)/omc
          enddo
        enddo
      enddo

# endif
# endif
# ifdef DAS_GOES_SST
      do j=JR_RANGE
        do i=IR_RANGE
          cft=sst_goes(i,j)-t_das(i,j,NDAS,itemp)
          if ((cft .gt. crt_hg) .or. (cft .lt. crt_lg))
     &                               goes_mask(i,j)=0.0
          sst_goes(i,j)=cft * goes_mask(i,j)
          goes_oin(i,j)=1.0/ogoes
        enddo
      enddo
# endif
!
#if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
# ifdef DAS_TMISST
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        sst_tmi(START_2D_ARRAY))
# endif
# ifdef DAS_MCSST
#  if !defined DAS_FDN_MCSST
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        sst_mc(START_2D_ARRAY))
#  else
      do ifndlev=1,nfndlev
       call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        sst_mc(START_2D_ARRAY,ifndlev))
      enddo
#  endif
# endif
# ifdef DAS_GOES_SST
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        sst_goes(START_2D_ARRAY))
# endif
#endif
!
# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
      return
      end
#else
      subroutine das_innov_sst_empty
      return
      end
#endif /* DAS_TMISST DAS_MCSST */
