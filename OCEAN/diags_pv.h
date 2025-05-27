! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_PV
# if defined DIAGNOSTICS_DIAPYCNAL
      real bflux(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE bflux(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real diffusivity(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE diffusivity(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      real Mrhs(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE Mrhs(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trhs(GLOBAL_2D_ARRAY,N,4+ntrc_pas)
!CSDISTRIBUTE_RESHAPE Trhs(BLOCK_PATTERN,*) BLOCK_CLAUSE


# ifdef AVERAGES
      real timediags_pv_avg
#  if defined DIAGNOSTICS_DIAPYCNAL
      real bflux_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE bflux_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real diffusivity_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE diffusivity_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
      real Mrhs_avg(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE Mrhs_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trhs_avg(GLOBAL_2D_ARRAY,N,4+ntrc_pas)
!CSDISTRIBUTE_RESHAPE Trhs_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE

# endif

# if defined DIAGNOSTICS_DIAPYCNAL
      common /diag_bflux/bflux
     &       /diag_diffusivity/diffusivity
# endif
      common /diag_Mrhs/Mrhs
     &       /diag_Trhs/Trhs

# ifdef AVERAGES
      common /diag_timediags_pv_avg/timediags_pv_avg
#  if defined DIAGNOSTICS_DIAPYCNAL
      common /diag_bflux_avg/bflux_avg
     &       /diag_diffusivity_avg/diffusivity_avg
#  endif
      common /diag_Mrhs_avg/Mrhs_avg
     &       /diag_Trhs_avg/Trhs_avg
# endif


# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_UV && ! defined DIAGNOSTICS_EK_FULL
      real MXadv(GLOBAL_2D_ARRAY,N,2)
      real MYadv(GLOBAL_2D_ARRAY,N,2)
      real MHdiff(GLOBAL_2D_ARRAY,N,2)
      real MHmix(GLOBAL_2D_ARRAY,N,2,2)
      real MVmix(GLOBAL_2D_ARRAY,N,2)
# if defined DIAGNOSTICS_BARO
      real MBaro(GLOBAL_2D_ARRAY,N,2)
# endif
# if defined M3FAST
      real Mfast(GLOBAL_2D_ARRAY,N,2)
# endif
      common /diag_MXadv/MXadv
     &       /diag_MYadv/MYadv
     &       /diag_MHdiff/MHdiff
     &       /diag_MHmix/MHmix
     &       /diag_MVmix/MVmix

# if defined DIAGNOSTICS_BARO
      common /diag_MBaro/MBaro
# endif
# if defined M3FAST
      common /diag_Mfast/Mfast
# endif
# endif

# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_UV
      real MVmix2(GLOBAL_2D_ARRAY,N,2)
      real Mrate(GLOBAL_2D_ARRAY,N,2)
      common /diag_MVmix2/MVmix2
     &       /diag_Mrate/Mrate
# endif

# if defined DIAGNOSTICS_PV && ! defined DIAGNOSTICS_TS

      real THmix(GLOBAL_2D_ARRAY,N,NT)
      real TVmix(GLOBAL_2D_ARRAY,N,NT)
      real TForc(GLOBAL_2D_ARRAY,N,NT)
      real Trate(GLOBAL_2D_ARRAY,N,NT)
      real TXadv(GLOBAL_2D_ARRAY,N,NT)
      real TYadv(GLOBAL_2D_ARRAY,N,NT)
      real TVadv(GLOBAL_2D_ARRAY,N,NT)

      common /diag_TForc/TForc
     &       /diag_THmix/THmix
     &       /diag_TVmix/TVmix
     &       /diag_Trate/Trate
     &       /diag_TXadv/TXadv
     &       /diag_TYadv/TYadv
     &       /diag_TVadv/TVadv

# endif

#  if defined DIAGNOSTICS_DIAPYCNAL && ! defined DIAGNOSTICS_TRACER_ISO

      real dbdx(GLOBAL_2D_ARRAY,N)
      real dbdy(GLOBAL_2D_ARRAY,N)
      real dbdz(GLOBAL_2D_ARRAY,0:N)

      common /diag_dbdx/dbdx
     &       /diag_dbdy/dbdy
     &       /diag_dbdz/dbdz

      real TF_xHmix(GLOBAL_2D_ARRAY,N,NTA)
      real TF_yHmix(GLOBAL_2D_ARRAY,N,NTA)
      real TF_zHmix(GLOBAL_2D_ARRAY,0:N,NTA)
      real TF_zVmix(GLOBAL_2D_ARRAY,0:N,NTA)
      real TF_Vadv(GLOBAL_2D_ARRAY,0:N,NTA)

      common /diag_TF_xHmix/TF_xHmix
     &       /diag_TF_yHmix/TF_yHmix
     &       /diag_TF_zHmix/TF_zHmix
     &       /diag_TF_zVmix/TF_zVmix
     &       /diag_TF_Vadv/TF_Vadv

# endif

#endif /* DIAGNOSTICS_PV */




