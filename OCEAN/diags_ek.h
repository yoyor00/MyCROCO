! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!
#ifdef DIAGNOSTICS_EK
!
!------------------------------------------------
!  Main arrays
!------------------------------------------------
!
! Depth integrated diags -------------
!
      real ekHadv(GLOBAL_2D_ARRAY)
      real ekHdiff(GLOBAL_2D_ARRAY)
      real ekVadv(GLOBAL_2D_ARRAY)
      real ekCor(GLOBAL_2D_ARRAY)
      real ekPrsgrd(GLOBAL_2D_ARRAY)
      real ekHmix(GLOBAL_2D_ARRAY)
      real ekVmix(GLOBAL_2D_ARRAY)
      real ekrate(GLOBAL_2D_ARRAY)
      real ekvol(GLOBAL_2D_ARRAY)
      real ekVmix2(GLOBAL_2D_ARRAY)
      real ekWind(GLOBAL_2D_ARRAY)
      real ekDrag(GLOBAL_2D_ARRAY)
# ifdef DIAGNOSTICS_BARO
      real ekBaro(GLOBAL_2D_ARRAY)
# endif
# ifdef M3FAST
      real ekfast(GLOBAL_2D_ARRAY)
# endif
# ifdef AVERAGES
      real timediags_ek_avg
      real ekHadv_avg(GLOBAL_2D_ARRAY)
      real ekHdiff_avg(GLOBAL_2D_ARRAY)
      real ekVadv_avg(GLOBAL_2D_ARRAY)
      real ekCor_avg(GLOBAL_2D_ARRAY)
      real ekPrsgrd_avg(GLOBAL_2D_ARRAY)
      real ekHmix_avg(GLOBAL_2D_ARRAY)
      real ekVmix_avg(GLOBAL_2D_ARRAY)
      real ekrate_avg(GLOBAL_2D_ARRAY)
      real ekvol_avg(GLOBAL_2D_ARRAY)
      real ekVmix2_avg(GLOBAL_2D_ARRAY)
      real ekWind_avg(GLOBAL_2D_ARRAY)
      real ekDrag_avg(GLOBAL_2D_ARRAY)
#  ifdef DIAGNOSTICS_BARO
      real ekBaro_avg(GLOBAL_2D_ARRAY)
#  endif
#  ifdef M3FAST
      real ekfast_avg(GLOBAL_2D_ARRAY)
#  endif
# endif /* AVERAGES */
      common /diag_ekHadv/ekHadv
     &       /diag_ekHdiff/ekHdiff
     &       /diag_ekVadv/ekVadv
     &       /diag_ekCor/ekCor
     &       /diag_ekPrsgrd/ekPrsgrd
     &       /diag_ekHmix/ekHmix
     &       /diag_ekVmix/ekVmix
     &       /diag_ekrate/ekrate
     &       /diag_ekvol/ekvol
     &       /diag_ekVmix2/ekVmix2
     &       /diag_ekWind/ekWind
     &       /diag_ekDrag/ekDrag
# ifdef DIAGNOSTICS_BARO
     &       /diag_ekBaro/ekBaro
# endif
# ifdef M3FAST
     &       /diag_ekfast/ekfast
# endif
# ifdef AVERAGES
      common /diag_timediags_ek_avg/timediags_ek_avg
      common /diag_ekHadv_avg/ekHadv_avg
     &       /diag_ekHdiff_avg/ekHdiff_avg
     &       /diag_ekVadv_avg/ekVadv_avg
     &       /diag_ekCor_avg/ekCor_avg
     &       /diag_ekPrsgrd_avg/ekPrsgrd_avg
     &       /diag_ekHmix_avg/ekHmix_avg
     &       /diag_ekVmix_avg/ekVmix_avg
     &       /diag_ekrate_avg/ekrate_avg
     &       /diag_ekvol_avg/ekvol_avg
     &       /diag_ekVmix2_avg/ekVmix2_avg
     &       /diag_ekWind_avg/ekWind_avg
     &       /diag_ekDrag_avg/ekDrag_avg
#  ifdef DIAGNOSTICS_BARO
     &       /diag_ekBaro_avg/ekBaro_avg
#  endif
#  ifdef M3FAST
     &       /diag_ekfast_avg/ekfast_avg
#  endif
# endif  /* AVERAGES */
!
! MLD integrated diags -------------
!
# ifdef DIAGNOSTICS_EK_MLD
      real ekHadv_mld(GLOBAL_2D_ARRAY)
      real ekHdiff_mld(GLOBAL_2D_ARRAY)
      real ekVadv_mld(GLOBAL_2D_ARRAY)
      real ekCor_mld(GLOBAL_2D_ARRAY)
      real ekPrsgrd_mld(GLOBAL_2D_ARRAY)
      real ekHmix_mld(GLOBAL_2D_ARRAY)
      real ekVmix_mld(GLOBAL_2D_ARRAY)
      real ekrate_mld(GLOBAL_2D_ARRAY)
      real ekvol_mld(GLOBAL_2D_ARRAY)
      real ekVmix2_mld(GLOBAL_2D_ARRAY)
#  ifdef DIAGNOSTICS_BARO
      real ekBaro_mld(GLOBAL_2D_ARRAY)
#  endif
#  ifdef AVERAGES
      real ekHadv_mld_avg(GLOBAL_2D_ARRAY)
      real ekHdiff_mld_avg(GLOBAL_2D_ARRAY)
      real ekVadv_mld_avg(GLOBAL_2D_ARRAY)
      real ekCor_mld_avg(GLOBAL_2D_ARRAY)
      real ekPrsgrd_mld_avg(GLOBAL_2D_ARRAY)
      real ekHmix_mld_avg(GLOBAL_2D_ARRAY)
      real ekVmix_mld_avg(GLOBAL_2D_ARRAY)
      real ekrate_mld_avg(GLOBAL_2D_ARRAY)
      real ekvol_mld_avg(GLOBAL_2D_ARRAY)
      real ekVmix2_mld_avg(GLOBAL_2D_ARRAY)
#   ifdef DIAGNOSTICS_BARO
      real ekBaro_mld_avg(GLOBAL_2D_ARRAY)
#   endif
#  endif /* AVERAGES */
      common /diag_ekHadv_mld/ekHadv_mld
     &       /diag_ekHdiff_mld/ekHdiff_mld
     &       /diag_ekVadv_mld/ekVadv_mld
     &       /diag_ekCor_mld/ekCor_mld
     &       /diag_ekPrsgrd_mld/ekPrsgrd_mld
     &       /diag_ekHmix_mld/ekHmix_mld
     &       /diag_ekVmix_mld/ekVmix_mld
     &       /diag_ekrate_mld/ekrate_mld
     &       /diag_ekvol_mld/ekvol_mld
     &       /diag_ekVmix2_mld/ekVmix2_mld
#  ifdef DIAGNOSTICS_BARO
     &       /diag_ekBaro_mld/ekBaro_mld
#  endif
#  ifdef AVERAGES
      common /diag_ekHadv_mld_avg/ekHadv_mld_avg
     &       /diag_ekHdiff_mld_avg/ekHdiff_mld_avg
     &       /diag_ekVadv_mld_avg/ekVadv_mld_avg
     &       /diag_ekCor_mld_avg/ekCor_mld_avg
     &       /diag_ekPrsgrd_mld_avg/ekPrsgrd_mld_avg
     &       /diag_ekHmix_mld_avg/ekHmix_mld_avg
     &       /diag_ekVmix_mld_avg/ekVmix_mld_avg
     &       /diag_ekrate_mld_avg/ekrate_mld_avg
     &       /diag_ekvol_mld_avg/ekvol_mld_avg
     &       /diag_ekVmix2_mld_avg/ekVmix2_mld_avg
#   ifdef DIAGNOSTICS_BARO
     &       /diag_ekBaro_mld_avg/ekBaro_mld_avg
#   endif
#  endif
# endif  /* DIAGNOSTICS_EK_MLD */
!
!------------------------------------------------
!  Work arrays
!------------------------------------------------
!
      real ekwrkHadv(GLOBAL_2D_ARRAY,2)
      real ekwrkHdiff(GLOBAL_2D_ARRAY,2)
      real ekwrkVadv(GLOBAL_2D_ARRAY,2)
      real ekwrkCor(GLOBAL_2D_ARRAY,2)
      real ekwrkPrsgrd(GLOBAL_2D_ARRAY,2)
      real ekwrkHmix(GLOBAL_2D_ARRAY,2,2)
      real ekwrkVmix(GLOBAL_2D_ARRAY,2)
      real ekwrkrate(GLOBAL_2D_ARRAY,2)
      real ekwrkvol(GLOBAL_2D_ARRAY,2)
      real ekwrkVmix2(GLOBAL_2D_ARRAY,2)
      real ekwrkwind(GLOBAL_2D_ARRAY,2)
      real ekwrkdrag(GLOBAL_2D_ARRAY,2)
# ifdef DIAGNOSTICS_BARO
      real ekwrkBaro(GLOBAL_2D_ARRAY,2)
# endif
# ifdef M3FAST
      real ekwrkfast(GLOBAL_2D_ARRAY,2)
# endif
      common /diag_ekwrkHadv/ekwrkHadv
     &       /diag_ekwrkHdiff/ekwrkHdiff
     &       /diag_ekwrkVadv/ekwrkVadv
     &       /diag_ekwrkCor/ekwrkCor
     &       /diag_ekwrkPrsgrd/ekwrkPrsgrd
     &       /diag_ekwrkHmix/ekwrkHmix
     &       /diag_ekwrkVmix/ekwrkVmix
     &       /diag_ekwrkrate/ekwrkrate
     &       /diag_ekwrkvol/ekwrkvol
     &       /diag_ekwrkVmix2/ekwrkVmix2
     &       /diag_ekwrkwind/ekwrkwind
     &       /diag_ekwrkdrag/ekwrkdrag
# ifdef DIAGNOSTICS_BARO
     &       /diag_ekwrkBaro/ekwrkBaro
# endif
# ifdef M3FAST
     &       /diag_ekwrkfast/ekwrkfast
# endif
# ifdef DIAGNOSTICS_EK_MLD
      real ekwrkHadv_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkHdiff_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkVadv_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkCor_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkPrsgrd_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkHmix_mld(GLOBAL_2D_ARRAY,2,2)
      real ekwrkVmix_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkrate_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkvol_mld(GLOBAL_2D_ARRAY,2)
      real ekwrkVmix2_mld(GLOBAL_2D_ARRAY,2)
#  if defined DIAGNOSTICS_BARO
      real ekwrkBaro_mld(GLOBAL_2D_ARRAY,2)
#  endif
      common /diag_ekwrkHadv_mld/ekwrkHadv_mld
     &       /diag_ekwrkHdiff_mld/ekwrkHdiff_mld
     &       /diag_ekwrkVadv_mld/ekwrkVadv_mld
     &       /diag_ekwrkCor_mld/ekwrkCor_mld
     &       /diag_ekwrkPrsgrd_mld/ekwrkPrsgrd_mld
     &       /diag_ekwrkHmix_mld/ekwrkHmix_mld
     &       /diag_ekwrkVmix_mld/ekwrkVmix_mld
     &       /diag_ekwrkrate_mld/ekwrkrate_mld
     &       /diag_ekwrkvol_mld/ekwrkvol_mld
     &       /diag_ekwrkVmix2_mld/ekwrkVmix2_mld
#  if defined DIAGNOSTICS_BARO
     &       /diag_ekwrkBaro_mld/ekwrkBaro_mld
#  endif
# endif /* DIAGNOSTICS_EK_MLD */
!
!------------------------------------------------
! KE diags from momentum terms
!------------------------------------------------
!
# if defined DIAGNOSTICS_EK_FULL && !defined DIAGNOSTICS_UV
      real MXadv(GLOBAL_2D_ARRAY,N,2)
      real MYadv(GLOBAL_2D_ARRAY,N,2)
      real MHdiff(GLOBAL_2D_ARRAY,N,2)
      real MVadv(GLOBAL_2D_ARRAY,N,2)
      real MCor(GLOBAL_2D_ARRAY,N,2)
      real MPrsgrd(GLOBAL_2D_ARRAY,N,2)
      real MHmix(GLOBAL_2D_ARRAY,N,2,2)
      real MVmix(GLOBAL_2D_ARRAY,N,2)
      real Mbody(GLOBAL_2D_ARRAY,N,2)
#  ifdef DIAGNOSTICS_BARO
      real MBaro(GLOBAL_2D_ARRAY,N,2)
#  endif
#  ifdef M3FAST
      real Mfast(GLOBAL_2D_ARRAY,N,2)
#  endif
      common /diag_MXadv/MXadv
     &       /diag_MYadv/MYadv
     &       /diag_MHdiff/MHdiff
     &       /diag_MVadv/MVadv
     &       /diag_MCor/MCor
     &       /diag_MPrsgrd/MPrsgrd
     &       /diag_MHmix/MHmix
     &       /diag_MVmix/MVmix
     &       /diag_Mbody/Mbody
#  ifdef DIAGNOSTICS_BARO
     &       /diag_MBaro/MBaro
#  endif
#  ifdef M3FAST
     &       /diag_Mfast/Mfast
#  endif
# endif /* DIAGNOSTICS_EK_FULL */

#endif /* DIAGNOSTICS_EK */




