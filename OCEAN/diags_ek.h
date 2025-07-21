! This is include file "diags_ek.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_EK
      real ekHadv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHdiff(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekvol(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekWind(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekWind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekDrag(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekDrag(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined DIAGNOSTICS_BARO
      real ekBaro(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekBaro(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
# if defined M3FAST
      real ekfast(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekfast(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif
      real ekwrkrate(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE ekwrkrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkvol(GLOBAL_2D_ARRAY,N,2)
!CSDISTRIBUTE_RESHAPE ekwrkvol(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkwind(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkwind(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekwrkdrag(GLOBAL_2D_ARRAY,2)
!CSDISTRIBUTE_RESHAPE ekwrkdrag(BLOCK_PATTERN,*) BLOCK_CLAUSE

# ifdef AVERAGES
      real timediags_ek_avg
      real ekHadv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHdiff_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHdiff_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVadv_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekCor_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekCor_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekPrsgrd_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekPrsgrd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekHmix_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekHmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekrate_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekrate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekvol_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekvol_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekVmix2_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekVmix2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekWind_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekWind_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real ekDrag_avg(GLOBAL_2D_ARRAY)
!CSDISTRIBUTE_RESHAPE ekDrag_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  if defined DIAGNOSTICS_BARO
      real ekBaro_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekBaro_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
#  endif
#  if defined M3FAST
      real ekfast_avg(GLOBAL_2D_ARRAY,N)
!CSDISTRIBUTE_RESHAPE ekfast_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
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
# if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro/ekBaro
# endif
# if defined M3FAST
     &       /diag_ekfast/ekfast
# endif
     &       /diag_ekwrkrate/ekwrkrate
     &       /diag_ekwrkvol/ekwrkvol
     &       /diag_ekwrkwind/ekwrkwind
     &       /diag_ekwrkdrag/ekwrkdrag
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
#  if defined DIAGNOSTICS_BARO
     &       /diag_ekBaro_avg/ekBaro_avg
#  endif
#  if defined M3FAST
     &       /diag_ekfast_avg/ekfast_avg
#  endif
# endif  /* AVERAGES */


#endif /* DIAGNOSTICS_EK */




