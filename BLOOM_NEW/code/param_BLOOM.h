! $Id: param_BLOOM_tracer.h  $

!!!***********************************************************************************
!  number of variables for module BLOOM depends on CPPkeys (modules bio add variables)
!!!***********************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  number of Bio advected variables without additional modules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_BLOOM_opt2
      INTEGER,PARAMETER :: ntrc_biol=17
      PARAMETER ( ntfix=4 )
#else
#    ifdef key_BLOOM_insed
      ! add PFe & ODU & NdetR & PdetR
      INTEGER,PARAMETER :: ntrc_biol=18
#    else
      INTEGER,PARAMETER :: ntrc_biol=14
#    endif
      PARAMETER ( ntfix=3 )
#endif


#if defined DIAGNOSTICS_BIO
        INTEGER,PARAMETER :: NumBLMdiag3d=34
        INTEGER,PARAMETER :: NumBLMdiag2d=21
        INTEGER,PARAMETER :: NumBLMdiag1d=0
#endif



