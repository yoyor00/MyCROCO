MODULE storng_testmpi

#include "cppdefs.h"
#if defined STOGEN

   !!======================================================================
   !!                       ***  MODULE  storng_test ***
   !! Random number generator, used only to test reproducibility
   !!        with different MPI domain decompositions
   !! For that purpose, a different seed is used
   !!        for each grid point of the global grid
   !! Do not use it for real simulations
   !!
   !!=====================================================================

   !!----------------------------------------------------------------------
   !! The module is based on (and includes) the
   !! 64-bit KISS (Keep It Simple Stupid) random number generator
   !! distributed by George Marsaglia :
   !! http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/a85bf5f2a97f5a55
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! kiss64_test      : 64-bit KISS random number generator (period ~ 2^250)
   !! kiss_seed_test   : Define seeds for KISS random number generator
   !! kiss_normal_test : Real random numbers with normal distribution N(0,1)
   !!----------------------------------------------------------------------
   USE stoexternal , only : wp, i4, i8, jpi, jpj, jpiglo, jpjglo, mig, mjg, nmember
   USE storng_kiss

   IMPLICIT NONE
   PRIVATE

   ! Public functions/subroutines
   PUBLIC :: kiss_seed_test, kiss_normal_test

   ! Default/initial seeds (64-bit kiss)
   INTEGER(KIND=i8), SAVE :: x8=1234567890987654321_i8
   INTEGER(KIND=i8), SAVE :: y8=362436362436362436_i8
   INTEGER(KIND=i8), SAVE :: z8=1066149217761810_i8
   INTEGER(KIND=i8), SAVE :: w8=123456123456123456_i8

   ! Parameters to generate real random variates
   REAL(KIND=wp), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0
   REAL(wp), PARAMETER    ::  scaling_uni1_8 = 0.5_wp/(real(huge(x8),wp)+1.0_wp)
   REAL(wp), PARAMETER    ::  scaling_uni2_8 = 1.0_wp/(real(huge(x8),wp)+1.0_wp)

   ! Variables to store 2 Gaussian random numbers with current index (ig)
   INTEGER(KIND=i8), SAVE :: ig=1
   REAL(KIND=wp), SAVE :: gran1, gran2

   ! Store state of the random number generator for each grid point
   INTEGER(KIND=i8), DIMENSION(:,:), ALLOCATABLE, SAVE :: x8_sav
   INTEGER(KIND=i8), DIMENSION(:,:), ALLOCATABLE, SAVE :: y8_sav
   INTEGER(KIND=i8), DIMENSION(:,:), ALLOCATABLE, SAVE :: z8_sav
   INTEGER(KIND=i8), DIMENSION(:,:), ALLOCATABLE, SAVE :: w8_sav

CONTAINS

   FUNCTION kiss64_test(ji,jj)
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION kiss64_test  ***
      !!
      !! ** Purpose :   64-bit KISS random number generator
      !!                (locally for each grid point)
      !!
      !! ** Method  :   combine several random number generators:
      !!                (1) Xorshift (XSH), period 2^64-1,
      !!                (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
      !!                (3) Congruential generator (CNG), period 2^64.
      !!
      !!                overall period:
      !!                (2^250+2^192+2^64-2^186-2^129)/6
      !!                            ~= 2^(247.42) or 10^(74.48)
      !!
      !!                set your own seeds with 'kiss_seed_test'
      ! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ji,jj ! index in grid of local domain
      INTEGER(KIND=i8) :: kiss64_test, t

      ! Use state of local random number generator
      x8 = x8_sav(ji,jj)
      y8 = y8_sav(ji,jj)
      z8 = z8_sav(ji,jj)
      w8 = w8_sav(ji,jj)

      ! Get next random number from KISS64
      t = ISHFT(x8,58) + w8
      IF (s(x8).eq.s(t)) THEN
         w8 = ISHFT(x8,-6) + s(x8)
      ELSE
         w8 = ISHFT(x8,-6) + 1 - s(x8+t)
      ENDIF
      x8 = t + x8
      y8 = m( m( m(y8,13_i8), -17_i8 ), 43_i8 )
      z8 = 6906969069_i8 * z8 + 1234567_i8

      kiss64_test = x8 + y8 + z8

      ! Save state of local random number generator
      x8_sav(ji,jj) = x8
      y8_sav(ji,jj) = y8
      z8_sav(ji,jj) = z8
      w8_sav(ji,jj) = w8

      CONTAINS

         FUNCTION s(k)
            INTEGER(KIND=i8) :: s, k
            s = ISHFT(k,-63)
         END FUNCTION s

         FUNCTION m(k, n)
            INTEGER(KIND=i8) :: m, k, n
            m =  IEOR(k, ISHFT(k, n) )
         END FUNCTION m

   END FUNCTION kiss64_test


   SUBROUTINE kiss_seed_test( )
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_seed_test  ***
      !!
      !! ** Purpose :   Define seeds for 'test_mpi' random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ji, jj, jseed, seedindex
      INTEGER(KIND=i8) :: zseed1_8, zseed2_8, zseed3_8, zseed4_8

      ! Allocate array to save state of random number generators
      ALLOCATE(x8_sav(jpi,jpj))
      ALLOCATE(y8_sav(jpi,jpj))
      ALLOCATE(z8_sav(jpi,jpj))
      ALLOCATE(w8_sav(jpi,jpj))

      ! Reset seed to default seed
      CALL kiss_reset( )
      ! Set unique initial seedindex for each ensemble member
      ! (with enough space inbetween to have one unique seed per grid point)
      seedindex = nmember * jpiglo * jpjglo
      ! Compute initial seed for this ensemble member using standard KISS code
      DO jseed = 0, seedindex
        zseed1_8 = kiss64() ; zseed2_8 = kiss64()
        zseed3_8 = kiss64() ; zseed4_8 = kiss64()
      END DO
      x8_sav(:,:) = zseed1_8
      y8_sav(:,:) = zseed2_8
      z8_sav(:,:) = zseed3_8
      w8_sav(:,:) = zseed4_8

      ! Seed every local random number generators
      ! with random seed, but reproducible for each global grid point
      DO ji = 1,jpi
      DO jj = 1,jpj
        ! Reset seed to the initial seed of this ensemble member
        CALL kiss_seed( x8_sav(ji,jj), y8_sav(ji,jj), z8_sav(ji,jj), w8_sav(ji,jj) )
        ! Set unique seedindex for each point of the global grid [stored in mig(ji), mjg(jj)]
        ! starting from pre-computed initial seed.
        ! For a given global grid point [mig(ji), mjg(jj)] and a given member,
        ! the seed will thus be the same, whatever the MPI domain decomposition.
        seedindex = ( mjg(jj) - 1 ) * jpiglo + mig(ji)
        ! Compute seed using standard KISS code
        DO jseed = 1, seedindex
          zseed1_8 = kiss64() ; zseed2_8 = kiss64()
          zseed3_8 = kiss64() ; zseed4_8 = kiss64()
        END DO
        ! Initialize each random chain with this seed
        x8_sav(ji,jj) = zseed1_8
        y8_sav(ji,jj) = zseed2_8
        z8_sav(ji,jj) = zseed3_8
        w8_sav(ji,jj) = zseed4_8
      ENDDO
      ENDDO

   END SUBROUTINE kiss_seed_test


   FUNCTION kiss_normal_test(ji,jj)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_normal_test  ***
      !!
      !! ** Purpose :   Real random numbers with Gaussian distribution N(0,1)
      !!
      !! ** Method  :   Generate 2 new Gaussian draws (gran1 and gran2)
      !!                from 2 uniform draws on [-1,1] (u1 and u2),
      !!                using the Marsaglia polar method
      !!                (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ji,jj ! index in grid of local domain

      REAL(KIND=wp) :: kiss_normal_test
      REAL(KIND=wp) :: u1, u2, rsq, fac

      IF (ig.EQ.1) THEN
         rsq = two
         DO WHILE ( (rsq.GE.one).OR. (rsq.EQ.zero) )
            u1 = REAL(kiss64_test(ji,jj),wp) * scaling_uni2_8
            u2 = REAL(kiss64_test(ji,jj),wp) * scaling_uni2_8
            rsq = u1*u1 + u2*u2
         ENDDO
         fac = SQRT(-two*LOG(rsq)/rsq)
         gran1 = u1 * fac
         gran2 = u2 * fac
      ENDIF

      ! Output one of the 2 draws
      IF (ig.EQ.1) THEN
         kiss_normal_test = gran1 ; ig = 2
      ELSE
         kiss_normal_test = gran2 ; ig = 1
      ENDIF

   END FUNCTION kiss_normal_test

#endif /* if defined STOGEN */

END MODULE storng_testmpi
