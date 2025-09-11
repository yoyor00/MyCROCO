!
! Create subroutines:
!
! - `exchange_2d_tile(...)`
! - `exchange_2d_1pts_tile(...)`
! - `exchange_2d_3pts_tile(...)`
! - `exchange_2d_lh_tile(...)`
!
! Duplication of code is avoided by self-including this file
! which can be found at the end.
!
! Note that these subroutines do not really exist, but are redefined
! with a precompiler in `exchange.F`
!

#include "latency_hiding_2d.h"

#if (!defined MP_3PTS) && (!defined MP_1PTS) && (!defined MP_LH_PTS)
! Use per default 2 points
      subroutine exchange_2d_tile (Istr,Iend,Jstr,Jend, A)
#elif defined MP_1PTS
      subroutine exchange_2d_1pts_tile (Istr,Iend,Jstr,Jend, A)
#elif defined MP_3PTS
      subroutine exchange_2d_3pts_tile (Istr,Iend,Jstr,Jend, A)
#elif defined MP_LH_PTS
      ! Latency hiding points
      subroutine exchange_2d_lh_tile (Istr,Iend,Jstr,Jend, A)
#else
#  error "Internal error"
#endif

!
! Set periodic boundary conditions (if any) for a two-dimensional
! field A of ZETA-, U-, V- or PSI-type. This file is designed to
! generate four different subroutines, by redefining (via CPP) the
! name of the subroutine exchange_2d_tile above and the starting
! indices ISTART = [Istr for U-,PSI-type; IstrR for V-,ZETA-type] 
! and JSTART = [Jstr for V-,PSI-type; JstrR for U-,ZETA-type]
! below. See also mounting file exchange.F  
!
      implicit none
#include "param.h"
#include "scalars.h"
      integer Npts,ipts,jpts
#if (!defined MP_1PTS) && (!defined MP_3PTS) && (!defined MP_LH_PTS)
      ! Use per default 2 points
      parameter (Npts=2)
#elif defined MP_1PTS
      parameter (Npts=1)
#elif defined MP_3PTS
      parameter (Npts=3)
#elif defined MP_LH_PTS
      ! For latency hiding we communicate more layers than necessary
      ! We do at least 2 layers
      parameter (Npts=2+MPI_OVERLAPPING_SCHWARZ_2D_ADD_LAYERS)
#else
#  error "Internal error"
#endif
      real A(GLOBAL_2D_ARRAY)
      integer Istr,Iend,Jstr,Jend, i,j
!
#include "compute_auxiliary_bounds.h"
!

!$acc kernels if(compute_on_device) default(present)  
#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define J_RANGE Jstr,Jend
# else
#  define J_RANGE JSTART,JendR
# endif
# ifdef MPI
      if (NP_XI.eq.1) then
# endif
        if (WESTERN_EDGE) then		    
          do j=J_RANGE
            do ipts=1,Npts
              A(Lm+ipts,j)=A(ipts,j)
	      enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_RANGE
            do ipts=1,Npts
              A(ipts-Npts,j)=A(Lm+ipts-Npts,j)
            enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef J_RANGE
#endif

#ifdef NS_PERIODIC
# ifdef EW_PERIODIC
#  define I_RANGE Istr,Iend
# else
#  define I_RANGE ISTART,IendR
# endif
# ifdef MPI
      if (NP_ETA.eq.1) then
# endif
        if (SOUTHERN_EDGE) then
          do i=I_RANGE
	      do jpts=1,Npts
              A(i,Mm+jpts)=A(i,jpts)
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_RANGE
	      do jpts=1,Npts
              A(i,jpts-Npts)=A(i,Mm+jpts-Npts)
	      enddo
          enddo
        endif
# ifdef MPI
      endif
# endif
# undef I_RANGE
#endif

#if defined EW_PERIODIC && defined NS_PERIODIC
# ifdef MPI
      if (NP_XI.eq.1 .and. NP_ETA.eq.1) then
# endif
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(Lm+ipts,Mm+jpts)=A(ipts,jpts)
	      enddo
	    enddo
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          do jpts=1,Npts
            do ipts=1,Npts
                  A(ipts-Npts,Mm+jpts)=A(Lm+ipts-Npts,jpts)
            enddo
          enddo
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(Lm+ipts,jpts-Npts)=A(ipts,Mm+jpts-Npts)
	      enddo
	    enddo
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
	    do jpts=1,Npts
	      do ipts=1,Npts
	        A(ipts-Npts,jpts-Npts)=A(Lm+ipts-Npts,Mm+jpts-Npts)
	      enddo
	    enddo
        endif
# ifdef MPI
      endif
# endif
#endif
!$acc end kernels

#ifdef MPI
#  if (!defined MP_3PTS) && (!defined MP_1PTS) && (!defined MP_LH_PTS)
      call MessPass2D_tile (Istr,Iend,Jstr,Jend,  A)
#  elif defined MP_1PTS
      call MessPass2D_1pts_tile (Istr,Iend,Jstr,Jend,  A)
#  elif defined MP_3PTS
      call MessPass2D_3pts_tile (Istr,Iend,Jstr,Jend,  A)
#  elif defined MP_LH_PTS
      call MessPass2D_lh_tile (Istr,Iend,Jstr,Jend,  A)
#  else
#    error "Internal error"
#  endif

#  ifdef  BAND_DEBUG          
      chkbandname='none'
#  endif

#endif

#if defined OPENMP && defined OPENACC
      if (.not.SOUTHERN_EDGE) then
!$acc update host(A(:,Jstr:Jstr+Npts-1))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update host(A(:,Jend-Npts+1:Jend))
      endif
C$OMP BARRIER
      if (.not.SOUTHERN_EDGE) then
!$acc update device(A(:,Jstr-Npts:Jstr-1))
      endif
      if (.not.NORTHERN_EDGE) then
!$acc update device(A(:,Jend+1:Jend+Npts))
      endif
#endif
           
      return
      end

!
! Create all variants of the exchange for tiles
! by including this file by itself.
!

#if (!defined MP_1PTS) && (!defined MP_3PTS) && (!defined MP_LH_PTS)

#  define MP_3PTS
#  include "exchange_2d_tile.h"
#  undef MP_3PTS

#  define MP_1PTS
#  include "exchange_2d_tile.h"
#  undef MP_1PTS

#  ifdef MPI_OVERLAPPING_SCHWARZ_2D
#    define MP_LH_PTS
#    include "exchange_2d_tile.h"
#    undef MP_LH_PTS
#  endif

#endif
