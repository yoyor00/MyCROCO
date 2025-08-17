!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if (defined DIAGNOSTICS_UV && defined AVERAGES)

      subroutine set_diagsM_avg(tile)

!      USE param

      implicit none
      integer tile
# include "param.h"
# ifdef SOLVE3D
#  include "work.h"
#  include "ncscrum.h"
# endif
# include "compute_tile_bounds.h"

      call set_diagsM_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end

      subroutine set_diagsM_avg_tile(Istr,Iend,Jstr,Jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Although this algorithm results in somewhat awkward controlling
! logic it has the advantage that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(ilc-1,nwrtdiaM_avg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(ilc-1,nwrtdiaM_avg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j, iflux, k, ilc
      real    cff, cff1
      real    cffux, cffuy, cffvx, cffvy, cffM
      real    cffwx, cffwy, cffuz, cffvz

# include "param.h"
# include "work.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "averages.h"
# include "diagnostics.h"
!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiaM_avg
! is a positive number
!

      if (ilc.gt.ntsdiaM_avg) then

        if (mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).eq.1 .OR. nwrtdiaM_avg.eq.1) then
          cff =1.0
          cff1=0.0
          if (ZEROTH_TILE) then
            timediaM_avg=time
!            MPI_master_only write(*,*) 'started averaging M',
!     &                          iic,ntsdiaM_avg,nwrtdiaM_avg
          endif
        elseif (mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).gt.1) then
          cff =1.0
          cff1=1.0
          if (ZEROTH_TILE) then
           timediaM_avg=timediaM_avg+time
          endif
        elseif (mod(ilc-ntsdiaM_avg,nwrtdiaM_avg).eq.0) then
          cff =1.0/float(nwrtdiaM_avg)
          cff1=1.0
          if (ZEROTH_TILE) then
            timediaM_avg=cff*(timediaM_avg+time)
!            MPI_master_only write(*,*) 'finished averaging M',
!     &                    iic,ntsdiaM_avg,nwrtdiaM_avg
          endif
        endif

        do iflux=1,2
          do k=1,N
            do j=JstrR,JendR
              do i=IstrR,IendR

                MXadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MXadv_avg(i,j,k,iflux) +
     &                            MXadv(i,j,k,iflux) )
                MYadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MYadv_avg(i,j,k,iflux) +
     &                            MYadv(i,j,k,iflux) )
                MVadv_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVadv_avg(i,j,k,iflux) +
     &                            MVadv(i,j,k,iflux) )
#  if defined DIAGNOSTICS_UV_VADV
                MVadvz_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVadvz_avg(i,j,k,iflux) +
     &                            MVadvz(i,j,k,iflux) )
#  endif
                MCor_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MCor_avg(i,j,k,iflux) +
     &                            MCor(i,j,k,iflux) )
                MPrsgrd_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MPrsgrd_avg(i,j,k,iflux) +
     &                            MPrsgrd(i,j,k,iflux) )
                MHmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHmix_avg(i,j,k,iflux) +
     &                            MHmix(i,j,k,iflux,nstp) )
                MHdiff_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MHdiff_avg(i,j,k,iflux) +
     &                            MHdiff(i,j,k,iflux) )
                MVmix_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix_avg(i,j,k,iflux) +
     &                            MVmix(i,j,k,iflux) )
                MVmix2_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MVmix2_avg(i,j,k,iflux) +
     &                            MVmix2(i,j,k,iflux) )
                Mrate_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mrate_avg(i,j,k,iflux) +
     &                            Mrate(i,j,k,iflux) )
#  if defined DIAGNOSTICS_BARO
                MBaro_avg(i,j,k,iflux) = cff *
     &                     ( cff1*MBaro_avg(i,j,k,iflux) +
     &                            MBaro(i,j,k,iflux) )
#  endif
#  if defined M3FAST
                Mfast_avg(i,j,k,iflux) = cff *
     &                     ( cff1*Mfast_avg(i,j,k,iflux) +
     &                            Mfast(i,j,k,iflux) )
#  endif

# if defined MRL_WCI
                Mvf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mvf_avg(i,j,k,iflux) +
     &                            Mvf(i,j,k,iflux) )
                Mbrk_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mbrk_avg(i,j,k,iflux) +
     &                            Mbrk(i,j,k,iflux) )
                MStCo_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MStCo_avg(i,j,k,iflux) +
     &                            MStCo(i,j,k,iflux) )
                MVvf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MVvf_avg(i,j,k,iflux) +
     &                            MVvf(i,j,k,iflux) )
                MPrscrt_avg(i,j,k,iflux)= cff *
     &                     ( cff1*MPrscrt_avg(i,j,k,iflux) +
     &                            MPrscrt(i,j,k,iflux) )
                Msbk_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Msbk_avg(i,j,k,iflux) +
     &                            Msbk(i,j,k,iflux) )
                Mbwf_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mbwf_avg(i,j,k,iflux) +
     &                            Mbwf(i,j,k,iflux) )
                Mfrc_avg(i,j,k,iflux)= cff *
     &                     ( cff1*Mfrc_avg(i,j,k,iflux) +
     &                            Mfrc(i,j,k,iflux) )
# endif
              enddo
            enddo
          enddo
        enddo

      endif

      return
      end
#else /* DIAGNOSTICS_UV && AVERAGES */
      subroutine set_diagsM_avg_empty
      end
#endif /* DIAGNOSTICS_UV && AVERAGES */
