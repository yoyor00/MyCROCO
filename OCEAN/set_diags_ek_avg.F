#include "cppdefs.h"
#if (defined DIAGNOSTICS_EK && defined AVERAGES)

      subroutine set_diags_ek_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_ek_avg_tile (istr,iend,jstr,jend)
      return
      end

      subroutine set_diags_ek_avg_tile (istr,iend,jstr,jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Thought this algorithm results in somewhat ackwad controlling
! logic it has the advantage that that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(iic-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(iic-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer istr,iend,jstr,jend, i,j, ilc, iflux
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# include "averages.h"
# include "diags_ek.h"

!
# include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntstart  ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_ek_avg
! is a positive number
!


      if (ilc.gt.ntsdiags_ek_avg) then
        if (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.1
     &                    .OR. nwrtdiags_ek_avg.eq.1) then

          if (ZEROTH_TILE) then
            timediags_ek_avg=time
!            MPI_master_only write(*,*) 'started averaging ek',iic,
!     &                            ntsdiags_ek_avg,nwrtdiags_ek_avg
          endif


          do j=JstrR,JendR
            do i=IstrR,IendR
              do k=1,N
                ekHadv_avg(i,j,k)=ekHadv(i,j,k)
                ekHdiff_avg(i,j,k)=ekHdiff(i,j,k)
                ekVadv_avg(i,j,k)=ekVadv(i,j,k)
                ekCor_avg(i,j,k)=ekCor(i,j,k)
                ekPrsgrd_avg(i,j,k)=ekPrsgrd(i,j,k)
                ekHmix_avg(i,j,k)=ekHmix(i,j,k)
                ekVmix_avg(i,j,k)=ekVmix(i,j,k)
                ekvol_avg(i,j,k)=ekvol(i,j,k)
                ekrate_avg(i,j,k)=ekrate(i,j,k)
                ekVmix2_avg(i,j,k)=ekVmix2(i,j,k)
# if defined DIAGNOSTICS_BARO
                ekBaro_avg(i,j,k)=ekBaro(i,j,k)
# endif
# if defined M3FAST
                ekfast_avg(i,j,k)=ekfast(i,j,k)
# endif
              enddo
              ekWind_avg(i,j)=ekWind(i,j)
              ekDrag_avg(i,j)=ekDrag(i,j)
            enddo
          enddo


        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).gt.1) then

          if (ZEROTH_TILE) timediags_ek_avg=timediags_ek_avg+time

            do j=JstrR,JendR
              do i=IstrR,IendR
                do k=1,N
                  ekHadv_avg(i,j,k) =
     &                       ekHadv_avg(i,j,k) +
     &                       ekHadv(i,j,k)
                  ekHdiff_avg(i,j,k) =
     &                       ekHdiff_avg(i,j,k) +
     &                       ekHdiff(i,j,k)
                  ekVadv_avg(i,j,k) =
     &                       ekVadv_avg(i,j,k) +
     &                       ekVadv(i,j,k)
                  ekCor_avg(i,j,k) =
     &                       ekCor_avg(i,j,k) +
     &                       ekCor(i,j,k)
                  ekPrsgrd_avg(i,j,k) =
     &                       ekPrsgrd_avg(i,j,k) +
     &                       ekPrsgrd(i,j,k)
                  ekHmix_avg(i,j,k) =
     &                       ekHmix_avg(i,j,k) +
     &                       ekHmix(i,j,k)
                  ekVmix_avg(i,j,k) =
     &                       ekVmix_avg(i,j,k) +
     &                       ekVmix(i,j,k)
                  ekVmix2_avg(i,j,k) =
     &                       ekVmix2_avg(i,j,k) +
     &                       ekVmix2(i,j,k)
# if defined DIAGNOSTICS_BARO
                  ekBaro_avg(i,j,k) =
     &                      ekBaro_avg(i,j,k) +
     &                      ekBaro(i,j,k)
# endif
# if defined M3FAST
                  ekfast_avg(i,j,k) =
     &                      ekfast_avg(i,j,k) +
     &                      ekfast(i,j,k)
# endif
                  ekrate_avg(i,j,k) =
     &                      ekrate_avg(i,j,k) +
     &                      ekrate(i,j,k)
                  ekvol_avg(i,j,k) =
     &                      ekvol_avg(i,j,k) +
     &                      ekvol(i,j,k)
                enddo

               ekWind_avg(i,j) =
     &                   ekWind_avg(i,j) +
     &                   ekWind(i,j)
               ekDrag_avg(i,j) =
     &                   ekDrag_avg(i,j) +
     &                   ekDrag(i,j)
           enddo
          enddo


        elseif (mod(ilc-ntsdiags_ek_avg,nwrtdiags_ek_avg).eq.0) then

          cff=1./float(nwrtdiags_ek_avg)
          if (ZEROTH_TILE) then
            timediags_ek_avg=cff*(timediags_ek_avg+time)
!            MPI_master_only write(*,*) 'finish aver. ek',iic,
!     &                       ntsdiags_ek_avg,nwrtdiags_ek_avg
          endif


            do j=JstrR,JendR
              do i=IstrR,IendR
                do k=1,N
                  ekHadv_avg(i,j,k) = cff *
     &                       ( ekHadv_avg(i,j,k) +
     &                         ekHadv(i,j,k) )
                  ekHdiff_avg(i,j,k) = cff *
     &                       ( ekHdiff_avg(i,j,k) +
     &                         ekHdiff(i,j,k) )
                  ekVadv_avg(i,j,k) = cff *
     &                       ( ekVadv_avg(i,j,k) +
     &                         ekVadv(i,j,k) )
                  ekCor_avg(i,j,k) = cff *
     &                       ( ekCor_avg(i,j,k) +
     &                         ekCor(i,j,k) )
                  ekPrsgrd_avg(i,j,k) = cff *
     &                       ( ekPrsgrd_avg(i,j,k) +
     &                         ekPrsgrd(i,j,k) )
                  ekHmix_avg(i,j,k) = cff *
     &                       ( ekHmix_avg(i,j,k) +
     &                         ekHmix(i,j,k) )
                  ekVmix_avg(i,j,k) = cff *
     &                       ( ekVmix_avg(i,j,k) +
     &                         ekVmix(i,j,k) )
                  ekVmix2_avg(i,j,k) = cff *
     &                       ( ekVmix2_avg(i,j,k) +
     &                         ekVmix2(i,j,k) )
                  ekrate_avg(i,j,k) = cff *
     &                       ( ekrate_avg(i,j,k) +
     &                         ekrate(i,j,k) )
                  ekvol_avg(i,j,k) = cff *
     &                      ( ekvol_avg(i,j,k) +
     &                        ekvol(i,j,k) )
# if defined DIAGNOSTICS_BARO
                  ekBaro_avg(i,j,k) = cff *
     &                       ( ekBaro_avg(i,j,k) +
     &                         ekBaro(i,j,k) )
# endif
# if defined M3FAST
                  ekfast_avg(i,j,k) = cff *
     &                       ( ekfast_avg(i,j,k) +
     &                         ekfast(i,j,k) )
# endif
                enddo
                ekWind_avg(i,j) = cff *
     &                     ( ekWind_avg(i,j) +
     &                       ekWind(i,j) )
                ekDrag_avg(i,j) = cff *
     &                     ( ekDrag_avg(i,j) +
     &                       ekDrag(i,j) )
              enddo
            enddo
        endif
      endif       !<-- iic.gt.ntsavg



      return
      end
#else /* DIAGNOSTICS_EK && AVERAGES */
      subroutine set_diags_ek_avg_empty
      end
#endif /* DIAGNOSTICS_EK && AVERAGES */
