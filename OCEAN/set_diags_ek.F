#include "cppdefs.h"
#if defined DIAGNOSTICS_EK


      subroutine set_diags_ek (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_ek_tile (istr,iend,jstr,jend,tile)
      return
      end

      subroutine set_diags_ek_tile (istr,iend,jstr,jend,tile)
!
!
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
# include "scalars.h"
# include "ncscrum.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "forces.h"
# include "work.h"
# include "diagnostics.h"
# include "diags_ek.h"
# include "compute_auxiliary_bounds.h"



!--------------------------------------------------

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE IstrR,IendR+1
#  define IV_RANGE IstrR,IendR+1
# endif
!
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR+1
#  define JV_RANGE JstrR,JendR+1
# endif

!--------------------------------------------------

      do k=1,N
        do j=JstrR,JendR+1
          do i=IstrR,IendR+1

            work2d(i,j)=om_u(i,j)* on_u(i,j)
     &           * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))

            work2d2(i,j)=om_v(i,j)*on_v(i,j)
     &           * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))

          enddo
        enddo

        do j=JstrR,JendR
          do i=IstrR,IendR

            ekHadv(i,j,k) = 0.5 * (
     &       (MXadv(i  ,j  ,k,1)+ MYadv(i  ,j  ,k,1))*work2d(i  ,j  )
     &      +(MXadv(i+1,j  ,k,1)+ MYadv(i+1,j  ,k,1))*work2d(i+1,j  )
     &      +(MXadv(i  ,j  ,k,2)+ MYadv(i  ,j  ,k,2))*work2d2(i  ,j  )
     &      +(MXadv(i  ,j+1,k,2)+ MYadv(i  ,j+1,k,2))*work2d2(i  ,j+1)
     &                      )

            ekHdiff(i,j,k) = 0.5 * (
     &                         MHdiff(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MHdiff(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MHdiff(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MHdiff(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )


            ekVadv(i,j,k) = 0.5 * (
     &                         MVadv(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MVadv(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MVadv(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MVadv(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )

#  if defined DIAGNOSTICS_UV_VADV
            ekVadvz(i,j,k) = 0.5 * (
     &                         MVadvz(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MVadvz(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MVadvz(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MVadvz(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )
#   endif
            ekCor(i,j,k) = 0.5 * (
     &                         MCor(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MCor(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MCor(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MCor(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )

            ekPrsgrd(i,j,k) = 0.5 * (
     &                        MPrsgrd(i  ,j  ,k,1) * work2d(i  ,j  )
     &                      + MPrsgrd(i+1,j  ,k,1) * work2d(i+1,j  )
     &                      + MPrsgrd(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                      + MPrsgrd(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )

            ekHmix(i,j,k) = 0.5 * (
     &                    MHmix(i  ,j  ,k,1,3-nstp) * work2d(i  ,j  )
     &                  + MHmix(i+1,j  ,k,1,3-nstp) * work2d(i+1,j  )
     &                  + MHmix(i  ,j  ,k,2,3-nstp) * work2d2(i  ,j  )
     &                  + MHmix(i  ,j+1,k,2,3-nstp) * work2d2(i  ,j+1)
     &                      )

            ekVmix(i,j,k) = 0.5 * (
     &                         MVmix(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MVmix(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MVmix(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MVmix(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )

!            ekrate(i,j,k) = 0.5 * (
!     &                         Mrate(i  ,j  ,k,1) * work2d(i  ,j  )
!     &                       + Mrate(i+1,j  ,k,1) * work2d(i+1,j  )
!     &                       + Mrate(i  ,j  ,k,2) * work2d2(i  ,j  )
!     &                       + Mrate(i  ,j+1,k,2) * work2d2(i  ,j+1)
!     &                      )

            ekrate(i,j,k) = 0.5 * (
     &                         ekwrkrate(i  ,j  ,k,1)
     &                       + ekwrkrate(i+1,j  ,k,1)
     &                       + ekwrkrate(i  ,j  ,k,2)
     &                       + ekwrkrate(i  ,j+1,k,2)
     &                      )

            ekvol(i,j,k) = 0.5 * (
     &                         ekwrkvol(i  ,j  ,k,1)
     &                       + ekwrkvol(i+1,j  ,k,1)
     &                       + ekwrkvol(i  ,j  ,k,2)
     &                       + ekwrkvol(i  ,j+1,k,2)
     &                      )

            ekVmix2(i,j,k) = 0.5 * (
     &                         MVmix2(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MVmix2(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MVmix2(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MVmix2(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )

# if defined DIAGNOSTICS_BARO
            ekBaro(i,j,k) = 0.5 * (
     &                         MBaro(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + MBaro(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + MBaro(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + MBaro(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )
# endif
# if defined M3FAST
            ekfast(i,j,k) = 0.5 * (
     &                         Mfast(i  ,j  ,k,1) * work2d(i  ,j  )
     &                       + Mfast(i+1,j  ,k,1) * work2d(i+1,j  )
     &                       + Mfast(i  ,j  ,k,2) * work2d2(i  ,j  )
     &                       + Mfast(i  ,j+1,k,2) * work2d2(i  ,j+1)
     &                      )
# endif

# ifdef BODYFORCE
          if (k.eq.1) then
            ekWind(i,j) = 0.5 * (ekwrkWind(i,j,1) + ekwrkWind(i+1,j,1))
     &                  + 0.5 * (ekwrkWind(i,j,2) + ekwrkWind(i,j+1,2))
            ekDrag(i,j) = 0.5 * (ekwrkDrag(i,j,1) + ekwrkDrag(i+1,j,1))
     &                  + 0.5 * (ekwrkDrag(i,j,2) + ekwrkDrag(i,j+1,2))
          endif
# else
          if (k.eq.1) then
            ekDrag(i,j) = 0.5 * (
     &                         bustr(i  ,j  ) * work2d(i  ,j  )
     &                       + bustr(i+1,j  ) * work2d(i+1,j  )
     &                       + bvstr(i  ,j  ) * work2d2(i  ,j  )
     &                       + bvstr(i  ,j+1) * work2d2(i  ,j+1)
     &                      )
          endif
          if (k.eq.N) then
#  ifdef STRESS_AT_RHO_POINTS
          ekWind(i,j) = 0.5 * (
     &               0.5*(sustr(i,j)+sustr(i-1,j)) * work2d(i  ,j  )
     &             + 0.5*(sustr(i,j)+sustr(i+1,j)) * work2d(i+1,j  )
     &             + 0.5*(svstr(i,j)+svstr(i,j-1)) * work2d2(i  ,j  )
     &             + 0.5*(svstr(i,j)+svstr(i,j+1)) * work2d2(i  ,j+1)
     &                      )
#  else
          ekWind(i,j) = 0.5 * (
     &                         sustr(i  ,j  ) * work2d(i  ,j  )
     &                       + sustr(i+1,j  ) * work2d(i+1,j  )
     &                       + svstr(i  ,j  ) * work2d2(i  ,j )
     &                       + svstr(i  ,j+1) * work2d2(i  ,j+1)
     &                      )
#  endif
          endif
# endif
          enddo
         enddo
       enddo


! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_EK*/
      subroutine set_diags_ek_empty
      end
#endif /* DIAGNOSTICS_EK */










