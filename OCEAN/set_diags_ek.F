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

      do j=JstrR,JendR+1
        do i=IstrR,IendR+1
          do k=1,N

            workr(i,j,k)=om_u(i,j)* on_u(i,j)
     &           * 0.5*(u(i,j,k,nstp)+u(i,j,k,nnew))


          enddo
        enddo
      enddo

      do j=JstrR,JendR+1
        do i=IstrR,IendR+1
          do k=1,N

            work(i,j,k)=om_v(i,j)*on_v(i,j)
     &           * 0.5*(v(i,j,k,nstp)+v(i,j,k,nnew))


          enddo
        enddo
      enddo

        do j=JstrR,JendR
          do i=IstrR,IendR
          do k=1,N

            ekHadv(i,j,k) = 0.5 * (
     &       (MXadv(i  ,j  ,k,1)+ MYadv(i  ,j  ,k,1))*workr(i  ,j  ,k)
     &      +(MXadv(i+1,j  ,k,1)+ MYadv(i+1,j  ,k,1))*workr(i+1,j  ,k)
     &      +(MXadv(i  ,j  ,k,2)+ MYadv(i  ,j  ,k,2))*work(i  ,j  ,k)
     &      +(MXadv(i  ,j+1,k,2)+ MYadv(i  ,j+1,k,2))*work(i  ,j+1,k)
     &                      )

            ekHdiff(i,j,k) = 0.5 * (
     &                         MHdiff(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MHdiff(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MHdiff(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MHdiff(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )


            ekVadv(i,j,k) = 0.5 * (
     &                         MVadv(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MVadv(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MVadv(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MVadv(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )

            ekCor(i,j,k) = 0.5 * (
     &                         MCor(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MCor(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MCor(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MCor(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )

            ekPrsgrd(i,j,k) = 0.5 * (
     &                        MPrsgrd(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                      + MPrsgrd(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                      + MPrsgrd(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                      + MPrsgrd(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )

            ekHmix(i,j,k) = 0.5 * (
     &                    MHmix(i  ,j  ,k,1,3-nstp) * workr(i  ,j  ,k)
     &                  + MHmix(i+1,j  ,k,1,3-nstp) * workr(i+1,j  ,k)
     &                  + MHmix(i  ,j  ,k,2,3-nstp) * work(i  ,j  ,k)
     &                  + MHmix(i  ,j+1,k,2,3-nstp) * work(i  ,j+1,k)
     &                      )

            ekVmix(i,j,k) = 0.5 * (
     &                         MVmix(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MVmix(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MVmix(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MVmix(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )

!            ekrate(i,j,k) = 0.5 * (
!     &                         Mrate(i  ,j  ,k,1) * workr(i  ,j  ,k)
!     &                       + Mrate(i+1,j  ,k,1) * workr(i+1,j  ,k)
!     &                       + Mrate(i  ,j  ,k,2) * work(i  ,j  ,k)
!     &                       + Mrate(i  ,j+1,k,2) * work(i  ,j+1,k)
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
     &                         MVmix2(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MVmix2(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MVmix2(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MVmix2(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )

# if defined DIAGNOSTICS_BARO
            ekBaro(i,j,k) = 0.5 * (
     &                         MBaro(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + MBaro(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + MBaro(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + MBaro(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )
# endif
# if defined M3FAST
            ekfast(i,j,k) = 0.5 * (
     &                         Mfast(i  ,j  ,k,1) * workr(i  ,j  ,k)
     &                       + Mfast(i+1,j  ,k,1) * workr(i+1,j  ,k)
     &                       + Mfast(i  ,j  ,k,2) * work(i  ,j  ,k)
     &                       + Mfast(i  ,j+1,k,2) * work(i  ,j+1,k)
     &                      )
# endif

          enddo
        enddo
      enddo

! ------- ------------- ------ ------ - -----

        do j=Jstr,JendR
          do i=Istr,IendR

          ekWind(i,j) = 0.5 * (ekwrkWind(i,j,1) + ekwrkWind(i+1,j,1))
     &                + 0.5 * (ekwrkWind(i,j,2) + ekwrkWind(i,j+1,2))


          ekDrag(i,j) = 0.5 * (ekwrkDrag(i,j,1) + ekwrkDrag(i+1,j,1))
     &                + 0.5 * (ekwrkDrag(i,j,2) + ekwrkDrag(i,j+1,2))


          enddo
         enddo



! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_EK*/
      subroutine set_diags_ek_empty
      end
#endif /* DIAGNOSTICS_EK */










