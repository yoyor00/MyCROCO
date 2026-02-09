#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine das_s2z_3d(tile)
      implicit none
!
! vertical interpolation from the model
! S-coordinate levels to physical heights
!
      integer tile
#include "param.h"
#include "compute_tile_bounds.h"
      call das_s2z_3d_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_s2z_3d_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc
#include "param.h"
#include "das_param.h"
#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "das_ocean.h"
       real s1d(N),  z1d(NDAS)
       real xs1d(N), xz1d(NDAS)
       integer KZ
       external das_spln1d
!
# include "compute_extended_bounds.h"
!
! sea surface height
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta_das(i,j)=zeta(i,j,1)
        enddo
      enddo
!
! interpolate tracers
!
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
             do k=1,N
               s1d(k)=z_r(i,j,k)
               xs1d(k)=t(i,j,k,1,itrc)
             enddo
             KZ=NDAS-nzr_das(i,j)+1
             do k=nzr_das(i,j),NDAS
               z1d(k-nzr_das(i,j)+1)=z_das(k) 
             enddo
             call das_spln1d(N,s1d,xs1d,KZ,z1d,xz1d)
             do k=nzr_das(i,j),NDAS
               t_das(i,j,k,itrc)=xz1d(k-nzr_das(i,j)+1)
             enddo
          enddo
        enddo
      enddo
!      do itrc=1,NT
!        do j=JstrR,JendR
!          do i=IstrR,IendR
!            t_das(i,j,NDAS,itrc)=t(i,j,N,itrc)
!          enddo
!        enddo
!      enddo
# ifdef EW_PERIODIC
#  define IR_RANGE Istr,Iend
#  define IU_RANGE Istr,Iend
# else
#  define IR_RANGE IstrR,IendR
#  define IU_RANGE  Istr,IendR
#  ifdef MPI             
      if (WEST_INTER) IstrR=Istr          ! computational boundary
      if (EAST_INTER) IendR=Iend          ! are filled during
#  endif                
# endif                
# ifdef NS_PERIODIC
#  define JR_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JR_RANGE JstrR,JendR
#  define JV_RANGE  Jstr,JendR
#  ifdef MPI
      if (SOUTH_INTER) JstrR=Jstr         ! same as above.
      if (NORTH_INTER) JendR=Jend         !
#  endif
# endif
!
! interpolate u
!
      do j=JR_RANGE
        do i=IU_RANGE
           do k=1,N
             s1d(k)=0.5*(z_r(i-1,j,k)+z_r(i,j,k))
             xs1d(k)=u(i,j,k,1)
           enddo
           KZ=NDAS-nzu_das(i,j)+1
           do k=nzu_das(i,j),NDAS
             z1d(k-nzu_das(i,j)+1)=z_das(k) 
           enddo
           call das_spln1d(N,s1d,xs1d,KZ,z1d,xz1d)
           do k=nzu_das(i,j),NDAS
             u_das(i,j,k)=xz1d(k-nzu_das(i,j)+1)
           enddo
        enddo
      enddo
      do j=JR_RANGE
        do i=IU_RANGE
          u_das(i,j,NDAS)=u(i,j,N,1)
        enddo
      enddo
!
! interpolate v
!
      do j=JV_RANGE
        do i=IR_RANGE
           do k=1,N
             s1d(k)=0.5*(z_r(i,j-1,k)+z_r(i,j,k))
             xs1d(k)=v(i,j,k,1)
           enddo
           KZ=NDAS-nzv_das(i,j)+1
           do k=nzv_das(i,j),NDAS
             z1d(k-nzv_das(i,j)+1)=z_das(k)
           enddo
           call das_spln1d(N,s1d,xs1d,KZ,z1d,xz1d)
           do k=nzv_das(i,j),NDAS
             v_das(i,j,k)=xz1d(k-nzv_das(i,j)+1)
           enddo
        enddo
      enddo
      do j=JR_RANGE
        do i=IU_RANGE
          v_das(i,j,NDAS)=v(i,j,N,1)
        enddo
      enddo
!
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE
!
! Exchange periodic boundaries, if so prescribed.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call das_exchange_u3d_tile (Istr,Iend,Jstr,Jend, u_das)
      call das_exchange_v3d_tile (Istr,Iend,Jstr,Jend, v_das)
# endif
      return
      end
#else
      subroutine das_s2z_empty
      return
      end
#endif /* SOLVE3D */
