#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine das_z2s(tile)
      implicit none
!
! vertical interpolation from the model
! S-coordinate levels to physical heights
!
      integer tile
#include "param.h"
#include "compute_tile_bounds.h"
      call das_z2s_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_z2s_tile (Istr,Iend,Jstr,Jend)
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
       integer KZ,KS
       external das_spln1d
!
# include "compute_extended_bounds.h"
!
! sea surface height
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,1)=zeta_das(i,j)
        enddo
      enddo
!
! interpolate tracers
!
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
             KS=N-nsr_das(i,j)+1
             do k=nsr_das(i,j),N
               s1d(k-nsr_das(i,j)+1)=z_r(i,j,k)
             enddo
             KZ=NDAS-nzr_das(i,j)+1
             do k=nzr_das(i,j),NDAS
               z1d(k-nzr_das(i,j)+1)=z_das(k) 
               xz1d(k-nzr_das(i,j)+1)=t_das(i,j,k,itrc)
             enddo
             call das_spln1d(KZ,z1d,xz1d,KS,s1d,xs1d)
             do k=nsr_das(i,j),N
               t(i,j,k,1,itrc)=xs1d(k-nsr_das(i,j)+1)
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
#ifdef EW_PERIODIC
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE  Istr,IendR
#endif
                                    ! see also below...
#ifdef NS_PERIODIC
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE  Jstr,JendR
#endif
!
! interpolate u
!
      do j=JR_RANGE
        do i=IU_RANGE
           KS=N-nsu_das(i,j)+1
           do k=nsu_das(i,j),N
             s1d(k-nsu_das(i,j)+1)=0.5*(z_r(i-1,j,k)+z_r(i,j,k))
           enddo
           KZ=NDAS-nzu_das(i,j)+1
           do k=nzu_das(i,j),NDAS
             z1d(k-nzu_das(i,j)+1)=z_das(k) 
             xz1d(k-nzu_das(i,j)+1)=u_das(i,j,k)
           enddo
           call das_spln1d(KZ,z1d,xz1d,KS,s1d,xs1d)
           do k=nsu_das(i,j),N
             u(i,j,k,1)=xs1d(k-nsu_das(i,j)+1)
           enddo
        enddo
      enddo
!
! interpolate v
!
      do j=JV_RANGE
        do i=IR_RANGE
           KS=N-nsv_das(i,j)+1
           do k=nsv_das(i,j),N
             s1d(k-nsv_das(i,j)+1)=0.5*(z_r(i,j-1,k)+z_r(i,j,k))
           enddo
           KZ=NDAS-nzv_das(i,j)+1
           do k=nzv_das(i,j),NDAS
             z1d(k-nzv_das(i,j)+1)=z_das(k)
             xz1d(k-nzv_das(i,j)+1)=v_das(i,j,k)
           enddo
           call das_spln1d(KZ,z1d,xz1d,KS,s1d,xs1d)
           do k=nsv_das(i,j),N
             v(i,j,k,1)=xs1d(k-nsv_das(i,j)+1)
           enddo
        enddo
      enddo
!
#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call das_exchange_u3d_tile (Istr,Iend,Jstr,Jend, u)
      call das_exchange_v3d_tile (Istr,Iend,Jstr,Jend, v)
# endif
      return
      end
#else
      subroutine das_z2s_empty
      return
      end
#endif /* SOLVE3D */
