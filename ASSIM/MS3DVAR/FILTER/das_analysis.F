#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine das_analysis(tile)
      implicit none
!
! vertical interpolation from the model
! S-coordinate levels to physical heights
!
      integer tile
#include "param.h"
#include "compute_tile_bounds.h"
      call das_analysis_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine das_analysis_tile (Istr,Iend,Jstr,Jend)
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j,k,itrc
#include "param.h"
#include "das_param.h"
#include "grid.h" 
#include "scalars.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "das_ocean.h"
      real s1d(N),  z1d(NDAS)
      real xs1d(N), xz1d(NDAS)

      external das_spln1d
!
# include "compute_extended_bounds.h"
!
! sea surface height
!
!      do j=JstrR,JendR
!        do i=IstrR,IendR
!          zeta(i,j)=zeta(i,j) + zeta_s(i,j)
!        enddo
!      enddo
!
! interpolate tracers
!
! Surface
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
            t(i,j,N,1,itrc)=t(i,j,N,1,itrc)+t_s(i,j,NDAS,itrc)
          enddo
        enddo
      enddo

      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
             do k=1,N
               s1d(k)=z_r(i,j,k)
             enddo
             do k=1,NDAS
               xz1d(k)=t_s(i,j,k,itrc)
             enddo
             call das_spln1d(NDAS,z_das,xz1d,N,s1d,xs1d)
             do k=1,N-1
               if (abs(s1d(k)) .lt. abs(z_das(1)) ) then
                 t(i,j,k,1,itrc)=t(i,j,k,1,itrc)+xs1d(k)
               endif
             enddo
          enddo
        enddo
      enddo
!
#ifdef EW_PERIODIC
# define IR_RANGE Istr,Iend
# define IU_RANGE Istr,Iend
#else
# define IR_RANGE IstrR,IendR
# define IU_RANGE  Istr,IendR
# ifdef MPI         
      if (WEST_INTER) IstrR=Istr          
      if (EAST_INTER) IendR=Iend         
# endif            
#endif            
#ifdef NS_PERIODIC
# define JR_RANGE Jstr,Jend
# define JV_RANGE Jstr,Jend
#else
# define JR_RANGE JstrR,JendR
# define JV_RANGE  Jstr,JendR
# ifdef MPI
      if (SOUTH_INTER) JstrR=Jstr        
      if (NORTH_INTER) JendR=Jend       
# endif
#endif
!
! interpolate u
!

      do j=JR_RANGE                                                                                          
        do i=IU_RANGE
          u(i,j,N,1)=u(i,j,N,1)+u_s(i,j,NDAS)
          ubar(i,j,1)=ubar(i,j,1) - u_s(i,j,NDAS)
     &                        *(Hz(i-1,j,N)+Hz(i,j,N))                                                  
     &                        /(z_w(i-1,j,0)+z_w(i,j,0))
        enddo
      enddo
      do j=JR_RANGE
        do i=IU_RANGE
           do k=1,N
             s1d(k)=0.5*(z_r(i-1,j,k)+z_r(i,j,k))
           enddo
           do k=1,NDAS
             xz1d(k)=u_s(i,j,k)
           enddo
           call das_spln1d(NDAS,z_das,xz1d,N,s1d,xs1d)
           do k=1,N-1
             if (abs(s1d(k)) .lt. abs(z_das(1)) ) then
               u(i,j,k,1)=u(i,j,k,1) + xs1d(k)
               ubar(i,j,1)=ubar(i,j,1) - xs1d(k)
     &                        *(Hz(i-1,j,k)+Hz(i,j,k))
     &                        /(z_w(i-1,j,0)+z_w(i,j,0))
             endif
           enddo   !k
        enddo
      enddo
!
! interpolate v
!
      do j=JV_RANGE        
        do i=IR_RANGE
          v(i,j,N,1)=v(i,j,N,1)+v_s(i,j,NDAS)
          vbar(i,j,1)=vbar(i,j,1) - v_s(i,j,NDAS)
     &                        *(Hz(i,j-1,N)+Hz(i,j,N))    
     &                        /(z_w(i,j-1,0)+z_w(i,j,0))
        enddo
      enddo
      do j=JV_RANGE
        do i=IR_RANGE
           do k=1,N
             s1d(k)=0.5*(z_r(i,j-1,k)+z_r(i,j,k))
           enddo
           do k=1,NDAS
             xz1d(k)=v_s(i,j,k)
           enddo
           call das_spln1d(NDAS,z_das,xz1d,N,s1d,xs1d)
           do k=1,N-1
             if (abs(s1d(k)) .lt. abs(z_das(1)) ) then
               v(i,j,k,1)=v(i,j,k,1) + xs1d(k)
               vbar(i,j,1)=vbar(i,j,1) - xs1d(k)
     &                        *(Hz(i,j-1,k)+Hz(i,j,k))
     &                        /(z_w(i,j-1,0)+z_w(i,j,0))
             endif
           enddo  !k
        enddo
      enddo    

!
! mask 2d
!
!      do j=JstrR,JendR
!        do i=IstrR,IendR
!          zeta(i,j)=zeta(i,j) * rmask(i,j)
!        enddo
!      enddo
      do j=JR_RANGE
        do i=IU_RANGE
          ubar(i,j,1)=ubar(i,j,1) *  umask(i,j)
        enddo
      enddo
      do j=JV_RANGE
        do i=IR_RANGE
          vbar(i,j,1)=vbar(i,j,1) *  vmask(i,j)
        enddo
      enddo
!
! mask 3d
!
      do k=1,N

        do itrc=1,NT                                     
          do j=JstrR,JendR                              
            do i=IstrR,IendR                           
              t(i,j,k,1,itrc)=t(i,j,k,1,itrc) * rmask(i,j)
            enddo                                     
          enddo                                      
        enddo                                       

        do j=JR_RANGE
          do i=IU_RANGE
            u(i,j,k,1)=u(i,j,k,1) *  umask(i,j)
          enddo
        enddo
        do j=JV_RANGE
          do i=IR_RANGE
            v(i,j,k,1)=v(i,j,k,1) *  vmask(i,j)
          enddo
        enddo
      enddo    !k
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
      subroutine das_analysis_empty
      return
      end
#endif /* SOLVE3D */
