!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#ifdef LERAY_TURB
!      
!Computation of the filter weights to apply to velocity components in
!advection operators (Petersen, 2008)
      subroutine calc_smoothing_weights()
      implicit none    
# include "param.h"
# include "leray_smoothing.h"

# ifdef LERAY_FILTER_9PTS     
      !j-4
      filter_weights(-4,-4) = evs*evs
      filter_weights(-3,-4) = evs*dvs
      filter_weights(-2,-4) = evs*cvs
      filter_weights(-1,-4) = evs*bvs
      filter_weights( 0,-4) = evs
      filter_weights( 1,-4) = evs*bvs
      filter_weights( 2,-4) = evs*cvs
      filter_weights( 3,-4) = evs*dvs
      filter_weights( 4,-4) = evs*evs
      !j-3
      filter_weights(-4,-3) = dvs*evs
      filter_weights( 4,-3) = dvs*evs
      !j-2
      filter_weights(-4,-2) = cvs*evs
      filter_weights( 4,-2) = cvs*evs
      !j-1
      filter_weights(-4,-1) = bvs*evs
      filter_weights( 4,-1) = bvs*evs
      !j
      filter_weights(-4, 0) = evs
      filter_weights( 4, 0) = evs
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
      !j-3
      filter_weights(-3,-3) = dvs*dvs
      filter_weights(-2,-3) = dvs*cvs
      filter_weights(-1,-3) = dvs*bvs
      filter_weights( 0,-3) = dvs
      filter_weights( 1,-3) = dvs*bvs
      filter_weights( 2,-3) = dvs*cvs
      filter_weights( 3,-3) = dvs*dvs
      !j-2
      filter_weights(-3,-2) = cvs*dvs
      filter_weights( 3,-2) = cvs*dvs
      !j-1
      filter_weights(-3,-1) = bvs*dvs
      filter_weights( 3,-1) = bvs*dvs
      !j
      filter_weights(-3, 0) = dvs
      filter_weights( 3, 0) = dvs
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
      !j-2
      filter_weights(-2,-2) = cvs*cvs
      filter_weights(-1,-2) = cvs*bvs
      filter_weights( 0,-2) = cvs
      filter_weights( 1,-2) = cvs*bvs
      filter_weights( 2,-2) = cvs*cvs
      !j-1
      filter_weights(-2,-1) = bvs*cvs
      filter_weights( 2,-1) = bvs*cvs
      !j
      filter_weights(-2, 0) = cvs
      filter_weights( 2, 0) = cvs
# endif      
      !j-1
      filter_weights(-1,-1) = bvs*bvs
      filter_weights( 0,-1) = bvs
      filter_weights( 1,-1) = bvs*bvs
      !j
      filter_weights(-1, 0) = bvs
      filter_weights( 0, 0) = avs
      filter_weights( 1, 0) = bvs
      !j+1
      filter_weights(:, 1) = filter_weights(:, -1)
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS      
      !j+2
      filter_weights(:, 2) = filter_weights(:, -2)
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS      
      !j+3
      filter_weights(:, 3) = filter_weights(:, -3)
# endif
# if defined LERAY_FILTER_9PTS      
      !j+4
      filter_weights(:, 4) = filter_weights(:, -4)
# endif      

      inv_weight_sum3 = 1./sum( filter_weights(-1:1,-1:1) )
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS      
      inv_weight_sum5 = 1./sum( filter_weights(-2:2,-2:2) )
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS      
      inv_weight_sum7 = 1./sum( filter_weights(-3:3,-3:3) )
# endif
# if defined LERAY_FILTER_9PTS      
      inv_weight_sum9 = 1./sum( filter_weights(-4:4,-4:4) )
# endif      

      return
      end subroutine
!
!Count the number of available points around each point to determine the
!size of the filter window:
!        - no available points (domain boundary or mask points)
      subroutine compute_fwidth_array(tile)
      implicit none
      integer  tile, trd, omp_get_thread_num
# include "param.h"      
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call compute_fwidth_array_tile (istr,iend,jstr,jend)
      return
      end

      subroutine compute_fwidth_array_tile(istr,iend,jstr,jend)
      implicit none
# include "param.h"
# include "ocean3d.h"
# include "leray_smoothing.h" 
# include "scalars.h"
# include "grid.h"      

      integer istr,iend,jstr,jend,i,j,k
      integer imin,imax,jmin,jmax
# include "compute_auxiliary_bounds.h"
!compute imin,jmin,imax,jmax according to periodic boundaries, INTER or .not.INTER
# ifdef NS_PERIODIC
      jmin=Jstr-filter_width/2
      jmax=Jend+filter_width/2
# else
#  ifdef MPI
      if (SOUTH_INTER) then
        jmin=Jstr-filter_width/2
      else
        jmin=Jstr
      endif
      if (NORTH_INTER) then
        jmax=Jend+filter_width/2
      else
        jmax=Jend
      endif
#  else
      jmin=Jstr
      jmax=Jend
#  endif
# endif      
# ifdef EW_PERIODIC
      imin=IstrU-filter_width/2
      imax=Iend+filter_width/2
# else
#  ifdef MPI
      if (WEST_INTER) then
        imin=IstrU-filter_width/2
      else
        imin=IstrU
      endif
      if (EAST_INTER) then
        imax=Iend+filter_width/2
      else
        imax=Iend
      endif
#  else
      imin=IstrU
      imax=Iend
#  endif
# endif      
      u_fwidth_array(IstrU:Iend,Jstr:Jend)=filter_width
      if (WESTERN_EDGE  .or. EASTERN_EDGE .or.
     &    SOUTHERN_EDGE .or. NORTHERN_EDGE) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            do k=filter_width/2,1,-1
              if(i-k.ge.imin .and. i+k.le.imax
     &     .and. j-k.ge.jmin .and. j+k.le.jmax
# ifdef MASKING 
     &     .and. count(umask(i-k:i+k,j-k:j+k).eq.1.0).eq.(2*k+1)**2
# endif
# ifdef WET_DRY
     &     .and. count(umask_wet(i-k:i+k,j-k:j+k).ne.0.0).eq.(2*k+1)**2
# endif           
     &                                                          ) then
                exit
              else
                u_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
# if defined MASKING | defined WET_DRY
      else
        do j=Jstr,Jend
          do i=IstrU,Iend
            do k=filter_width/2,1,-1
              if(
#  ifdef MASKING
     &          count(umask(i-k:i+k,j-k:j+k).eq.1.0).eq.(2*k+1)**2
#   ifdef WET_DRY
     &     .and.
#   endif              
#  endif
#  ifdef WET_DRY
     &          count(umask_wet(i-k:i+k,j-k:j+k).ne.0.0).eq.(2*k+1)**2
#  endif           
     &                                                          ) then
                exit
              else
                u_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
# endif        
      endif
            
            
# ifdef NS_PERIODIC
      jmin=JstrV-filter_width/2
      jmax=Jend+filter_width/2
# else
#  ifdef MPI
      if (SOUTH_INTER) then
        jmin=JstrV-filter_width/2
      else
        jmin=JstrV
      endif
      if (NORTH_INTER) then
        jmax=Jend+filter_width/2
      else
        jmax=Jend
      endif
#  else
      jmin=JstrV
      jmax=Jend
#  endif
# endif      
# ifdef EW_PERIODIC
      imin=Istr-filter_width/2
      imax=Iend+filter_width/2
# else
#  ifdef MPI
      if (WEST_INTER) then
        imin=Istr-filter_width/2
      else
        imin=Istr
      endif
      if (EAST_INTER) then
        imax=Iend+filter_width/2
      else
        imax=Iend
      endif
#  else
      imin=Istr
      imax=Iend
#  endif
# endif        
      v_fwidth_array(Istr:Iend,JstrV:Jend)=filter_width
      if (WESTERN_EDGE  .or. EASTERN_EDGE .or.
     &    SOUTHERN_EDGE .or. NORTHERN_EDGE) then
        do j=JstrV,Jend
          do i=Istr,Iend
            do k=filter_width/2,1,-1
              if(i-k.ge.imin .and. i+k.le.imax
     &     .and. j-k.ge.jmin .and. j+k.le.jmax
# ifdef MASKING 
     &     .and. count(vmask(i-k:i+k,j-k:j+k).eq.1.0).eq.(2*k+1)**2
# endif
# ifdef WET_DRY
     &     .and. count(vmask_wet(i-k:i+k,j-k:j+k).ne.0.0).eq.(2*k+1)**2
# endif           
     &                                                          ) then
                exit
              else
                v_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
# if defined MASKING | defined WET_DRY
      else
        do j=JstrV,Jend
          do i=Istr,Iend
            do k=filter_width/2,1,-1
              if(
#  ifdef MASKING
     &          count(vmask(i-k:i+k,j-k:j+k).eq.1.0).eq.(2*k+1)**2
#   ifdef WET_DRY
     &     .and.
#   endif              
#  endif
#  ifdef WET_DRY
     &          count(vmask_wet(i-k:i+k,j-k:j+k).ne.0.0).eq.(2*k+1)**2
#  endif           
     &                                                          ) then
                exit
              else
                v_fwidth_array(i,j)=2*(k-1)+1
              endif  
            enddo  
          enddo
        enddo
# endif        
      endif  

      return
      end subroutine

      subroutine leray_smoothing_3d(tile)
      implicit none
      integer  tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"      
# include "compute_tile_bounds.h"
      trd=omp_get_thread_num()
      call leray_smoothing_3d_tile (istr,iend,jstr,jend)
      return
      end

      subroutine leray_smoothing_3d_tile(istr,iend,jstr,jend)
      implicit none
# include "param.h"      
# include "ocean3d.h"
# include "leray_smoothing.h"
# include "scalars.h"      
# include "set_global_definitions.h" 

      integer istr,iend,jstr,jend
      integer i,j,k

# include "compute_auxiliary_bounds.h"     
#  ifdef EW_PERIODIC
#   define IU_RANGE Istr,Iend
#   define IV_RANGE Istr,Iend
#  else
#   define IU_RANGE Istr,IendR
#   define IV_RANGE IstrR,IendR
#  endif

#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr,Jend
#   define JV_RANGE Jstr,Jend
#  else
#   define JU_RANGE JstrR,JendR
#   define JV_RANGE Jstr,JendR
#  endif

      do k=1,N
        do j=JU_RANGE
          do i=IU_RANGE
# ifdef LERAY_FILTER_9PTS         
            if(u_fwidth_array(i,j).eq.9) then
              u_smooth(i,j,k) = sum( u(i-4:i+4,j-4:j+4,k,nrhs)
     &                          * filter_weights(-4:4,-4:4) )
     &                          * inv_weight_sum9
            endif 
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
            if(u_fwidth_array(i,j).eq.7) then 
              u_smooth(i,j,k) = sum( u(i-3:i+3,j-3:j+3,k,nrhs)
     &                          * filter_weights(-3:3,-3:3) )
     &                          * inv_weight_sum7
            endif  
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
            if(u_fwidth_array(i,j).eq.5) then 
              u_smooth(i,j,k) = sum( u(i-2:i+2,j-2:j+2,k,nrhs)
     &                          * filter_weights(-2:2,-2:2) )
     &                          * inv_weight_sum5
            endif  
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
            if(u_fwidth_array(i,j).eq.3) then 
              u_smooth(i,j,k) = sum( u(i-1:i+1,j-1:j+1,k,nrhs)
     &                          * filter_weights(-1:1,-1:1) )
     &                          * inv_weight_sum3
            endif
# endif
            if(u_fwidth_array(i,j).lt.3) then
              u_smooth(i,j,k) = u(i,j,k,nrhs)
            endif
          enddo
        enddo
      enddo  

      do k=1,N
        do j=JV_RANGE
          do i=IV_RANGE
# ifdef LERAY_FILTER_9PTS         
            if(v_fwidth_array(i,j).eq.9) then
              v_smooth(i,j,k) = sum( v(i-4:i+4,j-4:j+4,k,nrhs)
     &                          * filter_weights(-4:4,-4:4) )
     &                          * inv_weight_sum9
            endif 
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
            if(v_fwidth_array(i,j).eq.7) then 
              v_smooth(i,j,k) = sum( v(i-3:i+3,j-3:j+3,k,nrhs)
     &                          * filter_weights(-3:3,-3:3) )
     &                          * inv_weight_sum7
            endif  
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
            if(v_fwidth_array(i,j).eq.5) then 
              v_smooth(i,j,k) = sum( v(i-2:i+2,j-2:j+2,k,nrhs)
     &                          * filter_weights(-2:2,-2:2) )
     &                          * inv_weight_sum5
            endif  
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
            if(v_fwidth_array(i,j).eq.3) then 
              v_smooth(i,j,k) = sum( v(i-1:i+1,j-1:j+1,k,nrhs)
     &                          * filter_weights(-1:1,-1:1) )
     &                          * inv_weight_sum3
            endif
# endif
            if(v_fwidth_array(i,j).lt.3) then
              v_smooth(i,j,k) = v(i,j,k,nrhs)
            endif
          enddo
        enddo
      enddo  


# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef THREE_GHOST_POINTS_UV
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             u_smooth(START_2D_ARRAY,1))
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             v_smooth(START_2D_ARRAY,1))
#  else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend, 
     &                        u_smooth(START_2D_ARRAY,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend, 
     &                        v_smooth(START_2D_ARRAY,1))
#  endif
# endif 

      return
      end subroutine

#else
      subroutine leray_smoothing_empty
      return
      end
#endif
