! $Id: zoombc_3Dfast.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! This routine belongs to the specific CROCO package.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!
#include "cppdefs.h"

#if defined AGRIF && defined NBQ
!
!====================================================================
!                   subroutine unbq_bc_interp_tile
!====================================================================
!
      subroutine unbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
      implicit none
# include "param.h"
# include "boundary.h"  
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "zoom.h"
# include "ocean3d.h"
# include "nbq.h"

      integer Istr,Iend,Jstr,Jend, i,j,k
      real tinterp
      real :: ainterp(2)
      external qdmuinterp
      !
# include "compute_auxiliary_bounds.h"
!
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif

      if(Agrif_NbStepint() == 0) then
 !       write(*,*) 'WE INTERP iif = ', iif
 !       write(*,*)
# ifdef MASKING
        Agrif_UseSpecialValue = .true.
# endif
        Agrif_SpecialValue = 0.
         tinterp=1.
        Call Agrif_Set_bc(qdmunbqid,(/0,0/),
     &   InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(qdmunbqid,calledweight=tinterp,
     & procname=qdmuinterp)
        Agrif_UseSpecialValue=.false.
      endif

      ainterp(1)= (REAL(Agrif_Nbstepint()) + 1.) / Agrif_Rhot()
      ainterp(2)=1.- ainterp(1)
 !     write(*,*) 'U:' 
 !     write(*,*) maxval(Unbq_north)
 !     write(*,*) maxval(Unbq_south)
 !     write(*,*) maxval(Unbq_east)
 !     write(*,*) maxval(Unbq_west)
 !     write(*,*)


# ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N      
          do i=Istr,IendR
#  if defined NBQ_FRC_BRY     
            unbqbry_south(i,k)=
#  else
            unbqclm(i,Jstr-1,k)=
#  endif
     &  (ainterp(1)*Unbq_south(i,Jstr-1,k,1)+
     &   ainterp(2)*Unbq_south(i,Jstr-1,k,2))
#  ifdef MASKING     
     &    *umask(i,Jstr-1)
#  endif
          enddo
        enddo
      endif            
# endif

# ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=Istr,IendR
#  if defined NBQ_FRC_BRY     
            unbqbry_north(i,k)=
#  else   
            unbqclm(i,Jend+1,k)=
#  endif
     &  (ainterp(1)*Unbq_north(i,Jend+1,k,1)+
     &   ainterp(2)*Unbq_north(i,Jend+1,k,2))     
#  ifdef MASKING     
     &    *umask(i,Jend+1)
#  endif
          enddo
        enddo
      endif      
# endif
 
# ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  if defined NBQ_FRC_BRY     
           unbqbry_west(j,k)=
#  else     
            unbqclm(Istr,j,k)=
#  endif
     &  (ainterp(1)*Unbq_west(Istr,j,k,1)+
     &   ainterp(2)*Unbq_west(Istr,j,k,2))
#  ifdef MASKING     
     &    *umask(Istr,j)
#  endif
          enddo
        enddo
      endif      
# endif
      
# ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=1,N      
          do j=JstrR,JendR
#  if defined NBQ_FRC_BRY  
            unbqbry_east(j,k)=
#  else
            unbqclm(Iend+1,j,k)=
#  endif
     &  (ainterp(1)*Unbq_east(Iend+1,j,k,1)+
     &   ainterp(2)*Unbq_east(Iend+1,j,k,2))         
#  ifdef MASKING     
     &    *umask(Iend+1,j)
#  endif
          enddo
        enddo
      endif
# endif
    
      return
      end
!
!====================================================================
!                   subroutine vnbq_bc_interp_tile
!====================================================================
!
      subroutine vnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
!     implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"
# include "climat.h"
# include "boundary.h"
# include "zoom.h"
# include "coupling.h"
# include "nbq.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real tinterp
      real :: ainterp(2)
      external qdmvinterp   
  
# ifdef MPI
       include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   
!          
# include "compute_auxiliary_bounds.h"
!
      if(Agrif_NbStepint() == 0) then
# ifdef MASKING
        Agrif_UseSpecialValue = .true.
# endif
        Agrif_SpecialValue = 0.
         tinterp=1.
        Call Agrif_Set_bc(qdmvnbqid,(/0,0/),
     &   InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(qdmvnbqid,calledweight=tinterp,
     & procname=qdmvinterp)
        Agrif_UseSpecialValue=.false.
      endif

      ainterp(1)= (REAL(Agrif_Nbstepint()) + 1.) / Agrif_Rhot()
      ainterp(2)=1.- ainterp(1)
 !     write(*,*) 'V:' 
 !     write(*,*) maxval(Vnbq_north)
 !     write(*,*) maxval(Vnbq_south)
 !     write(*,*) maxval(Vnbq_east)
 !     write(*,*) maxval(Vnbq_west)
 !     write(*,*)


#  ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=1,N      
          do i=IstrR,IendR
#   if defined NBQ_FRC_BRY   
            vnbqbry_south(i,k)=
#   else
            vnbqclm(i,Jstr,k)=
#   endif
     &    (ainterp(1)*Vnbq_south(i,Jstr,k,1)+
     &     ainterp(2)*Vnbq_south(i,Jstr,k,2))
#   ifdef MASKING     
     &    *vmask(i,Jstr)
#   endif
          enddo
        enddo
      endif
#  endif

#  ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=1,N      
          do i=IstrR,IendR
#   if defined NBQ_FRC_BRY
            vnbqbry_north(i,k)=
#   else
            vnbqclm(i,Jend+1,k)=
#   endif
     &    (ainterp(1)*Vnbq_north(i,Jend+1,k,1)+
     &     ainterp(2)*Vnbq_north(i,Jend+1,k,2))     
#   ifdef MASKING     
     &    *vmask(i,Jend+1)
#   endif 
          enddo
        enddo
      endif      
#  endif
#  ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=1,N      
          do j=Jstr,JendR
#   if defined NBQ_FRC_BRY
            vnbqbry_west(j,k)=
#   else
            vnbqclm(Istr-1,j,k)=
#   endif
     &    (ainterp(1)*Vnbq_west(Istr-1,j,k,1)+
     &     ainterp(2)*Vnbq_west(Istr-1,j,k,2))  
#   ifdef MASKING     
     &    *vmask(Istr-1,j)
#   endif
          enddo
        enddo
      endif      
#  endif

#  ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=1,N      
          do j=Jstr,JendR
#   if defined M3_FRC_BRY  
            vnbqbry_east(j,k)=
#   else
            vnbqclm(Iend+1,j,k)=
#   endif
     &    (ainterp(1)*Vnbq_east(Iend+1,j,k,1)+
     &     ainterp(2)*Vnbq_east(Iend+1,j,k,2))      
#   ifdef MASKING     
     &    *vmask(Iend+1,j)
#   endif
          enddo
        enddo
      endif
#  endif    
       return
      end       
!
!====================================================================
!                   subroutine wnbq_bc_interp_tile
!====================================================================
!
      subroutine wnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
      use AGRIF_Util
!     implicit none
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"
# include "climat.h"
# include "boundary.h"
# include "zoom.h"
# include "coupling.h"
# include "nbq.h"
      integer Istr,Iend,Jstr,Jend, i,j,k
      real tinterp
      real :: ainterp(2)
      external qdmwinterp   
  
# ifdef MPI
       include 'mpif.h'
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif   
!          
# include "compute_auxiliary_bounds.h"
!
      if(Agrif_NbStepint() == 0) then
# ifdef MASKING
        Agrif_UseSpecialValue = .true.
# endif
        Agrif_SpecialValue = 0.
         tinterp=1.
        Call Agrif_Set_bc(qdmwnbqid,(/0,0/),
     &   InterpolationShouldbemade=.TRUE.)
        Call Agrif_Bc_variable(qdmwnbqid,calledweight=tinterp,
     & procname=qdmwinterp)
        Agrif_UseSpecialValue=.false.
      endif

      ainterp(1)= (REAL(Agrif_Nbstepint()) + 1.) / Agrif_Rhot()
      ainterp(2)=1.- ainterp(1)
 !     write(*,*) 'V:' 
 !     write(*,*) maxval(Vnbq_north)
 !     write(*,*) maxval(Vnbq_south)
 !     write(*,*) maxval(Vnbq_east)
 !     write(*,*) maxval(Vnbq_west)
 !     write(*,*)


#  ifdef AGRIF_OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do k=0,N      
          do i=IstrR,IendR
#   if defined NBQ_FRC_BRY   
            wnbqbry_south(i,k)=
#   else
            wnbqclm(i,Jstr-1,k)=
#   endif
     &    (ainterp(1)*Wnbq_south(i,Jstr-1,k,1)+
     &     ainterp(2)*Wnbq_south(i,Jstr-1,k,2))
#   ifdef MASKING     
     &    *rmask(i,Jstr-1)
#   endif
          enddo
        enddo
      endif
#  endif

#  ifdef AGRIF_OBC_NORTH
      if (NORTHERN_EDGE) then      
        do k=0,N      
          do i=IstrR,IendR
#   if defined NBQ_FRC_BRY
            wnbqbry_north(i,k)=
#   else
            wnbqclm(i,Jend+1,k)=
#   endif
     &    (ainterp(1)*Wnbq_north(i,Jend+1,k,1)+
     &     ainterp(2)*Wnbq_north(i,Jend+1,k,2))     
#   ifdef MASKING     
     &    *rmask(i,Jend+1)
#   endif 
          enddo
        enddo
      endif      
#  endif
#  ifdef AGRIF_OBC_WEST
      if (WESTERN_EDGE) then      
        do k=0,N      
          do j=JstrR,JendR
#   if defined NBQ_FRC_BRY
            wnbqbry_west(j,k)=
#   else
            wnbqclm(Istr-1,j,k)=
#   endif
     &    (ainterp(1)*Wnbq_west(Istr-1,j,k,1)+
     &     ainterp(2)*Wnbq_west(Istr-1,j,k,2))  
#   ifdef MASKING     
     &    *rmask(Istr-1,j)
#   endif
          enddo
        enddo
      endif      
#  endif

#  ifdef AGRIF_OBC_EAST
      if (EASTERN_EDGE) then      
        do k=0,N      
          do j=JstrR,JendR
#   if defined M3_FRC_BRY  
            wnbqbry_east(j,k)=
#   else
            wnbqclm(Iend+1,j,k)=
#   endif
     &    (ainterp(1)*Wnbq_east(Iend+1,j,k,1)+
     &     ainterp(2)*Wnbq_east(Iend+1,j,k,2))      
#   ifdef MASKING     
     &    *rmask(Iend+1,j)
#   endif
          enddo
        enddo
      endif
#  endif    
       return
      end         


      subroutine qdmuinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "nbq.h"
# include "zoom.h"

       integer i1,i2,j1,j2,k1,k2
       real tabres(i1:i2,j1:j2,k1:k2)
       logical before
       integer nb,ndir
       integer nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then
         tabres(i1:i2,j1:j2,k1:k2) = qdmu_nbq(i1:i2,j1:j2,k1:k2)
       else

#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      nsource=1
      ndest=2 

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
            if (.not.SOUTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Unbq_south(i1:i2,j1:j2,1:N,ndest)=
     &          Unbq_south(i1:i2,j1:j2,1:N,nsource)
             endif
                Unbq_south(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
            if (.not.NORTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Unbq_north(i1:i2,j1:j2,1:N,ndest)=
     &          Unbq_north(i1:i2,j1:j2,1:N,nsource)
       endif
                Unbq_north(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
            if (.not.WEST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Unbq_west(i1:i2,j1:j2,1:N,ndest)=
     &          Unbq_west(i1:i2,j1:j2,1:N,nsource)
        endif
                Unbq_west(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
            if (.not.EAST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Unbq_east(i1:i2,j1:j2,1:N,ndest)=
     &          Unbq_east(i1:i2,j1:j2,1:N,nsource)
       endif
                Unbq_east(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

       endif
      return
      end

      subroutine qdmvinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "nbq.h"
# include "zoom.h"

       integer i1,i2,j1,j2,k1,k2
       real tabres(i1:i2,j1:j2,k1:k2)
       logical before
       integer nb,ndir
       integer nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then
         tabres(i1:i2,j1:j2,k1:k2) = qdmv_nbq(i1:i2,j1:j2,k1:k2)
       else

#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      nsource=1
      ndest=2 

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
            if (.not.SOUTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then
                Vnbq_south(i1:i2,j1:j2,1:N,ndest)=
     &          Vnbq_south(i1:i2,j1:j2,1:N,nsource)
                endif

                Vnbq_south(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
            if (.not.NORTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then
               Vnbq_north(i1:i2,j1:j2,1:N,ndest)=
     &          Vnbq_north(i1:i2,j1:j2,1:N,nsource)
        endif
                Vnbq_north(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
            if (.not.WEST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Vnbq_west(i1:i2,j1:j2,1:N,ndest)=
     &          Vnbq_west(i1:i2,j1:j2,1:N,nsource)
             endif
                Vnbq_west(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
            if (.not.EAST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Vnbq_east(i1:i2,j1:j2,1:N,ndest)=
     &          Vnbq_east(i1:i2,j1:j2,1:N,nsource)
              endif
                Vnbq_east(i1:i2,j1:j2,1:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

       endif
      return
      end

      subroutine qdmwinterp(tabres,i1,i2,j1,j2,k1,k2,before,nb,ndir)
      implicit none
# include "param.h"
# include "ocean3d.h"
# include "scalars.h"
# include "nbq.h"
# include "zoom.h"

       integer i1,i2,j1,j2,k1,k2
       real tabres(i1:i2,j1:j2,k1:k2)
       logical before
       integer nb,ndir
       integer nsource, ndest
       logical :: western_side, eastern_side
       logical :: northern_side,southern_side

       if (before) then
         tabres(i1:i2,j1:j2,k1:k2) = qdmw_nbq(i1:i2,j1:j2,k1:k2)
       else

#  ifdef MPI
#   define LOCALLM Lmmpi
#   define LOCALMM Mmmpi
#  else
#   define LOCALLM Lm
#   define LOCALMM Mm
#  endif

         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      nsource=1
      ndest=2 

# ifdef AGRIF_OBC_SOUTH
        if (southern_side) then
#  ifdef MPI
            if (.not.SOUTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then
                Wnbq_south(i1:i2,j1:j2,0:N,ndest)=
     &          Wnbq_south(i1:i2,j1:j2,0:N,nsource)
                endif
                Wnbq_south(i1:i2,j1:j2,0:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_NORTH
        if (northern_side) then
#  ifdef MPI
            if (.not.NORTH_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then
                Wnbq_north(i1:i2,j1:j2,0:N,ndest)=
     &          Wnbq_north(i1:i2,j1:j2,0:N,nsource)
               endif
                Wnbq_north(i1:i2,j1:j2,0:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_WEST
        if (western_side) then
#  ifdef MPI
            if (.not.WEST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then
                Wnbq_west(i1:i2,j1:j2,0:N,ndest)=
     &          Wnbq_west(i1:i2,j1:j2,0:N,nsource)
                 endif
                Wnbq_west(i1:i2,j1:j2,0:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

# ifdef AGRIF_OBC_EAST
        if (eastern_side) then
#  ifdef MPI
            if (.not.EAST_INTER) then
#  endif
                   if (Agrif_NbStepint() == 0) then

                Wnbq_east(i1:i2,j1:j2,0:N,ndest)=
     &          Wnbq_east(i1:i2,j1:j2,0:N,nsource)
             endif
                Wnbq_east(i1:i2,j1:j2,0:N,nsource) = tabres
#  ifdef MPI
            endif
#  endif
        endif
# endif

       endif
      return
      end
!
!====================================================================
!                   subroutine rnbq_bc_interp_tile
!====================================================================
!
      subroutine rnbq_bc_interp_tile(Istr,Iend,Jstr,Jend)
!
! \TODO To be continued
!
!
      integer indx,itrc
      indx=indxR              ! \TODO check index
      itrc=1                  ! \TODO set index
      !call t3dbc_interp_tile(Istr,Iend,Jstr,Jend,indx,itrc)
      return
      end       
!
#else
      subroutine zoombc_3Dfast_empty()
      return
      end
#endif
