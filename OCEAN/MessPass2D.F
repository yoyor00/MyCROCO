!
! $Id: MessPass2D.F 1143 2013-05-17 08:17:23Z serena $
!
#include "cppdefs.h"
#ifdef MPI

# ifndef MP_3PTS
#  ifdef OPENACC_MPI
      subroutine MessPass2D_tile_device (Istr,Iend,Jstr,Jend, A)
#  else
      subroutine MessPass2D_tile (Istr,Iend,Jstr,Jend, A)
#  endif
# else
#  ifdef OPENACC_MPI
      subroutine MessPass2D_3pts_tile_device (Istr,Iend,Jstr,Jend, A)
#  else
      subroutine MessPass2D_3pts_tile (Istr,Iend,Jstr,Jend, A)
#  endif
# endif
!
!======================================================================
!
! This subroutine is designed for ROMS-MPI code. It exchanges domain
! boundary information, including 2 (or 3) ghost-cells in each direction.
!
! Ping Wang 9/15/99.
! Patrick Marchesiello 2012: generic number of ghost-cells Npts
!
!======================================================================
!     implicit none
# include "param.h"
# include "scalars.h"
# include "mpi_cpl.h"
      include 'mpif.h'
!
! Nb of boundary points involved in communication
!
      integer Npts,ipts,jpts
# ifndef MP_3PTS
      parameter (Npts=2)
# else
      parameter (Npts=3)
# endif
      real A(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer Istr,Iend,Jstr,Jend, i,j, isize,jsize,ksize,  iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj
      real buf_snd4(Npts*Npts),     buf_snd2(Npts*Npts), 
     &     buf_rev4(Npts*Npts),     buf_rev2(Npts*Npts),
     &     buf_snd1(Npts*Npts),     buf_snd3(Npts*Npts), 
     &     buf_rev1(Npts*Npts),     buf_rev3(Npts*Npts)
      integer sub_X,size_X, sub_E,size_E
#ifdef OPENACC_MPI
!$acc declare create(
!$acc&  buf_snd4,buf_snd2
!$acc& ,buf_rev4,buf_rev2
!$acc& ,buf_snd1,buf_snd3
!$acc& ,buf_rev1,buf_rev3
!$acc& )
#endif
#ifndef AGRIF      
      parameter (sub_X=Lm,  size_X=Npts*(sub_X+2*Npts)-1, 
     &           sub_E=Mm,  size_E=Npts*(sub_E+2*Npts)-1)
      real ibuf_sndN(0:size_X), ibuf_revN(0:size_X), 
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E)
#ifdef OPENACC_MPI
!$acc declare create(
!$acc&  ibuf_sndN,ibuf_revN
!$acc& ,ibuf_sndS,ibuf_revS
!$acc& ,jbuf_sndW,jbuf_sndE
!$acc& ,jbuf_revW,jbuf_revE
!$acc& )
#endif
#else
      real, dimension(:), allocatable :: 
     &     ibuf_sndN, ibuf_revN, 
     &     ibuf_sndS, ibuf_revS,
     &     jbuf_sndW, jbuf_sndE,
     &     jbuf_revW, jbuf_revE 
#endif
!
# include "compute_message_bounds.h"

# ifdef AGRIF
         sub_X=Lm
         size_X=Npts*(sub_X+2*Npts)-1  !7+Npts*sub_X
         sub_E=Mm
         size_E=Npts*(sub_E+2*Npts)-1  !7+Npts*sub_E
         Allocate(ibuf_sndN(0:size_X), ibuf_revN(0:size_X), 
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E))
# endif
!
      ksize=Npts*Npts
      isize=Npts*ishft                  ! sizes for side messages
      jsize=Npts*jshft                  ! in XI and ETA directions

#define write !
c*      write(*,'(2(6x,A3,I2,2x,A5,I3,2x,A5,I3))')
c*     &        'ii=',ii,'imin=',imin,'imax=',imax,
c*     &        'jj=',jj,'jmin=',jmin,'jmax=',jmax
                            !  Message passing split into two stages
                            !  in order to optimize Send-Recv pairing
                            !  in such a way that if one subdomain
      do iter=0,1           !  sends message to, say, its WESTERN
        mdii=mod(ii+iter,2) !  neighbor, that neighbor is preparing
        mdjj=mod(jj+iter,2) !  to receive this message first (i.e.
                            !  message coming from its EASTERN side),
                            !  rather than send his WEST
                            !  bound message, similarly to the first
                            !  subdomain.
!
! Prepare to receive and send: sides....
        if (mdii.eq.0) then
          if (WEST_INTER2) then
            write(*,*) 'MessPass2D: 1.1', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndW(j-jmin+(ipts-1)*jshft)=A(ipts,j)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( jbuf_revW, jbuf_sndW )    
#endif
           call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
          endif
        else
          if (EAST_INTER2) then
           write(*,*) 'MessPass2D: 1.2', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
           do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndE(j-jmin+(ipts-1)*jshft)=A(Lmmpi-Npts+ipts,j)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( jbuf_revE, jbuf_sndE )    
#endif
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
          endif
        endif

        if (mdjj.eq.0) then
          if (SOUTH_INTER2) then
            write(*,*) 'MessPass2D: 1.3', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndS(i-imin+(jpts-1)*ishft)=A(i,jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( ibuf_revS, ibuf_sndS )    
#endif
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
          endif
        else
          if (NORTH_INTER2) then
            write(*,*) 'MessPass2D: 1.4', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndN(i-imin+(jpts-1)*ishft)=A(i,Mmmpi-Npts+jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( ibuf_revN, ibuf_sndN )    
#endif
           call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
         endif
        endif
!
! ...corners:
!

        if (mdii.eq.0) then
          if (CORNER_SW) then
            write(*,*) 'MessPass2D: 1.5', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd1(ipts+Npts*(jpts-1))=A(ipts,jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( buf_rev1, buf_snd1 )    
#endif
            call MPI_Irecv (buf_rev1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               5, MPI_COMM_WORLD,        ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
         endif
        else
          if (CORNER_NE) then
            write(*,*) 'MessPass2D: 1.6', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd2(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( buf_rev2, buf_snd2 )    
#endif
            call MPI_Irecv (buf_rev2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               6, MPI_COMM_WORLD,        ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
          endif
        endif
        if (mdii.eq.1) then
          if (CORNER_SE) then
            write(*,*) 'MessPass2D: 1.7', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd3(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( buf_rev3, buf_snd3 )    
#endif
           call MPI_Irecv (buf_rev3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               8, MPI_COMM_WORLD, req(7), ierr)
            call MPI_Send  (buf_snd3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               7, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
          endif
       else
          if (CORNER_NW) then
            write(*,*) 'MessPass2D: 1.8', mynode
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd4(ipts+Npts*(jpts-1))=
     &                                A(ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
#ifdef OPENACC_MPI
!$acc end kernels          
!$acc host_data use_device( buf_rev4, buf_snd4 )    
#endif
           call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               7, MPI_COMM_WORLD, req(8), ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               8, MPI_COMM_WORLD,         ierr)
#ifdef OPENACC_MPI
!$acc end host_data    
#endif
         endif
        endif
      enddo   !<-- iter
!
! Wait for completion of receive and fill ghost points: sides...
!
      if (WEST_INTER2) then
        write(*,*) 'MessPass2D: 2.1', mynode
        call MPI_Wait (req(1),status(1,1),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do j=jmin,jmax
          do ipts=1,Npts
           A(ipts-Npts,j)=jbuf_revW(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if ( WEST_INTER .and. .not. WEST_INTER2 ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do j=jmin,jmax
          do ipts=1,Npts
           A(ipts-Npts,j)=A(ipts,j)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif

      if (EAST_INTER2) then
        write(*,*) 'MessPass2D: 2.2', mynode
        call MPI_Wait (req(2),status(1,2),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do j=jmin,jmax
          do ipts=1,Npts
           A(Lmmpi+ipts,j)=jbuf_revE(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if ( EAST_INTER .and. .not. EAST_INTER2 ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do j=jmin,jmax
          do ipts=1,Npts
           A(Lmmpi+ipts,j)=A(Lmmpi-Npts+ipts,j)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif

      if (SOUTH_INTER2) then
        write(*,*) 'MessPass2D: 2.3', mynode
        call MPI_Wait (req(3),status(1,3),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do i=imin,imax
          do jpts=1,Npts
           A(i,jpts-Npts)=ibuf_revS(i-imin+(jpts-1)*ishft)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if ( SOUTH_INTER .and. .not. SOUTH_INTER2 ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do i=imin,imax
          do jpts=1,Npts
           A(i,jpts-Npts)=A(i,jpts)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif

      if (NORTH_INTER2) then
        write(*,*) 'MessPass2D: 2.4', mynode
        call MPI_Wait (req(4),status(1,4),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do i=imin,imax
          do jpts=1,Npts
            A(i,Mmmpi+jpts)=ibuf_revN(i-imin+(jpts-1)*ishft)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if (NORTH_INTER .AND. .not. NORTH_INTER2) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do i=imin,imax
          do jpts=1,Npts
            A(i,Mmmpi+jpts)=A(i,Mmmpi-Npts+jpts)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif

!
! ...corners:
!
      if ( CORNER_SW) then
        write(*,*) 'MessPass2D: 2.5', mynode
        call MPI_Wait (req(5),status(1,5),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
          A(ipts-Npts,jpts-Npts)=buf_rev1(ipts+Npts*(jpts-1))
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if (.not.CORNER_SW .and.
     &  SOUTH_INTER .and.WEST_INTER ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
           A(ipts-Npts,jpts-Npts)=A(ipts,jpts)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
       endif  

      if ( CORNER_NE) then
        write(*,*) 'MessPass2D: 2.6', mynode
        call MPI_Wait (req(6),status(1,6),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,Mmmpi+jpts)=buf_rev2(ipts+Npts*(jpts-1))
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if (.not.CORNER_NE .and.
     &  NORTH_INTER .and.EAST_INTER ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,Mmmpi+jpts)=A(Lmmpi+ipts-NPts,Mmmpi+jpts-Npts)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
       endif  

      if (CORNER_SE) then
        write(*,*) 'MessPass2D: 2.7', mynode
        call MPI_Wait (req(7),status(1,7),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,jpts-Npts)=buf_rev3(ipts+Npts*(jpts-1))
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if (.not. CORNER_SE .and. 
     &  SOUTH_INTER .and. EAST_INTER ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,jpts-Npts)=A(Lmmpi+ipts-Npts,jpts)
          enddo
        enddo      
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif  


      if (CORNER_NW) then
        write(*,*) 'MessPass2D: 2.8', mynode
        call MPI_Wait (req(8),status(1,8),ierr)
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
        do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,Mmmpi+jpts)=buf_rev4(ipts+Npts*(jpts-1))
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
      endif
      if (.not. CORNER_NW .and.
     &   NORTH_INTER .and.  WEST_INTER ) then
#ifdef OPENACC_MPI
!$acc kernels default(present)         
!$acc loop independent collapse(2)
#endif
       do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,Mmmpi+jpts)=A(ipts,Mmmpi+jpts-NPts)
          enddo
        enddo
#ifdef OPENACC_MPI
!$acc end kernels          
#endif
       endif 


c      write(*,*) 'MessPass2D:    ', mynode,' exit'

#if defined AGRIF
         DeAllocate(ibuf_sndN, ibuf_revN, 
     &     ibuf_sndS, ibuf_revS,
     &     jbuf_sndW, jbuf_sndE,
     &     jbuf_revW, jbuf_revE)
#endif

      return
      end

# ifndef MP_3PTS
#  define MP_3PTS
#  include "MessPass2D.F"
#  undef MP_3PTS
# endif

# ifndef MP_3PTS
# ifndef OPENACC_MPI
#  define OPENACC_MPI
#  include "MessPass2D.F"
#  undef OPENACC_MPI
# endif
# endif


#else
      subroutine MessPass2D_empty
      return
      end
#endif
