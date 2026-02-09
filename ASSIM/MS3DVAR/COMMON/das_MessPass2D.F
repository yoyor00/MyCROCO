#include "cppdefs.h"
#ifdef MPI

      subroutine das_MessPass2D_tile (Istr,Iend,Jstr,Jend, A)
!
! This subroutine is designed for ROMS-MPI code. It exchanges domain
! boundary information, including 2 ghost-cells in each direction.
! Ping Wang 9/15/99.
!
      implicit none
# include "param.h"
# include "scalars.h"
      include 'mpif.h'

      integer A(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer Istr,Iend,Jstr,Jend, i,j, isize,jsize,   iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj

      integer sub_X,size_X, sub_E,size_E
      parameter (sub_X=(Lm+NSUB_X-1)/NSUB_X,  size_X=7+2*sub_X, 
     &           sub_E=(Mm+NSUB_E-1)/NSUB_E,  size_E=7+2*sub_E)

      integer buf_snd4(4),  ibuf_sndN(0:size_X),     buf_snd2(4), 
     &     buf_rev4(4),     ibuf_revN(0:size_X),     buf_rev2(4),

     &    jbuf_sndW(0:size_E),                jbuf_sndE(0:size_E),
     &    jbuf_revW(0:size_E),                jbuf_revE(0:size_E), 

     &     buf_snd1(4),     ibuf_sndS(0:size_X),     buf_snd3(4), 
     &     buf_rev1(4),     ibuf_revS(0:size_X),     buf_rev3(4)
c**
c      common /buffers_2D/
c     &     buf_snd4,     ibuf_sndN,     buf_snd2,
c     &     buf_rev4,     ibuf_revN,     buf_rev2,
c
c     &    jbuf_sndW,                    jbuf_sndE,
c     &    jbuf_revW,                    jbuf_revE,
c
c     &     buf_snd1,     ibuf_sndS,     buf_snd3,
c     &     buf_rev1,     ibuf_revS,     buf_rev3
c**
!
#include "compute_message_bounds.h"
!
      isize=2*ishft                  ! sizes for side messages
      jsize=2*jshft                  ! in XI and ETA directions

c*      write(*,'(2(6x,A3,I2,2x,A5,I3,2x,A5,I3))')
c*     &        'ii=',ii,'imin=',imin,'imax=',imax,
c*     &        'jj=',jj,'jmin=',jmin,'jmax=',jmax
#define write !
                            !  Message passing split into two stages
                            !  in order to optimize Send-Recv pairing
                            !  in such a way that if one subdomain
      do iter=0,1           !  sends message to, say, its WESTERN
       mdii=mod(ii+iter,2)  !  neighbor, that neighbor is preparing
       mdjj=mod(jj+iter,2)  !  to receive this message first (i.e.
                            !  message coming from its EASTERN side),
                            !  rather than! rather than send his WEST
                            !  bound message, similarly to the first
                            !  subdomain.
!
! Prepare to receive and send: sides....
!
        if (mdii.eq.0) then
          if (WEST_INTER) then
            write(*,*) 'MessPass2D: 1.1', mynode
            do j=jmin,jmax
              jbuf_sndW(j-jmin      )=A(1,j)
              jbuf_sndW(j-jmin+jshft)=A(2,j)
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_INTEGER,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_INTEGER,
     &                       p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            write(*,*) 'MessPass2D: 1.2', mynode
            do j=jmin,jmax
              jbuf_sndE(j-jmin      )=A(Lm-1,j)
              jbuf_sndE(j-jmin+jshft)=A(Lm  ,j)
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_INTEGER,
     &                         p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_INTEGER,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif

        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            write(*,*) 'MessPass2D: 1.3', mynode
            do i=imin,imax
              ibuf_sndS(i-imin      )=A(i,1)
              ibuf_sndS(i-imin+ishft)=A(i,2)
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_INTEGER,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_INTEGER,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            write(*,*) 'MessPass2D: 1.4', mynode
            do i=imin,imax
              ibuf_sndN(i-imin      )=A(i,Mm-1)
              ibuf_sndN(i-imin+ishft)=A(i,Mm  )
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_INTEGER,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_INTEGER,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
!
! ...corners:
!
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            write(*,*) 'MessPass2D: 1.5', mynode
            buf_snd1(1)=A(1,1)
            buf_snd1(2)=A(2,1)
            buf_snd1(3)=A(1,2)
            buf_snd1(4)=A(2,2)
            call MPI_Irecv (buf_rev1,4, MPI_INTEGER,  p_SW,
     &                               6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1,4, MPI_INTEGER,  p_SW,
     &                               5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            write(*,*) 'MessPass2D: 1.6', mynode
            buf_snd2(1)=A(Lm-1,Mm-1)
            buf_snd2(2)=A(Lm  ,Mm-1)
            buf_snd2(3)=A(Lm-1,Mm  )
            buf_snd2(4)=A(Lm  ,Mm  )
            call MPI_Irecv (buf_rev2,4, MPI_INTEGER,  p_NE,
     &                               5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2,4, MPI_INTEGER,  p_NE,
     &                               6, MPI_COMM_WORLD,        ierr)
          endif
        endif

        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            write(*,*) 'MessPass2D: 1.7', mynode
            buf_snd3(1)=A(Lm-1,1)
            buf_snd3(2)=A(Lm  ,1)
            buf_snd3(3)=A(Lm-1,2)
            buf_snd3(4)=A(Lm  ,2)
            call MPI_Irecv (buf_rev3,4, MPI_INTEGER,  p_SE,
     &                               8, MPI_COMM_WORLD, req(7), ierr)
            call MPI_Send  (buf_snd3,4, MPI_INTEGER,  p_SE,
     &                               7, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            write(*,*) 'MessPass2D: 1.8', mynode
            buf_snd4(1)=A(1,Mm-1)
            buf_snd4(2)=A(2,Mm-1)
            buf_snd4(3)=A(1,Mm  )
            buf_snd4(4)=A(2,Mm  )
            call MPI_Irecv (buf_rev4, 4, MPI_INTEGER, p_NW,
     &                               7, MPI_COMM_WORLD, req(8), ierr)
            call MPI_Send  (buf_snd4, 4, MPI_INTEGER, p_NW,
     &                               8, MPI_COMM_WORLD,         ierr)
          endif
        endif
      enddo   !<-- iter
!
! Wait for completion of receive and fill ghost points: sides...
!

      if (WEST_INTER) then
        write(*,*) 'MessPass2D: 2.1', mynode
        call MPI_Wait (req(1),status(1,1),ierr)
        do j=jmin,jmax
          A(-1,j)=jbuf_revW(j-jmin)
          A( 0,j)=jbuf_revW(j-jmin+jshft)
        enddo
      endif

      if (EAST_INTER) then
        write(*,*) 'MessPass2D: 2.2', mynode
        call MPI_Wait (req(2),status(1,2),ierr)
        do j=jmin,jmax
          A(Lm+1,j)=jbuf_revE(j-jmin)
          A(Lm+2,j)=jbuf_revE(j-jmin+jshft)
        enddo
      endif

      if (SOUTH_INTER) then
        write(*,*) 'MessPass2D: 2.3', mynode
        call MPI_Wait (req(3),status(1,3),ierr)
        do i=imin,imax
          A(i,-1)=ibuf_revS(i-imin )
          A(i, 0)=ibuf_revS(i-imin+ishft)
        enddo
      endif

      if (NORTH_INTER) then
        write(*,*) 'MessPass2D: 2.4', mynode
        call MPI_Wait (req(4),status(1,4),ierr)
        do i=imin,imax
          A(i,Mm+1)=ibuf_revN(i-imin)
          A(i,Mm+2)=ibuf_revN(i-imin+ishft)
        enddo
      endif
!
! ...corners:
!
      if (SOUTH_INTER .and. WEST_INTER) then
        write(*,*) 'MessPass2D: 2.5', mynode
        call MPI_Wait (req(5),status(1,5),ierr)
        A(-1,-1)=buf_rev1(1)
        A( 0,-1)=buf_rev1(2)
        A(-1, 0)=buf_rev1(3)
        A( 0, 0)=buf_rev1(4)
      endif

      if (NORTH_INTER .and. EAST_INTER) then
        write(*,*) 'MessPass2D: 2.6', mynode
        call MPI_Wait (req(6),status(1,6),ierr)
        A(Lm+1,Mm+1)=buf_rev2(1)
        A(Lm+2,Mm+1)=buf_rev2(2)
        A(Lm+1,Mm+2)=buf_rev2(3)
        A(Lm+2,Mm+2)=buf_rev2(4)
      endif

      if (SOUTH_INTER .and. EAST_INTER) then
        write(*,*) 'MessPass2D: 2.7', mynode
        call MPI_Wait (req(7),status(1,7),ierr)
        A(Lm+1,-1)=buf_rev3(1)
        A(Lm+2,-1)=buf_rev3(2)
        A(Lm+1, 0)=buf_rev3(3)
        A(Lm+2, 0)=buf_rev3(4)
      endif

      if (NORTH_INTER .and. WEST_INTER) then
        write(*,*) 'MessPass2D: 2.8', mynode
        call MPI_Wait (req(8),status(1,8),ierr)
        A(-1,Mm+1)=buf_rev4(1)
        A( 0,Mm+1)=buf_rev4(2)
        A(-1,Mm+2)=buf_rev4(3)
        A( 0,Mm+2)=buf_rev4(4)
      endif
      write(*,*) 'MessPass2D:    ', mynode,' exit'
      return
      end
#else
      subroutine das_MessPass2D_empty
      return
      end
#endif
