c     this file is not differentiated

# include "cppdefs.h"

      subroutine read_obs()
      implicit none
# include "netcdf.inc"
# include "adinc.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range,k
      character*3 mynode_str
      
#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      write(mynode_str, 910) mynode
      open(unit=333, file='obs.'//mynode_str//'.dat', form='formatted',
     &     access='stream')

      do tile=my_first,my_last
         call read_obs_2d(tile,333)
      end do

      close (333)
 910  format(I3.3)
      end subroutine read_obs

      subroutine read_obs_2d(tile,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile,fid
# include "compute_tile_bounds.h"

      call read_obs_2d_tile(Istr,Iend,Jstr,Jend,fid)

      end subroutine read_obs_2d

      subroutine read_obs_2d_tile(Istr,Iend,Jstr,Jend,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k,mk,ik
      integer Istr,Iend,Jstr,Jend,fid,i,j,m,l

      mk = ad_ast + ad_nt*ad_ns+3
      do k=1,mk
         do l = 1,ad_nobs
            read(fid,*) ik,i,j,ad_obs(l,ik)
         end do
      end do
      
      end subroutine read_obs_2d_tile

      subroutine write_obs(k,fid)
      implicit none
# include "netcdf.inc"
# include "adinc.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range,k,fid

#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      do tile=my_first,my_last
         call write_obs_2d(tile,k,fid)
      end do


      end subroutine write_obs

      subroutine write_obs_2d(tile,k,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile,k,fid
# include "compute_tile_bounds.h"

      call write_obs_2d_tile(Istr,Iend,Jstr,Jend,k,fid)

      end subroutine write_obs_2d

      subroutine write_obs_2d_tile(Istr,Iend,Jstr,Jend,k,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,fid,i,j,k,l

      do l=1,ad_nobs
         i=ad_obs_i(l)
         j=ad_obs_j(l)

         if (i.eq.0) then
            write(*,*) 'WWWWWWWWWWWWWWWWWWWWWWW',i,l,ad_nobs
         end if
         write(fid,*) k,i,j,zeta(i,j,knew)
      end do

      end subroutine write_obs_2d_tile


      subroutine write_grid()
      implicit none
# include "netcdf.inc"
# include "adinc.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range,k,fid

#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      do tile=my_first,my_last
         call write_grid_2d(tile)
      end do


      end subroutine write_grid

      subroutine write_grid_2d(tile)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile,k,fid
# include "compute_tile_bounds.h"

      call write_grid_2d_tile(Istr,Iend,Jstr,Jend,tile)

      end subroutine write_grid_2d

      subroutine write_grid_2d_tile(Istr,Iend,Jstr,Jend,tile)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend,tile,i,j,k
      character*6 tile_str

      write(tile_str,999) mynode,tile
      open(unit=799, file='grid.'//tile_str)
      write(799, *) Iend-Istr+1,Jend-Jstr+1,0
      do j=Jstr,Jend
         do i=Istr,Iend
#ifdef SPHERICAL
            write(799, *) i,j,lonr(i,j),latr(i,j)
#else
            write(799, *) i,j,xr(i,j),yr(i,j)
#endif
         end do
      end do
      close(799)
 999  format(I2.2,'-',I2.2)
      end subroutine write_grid_2d_tile



c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_pos()
      implicit none
# include "netcdf.inc"
# include "adinc.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range,k
      character*2 mynode_str
      
#ifdef OPENMP
      ntrds=omp_get_num_threads()
      trd=omp_get_thread_num()
#else
      ntrds=1
      trd=0
#endif
      range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
      my_first=trd*range
      my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)

      write(mynode_str, 910) mynode
      open(unit=333, file='tideg-'//mynode_str//'.txt', form='formatted',
     &     access='sequential',action='read')

      do tile=my_first,my_last
         call read_pos_2d(tile,333)
      end do

      close (333)
 910  format(I2.2)
      end subroutine read_pos

      subroutine read_pos_2d(tile,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile,fid
# include "compute_tile_bounds.h"

      call read_pos_2d_tile(Istr,Iend,Jstr,Jend,fid)

      end subroutine read_pos_2d

      subroutine read_pos_2d_tile(Istr,Iend,Jstr,Jend,fid)
      implicit none
# include "adinc.h"
# include "adparam.h"

      integer k,mk
      integer Istr,Iend,Jstr,Jend,fid,i,j,m,l

      l=0
      do
         l = l+1
         read(fid,*,END=99) ad_obs_i(l), ad_obs_j(l)
      end do
 99   ad_nobs = l-1
      if (ad_obs_i(ad_nobs).eq.0) then
         write(*,*) 'RRRRRRRRRRRRRRRRRRR', ad_obs_i(l), l, ad_nobs
      end if
      end subroutine read_pos_2d_tile
