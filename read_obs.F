c     this file is not differentiated

# include "cppdefs.h"

      subroutine read_obs()
      implicit none
# include "netcdf.inc"
# include "adinc.h"
# include "adparam.h"

      integer ntrds,trd,my_first,my_last,tile,range

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

      open(unit=333, file='obs.dat', form='unformatted', access='stream')

      do tile=my_first,my_last
         call read_obs_2d(tile,333)
      end do

      close (333)

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

      integer Istr,Iend,Jstr,Jend,fid,i,j,k

      do k=ad_ast,ad_ast+ad_nobs-1
         do j=Jstr,Jend
            do i=Istr,Iend
               read(fid) ad_obs(i,j,k-ad_ast+1)
               write (*,*) 'ad_obs:',i,j,k,ad_obs(i,j,k-ad_ast+1)
            end do
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

      integer Istr,Iend,Jstr,Jend,fid,i,j,k

      do j=Jstr,Jend
         do i=Istr,Iend
            write(fid) zeta(i,j,knew)
            write (*,*) 'zeta:',i,j,k,zeta(i,j,knew)
         end do
      end do

      end subroutine write_obs_2d_tile
