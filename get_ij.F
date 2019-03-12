      subroutine get_ij(lon,lat,i,j,out)
      implicit none
# include "adinc.h"
# include "adparam.h"
      double precision lon,lat
      logical out
      integer ntrds,trd,my_first,my_last,tile,range,i,j

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
         call get_ij_2d(tile,lon,lat,i,j,out)
      end do

      end subroutine get_ij

      subroutine get_ij_2d(tile,lon,lat,i,j,out)
      implicit none
# include "adinc.h"
# include "adparam.h"
      integer tile,i,j
      logical out
      double precision lon,lat
# include "compute_tile_bounds.h"

      call get_ij_2d_tile(Istr,Iend,Jstr,Jend,lon,lat,i,j,out)

      end subroutine get_ij_2d

      subroutine get_ij_2d_tile(Istr,Iend,Jstr,Jend,lon,lat,i,j,out)
      implicit none

# include "adinc.h"
# include "adparam.h"      

      integer Istr, Iend, Jstr, Jend, i, j, istep, jstep
      double precision lon,lat
      double precision distance
      double precision dlon,dlat
      logical out

c      write(*,*) '==========getij=========='
      
      if ((lon.gt.lonr(Iend,Jend)).or.
     &     (lon.lt.lonr(Istr,Jstr)).or.
     &     (lat.gt.latr(Iend,Jend)).or.
     &     (lat.lt.latr(Istr,Jstr))) then
         out = .true.
c         write(*,*) '|||',mynode,lon,lat,lonr(Iend,Jend),lonr(Istr,Jstr)
c     &        ,latr(Iend,Jend),latr(Istr,Jstr)
      else
c         write(*,*) '+++',mynode,lon,lat,lonr(Iend,Jend),lonr(Istr,Jstr)
c     &        ,latr(Iend,Jend),latr(Istr,Jstr)
         out = .false.
      
         i = (Iend-Istr)/2
         j = (Jend-Jstr)/2
         istep = i/2
         jstep = j/2
         
         do while (max(istep,jstep).gt.0)
            dlon = max(abs(lonr(i,j)-lonr(i+1,j)),
     &           abs(lonr(i,j)-lonr(i-1,j)))
            dlat = max(abs(latr(i,j)-latr(i,j+1)),
     &           abs(latr(i,j)-latr(i,j-1)))
            
            
            if (abs(lonr(i,j)-lon).gt.dlon) then
            
               if (lonr(i,j).lt.lon) then
c                  write(*,*) '|lon+:',lon,lonr(i,j),lonr(i+istep,j)
                  i = i+istep
               else
c                  write(*,*) '|lon-:',lon,lonr(i,j),lonr(i-istep,j)
                  i = i-istep
               end if
               istep = max(istep/2,1)
            else
               istep=0
            end if

            if (abs(latr(i,j)-lat).gt.dlat) then
               if (latr(i,j).lt.lat) then
c                  write(*,*) '|lat+:',lat,latr(i,j),latr(i,j+jstep)
                  j = j+jstep
               else
c                  write(*,*) '|lat-:',lat,latr(i,j),latr(i,j-jstep)
                  j = j-jstep
               end if
               jstep = max(jstep/2,1)
            else
               jstep = 0
            end if

c            write(*,*) 'distance after:',distance(lon,lat,lonr(i,j),
c     &           latr(i,j))

         end do

      end if
      end subroutine get_ij_2d_tile
