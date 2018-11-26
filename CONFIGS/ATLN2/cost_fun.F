
      subroutine cost_fun(ad_x, cost)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      double precision cost
      double precision xcost
      integer ta
      
      call set_state(ad_x)

      cost = 0

C$AD BINOMIAL-CKP ad_nt 40 1
      do ta=1,ad_nt-1
         call step()
         call cost_fun_step(cost, ta)
      end do

      write (*,*) 'Zob=',Zobt

      write (*,*) 'cost=',cost

      return
      end subroutine cost_fun

C$AD NOCHECKPOINT
      subroutine cost_fun_step(icost, ta)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      double precision icost
      integer ta
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

      do tile=my_first,my_last
         call cost_fun_step_2d(tile, icost, ta)
      end do

      end subroutine cost_fun_step

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d(tile, icost, ta)

      implicit none
# include "cppdefs.h"
# include "param.h"
      integer ta
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call cost_fun_step_2d_tile(Istr, Iend, Jstr, Jend, icost, ta)

      end subroutine cost_fun_step_2d

C$AD NOCHECKPOINT
      subroutine cost_fun_step_2d_tile(Istr,Iend,Jstr,Jend,icost,ta)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"      

      integer Istr,Iend,Jstr,Jend,i,j
      integer ta, k
      double precision delta_z, delta_zob, icost

      do k=1,ad_array_size/nnodes
         i=ad_i(k)
         j=ad_j(k)


         delta_z = (zeta(i,j,knew)-ad_obs(i,j,ta+1))
         delta_zob = (zob(i,j) - zob_bck(i,j))
         write(*,*) 'AAAAA',mynode,i,j,knew,sim_iicroot,k,sim_iicroot+k,
     &        zeta(i,j,knew),ad_obs(i,j,ta+1),zob(i,j),zob_bck(i,j)
         write(*,*) 'BBBBB',mynode,delta_z,delta_zob
         icost = icost + delta_z*delta_z + 1.D-5*delta_zob*delta_zob
      end do
      
      end subroutine cost_fun_step_2d_tile

C$AD NOCHECKPOINT
      subroutine set_state(ad_x)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
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

      do tile=my_first,my_last
         call set_state_2d(tile, ad_x)
      end do

      end subroutine set_state


C$AD NOCHECKPOINT
      subroutine set_state_2d(tile, ad_x)

      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"
      integer tile
# include "compute_tile_bounds.h"

      call set_state_2d_tile(Istr, Iend, Jstr, Jend, ad_x)

      end subroutine set_state_2d

C$AD NOCHECKPOINT
      subroutine set_state_2d_tile(Istr,Iend,Jstr,Jend,ad_x)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      integer Istr,Iend,Jstr,Jend, i,j,k
      double precision lzob

      integer ad_anode_size

#include "compute_extended_bounds.h"
      ad_anode_size=ad_array_size/nnodes

      do j=Jstr,Jend
         do i=Istr,Iend
c            k=ad_colloc(i,j)
c            if (k.lt.1.or.k.gt.ad_array_size) then
c               write (*,*) 'AD BIG PROBLEM k=',k
c               stop
c            endif
            zob(i,j)=0

            do k=1,ad_anode_size
               lzob=1.D-7 + 1/(1+exp(ad_x(k)))
               zob(i,j) = zob(i,j) + lzob*W(i,j,k)/SkW(i,j)
            end do
c            zob(i,j)=1.D-7 + 1/(1+exp(ad_x(k)))
         end do
      end do
      
      write(*, *) '**** set_state ****'

#ifdef MPI
      write(*, *) 'mynode=',mynode
#endif
      write(*, *) 'min zob=',minval(zob)
      write(*, *) 'max zob=',maxval(zob)
      write(*, *) '*******************'

 900  format(1PE16.9)

      return
      end subroutine set_state_2d_tile


      subroutine init_control()
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"      
# ifdef MPI
# include "mpi_roms.h"
# ifdef AMPI
# include "ampi/ampif.h"
# else
      include 'mpif.h'
# endif
# endif
      integer ierr

      
      integer k
      
      do k=1,ad_array_size/nnodes
         ad_x(k)=log(1./Zobt - 1.)
         ad_g(k)=0.
      end do
      
      call init_collocation_points()
      
      end subroutine init_control

      subroutine save_croco_state()
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      ubar_bck(:,:,:) = ubar(:,:,:)
      vbar_bck(:,:,:) = vbar(:,:,:)
      zeta_bck(:,:,:) = zeta(:,:,:)
      zob_bck(:,:) = zob(:,:)
      zobt_bck = zobt

      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) 'minzob',minval(zob_bck)
      write(*,*) 'maxzob',maxval(zob_bck)
      write(*,*) '----------'

      end subroutine save_croco_state
      
      subroutine restore_croco_state()
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"

      ubar(:,:,:) = ubar_bck(:,:,:)
      vbar(:,:,:) = vbar_bck(:,:,:)
      zeta(:,:,:) = zeta_bck(:,:,:)
      zob(:,:) = zob_bck(:,:)
      Zobt = Zobt_bck

      kstp = kstp_bck
      krhs = krhs_bck
      knew = knew_bck
      iic = iic_bck

      write(*,*) 'restore+++'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '++++++++++'

      end subroutine restore_croco_state

c from rosettacode.org 
      function to_radian(degree) result(rad)
          ! degrees to radians
          real,intent(in) :: degree
          real, parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
          real :: rad
 
          rad = degree*deg_to_rad
      end function to_radian
 
      function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
          ! great circle distance -- adapted from Matlab 
          real,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real :: a,c,dist,dlat,dlon,lat1,lat2
          real,parameter :: radius = 6372.8 
 
          dlat = to_radian(deglat2-deglat1)
          dlon = to_radian(deglon2-deglon1)
          lat1 = to_radian(deglat1)
          lat2 = to_radian(deglat2)
          a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
          c = 2*asin(sqrt(a))
          dist = radius*c
      end function haversine

      subroutine init_collocation_points()
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
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

      do tile=my_first,my_last
         call init_collocation_points_2d(tile)
      end do      

      end subroutine init_collocation_points

      subroutine init_collocation_points_2d(tile)

      implicit none
# include "cppdefs.h"
# include "param.h"
      integer k
      integer tile
      double precision icost
# include "compute_tile_bounds.h"

      call init_collocation_points_2d_tile(Istr, Iend, Jstr, Jend)

      end subroutine init_collocation_points_2d

      subroutine init_collocation_points_2d_tile(Istr,Iend,Jstr,Jend)

      implicit none
# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# include "adparam.h"      

      integer Istr,Iend,Jstr,Jend,i,j,k,istep,jstep
      double precision mindist,dist,rk,haversine

      integer ad_anode_size

      ad_anode_size = ad_array_size/nnodes

      
      istep=(Iend+1-Istr)/(2*sqrt(real(ad_array_size)))
      jstep=(Iend+1-Istr)/(2*sqrt(real(ad_array_size)))
      k=1
      jloop: do j=Jstr,Jend,istep
      do i=Istr,Iend,jstep
c         call random_number(r)
c         i=Istr+floor((Iend+1-Istr)*r)
c     j=Jstr+floor((Jend+1-Jstr)*r)
         if (h(i,j).gt.0.1.and.ad_obs(i,j,1).ne.spval) then
               write(*,*) 'new collocation point:',mynode,i,j,h(i,j),
     &           ad_obs(i,j,1)
               ad_i(k)=i
               ad_j(k)=j
               k=k+1
               if (k.gt.ad_anode_size) then
                  exit jloop
               endif
            endif
         enddo
      enddo jloop
      
      
      do j=Jstr,Jend
         do i=Istr,Iend
            SkW(i,j)=0.
            do k=1,ad_anode_size
               dist=haversine(latr(i,j),lonr(i,j),latr(ad_i(k),ad_j(k)),
     &              lonr(ad_i(k),ad_j(k)))
               rk=max(0.01*sqrt(g*h(i,j))*TM2,10000.)
               W(i,j,k)=exp(-dist/rk)
               SkW(i,j)=SkW(i,j)+W(i,j,k)
            end do
         end do
      end do

      do j=Jstr,Jend
         do i=Istr,Iend
            mindist=1.D230
            do k=1,ad_anode_size
               dist=haversine(latr(i,j),lonr(i,j),latr(ad_i(k),ad_j(k)),
     &              lonr(ad_i(k),ad_j(k)))
               if (dist.lt.mindist) then
                  ad_colloc(i,j)=k
               endif
            enddo
         enddo
      enddo
      
      end subroutine init_collocation_points_2d_tile

