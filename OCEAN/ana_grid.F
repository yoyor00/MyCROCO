! $Id: ana_grid.F 1620 2015-01-08 10:47:13Z marchesiello $
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
#ifdef ANA_GRID
      subroutine ana_grid (tile)
      implicit none
# include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      call ana_grid_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine ana_grid_tile (Istr,Iend,Jstr,Jend)
!
!  Set up model grid using analytical expressions:
!-------------------------------------------------------------------
!
!  INPUT: Grid configuration parameters
!  ------------------------------------
! Length_XI      are the physical dimensions of the computational
! Length_ETA     domain [usually measured in  meters];
!
! depth          is the maximum depth [meters, positive];
! f0,beta        are coriolis parameters which set up a beta-plane
!                                                  [1/s, 1/(m*s)].
!
!  OUTPUT: stored in common blocks, see files "scalars.h" "grid.h"
!  ---------------------------------------------------------------
!  xl,el      Physical dimensions of the computational domain 
!                                    [usually measured in  meters];  
!  h          Model bathymetry [meters, positive] at RHO-points.
!  hmin,hmax  Minimum and maximum values of depth of bathymetry [m].
!  f          Coriolis parameter (1/seconds) at RHO-points.
!  pm,pn      Coordinate transformation metric "m" [1/meters]
!                    associated with the differential distances in 
!                    XI- and ETA-directions, both are at RHO-points.
!  xp,xr      XI-coordinates [m] at PSI- and RHO-points.
!  yp,yr      ETA-coordinates [m] at PSI- and RHO-points.
!-------------------------------------------------------------------
!
# ifdef AGRIF
      use Agrif_Util
# endif
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
# include "param.h"
# include "grid.h"
# include "scalars.h"
!
      real Length_XI, Length_ETA, depth, f0, beta, 
     &     cff, y, x0, y0, dx, dy, NSUB
!
# ifdef DUNE 
#   ifdef DUNE3D
      real Rij,cffi,cffj,cff1
#   elif defined ANA_DUNE
      real cff2
#   endif
# endif
# ifdef TIDAL_FLAT
      real depth_shallow,depth_deep
# endif
# ifdef INTERNAL
      real cff1, ridge_width, ridge_height
# endif
# ifdef RIP
      real h0,per,xx,yy,xs,db,alpha,lambda,eps
# endif
# ifdef SANDBAR
      real p1,p2,p3,p4,p5,p6,p7,xx
# endif
# ifdef RIVER
      integer is,i0,j0
#  include "sources.h"
# endif
#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
!
# include "compute_extended_bounds.h"
!
!----------------------------------------------------------------------
! Configuration parameters (h,f,beta)
!----------------------------------------------------------------------
!
# if defined BASIN
                     depth=5000.
                     f0=1.E-4 
                     beta=2.E-11
# elif defined SINGLE_COLUMN
#  if defined KATO_PHILIPS || defined WILLIS_DEARDORFF 
                     depth=50.
                     f0=0.
                     beta=0.
#  elif defined DIURNAL_CYCLE
                     depth=150.
                     f0=0.
                     beta=0.
#  elif defined FORCED_EKBBL
                     depth=1500.
                     f0=1.E-4
                     beta=0.
#  elif defined FORCED_DBLEEK
                     depth=25.0
                     f0=1.E-4
                     beta=0.
#  elif defined FORCED_NONROTBBL || defined FORCED_OSCNONROTBBL
                     depth=5.
                     f0=0.0
                     beta=0.
#  endif 
# elif defined CANYON
                     depth=4000.
                     f0=1.E-4
                     beta=0.
# elif defined PLUME
                     depth=1300.
                     beta=0.
                     f0=0.5E-4 
# elif defined EQUATOR
                     depth=5000.
                     f0=0.
                     beta=2.2829e-11
# elif defined KH_INST 
                     depth=256.     
                     f0=0.
                     beta=0.
# elif defined GRAV_ADJ
#  ifdef NBQ
                     depth=0.3
#  else
                     depth=20.
#  endif
                     f0=0.
                     beta=0.
# elif defined ISOLITON
                     depth=0.29 ! Horn et al. (2001)
                     f0=0.
                     beta=0.
# elif defined ACOUSTIC
                     depth=128.
                     f0=0.
                     beta=0.
# elif defined INNERSHELF
                     depth=0.
                     f0=4*pi/86400*sin(-21*pi/180)
                     beta=0.
# elif defined INTERNAL
                     depth=2000.
                     f0=1.e-4
                     beta=0.
                     ridge_width=30.e3
                     ridge_height=1600.
# elif defined TS_HADV_TEST
                     depth=10.
                     f0=0.
                     beta=0.                     
# elif defined OVERFLOW
                     depth=40.
                     f0=0.
                     beta=0.
# elif defined SEAMOUNT
                     depth=4500.
                     beta=0.
                     f0=1.E-4
# elif defined SHELFRONT
                     depth=1660.
                     f0=1.E-4
                     beta=0.
# elif defined SOLITON
                     depth=1.
                     f0=0.
                     beta=1.
# elif defined UPWELLING
                     depth=150.
                     f0=-8.26E-5
                     beta=0.
# elif defined RIVER
                     depth=150.
                     f0=8.26E-5
                     beta=0.
# elif defined JET
                     depth=4000.
                     f0=1.E-4
                     beta=1.6E-11
# elif defined SHOREFACE
                     depth=15.
                     f0=0.   
                     beta=0.
# elif defined SANDBAR
                     depth=4.1
                     f0=0.   
                     beta=0.
# elif defined SWASH
                     depth=1.
                     f0=0.   
                     beta=0.
# elif defined RIP
                     depth=15.
                     f0=0.   
                     beta=0.
# elif defined THACKER
                     depth=10.
                     beta=0.
#  ifdef THACKER_2DV
                     f0=0.
#  else
                     f0=1.E-4
#  endif
# elif defined TANK
                     depth=10.
                     f0=0.
                     beta=0.
# elif defined MOVING_BATHY
                     depth=0.394
                     f0=0.
                     beta=0.
# elif defined DUNE
#  ifdef ANA_DUNE
                     depth=6.                 
#  else
                     depth=4.
#  endif
                     f0=0.
                     beta=0.
# elif defined SED_TOY
                     depth=5.
                     f0=1.E-4
                     beta=0.
# elif defined TIDAL_FLAT
                     depth_shallow=-4.
                     depth_deep=16.
                     f0=0.
                     beta=0.
# else
                     depth=???
                     f0=???
                     beta=???
# endif
!
!----------------------------------------------------------------------
! Grid dimensions (Length_XI, Length_ETA)
!----------------------------------------------------------------------
!
# define Length_XI_ISO_GRID Length_ETA*float(LLm0)/float(MMm0)
# define Length_ETA_ISO_GRID Length_XI*float(MMm0)/float(LLm0)
!
# ifdef AGRIF
      if (Agrif_Root()) then
# endif
# if defined BASIN
                     Length_XI =3600.0e+3
                     Length_ETA=2800.0e+3
# elif defined SINGLE_COLUMN
                     Length_XI = 10.0e+3
                     Length_ETA= 10.0e+3
# elif defined CANYON
                     Length_XI =128.0e+3
                     Length_ETA= 96.0e+3
# elif defined EQUATOR
                     Length_ETA=3200.0e+3
                     Length_XI =4000.0e+3
# elif defined KH_INST 
                     Length_XI =256.0
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined ACOUSTIC
                     Length_XI =128.0
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined GRAV_ADJ
#  ifdef NBQ
                     Length_XI =3.0
#  else
                     Length_XI =64.0e+3
#  endif
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined ISOLITON
                     Length_XI =6.0
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined INNERSHELF
                     Length_XI =200.0e+3
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined INTERNAL
                     Length_XI =1200.0e+3
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined TS_HADV_TEST
                     Length_XI = 1.0e+3
                     Length_ETA= 1.0e+3
# elif defined OVERFLOW
                     Length_ETA=64.0E+3
                     Length_XI =Length_XI_ISO_GRID
# elif defined SEAMOUNT
                     Length_XI =512.0e+3
                     Length_ETA=512.0e+3
# elif defined SHELFRONT
                     Length_ETA=200.0e+3
                     Length_XI =Length_XI_ISO_GRID
# elif defined SOLITON
                     Length_XI =48.
                     Length_ETA=16.
# elif defined UPWELLING
                     Length_ETA=80.e+3
                     Length_XI =Length_XI_ISO_GRID
# elif defined RIVER
                     Length_ETA=80.e+3
                     Length_XI =Length_XI_ISO_GRID
# elif defined JET
                     Length_XI = 500.e+3 
                     Length_ETA=2000.e+3 
# elif defined SHOREFACE
                     Length_XI =1180.0 
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined SANDBAR
                     Length_XI =200.
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined SWASH
                     Length_XI =110.0 
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined RIP
                     Length_XI =768.0 
                     Length_ETA=768.0
# elif defined THACKER
                     Length_XI =200.e+3 
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined TANK
#  ifndef TANKY
                     Length_XI =10.0
                     Length_ETA=Length_ETA_ISO_GRID
#  else
                     Length_ETA=10.0
                     Length_XI =Length_XI_ISO_GRID
#  endif
# elif defined MOVING_BATHY
                     Length_XI =4.0
                     Length_ETA=Length_ETA_ISO_GRID
# elif defined DUNE
#  ifdef ANA_DUNE
                     Length_XI = 300.0
#  else
                     Length_XI = 100.0
#  endif
#  ifdef MUSTANG
                     Length_ETA=Length_XI
#  else                     
                     Length_ETA=Length_ETA_ISO_GRID
#  endif
# elif defined SED_TOY
                     Length_XI = 10.e+3
                     Length_ETA= 10.e+3
# elif defined TIDAL_FLAT
                     Length_ETA= 1500.0
                     Length_XI = 100000.0
# else
                     Length_XI =???
                     Length_ETA=??? 
# endif

# undef Length_XI_ISO_GRID
# undef Length_ETA_ISO_GRID

!----------------------------------------------------------------------
! Copy physical dimensions of the grid into globally vizible variables
!----------------------------------------------------------------------
!
        xl=Length_XI
        el=Length_ETA

# ifdef AGRIF
      else      
        xl=Agrif_Parent(xl)/Agrif_Parent(LLm)*LLm/Agrif_Rhox()
        el=Agrif_Parent(el)/Agrif_Parent(MMm)*MMm/Agrif_Rhoy()
        Length_XI=xl
        Length_ETA=el
      endif
# endif /* AGRIF */
!
!
!----------------------------------------------------------------------
! Set grid spacings for rectangular grids
!----------------------------------------------------------------------
!
      dx=Length_XI/float(LLm)    
      dy=Length_ETA/float(MMm)   
!
!----------------------------------------------------------------------
! Set reference point location
!----------------------------------------------------------------------
!
      x0=0.
      y0=0.
# ifdef EQUATOR
      y0=y0-0.5*Length_ETA  
# endif
# if defined INTERNAL || defined THACKER || defined MOVING_BATHY
      x0=x0-0.5*Length_XI
      y0=y0-0.5*Length_ETA  
# endif

CR      write(*,'(4(A3,pe15.9,1x),I3)') 'dx=',dx, 'dy=',dy,
CR     &                                'x0=',x0, 'y0=',y0, mynode
!
!----------------------------------------------------------------------
! Setup rectangular grid:  coordinates (XI,ETA) at PSI- and RHO-points.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
         xp(i,j) = x0 + dx*float(i+iminmpi-1-1)
         xr(i,j) = x0 + dx*(float(i+iminmpi-1-1)+0.5)
         yp(i,j) = y0 + dy*float(j+jminmpi-1-1)
         yr(i,j) = y0 + dy*(float(j+jminmpi-1-1)+0.5)
        enddo
      enddo
!
!----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          pm(i,j)=1./dx
          pn(i,j)=1./dy
        enddo
      enddo
!
!----------------------------------------------------------------------
!  Set Coriolis parameter [1/s] at RHO-points.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef EQUATOR
          f(i,j)=f0+beta*yr(i,j)
#  else
          f(i,j)=f0+beta*(yr(i,j)-(0.5*el))
#  endif
        enddo
      enddo

# if defined UV_COR_NT || defined CROCO_QH
!
!  Horizontal (non-traditional) Coriolis parameter
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          if (beta.eq.0) then     ! f-plan
            e(i,j)=2.*Erotation*cos(asin(f0/(2.*Erotation)))
          else
            e(i,j)=beta*Eradius   ! beta plan
          endif
        enddo
      enddo
# endif
!
!----------------------------------------------------------------------
!  Set bathymetry [meters; positive] at RHO-points.
!----------------------------------------------------------------------
!
# if defined BASIN || defined EQUATOR || defined GRAV_ADJ \
                   || defined SOLITON || defined ISOLITON \
                   || defined KH_INST || defined TS_HADV_TEST \
                   || defined ACOUSTIC || defined SINGLE_COLUMN \
                   || defined JET
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth
        enddo
      enddo

# elif defined CANYON
      do j=JstrR,JendR
        do i=IstrR,IendR
          cff=32000.-16000.*(sin(pi*xr(i,j)/Length_XI))**24
          h(i,j)=20.+0.5*(depth-20.)*(1.+tanh((yr(i,j)-cff)/10000.))
        enddo
      enddo

# elif defined INNERSHELF
      do j=JstrR,JendR
        do i=IstrR,IendR
!         --- constant slope ---
          h(i,j)=1.e-3*(xl+0.5*dx-xr(i,j))+4.
!         --- shelf/slope tanh ---
!          h(i,j)=30.+500.*(1+sinh((150.e3-xr(i,j))/40.e3)
!     &                      /cosh((150.e3-xr(i,j))/40.e3))
        enddo
      enddo

# elif defined INTERNAL
      cff1=1./(ridge_width*ridge_width)
      do j=JstrR,JendR
        do i=IstrR,IendR
          if ((xr(i,j).gt.ridge_width).or.
     &       (xr(i,j).lt.(-ridge_width))) then
            cff=0.
          else
            cff=1-(xr(i,j)*xr(i,j)*cff1);
          endif
          h(i,j)=depth-ridge_height*cff*cff;
#   ifdef INTERNALSHELF
          if (xr(i,j).ge.0.) then
            h(i,j)=depth-ridge_height
          endif
#   endif
        enddo
      enddo

# elif defined OVERFLOW
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=20.+0.5*(depth-20.)*(1.+
     &                         tanh((yr(i,j)-40000.)/5000.))
        enddo
      enddo

# elif defined SEAMOUNT
      cff=(1./50.e3)**2
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth*(1-0.6*exp(-cff*((xr(i,j)-0.5*xl)**2+
     &                                  (yr(i,j)-0.5*el)**2))) 
        enddo
      enddo

# elif defined SHELFRONT
      do j=JstrR,JendR
        do i=IstrR,IendR
          cff=yr(i,j)/1000.
          if (cff.lt.50.) then
            h(i,j)=50.+2.*cff
          elseif (cff.lt.60.) then
            h(i,j)=160.+1.5*(cff-50.)**2-0.1*(cff-60.0)**2
          elseif (cff.lt.100.) then
            h(i,j)=310.+30.*(cff-60.)
          elseif (cff.lt.110.) then
            h(i,j)=1660.-1.5*(cff-110.)**2
          else
            h(i,j)=1660.
          endif
        enddo
      enddo

# elif defined UPWELLING
      do j=JstrR,JendR
#  ifdef MPI
        y=dy*float(j+jj*Mm)
#  else
        y=dy*float(j)
#  endif
        if (y.gt.Length_ETA/2.) y=Length_ETA-y+dy
        cff=min(depth,84.5+66.526*tanh(0.00015*(y-0.125*Length_ETA)))
        do i=IstrR,IendR
          h(i,j)=cff
        enddo
      enddo

# elif defined RIVER
      do i=IstrR,IendR
#  ifdef MPI
        cff=(float(i +ii*Lm)-0.5)/float(LLm)
#  else
        cff=(float(i       )-0.5)/float(LLm)
#  endif
        if (cff.lt.0.05) then
          h(i,JstrR)=15.
        elseif (cff.lt.0.15) then
          h(i,JstrR)=15.+843.75*(cff-0.05)**2
        elseif (cff.lt.0.85) then
          h(i,JstrR)=15.+168.75*(cff-0.1)
        elseif (cff.lt.0.95) then
          h(i,JstrR)=150.-843.75*(cff-0.95)**2
        else
          h(i,JstrR)=150.
        endif
      enddo
      do j=JstrR+1,JendR
        do i=IstrR,IendR
          h(i,j)=h(i,JstrR)
        enddo
      enddo

# elif defined SHOREFACE
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
#  ifdef WET_DRY
          h(i,j)=12.0-0.0125*Length_XI/float(LLm+1)*cff ! original fonction
#  else
          h(i,j)=12.0-0.0125*Length_XI/float(LLm+1)*cff + 3.
#  endif
#  ifdef MASKING
          rmask(i,j)=1.0
          if (h(i,j).lt.0.01) then
            h(i,j)=0.01
            rmask(i,j)=0.0
            rmask(48,j)=1.0
          endif
#  endif
        enddo
      enddo

# elif defined SANDBAR
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
#  ifdef SANDBAR_OFFSHORE
! Roelvink & Reniers 1995: LIP-1B
          xx=xr(i,j)
          if (xx.lt.19) then
            h(i,j)=-depth
          elseif (xx.lt.52) then
            p1=+1.520229705740388e-08
            p2=-3.423212989640106e-06
            p3=+3.159005532947043e-04
            p4=-1.526398548035649e-02
            p5=+4.066341504316825e-01
            p6=-5.608106604351117e+00
            p7=+2.702732753680607e+01
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.120) then
            p1=-4.642289413679605e-11
            p2=+2.870629635648417e-08
            p3=-7.130431368316412e-06
            p4=+9.126154058771732e-04
            p5=-6.357269253701968e-02
            p6=+2.303201222063785e+00
            p7=-3.657085297281925e+01
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.140) then
            p1=+3.099533433018596e-08
            p2=-2.425187846848144e-05
            p3=+7.897142871067360e-03
            p4=-1.369896736259833e+00
            p5=+1.335157006848632e+02
            p6=-6.932492254089954e+03
            p7=+1.498142545396247e+05
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.160) then
            p1=-2.524608914887271e-08
            p2=+2.387297325503284e-05
            p3=-9.398023259277126e-03
            p4=+1.971403061349242e+00
            p5=-2.323952920877027e+02
            p6=+1.459649630832179e+04
            p7=-3.816025205750951e+05
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          else
            p1=+3.693458909496914e-09
            p2=-4.051059061933883e-06
            p3=+1.847329725370577e-03
            p4=-4.483034673182216e-01
            p5=+6.106360892732523e+01
            p6=-4.426507173590056e+03
            p7=+1.334155898180324e+05
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          endif
#  elif defined SANDBAR_ONSHORE
! Roelvink & Reniers 1995: LIP-1C
          xx=xr(i,j)
          if (xx.lt.19) then
            h(i,j)=-depth
          elseif (xx.lt.52) then
            p1=-1.550500196449383e-09
            p2=+3.095246505454406e-07
            p3=-2.524513106364847e-05
            p4=+1.072358697579965e-03
            p5=-2.475915084067423e-02
            p6=+3.341039769921187e-01
            p7=-6.243997999730277e+00
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.120) then
            p1=+2.879040803660182e-11
            p2=-3.112995372255505e-09
            p3=-2.128994193906542e-06
            p4=+5.806818410041643e-04
            p5=-5.889240225358890e-02
            p6=+2.719179108579184e+00
            p7=-5.021582880671649e+01
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.140) then
            p1=+6.722443006700593e-07
            p2=-5.191501807411085e-04
            p3=+1.669215462536683e-01
            p4=-2.860227732993812e+01
            p5=+2.754780791103959e+03
            p6=-1.414013613750512e+05
            p7=+3.022008907301721e+06
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          elseif (xx.lt.160) then
            p1=-3.300743083271716e-07
            p2=+2.941386493481447e-04
            p3=-1.091628962526408e-01
            p4=+2.159706871907915e+01
            p5=-2.402365333342808e+03
            p6=+1.424582143575532e+05
            p7=-3.518304916201898e+06
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          else
            p1=+3.341378689344706e-09
            p2=-3.683115456141188e-06
            p3=+1.686167358885289e-03
            p4=-4.104224500606568e-01
            p5=+5.602413257412504e+01
            p6=-4.066824640330583e+03
            p7=+1.226609757884486e+05
            h(i,j) = p1*xx**6 + p2*xx**5 + p3*xx**4 +
     &               p4*xx**3 + p5*xx**2 + p6*xx    + p7
          endif
#  endif
          h(i,j) = min(-h(i,j),depth)
        enddo
      enddo

# elif defined SWASH
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
          h(i,j)=1.0571-0.0125*Length_XI/float(LLm+1)*cff  ! GLOBEX profil
#  ifndef WET_DRY
     &                                              + 0.6
#  endif
          if (xr(i,j).lt.16.57) h(i,j)=0.85               ! GLOBEX
        enddo
      enddo

# elif defined RIP
#  ifdef GRANDPOPO
!
!  Idealized Grand Popo Beach in Benin
!
      xs=85;        ! inner surf zone 
      db=50;        ! distance from xs to sand bar
      alpha=0.01;
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef OBC_EAST
          xx=xr(i,j)
#  else
          xx=Length_XI-xr(i,j)
#  endif
          h(i,j)=-4.5 
     &       -1.7*exp(-6*(((xx-xs-db)/db)**2))
     &       +3.1*(1+tanh(0.025*(xx-xs)))
     &       +0.014*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha)
        enddo
      enddo
#  else
!
! Weir et al. 2011 after Lippmann et al. 1999 (idealization of Duck beach)
!
      xs=150      ! inner surf zone 
      db=80       ! distance from xs to sand bar
      alpha=0.02
      lambda=256  ! sand bar perturbation wavelength
#  ifdef RIP_TOPO_2D
      eps=0.0     ! no sand bar perturbation
#  else
      eps=0.1     ! sand bar perturbation magnitude
#  endif
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef OBC_EAST
          xx=xr(i,j)
#  else
          xx=Length_XI-xr(i,j)
#  endif
          yy=yr(i,j);
          h0=-1.5*exp(-5*(((xx-xs-db)/db)**2))
     &       +1.35*(1+tanh(0.0253*(xx-xs)))
     &       +0.0032*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha)
          per=eps*cos(2*pi*yy/lambda)*exp(-5*(((xx-xs-db)/db)**2))
          h(i,j)=(1+per)*h0!-0.8;
        enddo
      enddo
#  endif

# elif defined THACKER
      cff=80.e3
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef THACKER_2DV
          h(i,j)=depth*(1-(xr(i,j)/cff)**2)
#  else
          h(i,j)=depth*(1-(xr(i,j)/cff)**2 -
     &                    (yr(i,j)/cff)**2 )
#  endif
        enddo
      enddo

# elif defined TANK
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth
        enddo
      enddo
# elif defined MOVING_BATHY
      cff=1./(0.03686)**2
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth-0.1*exp(-cff*xr(i,j)**2)
        enddo
      enddo

# elif defined DUNE
#  ifdef ANA_DUNE
      cff2=0.01
#  elif defined DUNE3D
      cff1=0.5*(1./10)**2
#  endif
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
         cff=(float(i +ii*Lm) )
#  else  
         cff=(float(i       ) )
#  endif
#  ifdef ANA_DUNE
         h(i,j)=depth-2.*exp(-cff2*(xr(i,j)-0.5*xl)**2)
#  elif defined DUNE3D
#    ifdef MPI
          cffi=(float(i +ii*Lm) )
          cffj=(float(j +jj*Mm) )
#    else   
          cffi=(float(i       ) )
          cffj=(float(j       ) )
#    endif 
         if (xr(i,j).lt.4.) then
          h(i,j)=depth
         else 
          Rij = sqrt( (cffi-28.)**2.+
     &                (cffj-25.)**2. )
          h(i,j)=depth*(1.-0.5*exp(-cff1*(2.*(Rij))**2.0))
         endif
#  else 
         h(i,j)=depth*(1.-0.5*exp(-5.e-3*(2.*(cff-25.))**2))
#  endif
        enddo
      enddo

# elif defined SED_TOY
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth
!          rmask(i,j)=1.0
        enddo
      enddo

# elif defined TIDAL_FLAT
      do j=JstrR,JendR
        h(IstrR,j)=depth_deep
        h(IstrR+1,j)=depth_deep
        rmask(IstrR,j)=1.0
        rmask(IstrR+1,j)=1.0
        do i=IstrR+2,IendR
          h(i,j)=depth_deep-REAL(i+iminmpi-2)/REAL(LLm)
     &                 *(depth_deep-depth_shallow)
          rmask(i,j)=1.0
        enddo
      enddo
      rmask(IstrR:IendR,JstrR:1)=0.0
      rmask(IstrR:IendR,3:JendR)=0.0
      rmask(IendR-1:IendR,2)=0.0

# else
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=???
        enddo
      enddo
# endif

!----------------------------------------------------------------------
! initialisation of z0b in case of analytical cases
!----------------------------------------------------------------------
      do j=JstrR,JendR
        do i=IstrR,IendR
          zob(i,j)=zobt
        end do
      end do
!
!----------------------------------------------------------------------
!  Set masking at RHO-points.
!----------------------------------------------------------------------
!
# ifdef MASKING
#  ifdef RIVER
      do j=JstrR,JendR             ! Set mask to all-water status,
        do i=IstrR,IendR           ! then mask out 8-point wide strip
          rmask(i,j)=1.            ! of land on the west, and finally
        enddo                      ! carve two 1-point wide channels
      enddo                        ! through that strip.
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          do i=IstrR,8
            rmask(i,j)=0.            ! <-- strip of land
          enddo
        enddo
        do is=1,Nsrc
#  ifdef MPI
          i0=Isrc_mpi(is,mynode)
          j0=Jsrc_mpi(is,mynode)
#  else
          i0=Isrc(is)
          j0=Jsrc(is)
#  endif
          if (IstrR.le.i0 .and. i0.le.IendR .and.
     &        JstrR.le.j0 .and. j0.le.JendR) then
            if (is.eq.1) then
              do j=j0,j0+6
                rmask(3,j)=1.      ! <-- upper channel with corner
              enddo
              do i=IstrR+3,8
                rmask(i,j0+6)=1.
              enddo
            else
              do i=IstrR+3,8
                rmask(i,j0)=1.     ! <-- lower channel along xi
              enddo
            endif
          endif
        enddo
      endif ! WESTERN_EDGE
#  endif
# endif
!
#else
      subroutine ana_grid_empty
#endif /* ANA_GRID */
      return
      end

