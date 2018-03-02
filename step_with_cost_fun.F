c     this file is not differentiated

# include "cppdefs.h"

      subroutine read_obs()
            implicit none
# include "netcdf.inc"
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# include "adparam.h"
      
      character*(*) FILE_NAME
      parameter (FILE_NAME='obs.nc')

      integer ncid, varid, retval, nf_fread

      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, 'zeta', varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_fread(obs(START_2D_ARRAY,1), ncid, varid, 0, r2dvar)

      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      end subroutine read_obs
      
      subroutine state_control(iicroot)
    
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "adparam.h"
# include "scalars.h"

      integer iicroot
      
C     params m1qn3
      
      integer indic,Sn,imp,io,imode(3),omode,niter,nsim,iz(5),ndz,
     &     reverse
      external simul_rc,euclid,ctonbe,ctcabe
      character*3 normtype

      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision dx, df1,epsrel,
     &     z(5*ad_array_size+1)

      parameter (Sn = ad_array_size,
     &     ndz=10*ad_array_size*ad_array_size+1)

      double precision f

      double precision gg(Sn)

      double precision dz(ndz)

      integer i

      sim_iicroot = iicroot
      
c     initial state: should call some user subroutine here
      do i=2,Sn
         x(i) = 60
      end do

      x(1) = rdrg2

      write(*,*) 'start state_control'
      open(unit=777, file="rdrg2.dat", action="write",
     &     form='FORMATTED')
      open(unit=778, file="h.dat", action="write",
     &     form='FORMATTED')

      call simul (indic,Sn,x,f,gg,izs,rzs,dzs)

      dx=1.e-10
      df1=f
      epsrel=1.e-5
      niter=200
      nsim=200
      normtype = 'dfn'
      io=6
      imp=5
      imode(1)=0
      imode(2)=0
      imode(3)=0
      reverse=1

      indic = 4

c     reverse mode 
      do while (reverse.ge.1)
         call simul(indic,Sn,x,f,gg,izs,rzs,dzs)
         MPI_master_only call m1qn3 (simul_rc,euclid,ctonbe,ctcabe,Sn,x,
     &        f,gg,dx,df1,epsrel,
     &        normtype,imp,io,imode,omode,niter,nsim,iz,dz,ndz,
     &        reverse,indic,izs,rzs,dzs)
      end do
C ...      
      call set_state(x)

      close (unit=777)
      close (unit=778)
      write(*,*) 'end state_control'
      
      return
      end subroutine state_control

      subroutine simul(indic,Sn,x,f,gg,izs,rzs,dzs)
      implicit none
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "adparam.h"
      
      integer indic, Sn
      double precision f, gg(Sn)
      real izs(1)
      real rzs(1)
      real dzs(1)

      double precision cost, costb

      if (indic.eq.1) then
         write(*,*) 'indic = 1'
      end if

      if (indic.eq.4) then
         write(*,*) 'indic = 4'

         call save_croco_state()

c     It is no trivial to get direct result from reverse mode with tapenade
c     cf Original statements and results disappearing from the Reverse code:
c     https://www-sop.inria.fr/tropics/tapenade/faq.html#adjlive
c     this need -nooptim adjointliveness but it is not sufficient,
c     and it needs some modifications in generated code.
c     Here we call directly the computation of cost_fun.
         call cost_fun(x, cost)
         f = cost
         call restore_croco_state()
         call save_croco_state()
         costb = 1
         call cost_fun_b(x, gg, cost, costb)
         call restore_croco_state()
         
      end if
      return
      end

      subroutine fun_val(iicroot)
# include "cppdefs.h"
# include "param.h"
# include "scalars.h"
# include "adparam.h"

      integer iicroot
      integer i, NVAL
      parameter (NVAL=60)
      real dx, cost, xv(NVAL), f(NVAL), gg(NVAL), costb

      sim_iicroot = iicroot
      
      dx = rdrg2 / NVAL
      
      x(1) = 0
      do i = 2,ad_array_size
         x(i)  = 60
      end do
      
      do i=1,NVAL
         xv(i) = 0.5*rdrg2 + i*dx
         x(1) = xv(i)
         call save_croco_state()
         call cost_fun(x, cost)
         f(i) = cost
         call restore_croco_state()
         call save_croco_state()
         costb = 1
         x(1) = xv(i)
         call cost_fun_b(x, gg(i), cost, costb)
         call restore_croco_state()
      end do

      do i=1,NVAL
         write(*,*) xv(i), f(i), gg(i)
      end do
      
      end subroutine fun_val

      subroutine save_croco_state()
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"

      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real h_bck(GLOBAL_2D_ARRAY)
      real rdrg2_bck

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

      common /backup/ ubar_bck, vbar_bck, zeta_bck, h_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, rdrg2_bck

      ubar_bck(:,:,:) = ubar(:,:,:)
      vbar_bck(:,:,:) = vbar(:,:,:)
      zeta_bck(:,:,:) = zeta(:,:,:)
      h_bck(:,:) = h(:,:)
      rdrg2_bck = rdrg2

      kstp_bck = kstp
      krhs_bck = krhs
      knew_bck = knew
      iic_bck = iic

      write(*,*) 'save------'
      write(*,*) 'knew=',knew
      write(*,*) 'kstp=',kstp
      write(*,*) 'krhs=',krhs
      write(*,*) 'iic=',iic
      write(*,*) '----------'

      end subroutine save_croco_state

      subroutine restore_croco_state()
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "adparam.h"

      real ubar_bck(GLOBAL_2D_ARRAY,4)
      real vbar_bck(GLOBAL_2D_ARRAY,4)
      real zeta_bck(GLOBAL_2D_ARRAY,4)
      real h_bck(GLOBAL_2D_ARRAY)
      real rdrg2_bck

      integer kstp_bck
      integer krhs_bck
      integer knew_bck
      integer iic_bck

      common /backup/ ubar_bck, vbar_bck, zeta_bck, h_bck,
     &     kstp_bck, krhs_bck, knew_bck, iic_bck, rdrg2_bck

      ubar(:,:,:) = ubar_bck(:,:,:)
      vbar(:,:,:) = vbar_bck(:,:,:)
      zeta(:,:,:) = zeta_bck(:,:,:)
      h(:,:) = h_bck(:,:)
      rdrg2 = rdrg2_bck

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
      
      subroutine handle_err(errcode)
      implicit none
# include "netcdf.inc"
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
