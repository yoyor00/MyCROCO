! $Id: get_tides.F 1458 2014-02-03 15:01:25Z gcambon $
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
                              ! Read gridded tidal parameters
      subroutine get_tides    ! from forcing netCDF file. 
                              ! Scale them. 
#if defined SSH_TIDES || defined UV_TIDES
      USE MOD_TIDES_MAS
      implicit none

#ifndef TIDES_MAS
# include "param.h"
#endif
# include "scalars.h"
# include "ncscrum.h"
# include "tides.h"
# include "grid.h"
      real cff
      integer i,j, itide, ierr, lstr,lvar, lenstr, nf_fread, 
     &        status, varid,ind
#ifdef TIDES_MAS
      integer :: nn,kount,k,l
      integer :: xt
#endif

#ifdef MPI
#define LOCALLM Lmmpi
#define LOCALMM Mmmpi
#else
#define LOCALLM Lm
#define LOCALMM Mm
#endif    
#include "netcdf.inc"

# ifdef TIDES_MAS
! initialisation du tableau des constituants harmoniques
!-------------------------------------------------
      call tide_data
       nu(1)=im2
       nu(2)=ik1
       nu(3)=ik2
       nu(4)=i2n22nm2
       nu(5)=im4
       nu(6)=imf
       nu(7)=imm
       nu(8)=in2
       nu(9)=io1
       nu(10)=ip1
       nu(11)=iq1
       nu(12)=is2
       nu(13)=im6
       nu(14)=imn4
       nu(15)=ims4
       nu(16)=imk4


! on reordonne les ondes pour decrire le spectre de facon compatible
! avec la liste SIMON
!-------------------------------------------------
!
      kount=1
      do while(kount/=0)
       kount=0
       do nn=1,Ntides-1
        if(nu(nn) > nu(nn+1)) then
         kount=kount+1
         xt=nu(nn)
         nu(nn)=nu(nn+1)
         nu(nn+1)=xt
        end if
       end do
      end do
      l_presence(:)=.false.
      l_presence(nu(1:Ntides))=.true.

      ! on prend que M2
      !l_presence(nu(1))=.true.
      ! on prend que M2
!      l_presence(im2)=.true.
!      l_presence(is2)=.true.
!      l_presence(ik2)=.true.
!      l_presence(in2)=.true.
!      l_presence(io1)=.true.
!      l_presence(is1)=.true.
!      l_presence(ip1)=.true.
!      l_presence(ik1)=.true.
     

      rr0(53)=-1000
      rr1(54)=32


# endif
!
! Inquire about the contents of forcing netCDF file:
! variables and dimensions. Check for consistency.
!
      if (may_day_flag.ne.0) return      !-->  EXIT
      lstr=lenstr(frcname)
!
! If not opened yet, open forcing netCDF file for reading. Find and
! save IDs for relevant variables, determine whether momentum stress
! components exist as fields or scalars.
!
      if (ncidfrc.eq.-1) then
        ierr=nf_open (frcname(1:lstr), nf_nowrite, ncidfrc)
        if (ierr .ne. nf_noerr) goto 4                 !--> ERROR 
      endif

#ifdef TIDES_MAS
!test
        status=nf_inq_varid(ncidfrc,'tide_period',varid)
        if (status.ne.NF_NOERR) then
          write(6,3) 'id_tide_per',frcname(1:lstr)
          STOP
        endif
        status=nf_get_var(ncidfrc,varid,Tperiod)
        MPI_master_only write(stdout,*)'TIDE PERIODS',Tperiod
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_per',frcname(1:lstr)
          STOP
        endif
#endif

      do itide=1,Ntides
#ifdef TIDES_MAS
      ! search for the corresponding tide
      Do l=1,Ntides
       MPI_master_only write(stdout,*)'tide ',trim(nommar(nu(itide))),
     &        itide,l,nu(itide),Tperiod(l),360./frequence(nu(itide))
        if (abs(Tperiod(l)-360./frequence(nu(itide))) .le. 1e-4  ) then
          ind=l
          exit
       endif
      end do  
#else
      ind=itide
#endif
# ifndef TIDES_MAS
!
!  Tidal Period.
!
        status=nf_inq_varid(ncidfrc,'tide_period',varid)
        status=nf_get_var1_FTYPE(ncidfrc,varid,itide,Tperiod(itide))
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_period',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        Tperiod(itide)=Tperiod(itide)*3600.
#endif

# ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase.
!
        status=nf_inq_varid(ncidfrc,'tide_Eamp',varid)
        status=nf_fread(SSH_Tamp(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Eamp',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

#ifdef TIDES_MAS_DBG
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
        if (l_presence(nu(itide))
     &    .and. i+iminmpi-1 .eq. 1 .and. j+jminmpi-1 .eq. 149)
     &         write(stdout,*),'SSH_TAMP',h(i,j),SSH_Tamp(i,j,itide)
          enddo
        enddo
#endif

        status=nf_inq_varid(ncidfrc,'tide_Ephase',varid)
        status=nf_fread(SSH_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Ephase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
# ifndef TIDES_MAS
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
            SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
#else
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
           if (SSH_Tphase(i,j,itide)<0.0) 
     &        SSH_Tphase(i,j,itide)=SSH_Tphase(i,j,itide)+360.
          enddo
        enddo
# endif

# endif
# ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
        MPI_master_only write(stdout,*)'write Cangle'
        status=nf_inq_varid(ncidfrc,'tide_Cangle',varid)
        status=nf_fread(UV_Tangle(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cangle',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        MPI_master_only write(stdout,*)'write Cphase'
        status=nf_inq_varid(ncidfrc,'tide_Cphase',varid)
        status=nf_fread(UV_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cphase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        MPI_master_only write(stdout,*)'write Cmax'
        status=nf_inq_varid(ncidfrc,'tide_Cmax',varid)
        status=nf_fread(UV_Tmajor(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cmax',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif

        MPI_master_only write(stdout,*)'write Cmin'
        status=nf_inq_varid(ncidfrc,'tide_Cmin',varid)
        status=nf_fread(UV_Tminor(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,ind,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Cmin',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1          
            UV_Tangle(i,j,itide)=UV_Tangle(i,j,itide)*deg2rad
            UV_Tphase(i,j,itide)=UV_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
# endif
# ifdef POT_TIDES
!
!  Tidal potential amplitude and phase.
!
        status=nf_inq_varid(ncidfrc,'tide_Pamp',varid)
        status=nf_fread(POT_Tamp(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,itide,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Pamp',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        status=nf_inq_varid(ncidfrc,'tide_Pphase',varid)
        status=nf_fread(POT_Tphase(START_2D_ARRAY,itide),
     &                              ncidfrc,varid,itide,r2dvar)
        if (status.ne.NF_NOERR) then
          write(6,3) 'tide_Pphase',itide,frcname(1:lstr)
          stop                           !-->  EXIT
        endif
        do j=0,LOCALMM+1
          do i=0,LOCALLM+1
            POT_Tphase(i,j,itide)=POT_Tphase(i,j,itide)*deg2rad
          enddo
        enddo
# endif

      enddo !Ntides
#ifdef TIDES_MAS
      Tperiod=Tperiod*3600.0
#endif
      return

!
! Sort out error messages: The following portion of the code is 
!===== === ===== ========= not accessed unless something goes wrong.
!
  3   format(/,' GET_TIDES - unable to find forcing variable: ',A,
     &                          /,14x,'in forcing netCDF file: ',A)
      goto 99
  4   write(stdout,5) frcname(1:lstr)
  5   format(/,' GET_TIDES - unable to open forcing netCDF file:',
     &                                                        1x,A)
  6   format(/,' GET_TIDES - error while reading variable: ',A,2x,
     &                                      ' at TIME index = ',i4)

 99   may_day_flag=2

#endif /* SSH_TIDES || UV_TIDES */
      return  
      end
