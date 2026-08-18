!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      subroutine clm_tides (tile)
!
!================================================== Robert Hetland ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to !
!  sea surface height and 2D momentum climatologies, respectively.   !
!                                                                    !
!=====================================================================
!
      implicit none
      integer tile, trd, omp_get_thread_num
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
!
      trd=omp_get_thread_num()
      call clm_tides_tile (Istr,Iend,Jstr,Jend,
     &                     A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd),
     &                     A2d(1,4,trd),A2d(1,5,trd),A2d(1,6,trd),
     &                                               A2d(1,7,trd))
      return
      end
!
!*********************************************************************
      subroutine clm_tides_tile (Istr,Iend,Jstr,Jend,Cangle,Sangle,
     &                           Cphase,Sphase,Etide,Utide,Vtide)
!*********************************************************************
!
# ifdef TIDES_MAS
      USE MOD_TIDES_MAS
# endif
      USE croco_namelist, only:dt
      implicit none
# ifndef TIDES_MAS
#  include "param.h"
# endif
# include "climat.h"
# include "grid.h"
# include "scalars.h"
# include "tides.h"
# include "boundary.h"
!
# ifdef TIDES_MAS
      CHARACTER (len=19)    :: cdate,tool_sectodat
      logical :: flag=.false.
      INTEGER               :: ian,imois,ijour,iheure,iminu,isec,jjulien
      integer :: nn,k,kk,ht
      INTEGER, dimension(0:10) :: esp_corres=(/1,2,3,4,5,0,6,0,7,0,8/)
      complex :: u_interp, v_interp
# endif
      integer
     &     Iend, Istr, Jend, Jstr, i, j, itide
      real
     &     angle, phase, omega_loc, ramp, cff
      real
     &        Cangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Cphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sangle(PRIVATE_2D_SCRATCH_ARRAY),
     &        Sphase(PRIVATE_2D_SCRATCH_ARRAY),
     &        Etide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Utide(PRIVATE_2D_SCRATCH_ARRAY),
     &        Vtide(PRIVATE_2D_SCRATCH_ARRAY)
!
# include "compute_auxiliary_bounds.h"
# ifdef MPI
#  define IRANGE Istr-1,Iend+1
#  define JRANGE Jstr-1,Jend+1
# else
#  define IRANGE IstrR,IendR
#  define JRANGE JstrR,JendR
# endif
!
      angle=0.
      phase=0.
      ramp=1.
# ifdef TIDERAMP
      ramp=TANH(dt*sec2day*float(iic-1))
# endif
      cff=2.*pi*time
!
# ifdef TIDES_MAS
!-------------
! MAS INIT
!-------------
      IF (FIRST_RST_TIME_STEP) THEN
!$acc kernels if(compute_on_device) default(present)

! lecture des conditions aux limites
!-----------------------------------
! a faire une bonne fois pour toute
!----------------------------------
!

! Trouve les ondes presente et stocke les infos
        num(:) = 0
        nomd(:,:) = 0
        fr1(:,:) = 0.0
        fr2(:,:) = 0.0

        DO ht = 1,tide_nbharmcp
          IF (l_presence(nu(ht))) THEN
            nol=ide(nu(ht))
            ir0=rr0(nu(ht))
            ir1=rr1(nu(ht))
            ir2=rr2(nu(ht))
            k=esp_corres(int(nol/100000))
            num(k)=num(k)+1
            nn=num(k)

            k_of(ht) = k
            nn_of(ht) = nn

            nomd(nn,k)=nol
            fr1(nn,k)=ir1/float(ir0)
            fr2(nn,k)=ir2/float(ir0)
          ENDIF
        ENDDO

#  if defined SSH_TIDES
#   if defined OBC_WEST
        rouest(:,:,:) = 0.0
        gouest(:,:,:) = 0.0
#   endif
#   if defined OBC_EAST
        rest(:,:,:) = 0.0
        gest(:,:,:) = 0.0
#   endif
#   if defined OBC_NORTH
        rnord(:,:,:) = 0.0
        gnord(:,:,:) = 0.0
#   endif
#   if defined OBC_SOUTH
        rsud(:,:,:) = 0.0
        gsud(:,:,:) = 0.0
#   endif
        do ht=1,tide_nbharmcp
          if (l_presence(nu(ht))) then
            ir0=rr0(nu(ht))
            k  = k_of(ht)
            nn = nn_of(ht)
#   if defined OBC_WEST
            IF (WESTERN_EDGE) THEN
              do j=JstrR,JendR
                rouest(nn,k,j)=SSH_Tamp(IstrR,j,ht)*100.0
                gouest(nn,k,j)=SSH_Tphase(IstrR,j,ht)
                if(ir0 < 0) gouest(nn,k,j)=gouest(nn,k,j)+180.0
                hnouest(129,j)=-999.
              enddo
            endif
#   endif
#   if defined OBC_EAST
            IF (EASTERN_EDGE) THEN
              do j=JstrR,JendR
                rest(nn,k,j)=SSH_Tamp(IendR,j,ht)*100.0
                gest(nn,k,j)=SSH_Tphase(IendR,j,ht)
                if(ir0 < 0) gest(nn,k,j)=gest(nn,k,j)+180.0
                hnest(129,j)=-999.
              enddo
            endif
#   endif
#   if defined OBC_NORTH
            IF (NORTHERN_EDGE) THEN
              do i=IstrR,IendR
                rnord(nn,k,i)=SSH_Tamp(i,JendR,ht)*100.0
                gnord(nn,k,i)=SSH_Tphase(i,JendR,ht)
                if(ir0 < 0) gnord(nn,k,i)=gnord(nn,k,i)+180.0
                hnnord(129,i)=-999.
              enddo
            endif
#   endif
#   if defined OBC_SOUTH
            IF (SOUTHERN_EDGE) THEN
              do i=IstrR,IendR
                rsud(nn,k,i)=SSH_Tamp(i,JstrR,ht)*100.0
                gsud(nn,k,i)=SSH_Tphase(i,JstrR,ht)
                if(ir0 < 0) gsud(nn,k,i)=gsud(nn,k,i)+180.0
                hnsud(129,i)=-999.
              enddo
            endif
#   endif
          endif
        enddo
#  endif

#  if defined UV_TIDES
#   if defined OBC_WEST
        urouest(:,:,:) = 0.0
        ugouest(:,:,:) = 0.0
        vrouest(:,:,:) = 0.0
        vgouest(:,:,:) = 0.0
#   endif
#   if defined OBC_EAST
        urest(:,:,:) = 0.0
        ugest(:,:,:) = 0.0
        vrest(:,:,:) = 0.0
        vgest(:,:,:) = 0.0
#   endif
#   if defined OBC_NORTH
        urnord(:,:,:) = 0.0
        ugnord(:,:,:) = 0.0
        vrnord(:,:,:) = 0.0
        vgnord(:,:,:) = 0.0
#   endif
#   if defined OBC_SOUTH
        ursud(:,:,:) = 0.0
        ugsud(:,:,:) = 0.0
        vrsud(:,:,:) = 0.0
        vgsud(:,:,:) = 0.0
#   endif

        do ht=1,tide_nbharmcp
          if (l_presence(nu(ht))) then
            ir0=rr0(nu(ht))
            k  = k_of(ht)
            nn = nn_of(ht)
            do j=JstrR,JendR
              do i=IstrR,IendR
                angle=UV_Tangle(i,j,ht)
#   ifdef CURVGRID
     &               -angler(i,j)*rad2deg
#   endif
                call ep2ap(UV_Tmajor(i,j,ht),UV_Tminor(i,j,ht),
     &                     angle,UV_Tphase(i,j,ht),
     &                     U_cmplx(i,j,ht), V_cmplx(i,j,ht))
              enddo
            enddo
#   if defined OBC_WEST
            IF (WESTERN_EDGE) THEN
              do j=JstrR,JendR
                U_interp = 0.5 * (U_cmplx(IstrR, j, ht) +
     &                            U_cmplx(IstrR+1, j, ht))
                urouest(nn,k,j) = ABS(U_interp) * 100.0
                ugouest(nn,k,j) = ATAN2(AIMAG(U_interp),
     &                                  REAL(U_interp)) * rad2deg
                if (ugouest(nn,k,j) < 0.0d0)
     &            ugouest(nn,k,j) = ugouest(nn,k,j) + 360.0d0
                if(ir0 < 0) ugouest(nn,k,j)=ugouest(nn,k,j)+180.0
                unouest(129,j)=-999.

                if (j>JstrR) then
                  v_interp = 0.5 * (V_cmplx(IstrR, j, ht) +
     &                              V_cmplx(IstrR, j-1, ht))

                  vrouest(nn,k,j) = ABS(v_interp) * 100.0
                  vgouest(nn,k,j) = ATAN2(AIMAG(v_interp),
     &                                    REAL(v_interp)) * rad2deg
                  if (vgouest(nn,k,j) < 0.0d0)
     &              vgouest(nn,k,j) = vgouest(nn,k,j) + 360.0d0
                  if(ir0 < 0) vgouest(nn,k,j)=vgouest(nn,k,j)+180.0
                  vnouest(129,j)=-999.
                endif
              enddo
            endif
#   endif
#   if defined OBC_EAST
            IF (EASTERN_EDGE) THEN
              do j=JstrR,JendR
                U_interp = 0.5 * (U_cmplx(IendR, j, ht) +
     &                            U_cmplx(IendR-1, j, ht)) 
                urest(nn,k,j) = ABS(U_interp) * 100.0
                ugest(nn,k,j) = ATAN2(AIMAG(U_interp),
     &                                REAL(U_interp)) * rad2deg
                if (ugest(nn,k,j) < 0.0d0)
     &            ugest(nn,k,j) = ugest(nn,k,j) + 360.0d0
                if(ir0 < 0) ugest(nn,k,j)=ugest(nn,k,j)+180.0
                unest(129,j)=-999.

                if (j>JstrR) then
                  v_interp = 0.5 * (V_cmplx(IendR, j, ht) +
     &                              V_cmplx(IendR, j-1, ht))

                  vrest(nn,k,j) = ABS(v_interp) * 100.0
                  vgest(nn,k,j) = ATAN2(AIMAG(v_interp),
     &                                  REAL(v_interp)) * rad2deg
                  if (vgest(nn,k,j) < 0.0d0)
     &              vgest(nn,k,j) = vgest(nn,k,j) + 360.0d0
                  if(ir0 < 0) vgest(nn,k,j)=vgest(nn,k,j)+180.0
                  vnest(129,j)=-999.
                endif
              enddo
            endif
#   endif
#   if defined OBC_NORTH
            IF (NORTHERN_EDGE) THEN
              do i=IstrR,IendR
                if (i > IstrR) then
                  U_interp = 0.5 * (U_cmplx(i, JendR, ht) +
     &                              U_cmplx(i-1, JendR, ht))
                  urnord(nn,k,i) = ABS(U_interp) * 100.0
                  ugnord(nn,k,i) = ATAN2(AIMAG(U_interp),
     &                                    REAL(U_interp)) * rad2deg
                  if (ugnord(nn,k,i) < 0.0d0)
     &              ugnord(nn,k,i) = ugnord(nn,k,i) + 360.0d0
                  if(ir0 < 0) ugnord(nn,k,i)=ugnord(nn,k,i)+180.0
                  unnord(129,i)=-999.
                endif

                V_interp = 0.5 * (V_cmplx(i, JendR, ht) +
     &                            V_cmplx(i, JendR-1, ht))

                vrnord(nn,k,i) = ABS(v_interp) * 100.0
                vgnord(nn,k,i) = ATAN2(AIMAG(v_interp),
     &                                REAL(v_interp)) * rad2deg
                if (vgnord(nn,k,i) < 0.0d0)
     &            vgnord(nn,k,i) = vgnord(nn,k,i) + 360.0d0
                if(ir0 < 0) vgnord(nn,k,i)=vgnord(nn,k,i)+180.0
                vnnord(129,i)=-999.
              enddo
            endif
#   endif
#   if defined OBC_SOUTH
            IF (SOUTHERN_EDGE) THEN
              do i=IstrR,IendR
                if (i>IstrR) then
                  U_interp = 0.5 * (U_cmplx(i, JstrR, ht) +
     &                              U_cmplx(i-1, JstrR, ht))
                  ursud(nn,k,i) = ABS(U_interp) * 100.0
                  ugsud(nn,k,i) = ATAN2(AIMAG(U_interp),
     &                                  REAL(U_interp)) * rad2deg
                  if (ugsud(nn,k,i) < 0.0d0)
     &              ugsud(nn,k,i) = ugsud(nn,k,i) + 360.0d0
                  if(ir0 < 0) ugsud(nn,k,i)=ugsud(nn,k,i)+180.0
                  unsud(129,i)=-999.
                endif

                v_interp = 0.5 * (V_cmplx(i, JstrR+1, ht) +
     &                            V_cmplx(i, JstrR, ht))

                vrsud(nn,k,i) = ABS(v_interp) * 100.0
                vgsud(nn,k,i) = ATAN2(AIMAG(v_interp),
     &                                    REAL(v_interp)) * rad2deg
                if (vgsud(nn,k,i) < 0.0d0)
     &            vgsud(nn,k,i) = vgsud(nn,k,i) + 360.0d0
                if(ir0 < 0) vgsud(nn,k,i)=vgsud(nn,k,i)+180.0
                vnsud(129,i)=-999.
              enddo
            endif
#   endif
          endif
        enddo
#  endif 
      ENDIF ! FIRST TIME STEP
# endif /* MAS INIT */

# if defined TIDES_MAS
      cdate = tool_sectodat(time)
      CALL tool_decompdate(cdate,ijour,imois,ian,iheure,iminu,isec)
      tempis =iheure*3600.0+iminu*60.0+isec
# endif

# if defined SSH_TIDES
#  ifndef TIDES_MAS
!$acc kernels if(compute_on_device) default(present)
!---------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!---------------------------------------------------------------------
!
      do j=JRANGE
        do i=IRANGE
          Etide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JRANGE
            do i=IRANGE
              Etide(i,j)=Etide(i,j)+ramp*
     &                   SSH_Tamp(i,j,itide)*
     &                   COS(omega_loc-SSH_Tphase(i,j,itide))
#   ifdef MASKING
              Etide(i,j)=Etide(i,j)*rmask(i,j)
#   endif
            enddo
          enddo
        endif
      enddo
!$acc end kernels     
#  else
! calcul des 128 hauteurs par jour
!---------------------------------
      t_0=0.0
      k=0
      if(ijour /= jourmem) then
#   if defined OBC_WEST
       if (WESTERN_EDGE) then
        do j=JstrR,JendR
          call ma1(flag,mynode,j,1,
     &         t_0,ijour,imois,ian,hnouest(:,j),rouest(:,:,j),
     &         gouest(:,:,j),nomd,fr1,fr2,num)
        enddo
       endif
#   endif
#   if defined OBC_EAST
       if (EASTERN_EDGE) then
        do j=JstrR,JendR
          call ma1(flag,mynode,j,1,
     &         t_0,ijour,imois,ian,hnest(:,j),rest(:,:,j),
     &         gest(:,:,j),nomd,fr1,fr2,num)
        enddo
       endif
#   endif
#   if defined OBC_NORTH
       if (NORTHERN_EDGE) then
        do i=IstrR,IendR
          call ma1(flag,mynode,i,1,
     &         t_0,ijour,imois,ian,hnnord(:,i),rnord(:,:,i),
     &         gnord(:,:,i),nomd,fr1,fr2,num)
        enddo
       endif
#   endif
#   if defined OBC_SOUTH
       if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
          call ma1(flag,mynode,i,1,
     &         t_0,ijour,imois,ian,hnsud(:,i),rsud(:,:,i),
     &         gsud(:,:,i),nomd,fr1,fr2,num)
        enddo
       endif
#   endif
      endif

! interpolation parmi les 128 hauteurs par jour
!----------------------------------------------
#   if defined OBC_WEST
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          call ma7(flag,tempis,hnouest(:,j),Etide(IstrR,j),k)
        enddo
      endif
#   endif
#   if defined OBC_EAST
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          call ma7(flag,tempis,hnest(:,j),Etide(IendR,j),k)
        enddo
      endif
#   endif
#   if defined OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=IstrR,IendR
          call ma7(flag,tempis,hnnord(:,i),Etide(i,JendR),k)
        enddo
      endif
#   endif
#   if defined OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
          call ma7(flag,tempis,hnsud(:,i),Etide(i,JstrR),k)
        enddo
      endif
#   endif
#   ifdef MASKING
       do j=JRANGE
         do i=IRANGE
           Etide(i,j)=Etide(i,j)*ramp*rmask(i,j)
         enddo
       enddo
#   endif
!$acc end kernels
#  endif /* MAS COMPOSITION */
!$acc kernels if(compute_on_device) default(present)
!----------------------------------------------------------------
!  Add sub-tidal forcing and adjust climatology to include tides.
!----------------------------------------------------------------
#  if defined ZCLIMATOLOGY
      do j=JRANGE
        do i=IRANGE
          ssh(i,j)=ssh(i,j)+Etide(i,j)
#   ifdef WET_DRY
          if (ssh(i,j) .lt. Dcrit(i,j)-h(i,j)) then
            ssh(i,j)=Dcrit(i,j)-h(i,j)
          endif
#   endif
        enddo
      enddo
#  endif /* ZCLIMATOLOGY */
!----------------------------------------------------------------
!  If appropriate, load tidal forcing into boundary arrays.
!----------------------------------------------------------------
#  ifdef Z_FRC_BRY
#   if defined OBC_WEST
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_west(j)=zetabry_west(j)+Etide(Istr-1,j)
#    ifdef WET_DRY
          if (zetabry_west(j) .lt. Dcrit(Istr-1,j)-h(Istr-1,j)) then
            zetabry_west(j)=Dcrit(Istr-1,j)-h(Istr-1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_EAST
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          zetabry_east(j)=zetabry_east(j)+Etide(Iend+1,j)
#    ifdef WET_DRY
          if (zetabry_east(j) .lt. Dcrit(Iend+1,j)-h(Iend+1,j)) then
            zetabry_east(j)=Dcrit(Iend+1,j)-h(Iend+1,j)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_south(i)=zetabry_south(i)+Etide(i,Jstr-1)
#    ifdef WET_DRY
          if (zetabry_south(i) .lt. Dcrit(i,Jstr-1)-h(i,Jstr-1)) then
            zetabry_south(i)=Dcrit(i,Jstr-1)-h(i,Jstr-1)
          endif
#    endif
        enddo
      endif
#   endif
#   if defined OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=IstrR,IendR
          zetabry_north(i)=zetabry_north(i)+Etide(i,Jend+1)
#    ifdef WET_DRY
          if (zetabry_north(i) .lt. Dcrit(i,Jend+1)-h(i,Jend+1)) then
            zetabry_north(i)=Dcrit(i,Jend+1)-h(i,Jend+1)
          endif
#    endif
        enddo
      endif
#   endif
#  endif /* Z_FRC_BRY */
!$acc end kernels
# endif /* SSH_TIDES */
!
# if defined UV_TIDES
!$acc kernels if(compute_on_device) default(present)
!
!---------------------------------------------------------------------
!  Add tidal currents (m/s) to 2D momentum climatologies.
!---------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=Istr,IendR
          Utide(i,j)=0.
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          Vtide(i,j)=0.
        enddo
      enddo
#  ifndef TIDES_MAS
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=Jstr-1,JendR
            do i=Istr-1,IendR
              angle=UV_Tangle(i,j,itide)
#  ifdef CURVGRID
     &                      -angler(i,j)
#  endif
              phase=omega_loc-UV_Tphase(i,j,itide)
              Cangle(i,j)=COS(angle)
              Cphase(i,j)=COS(phase)
              Sangle(i,j)=SIN(angle)
              Sphase(i,j)=SIN(phase)
            enddo
          enddo
          do j=JstrR,JendR
            do i=Istr,IendR
              Utide(i,j)=Utide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i-1,j,itide)+
     &                                UV_Tmajor(i  ,j,itide))*
     &                               (Cangle(i-1,j)+Cangle(i,j))*
     &                               (Cphase(i-1,j)+Cphase(i,j))-
     &                               (UV_Tminor(i-1,j,itide)+
     &                                UV_Tminor(i  ,j,itide))*
     &                               (Sangle(i-1,j)+Sangle(i,j))*
     &                               (Sphase(i-1,j)+Sphase(i,j)))
#  ifdef MASKING
             Utide(i,j)=Utide(i,j)*umask(i,j)
#  endif
            enddo
          enddo
          do j=Jstr,JendR
            do i=IstrR,IendR
              Vtide(i,j)=Vtide(i,j)+ramp*
     &                        0.125*((UV_Tmajor(i,j-1,itide)+
     &                                UV_Tmajor(i,j  ,itide))*
     &                               (Sangle(i,j-1)+Sangle(i,j))*
     &                               (Cphase(i,j-1)+Cphase(i,j))+
     &                               (UV_Tminor(i,j-1,itide)+
     &                                UV_Tminor(i,j  ,itide))*
     &                               (Cangle(i,j-1)+Cangle(i,j))*
     &                               (Sphase(i,j-1)+Sphase(i,j)))
#  ifdef MASKING
              Vtide(i,j)=Vtide(i,j)*vmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo
# else 
! calcul des 128 hauteurs par jour
!---------------------------------
!
      t_0=0.0
      k=0
      if(ijour /= jourmem) then
#  if defined OBC_WEST
        if (WESTERN_EDGE) then
          do j=JstrR,JendR
            call ma1(flag,mynode,j,1,
     &          t_0,ijour,imois,ian,unouest(:,j),urouest(:,:,j),
     &          ugouest(:,:,j),nomd,fr1,fr2,num)
            call ma1(flag,mynode,j,1,
     &          t_0,ijour,imois,ian,vnouest(:,j),vrouest(:,:,j),
     &          vgouest(:,:,j),nomd,fr1,fr2,num)
          end do
        end if
#  endif
#  if defined OBC_EAST
        if (EASTERN_EDGE) then
          do j=JstrR,JendR
            call ma1(flag,mynode,j,1,
     &          t_0,ijour,imois,ian,unest(:,j),urest(:,:,j),
     &          ugest(:,:,j),nomd,fr1,fr2,num)
          enddo
          do j=Jstr,JendR
            call ma1(flag,mynode,j,1,
     &          t_0,ijour,imois,ian,vnest(:,j),vrest(:,:,j),
     &          vgest(:,:,j),nomd,fr1,fr2,num)
          end do
        end if
#  endif
#  if defined OBC_NORTH
        if (NORTHERN_EDGE) then
          do i=IstrR,IendR
            call ma1(flag,mynode,i,1,
     &          t_0,ijour,imois,ian,unnord(:,i),urnord(:,:,i),
     &          ugnord(:,:,i),nomd,fr1,fr2,num)
            call ma1(flag,mynode,i,1,
     &          t_0,ijour,imois,ian,vnnord(:,i),vrnord(:,:,i),
     &          vgnord(:,:,i),nomd,fr1,fr2,num)
          end do
        end if
#  endif
#  if defined OBC_SOUTH
        if (SOUTHERN_EDGE) then
          do i=IstrR,IendR
            call ma1(flag,mynode,i,1,
     &          t_0,ijour,imois,ian,unsud(:,i),ursud(:,:,i),
     &          ugsud(:,:,i),nomd,fr1,fr2,num)
            call ma1(flag,mynode,i,1,
     &          t_0,ijour,imois,ian,vnsud(:,i),vrsud(:,:,i),
     &          vgsud(:,:,i),nomd,fr1,fr2,num)
          end do
        end if
#  endif
      end if

! interpolation parmi les 128 hauteurs par jour
!----------------------------------------------
#  if defined OBC_WEST
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
            call ma7(flag,tempis,unouest(:,j),Utide(Istr,j),k)
            call ma7(flag,tempis,vnouest(:,j),Vtide(IstrR,j),k)
        end do
      end if
#  endif
#  if defined OBC_EAST
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
            call ma7(flag,tempis,unest(:,j),Utide(IendR,j),k)
            call ma7(flag,tempis,vnest(:,j),Vtide(IendR,j),k)
        end do
      end if
#  endif
#  if defined OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=IstrR,IendR
            call ma7(flag,tempis,unnord(:,i),Utide(i,JendR),k)
            call ma7(flag,tempis,vnnord(:,i),Vtide(i,JendR),k)
        end do
      end if
#  endif
#  if defined OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=IstrR,IendR
            call ma7(flag,tempis,unsud(:,j),Utide(i,Jstr),k)
            call ma7(flag,tempis,vnsud(:,j),Vtide(i,JstrR),k)
        end do
      end if
#  endif


#   ifdef MASKING
       do j=JstrR,JendR
         do i=Istr,IendR
           Utide(i,j)= Utide(i,j)*ramp*umask(i,j)
         enddo
       enddo
       do j=Jstr,JendR
         do i=IstrR,IendR
           Vtide(i,j)= Vtide(i,j)*ramp*vmask(i,j)
         end do
       end do
#   endif
#  endif 

# if defined TIDES_MAS
      jourmem = ijour
# endif

#  if defined M2CLIMATOLOGY
!
!  Add sub-tidal forcing and adjust climatology to include tides.
!
      do j=JstrR,JendR
        do i=Istr,IendR
          ubclm(i,j)=ubclm(i,j)+Utide(i,j)
        enddo
      enddo
      do j=Jstr,JendR
        do i=IstrR,IendR
          vbclm(i,j)=vbclm(i,j)+Vtide(i,j)
        enddo
      enddo
#  endif
!
!  If appropriate, load tidal forcing into boundary arrays.
!
#  ifdef M2_FRC_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_west(j)=ubarbry_west(j)+Utide(Istr,j)
        enddo
        do j=Jstr,JendR
          vbarbry_west(j)=vbarbry_west(j)+Vtide(Istr-1,j)
        enddo
      endif
#   endif
#   ifdef OBC_EAST
      if (EASTERN_EDGE) THEN
        do j=JstrR,JendR
          ubarbry_east(j)=ubarbry_east(j)+Utide(Iend+1,j)
        enddo
        do j=Jstr,JendR
          vbarbry_east(j)=vbarbry_east(j)+Vtide(Iend+1,j)
        enddo
      endif
#   endif
#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_south(i)=ubarbry_south(i)+Utide(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          vbarbry_south(i)=vbarbry_south(i)+Vtide(i,Jstr)
        enddo
      endif
#   endif
#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) THEN
        do i=Istr,IendR
          ubarbry_north(i)=ubarbry_north(i)+Utide(i,Jend+1)
        enddo
        do i=IstrR,IendR
          vbarbry_north(i)=vbarbry_north(i)+Vtide(i,Jend+1)
        enddo
      endif
#   endif
#  endif /* M2_FRC_BRY */
!$acc end kernels
# endif /* UV_TIDES */

# if defined POT_TIDES
!$acc kernels if(compute_on_device) default(present)
!
!-----------------------------------------------------------------------
!  Compute tidal potential (m)
!  -- > to be applied as a pressure-gradient force
!       in prsgrd (if SOLVE3D) or step2d (if not SOLVE3D)
!-----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          Ptide(i,j)=0.
        enddo
      enddo
      do itide=1,Ntides
        if (Tperiod(itide).gt.0.) then
          omega_loc=cff/Tperiod(itide)
          do j=JstrR,JendR
            do i=IstrR,IendR
              Ptide(i,j)=Ptide(i,j)+
     &                   ramp*POT_Tamp(i,j,itide)*
     &                   COS(omega_loc-POT_Tphase(i,j,itide))
#  ifdef MASKING
              Ptide(i,j)=Ptide(i,j)*rmask(i,j)
#  endif
            enddo
          enddo
        endif
      enddo

!$acc end kernels
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        Ptide(START_2D_ARRAY))
#  endif

# endif /* POT_TIDES */

#else
      subroutine clm_tides_empty
#endif /* SSH_TIDES || UV_TIDES */
      return
      end
