!---------------------------------------------------------------------------
 MODULE sed_MUSTANG_CROCO
!---------------------------------------------------------------------------

#include "cppdefs.h"
#if defined MUSTANG 

   !&E==========================================================================
   !&E                   ***  MODULE  sed_MUSTANG_CROCO  ***
   !&E
   !&E
   !&E ** Purpose : concerns subroutines related to sediment dynamics link to hydrodynamic model
   !&E              to be used in CROCO system
   !&E 
   !&E ** Description :
   !&E     subroutine sed_MUSTANG_settlveloc     ! settling velocity in the water column
   !&E     subroutine sed_gradvit         ! calcul gradient de vitesse, u*
   !&E     subroutine sed_skinstress     ! computes the skin stress
   !&E     subroutine sedinit_fromfile  ! reads filrepsed where all results of 
   !&E                                   ! sediment dyn. are stored (depends on hydro model)
   !&E     subroutine sed_obc_corflu      ! echange MPI des corrections de flux horizontaux pour sables
   !&E     subroutine sed_exchange_w2s_MARS  ! echange MPI de flux verticaux a l'interface WS
   !&E     subroutine sed_exchange_s2w_MARS  ! echange MPI de flux verticaux a l'interface WS
   !&E     subroutine sed_outsaverestart ! save related field for future runs
   !&E
   !&E==========================================================================

#include "coupler_define_MUSTANG.h"

   !! * Modules used
   USE comMUSTANG
   USE comsubstance
   USE module_MUSTANG
   USE module_substance
   
   IMPLICIT NONE

   !! * Accessibility 
   PUBLIC sedinit_fromfile, sed_skinstress, sed_gradvit, sed_MUSTANG_settlveloc
# ifdef key_MUSTANG_bedload
   PUBLIC sed_bottom_slope
# endif

   PRIVATE
   
 CONTAINS
 
 !!===========================================================================================
 
  SUBROUTINE sed_MUSTANG_settlveloc(ifirst, ilast, jfirst, jlast,   &
                     WATER_CONCENTRATION) 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_settlveloc  ***
   !&E
   !&E ** Purpose : settling velocity computation
   !&E
   !&E ** Description : use arguments and common variable 
   !&E  arguments IN : 
   !&E         WATER_CONCENTRATION=t : WATER_CONCENTRATION 
   !&E  arguments OUT:
   !&E         WAT_SETTL: settling velocities for CROCO
   !&E         ws3_bottom_MUSTANG: WAT_SETTL_MUSTANG in  bottom cell
   !&E
   !&E  need to be know by hydrodynamic code:
   !&E          GRAVITY
   !&E         kmax=NB_LAYER_WAT  : connu via coupleur_dimhydro_MUSTANG.h
   !&E          
   !&E  need to be know by code treated substance (if not ==> coupler_MUSTANG.F90)
   !&E         imud1, nvpc, nvp, nv_adv, isand1,isand2
   !&E         f_ws(iv) (if key_MUSTANG_flocmod)
   !&E         ws_free_opt,ws_free_para,ws_free_min,ws_free_max,ws_hind_opt,ws_hind_para :    
   !&E     
   !&E  use module MUSTANG variables  :
   !&E         ros(iv)
   !&E         ws_sand(iv)
   !&E        
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E
   !&E--------------------------------------------------------------------------

   !! * Arguments
   INTEGER, INTENT(IN)                               :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh), DIMENSION(ARRAY_WATER_CONC), INTENT(IN)  :: WATER_CONCENTRATION  ! CROCO : directly t 
   
   !! * Local declarations
   INTEGER                            :: iv, k, ivpc, i, j
   REAL(KIND=rsh)                     :: cmes, phi, phiv, De
   REAL(KIND=rsh), PARAMETER           :: nuw = 0.00000102_rsh

   !!---------------------------------------------------------------------------
   !! * Executable part

      DO j = jfirst, jlast
      DO i = ifirst, ilast
           
         IF(htot(i,j) > h0fond)  THEN
           DO k=1,NB_LAYER_WAT
              cmes=0.0_rsh
              DO ivpc=imud1,nvpc
                !cmes=cmes+WATER_CONCENTRATION(ivpc,k,i,j)
                cmes=cmes+t(i,j,k,1,itemp+ntrc_salt+ivpc)
              ENDDO
              cmes=MAX(0.0_rsh,cmes)

! first calculating sand settling velocity
              DO iv=isand1,isand2
                WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=ws_sand(iv)
              ENDDO
              
! next mud settling velocity 
#ifdef key_MUSTANG_flocmod
              DO iv=imud1,nvp         
                 WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=f_ws(iv)  
              ENDDO
      
#else
      
              DO iv=imud1,nvp
! free settling velocity - flocculation
    
                IF(ws_free_opt(iv) == 0) THEN ! constant settling velocity
                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=ws_free_min(iv)
                ELSEIF (ws_free_opt(iv) == 1) THEN ! Van Leussen 1994

                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=ws_free_para(1,iv)*cmes**ws_free_para(2,iv) &
                              *(1._rsh+ws_free_para(3,iv)*gradvit(k,i,j))/  &
                               (1._rsh+ws_free_para(4,iv)*gradvit(k,i,j)**2._rsh)

                ELSEIF (ws_free_opt(iv) == 2) THEN ! Winterwerp 1999
                  De=ws_free_para(1,iv)+ws_free_para(2,iv)*cmes/(ws_free_para(3,iv)*sqrt(gradvit(k,i,j)))
                  IF (De.GT.sqrt(nuw/gradvit(k,i,j))) THEN 
                    De=sqrt(nuw/gradvit(k,i,j)) ! in case of large C/low G limit floc size to kolmogorov microscale
                  ENDIF
                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=(ros(iv)-RHOREF)*GRAVITY/(18._rsh*RHOREF*nuw)  &
                                *ws_free_para(1,iv)**(3._rsh-ws_free_para(4,iv))  &
                                *De**(ws_free_para(4,iv)-1._rsh)
                
                ELSEIF (ws_free_opt(iv) == 3) THEN ! Wolanski et al., 1989
                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=ws_free_para(1,iv)*cmes**ws_free_para(2,iv)

                ENDIF


! Hindered settling
! if ws_hind_opt.EQ.0 : no hindered settling...WAT_SETTL unchanged
                IF (ws_hind_opt(iv) == 1) THEN ! Scott, 1984
                  phi=MIN(1.0_rsh,cmes/ws_hind_para(1,iv))
                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=  &
                         WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)*(1._rsh-phi)**ws_hind_para(2,iv)
                ELSEIF (ws_hind_opt(iv) == 2) THEN ! Winterwerp, 2002 - ATTENTION : ros(iv) must be the same for all MUDS variables
                  phi=cmes/ros(iv)
                  IF (ws_free_opt(iv) == 2) THEN
                    phi=cmes/ros(iv)
                    phiv=phi*(De/ws_free_para(1,iv))**(3._rsh-ws_free_para(4,iv))
                  ELSE
                    phiv=cmes/ws_hind_para(1,iv)
                  ENDIF

                  WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)* &
                       (1._rsh-phiv)**ws_hind_para(2,iv)*(1._rsh-phi) /(1._rsh+2.5_rsh*phiv)
                ELSEIF (ws_hind_opt(iv) == 3) THEN ! wolanski et al., 1989
                  IF (ws_free_opt(iv) == 3) THEN
                    WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=WAT_SETTL(i,j,k,iv)/ &
                       (cmes**2._rsh+ws_hind_para(1,iv)**2._rsh)**ws_hind_para(2,iv)
                  ELSE
                    
                  ENDIF
                ENDIF

! limiting ws with min/max values...
! necessary if ndt_part not updated during the simulation, ndt_part calculated in subreaddat from a maximum settling velocity
                WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=max(ws_free_min(iv), &
                            min(ws_free_max(iv),WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)))

              ENDDO
              

            ! rverney : attention update parsub_newdt with flocculation! 

#endif
              DO iv=nvpc+1,nvp
                IF(irkm_var_assoc(iv) < imud1 .AND. irkm_var_assoc(iv) > 0) THEN    ! sorbed substances on sands
                   WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=WAT_SETTL(i,j,k,itemp+ntrc_salt+irkm_var_assoc(iv))
                ENDIF
              ENDDO
              DO iv=nvp+1,nv_adv
                WAT_SETTL(i,j,k,itemp+ntrc_salt+iv)=0.0_rsh
              ENDDO
              
           ENDDO
           ws3_bottom_MUSTANG(1:nvp,i,j)=WAT_SETTL(i,j,1,itemp+ntrc_salt+1:nvp+itemp+ntrc_salt)
         ELSE
           WAT_SETTL(i,j,:,:)=0.0_rsh
           ws3_bottom_MUSTANG(:,i,j)=0.0_rsh
         ENDIF

       ENDDO
       ENDDO

  END SUBROUTINE sed_MUSTANG_settlveloc     

!!===========================================================================================

  SUBROUTINE sed_gradvit(ifirst, ilast, jfirst, jlast)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_gradvit  ***
   !&E
   !&E ** Purpose : calculation of the turbulence energy G  depending on hydro code ?
   !&E
   !&E ** Description : G= sqrt(turbulence dissipation/viscosity)
   !&E                 to be programmed using hydrodynamic knowledge
   !&E           using htot, RHOREF, sig, epn, nz ..
   !&E
   !&E     output : gradvit (in comMUSTANG)
   !&E
   !&E ** Called by :  step
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#  include "mixing.h"
#  include "ocean3d.h"

   !! * Arguments 
   INTEGER, INTENT(IN)       :: ifirst, ilast, jfirst, jlast

   !! * Local declarations
   INTEGER        :: i, j, k
   REAL(KIND=rsh) :: dist_surf_on_bottom, nuw

      nuw = 1.0e-6
      DO j = jfirst, jlast
      DO i = ifirst, ilast
          IF(htot(i, j) .GT. h0fond)  THEN
           DO k=1, N
            dist_surf_on_bottom = ((z_w(i, j, N) - z_r(i, j, k)) / (z_r(i, j, k) - z_w(i, j, 0)))
            gradvit(k, i, j) = sqrt(ustarbot(i, j)**3._rsh / 0.4_rsh / htot(i, j)/(nuw + epsilon_MUSTANG) * dist_surf_on_bottom) 
           END DO
           ! gradvit : G=sqrt( turbulence dissipation rate/ vertical viscosity coefficient)
           !  if  turbulence dissipation rate has not been already evaluated: 
           ! use empirical formula from   Nezu and Nakawaga (1993)
           !  turbulence dissipation_rate= ustarbot**3 /Karman/Htot * (distance from surface/distance from bottom)
          ENDIF
      ENDDO
      ENDDO

  END SUBROUTINE sed_gradvit

   !!==============================================================================
   SUBROUTINE sed_skinstress(ifirst, ilast, jfirst, jlast)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_skinstress  ***
   !&E
   !&E ** Purpose : computes  bottom shear stress
   !&E
   !&E ** Description :  
   !&E 
   !&E Compute bottom skin-friction stress due to combined maximum wave and current 
   !&E interaction
   !&E 
   !&E Available options to compute tauskin (combine current + wave (if WAVE_OFFLINE 
   !&E is defined)) and tauskin_x/tauskin_y (components / rho) :
   !&E - default : Soulsby formulation (with z0sed)
   !&E - BBL : d50 is constant (160microns) in the bustrw/bvstrw computation see 
   !&E bbl.F (and optionnal key_tauskin* do not apply)
   !&E 
   !&E Available options to compute tauskin_c
   !&E - default : compute at u,v points and then compute at rho point 
   !&E     with ubar if key_tauskin_c_ubar is defined, with bottom u(1) if not, 
   !&E     this case use 12 u,v points to compute tauc at rho point
   !&E     if key_tauskin_c_upwind is defined, x/y component used are tauskin 
   !&E     computed at uv point upwind from current
   !&E - key_tauskin_c_rho : compute at rho point from immediate u,v value 
   !&E     with ubar if key_tauskin_c_ubar is definde, with bottom u(1) if not
   !&E
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE module_substance, ONLY : bustr, bvstr
# ifdef WAVE_OFFLINE
   USE module_substance, ONLY : Uwave, Dwave, Pwave
# endif

#  ifdef BBL
#  include "bbl.h"
#  endif

#  include "grid.h"
#  include "ocean3d.h"
#  include "ocean2d.h"


   INTEGER, INTENT(IN)  :: ifirst, ilast, jfirst, jlast
   INTEGER i,j,k
   REAL    speed ! current speed
   REAL(KIND=rsh)    :: urho, vrho

   REAL(KIND=rsh),DIMENSION(ifirst-1: ilast+1, jfirst-1: jlast+1)  :: Zr
   REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY) :: tauskin_c_u ! bottom stress due to current / rho, compute at u point
   REAL(KIND=rsh),DIMENSION(GLOBAL_2D_ARRAY) :: tauskin_c_v ! bottom stress due to current / rho, compute at v point

# ifdef WAVE_OFFLINE
   REAL(KIND=rsh)    :: fws2ij, speedbar, alpha, beta, cosamb, sinamb, tauskin_cw
# endif
# ifdef key_tauskin_c_upwind   
   REAL(KIND=rsh)    :: cff1, cff2, cff3, cff4
# endif  

!-----------------------------------------------------------------------
!  Compute bottom skin-friction stress due to combined maximum wave
!  and current interaction
!-----------------------------------------------------------------------
!
#  ifdef BBL /*warning, d50 is constant (160microns) in the bustrw/bvstrw computation see bbl.F */
      do j = jfirst, jlast
        do i = ifirst, ilast
          tauskin(i, j) = sqrt( bustrw(i, j)**2 + bvstrw(i, j)**2) * RHOREF
#  ifdef WET_DRY AND MASKING
          tauskin(i, j) = tauskin(i, j) * rmask_wet(i, j)
#  endif
#  ifdef key_MUSTANG_bedload
          urho = 0.5 * (u(i, j, 1, nnew) + u(i+1, j, 1, nnew))
          vrho = 0.5 * (v(i, j, 1, nnew) + v(i, j+1, 1, nnew))
          speed = SQRT(urho**2 + vrho**2)
          tauskin_c(i, j) = tauskin(i, j) ! used in eval_bedload
          tauskin_x(i, j) = urho / (speed + epsilon_MUSTANG) * tauskin_c(i, j) / (rho(i, j, 1) + rho0)
          tauskin_y(i, j) = vrho / (speed + epsilon_MUSTANG) * tauskin_c(i, j) / (rho(i, j, 1) + rho0)
#  endif


#  else /* else on #ifdef BBL */

      do j = jfirst-1, jlast+1
        do i = ifirst-1, ilast+1
          Zr(i, j) =  max(z_r(i, j, 1) - z_w(i, j, 0), z0sed(i, j) + 1.E-4)
        enddo
      enddo

# ifdef key_tauskin_c_rho
  do j = jfirst, jlast
    do i = ifirst-1, ilast
      raphbx(i, j) = ABS(u(i+1, j, 1, nnew)) / (ABS(u(i+1, j, 1, nnew)) + epsilon_MUSTANG)
    enddo
  enddo

  do j = jfirst-1, jlast
    do i = ifirst, ilast
      raphby(i, j) = ABS(v(i, j+1, 1, nnew)) / (ABS(v(i, j+1, 1, nnew)) + epsilon_MUSTANG)
    enddo
  enddo

  do j = jfirst, jlast
    do i = ifirst, ilast
#ifdef key_tauskin_c_ubar
#ifdef MPI
    if ((float(j + jj * Mm) .GE. 0) .and. (float(i + ii * Lm) .GE. 0) )then
#else
    if ((float(j       ) .GE. 0) .and. (float(i       ) .GE. 0)) then
#endif
      urho = 0.5 * (ubar(i, j, nnew) + ubar(i+1, j, nnew))
      vrho = 0.5 * (vbar(i, j, nnew) + vbar(i, j+1, nnew))  
    else
      urho = 0.
      vrho = 0.
    endif  
    speed = SQRT(urho**2 + vrho**2)
    tauskin_c(i, j) = 0.16_rsh * (LOG( ( z_w(i, j, N) - z_w(i, j, 0)) /   &
                  (z0sed(i, j) * 2.718)))**(-2) * speed**2                &
                  * (rho(i, j, 1) + rho0)
# else
    urho = 0.5 * (u(i, j, 1, nnew) + u(i+1, j, 1, nnew))
    vrho = 0.5 * (v(i, j, 1, nnew) + v(i, j+1, 1, nnew))
    speed = SQRT(urho**2 + vrho**2)
    tauskin_c(i, j) = 0.16_rsh * (LOG( ( Zr(i, j) ) /   &
                  (z0sed(i, j))))**(-2) * speed**2                &
                  * (rho(i, j, 1) + rho0)
# endif
      
      tauskin_x(i, j) = urho / (speed + epsilon_MUSTANG) * tauskin_c(i, j) / (rho(i, j, 1) + rho0)
      tauskin_y(i, j) = vrho / (speed + epsilon_MUSTANG) * tauskin_c(i, j) / (rho(i, j, 1) + rho0)
# else /* else on #ifdef key_tauskin_c_rho */

      do j = jfirst, jlast
        do i = ifirst-1, ilast
          raphbx(i,j)=ABS(u(i+1,j,1,nnew))/(ABS(u(i+1,j,1,nnew))+epsilon_MUSTANG)
#ifdef key_tauskin_c_ubar
          speed=SQRT(0.0625_rsh*(vbar(i,j+1,nnew)+vbar(i+1,j+1,nnew)+  &
                  vbar(i,j,nnew)+vbar(i+1,j,nnew))**2         &
                  +ubar(i+1,j,nnew)**2) *ubar(i+1,j,nnew)
#ifdef MPI
          if (float(i +ii*Lm) .GE. 0) then
#else
          if (float(i       ) .GE. 0) then
#endif
            tauskin_c_u(i,j)=0.16_rsh*(LOG( ( z_w(i,j,N)-z_w(i,j,0))/   &
                  (z0sed(i,j)*2.718)))**(-2)*speed
          else
            tauskin_c_u(i,j)=0.
          endif

#else
          speed=SQRT(0.0625_rsh*(v(i,j+1,1,nnew)+v(i+1,j+1,1,nnew)+  &
                  v(i,j,1,nnew)+v(i+1,j,1,nnew))**2         &
                  +u(i+1,j,1,nnew)**2) *u(i+1,j,1,nnew)   
          tauskin_c_u(i,j)=0.16_rsh*(LOG(0.5*(Zr(i+1,j)+Zr(i,j))/   &
                  z0sed(i,j)))**(-2)*speed

#endif
       enddo
      enddo  
      DO j=jfirst-1,jlast
        DO i=ifirst,ilast
          raphby(i,j)=ABS(v(i,j+1,1,nnew))/(ABS(v(i,j+1,1,nnew))+epsilon_MUSTANG)
#ifdef key_tauskin_c_ubar
          speed=SQRT(0.0625_rsh*(ubar(i+1,j,nnew)+ubar(i+1,j+1,nnew)+  &
                  ubar(i,j+1,nnew)+ubar(i,j,nnew))**2         &
                  +vbar(i,j+1,nnew)**2) *vbar(i,j+1,nnew)
#ifdef MPI
          if (float(j +jj*Mm) .GE. 0) then
#else
          if (float(j       ) .GE. 0) then
#endif
            tauskin_c_v(i,j)=0.16_rsh*(LOG( (z_w(i,j,N)-z_w(i,j,0))  /   &
                  (z0sed(i,j)*2.718)))**(-2)*speed
          else
            tauskin_c_v(i,j)=0.
          endif
#else
          speed=SQRT(0.0625_rsh*(u(i+1,j,1,nnew)+u(i+1,j+1,1,nnew)+  &
                  u(i,j+1,1,nnew)+u(i,j,1,nnew))**2         &
                  +v(i,j+1,1,nnew)**2) *v(i,j+1,1,nnew)
          tauskin_c_v(i,j)=0.16_rsh*(LOG(0.5*(Zr(i,j+1)+Zr(i,j))/   &
                  z0sed(i,j)))**(-2)*speed
#endif
       enddo
      enddo
  
      DO j=jfirst,jlast
        DO i=ifirst,ilast
# ifdef key_tauskin_c_upwind
            cff1=0.5*(1.0+SIGN(1.0,tauskin_c_u(i,j)))
            cff2=0.5*(1.0-SIGN(1.0,tauskin_c_u(i,j)))
            cff3=0.5*(1.0+SIGN(1.0,tauskin_c_u(i-1  ,j)))
            cff4=0.5*(1.0-SIGN(1.0,tauskin_c_u(i-1 ,j)))
            tauskin_x(i,j)=cff3*(cff1*tauskin_c_u(i-1,j)+                     &
                       cff2*0.5*(tauskin_c_u(i-1,j)+tauskin_c_u(i,j)))+    &
                       cff4*(cff2*tauskin_c_u(i,j)+                   &
                       cff1*0.5*(tauskin_c_u(i-1,j)+tauskin_c_u(i,j)))

            cff1=0.5*(1.0+SIGN(1.0,tauskin_c_v(i,j)))
            cff2=0.5*(1.0-SIGN(1.0,tauskin_c_v(i,j)))
            cff3=0.5*(1.0+SIGN(1.0,tauskin_c_v(i,j-1)))
            cff4=0.5*(1.0-SIGN(1.0,tauskin_c_v(i,j-1)))
            tauskin_y(i,j)=cff3*(cff1*tauskin_c_v(i,j-1)+                      &
                      cff2*0.5*(tauskin_c_v(i,j-1)+tauskin_c_v(i,j)))+      &
                       cff4*(cff2*tauskin_c_v(i,j)+                     &
                       cff1*0.5*(tauskin_c_v(i,j-1)+tauskin_c_v(i,j)))
# else
            tauskin_x(i,j) = (tauskin_c_u(i,j)*raphbx(i,j)+tauskin_c_u(i-1,j)   &
                   *raphbx(i-1,j))/(raphbx(i,j)                 &
                    +raphbx(i-1,j)+epsilon_MUSTANG)

            tauskin_y(i,j) = (tauskin_c_v(i,j)*raphby(i,j)+tauskin_c_v(i,j-1)   &
                    *raphby(i,j-1))/(raphby(i,j)                 &
                    +raphby(i,j-1)+epsilon_MUSTANG)
# endif
   tauskin_c(i,j)=SQRT(tauskin_x(i,j)**2+tauskin_y(i,j)**2)*(rho(i,j,1)+rho0)

# endif /* end of #ifdef key_tauskin_c_rho */

# ifdef WAVE_OFFLINE
            fws2ij=fws2
            IF (l_fricwave .AND. Uwave(i,j)*Pwave(i,j) > 0.0_rsh .AND.  &
                     Pwave(i,j) > 0.001_rsh.AND.Uwave(i,j)>0.001_rsh) THEN
              fws2ij=0.5_rsh*1.39_rsh*(Uwave(i,j)*Pwave(i,j)/REAL(2.0_rlg*pi*z0sed(i,j),rsh))**(-0.52_rsh)
            ENDIF

            ! calculation of shear stress due to waves tauskin_w 
            ! --------------------------
            tauskin_w(i,j) = ((rho(i,j,1)+rho0) * fws2ij * Uwave(i,j)**2)

            speedbar = SQRT(ubar(i,j,nnew)**2+vbar(i,j,nnew)**2)
            IF(tauskin_c(i,j) > 0.0_rsh .AND. tauskin_w(i,j) > 0.0_rsh .AND. speedbar> 0.0_rsh) THEN

            ! calculation of shear stress with the formula of Soulsby (1995)
            ! ======================================================
            ! calculation of  tauskin_c (current) influenced by the waves
            ! ----------------------------------------------------
              tauskin_cw = tauskin_c(i,j)*(1+(1.2*(tauskin_w(i,j)/(tauskin_w(i,j)+tauskin_c(i,j)))**3.2))

            ! calculating the difference in angle between the direction of the waves and the current
            ! ---------------------------------------------------------------------------
            ! calculating the direction of the current relative to the north
              alpha = ACOS(vbar(i,j,nnew)/speedbar)   ! en radians
            ! calculation of wave orientation relative to north
              beta = Dwave(i,j)   ! beta and Dwave in radians
            ! calculation of cos(alpha-beta) and sin(alpha-beta)
              cosamb = ABS(COS(alpha-beta))
              sinamb = ABS(SIN(alpha-beta))

            ! calculation of tauskin (waves + current)
            ! -----------------------------------
              tauskin(i,j) = SQRT( (tauskin_cw + tauskin_w(i,j) * cosamb)**2 + (tauskin_w(i,j) * sinamb)**2 )

          ELSE
            tauskin(i,j) = tauskin_w(i,j) + tauskin_c(i,j)
          ENDIF


# else
    tauskin(i,j) = tauskin_c(i,j)
# endif

# if defined WET_DRY && defined MASKING 
    tauskin(i,j) = tauskin(i,j) * rmask_wet(i,j)
#  endif


!!!!!!!!!!!!!!!! FORCING :     u*=ustarbot    
    if (htot(i, j) .GT. h0fond) then
      if (tauskin(i, j) < 0.) then
        ustarbot(i, j) = 0.0_rsh
      else
        ustarbot(i, j) = (tauskin(i, j) / RHOREF)**0.5_rsh
      endif
    endif

#  endif  /* end of #ifdef BBL */
        enddo
      enddo



  END SUBROUTINE sed_skinstress
   !!==============================================================================
#ifdef key_MUSTANG_bedload
  SUBROUTINE sed_bottom_slope(ifirst, ilast, jfirst, jlast, bathy)
   !&E--------------------------------------------------------------------------                         
   !&E                 ***  ROUTINE sed_bottom_slope  ***
   !&E
   !&E ** Purpose : evaluation of bottom slope
   !&E
   !&E ** Description : depend on host model grid, compute slope_dhdx and 
   !&E    slope_dhdy from bathy of neigbour cells if they are not masked
   !&E
   !&E ** Called by :  
   !&E
   !&E ** External calls :
   !&E
   !&E--------------------------------------------------------------------------
   !! * Arguments
   INTEGER, INTENT(IN)                                  :: ifirst, ilast, jfirst, jlast
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)  :: bathy  ! bathymetry (m)

   !! * Local declarations
   INTEGER          :: i, j

      DO j = jfirst, jlast
        DO i = ifirst, ilast
          IF (bathy(i+1, j).LE. -valmanq .OR. bathy(i-1, j).LE. -valmanq) then
             slope_dhdx(i, j) = 0.0_rsh
          ELSE
             slope_dhdx(i, j) = -1.0_rsh*(-bathy(i+1, j)+bathy(i-1, j)) / (2.0_rsh * CELL_DX(i, j))
          ENDIF
          IF (bathy(i, j+1).LE. -valmanq .OR. bathy(i, j-1).LE. -valmanq) then
             slope_dhdy(i, j) = 0.0_rsh
          ELSE
             slope_dhdy(i, j) = -1.0_rsh*(-bathy(i, j+1)+bathy(i, j-1)) / (2.0_rsh * CELL_DY(i, j))
          ENDIF
        ENDDO
      ENDDO

  END SUBROUTINE sed_bottom_slope
#endif

    !!==============================================================================
  SUBROUTINE sedinit_fromfile(BATHY_H0)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sedinit_fromfile  ***
   !&E
   !&E ** Purpose : manages the fields to be re-initialized in case of the 
   !&E              continuation of a previous run
   !&E
   !&E ** Description : open and read a netcdf file, written during a previous run 
   !&E                          (save file created by sed_outsaverestart) 
   !&E
   !&E ** Called by :  MUSTANG_init_sediment 
   !&E
   !&E ** External calls : ionc4_openr, ionc4_read_time,
   !&E                     ionc4_read_subxyt, ionc4_read_subzxyt
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
      implicit none

      !! * Arguments
      REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)       :: BATHY_H0                         

# include "netcdf.inc"
      real time_scale
      integer iv,k,itrc,indWrk
      integer ncid, indx, varid,  ierr, lstr, lvar, latt, lenstr,       &
      start(2), count(2), ibuff(6),   nf_fread, checkdims
      character units*180,nomcv*30
      character inised_name*180
      real tmp(GLOBAL_2D_ARRAY)
      real tmp3d(GLOBAL_2D_ARRAY,ksdmin:ksdmax)

     tmp(PROC_IN_ARRAY)=0
     tmp3d(PROC_IN_ARRAY,ksdmin:ksdmax)=0
     z0sed(PROC_IN_ARRAY)=z0seduni
     dzsmax(PROC_IN_ARRAY)=dzsmaxuni
     ksmi(PROC_IN_ARRAY)=ksmiuni
     ksma(PROC_IN_ARRAY)=0
     hsed(PROC_IN_ARRAY)=-valmanq
     dzs(ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq
     cv_sed(-1:nv_tot,ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq
     c_sedtot(ksdmin:ksdmax,PROC_IN_ARRAY)=-valmanq
!
! Open initial conditions netCDF file for reading. Check that all
! spatial dimensions in that file are consistent with the model
! arrays, determine how many time records are available in the file
! and set record from which the dada will be read.
!
! The record is set as follows: (1) if only one time record is
! available in the file, then that record is used REGARDLESS of
! value of nrrec supplied from the parameter file; (2) if the
! file has multiple records and nrrec is positive, then nrrec is
! used, provided that nrrec is within the available records; 
! (3) if the file has multiple records and nrrec<0, then THE LAST 
! available record is used.

!      if (may_day_flag.ne.0) return      !-->  EXIT
      inised_name=filrepsed
      lstr=lenstr(inised_name)
      ierr=nf_open(inised_name(1:lstr), nf_nowrite, ncid)
      if (ierr.eq.nf_noerr) then
        ierr=checkdims (ncid, inised_name, lstr, indx)

        if (ierr.ne.nf_noerr) then
         goto 99
        elseif (indx.eq.0) then
          indx=1
        elseif (indx.gt.0 .and. nrrec.gt.0 .and. nrrec.le.indx) then
          indx=nrrec
        elseif (indx.gt.0 .and. nrrec.gt.indx) then
          write(stdout,'(/1x,A,I4,A/16x,A,I4,A/16x,3A/)')   &
                 'SEDINIT_FROMFILE ERROR: requested restart time record', &
                  nrrec, ' exceeds',  'number of available records', &
                  indx,'in netCDF file', '''',ininame(1:lstr),'''.'
          goto 99                                        !--> ERROR
        endif
      else
        write(stdout,'(/1x,2A/15x,3A)') 'SEDINIT_FROMFILE ERROR: Cannot ', &
                    'open netCDF file', '''', ininame(1:lstr) ,'''.'
        goto 99                                           !--> ERROR
      endif


! ksmi, ksma
      ierr=nf_inq_varid (ncid,'ksmi', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp, ncid, varid, indx, 0)
        ksmi(PROC_IN_ARRAY)=INT(tmp(PROC_IN_ARRAY))
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'ksmi', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'ksmi', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif
      
       ierr=nf_inq_varid (ncid,'ksma', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp, ncid, varid, indx, 0)
        ksma(PROC_IN_ARRAY)=INT(tmp(PROC_IN_ARRAY))

        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'ksma', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'ksma', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif

     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) ksmi(PROC_IN_ARRAY)=1
     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) ksma(PROC_IN_ARRAY)=0
!  DZS
      ierr=nf_inq_varid (ncid,'DZS', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp3d, ncid, varid, indx, 12)
         do k=ksdmin,ksdmax
            dzs(k,:,:)=tmp3d(:,:,k)
         enddo
        
        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) 'DZS', indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        MPI_master_only  write(stdout,1) 'DZS', inised_name(1:lstr)
        goto 99                                           !--> ERROR
      endif


! CVSED
 
      c_sedtot(:,:,:)=0.0_rsh
      do iv=-1,nv_tot 
       
       indWrk=indxT+ntrc_salt+ntrc_substot+iv+4+2
       nomcv=vname(1,indWrk)

       ierr=nf_inq_varid (ncid,nomcv, varid)
       if (ierr .eq. nf_noerr) then
        ierr=nf_fread (tmp3d, ncid, varid, indx, 12)

        do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=tmp3d(:,:,k)
            c_sedtot(k,:,:)=c_sedtot(k,:,:)+cv_sed(iv,k,:,:)*typart(iv)  
        enddo

        if (ierr .ne. nf_noerr) then
          MPI_master_only write(stdout,2) nomcv, indx, inised_name(1:lstr)
          goto 99                                         !--> ERROR
        endif
       
       else
        if (iv > nvpc .and. iv < nvp+1 ) then
         ! not constitutive particulate  variables (Initial in M/Msed converted to M/m3 sed)
         IF (irkm_var_assoc(iv) >0) THEN
           do k=ksdmin,ksdmax
             cv_sed(iv,k,:,:)=cini_sed(iv)*cv_sed(irkm_var_assoc(iv),k,:,:)
           enddo
         ELSE
           do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=cini_sed(iv)*c_sedtot(k,:,:)
           enddo
         END IF
         MPI_master_only  write(stdout,3) nomcv, inised_name(1:lstr)
        else if (iv > nvp) then
          ! dissolved variables (M/m3 EI)
          do k=ksdmin,ksdmax
            cv_sed(iv,k,:,:)=cini_sed(iv)             
          enddo
         MPI_master_only  write(stdout,3) nomcv,vname(1,indxT+ntrc_salt+iv), &
                                          inised_name(1:lstr)
        else
         MPI_master_only  write(stdout,1) nomcv, inised_name(1:lstr)
         goto 99                          !--> ERROR
        endif
       endif
      enddo
     
     do k=ksdmin,ksdmax
     WHERE (BATHY_H0(PROC_IN_ARRAY) == -valmanq) c_sedtot(k,PROC_IN_ARRAY)=-valmanq
     enddo

! Close input NetCDF file.
!
      ierr=nf_close(ncid)

  1   format(/1x,'SEDINIT_FROMFILE - unable to find variable:',    1x,A,  &
                                 /15x,'in input NetCDF file:',1x,A/)
  2   format(/1x,'SEDINIT_FROMFILE - error while reading variable:',1x, A, &
         2x,'at time record =',i4/15x,'in input NetCDF file:',1x,A/)
  3   format(/1x,'SEDINIT_FROMFILE - unable to find variable:',    1x,A,/10x,A, &
     &                            /15x,'in input NetCDF file:',1x,A,  &
         1x,'-> analytical value (cini_sed)'/)
      return
  99  may_day_flag=2
      return
      
  END SUBROUTINE sedinit_fromfile
    !!==============================================================================


  
  
! **TODO** code for CROCO
! SUBROUTINE bathy_actu_fromfile(h0)
! SUBROUTINE sed_obc_corflu
! SUBROUTINE sed_exchange_corflu_MARS
! SUBROUTINE sed_exchange_w2s_MARS
! SUBROUTINE sed_exchange_s2w_MARS 
! subroutine sed_exchange_hxe_MARS(iwhat,xh0,xssh) 
! SUBROUTINE sed_exchange_maskbedload_MARS
! SUBROUTINE sed_exchange_flxbedload_MARS
! SUBROUTINE sed_outsaverestart(h0)

#endif

END MODULE sed_MUSTANG_CROCO
