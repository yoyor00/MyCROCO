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
   !! **TODO** Code and clean subroutine not yet integrated (key_MARS) !  les routines pas encore programmes pour CROCO sont mises sous cle key_MARS  
#ifdef key_MARS
   PUBLIC sed_obc_corflu, bathy_actu_fromfile
   PUBLIC sed_outsaverestart,sed_exchange_s2w_MARS,sed_exchange_w2s_MARS
#endif
   PUBLIC sedinit_fromfile, sed_skinstress, sed_gradvit, sed_MUSTANG_settlveloc
# ifdef key_MUSTANG_bedload
   PUBLIC sed_bottom_slope
# endif

   PRIVATE

   
 CONTAINS
 
 !!===========================================================================================
 
  SUBROUTINE sed_MUSTANG_settlveloc(ifirst, ilast, jfirst, jlast,   &
                     h0fond, WATER_CONCENTRATION) 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_MUSTANG_settlveloc  ***
   !&E
   !&E ** Purpose : settling velocity computation
   !&E
   !&E ** Description : use arguments and common variable 
   !&E  arguments IN : 
   !&E         h0fond, : RESIDUAL_THICKNESS_WAT (m) 
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
   REAL(KIND=rsh), INTENT(IN)                         :: h0fond
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

  SUBROUTINE sed_gradvit(ifirst,ilast,jfirst,jlast,h0fond)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_gradvit  ***
   !&E
   !&E ** Purpose : calculation of the turbulence energy G  depending on hydro code ?
   !&E
   !&E ** Description : G= sqrt(turbulence dissipation/viscosity)
   !&E                 to be programmed using hydrodynamic knowledge
   !&E           using htot, RESIDUAL_THICKNESS_WAT (h0fond), RHOREF, sig, epn, nz ..
   !&E
   !&E     output : gradvit (in comMUSTANG)
   !&E
   !&E ** Called by :  step
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#ifdef key_MARS
   USE comvars3d,    ONLY : sig,nz,turb_nbeq, turb_2eq_option
   USE comturb,      ONLY : eps
#endif
#  include "mixing.h"
#  include "ocean3d.h"

   !! * Arguments 
   INTEGER, INTENT(IN)                        :: ifirst,ilast,jfirst,jlast
   REAL(KIND=rsh),INTENT(IN)                  :: h0fond

   !! * Local declarations
   INTEGER                            :: i,j,k
   REAL(KIND=rsh) :: dist_surf_on_bottom,nuw

      nuw=1.0e-6
      DO j=jfirst,jlast
#ifdef key_MARS
      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else
      DO i=ifirst,ilast
#endif

          IF(htot(i,j).GT.h0fond)  THEN
#ifdef key_MARS
            DO k=1,kmax
              IF (turb_nbeq == 2 .AND. turb_2eq_option == 2) THEN ! using k-epsilon closure scheme
                gradvit(k,i,j)=sqrt(eps(k,i,j)/nz(k,i,j))
              ELSE ! using nezu and nakagawa
                fzu=f_zu(h0(i,j),ssh(i,j),k,i,j)
                dist_surf_on_bottom=(-fzu+ssh(i,j))/(htot(i,j)+fzu-ssh(i,j))
                gradvit(k,i,j)=sqrt(ustarbot(i,j)**3._rsh/0.4_rsh/htot(i,j)/nz(k,i,j)*dist_surf_on_bottom)
              ENDIF
            ENDDO
#else

           DO k=1,N
            dist_surf_on_bottom=((z_w(i,j,N)-z_r(i,j,k))/(z_r(i,j,k)-z_w(i,j,0)))
            gradvit(k,i,j)=sqrt(ustarbot(i,j)**3._rsh/0.4_rsh/htot(i,j)/(nuw+epsilon_MUSTANG )*dist_surf_on_bottom) 
           END DO

           ! gradvit : G=sqrt( turbulence dissipation rate/ vertical viscosity coefficient)
           !  if  turbulence dissipation rate has not been already evaluated: 
           ! use empirical formula from   Nezu and Nakawaga (1993)
           !  turbulence dissipation_rate= ustarbot**3 /Karman/Htot * (distance from surface/distance from bottom)
#endif

          ENDIF
#ifdef key_MARS
        ENDIF
#endif
      ENDDO
      ENDDO

  END SUBROUTINE sed_gradvit

   !!==============================================================================
   SUBROUTINE sed_skinstress(ifirst,ilast,jfirst,jlast)
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_skinstress  ***
   !&E
   !&E ** Purpose : computes  bottom shear stress
   !&E
   !&E ** Description :  bottom shear stress related to current and also to waves
   !&E            to be programmed from hydrodynamic variables
   !&E                   RHOREF,  hx, hy, ssh, h0fond, u,uz, v, vz, hminfrot,hmaxfrot,
   !&E         using some parameters and variables of MUSTANG :
   !&E                    z0sed,  htncrit
   !&E               and Wave coupling
   !&E
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


   INTEGER, INTENT(IN)                        :: ifirst,ilast,jfirst,jlast
   INTEGER i,j,k
   REAL    vit2

# ifdef MPI
   REAL(KIND=rsh),DIMENSION(ifirst-2:ilast+1,jfirst-2:jlast+1)  :: Zr
   REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY) :: workexch  
# else
   REAL(KIND=rsh),DIMENSION(ifirst-1:ilast+1,jfirst-1:jlast+1)  :: Zr
# endif

# ifdef WAVE_OFFLINE
   REAL(KIND=rsh)                                   :: fws2ij,courant,alpha,beta,cosamb,sinamb
# endif
# ifdef key_tenfon_upwind   
   REAL(KIND=rsh)                                   :: cff4,cff1,cff2,cff3
# endif  
   
   !!--------------------------------------------------------------------------
   !! * Executable part

!-----------------------------------------------------------------------
!  Compute bottom skin-friction stress due to combined maximum wave
!  and current interaction
!-----------------------------------------------------------------------
!
#  ifdef BBL /*warning, d50 is constant (160microns) in the bustrw/bvstrw computation see bbl.F */
      do j=jfirst,jlast
        do i=ifirst,ilast
          tenfon(i,j)=sqrt( bustrw(i,j)**2 &
                        +bvstrw(i,j)**2)*RHOREF
#  ifdef WET_DRY AND MASKING
          tenfon(i,j)=tenfon(i,j)*rmask_wet(i,j)
#  endif
        enddo
      enddo

#  else /* else on #ifdef BBL */

# ifdef MPI
      DO j=jfirst-2,jlast+1
        DO i=ifirst-2,ilast+1
# else
      DO j=jfirst-1,jlast+1
        DO i=ifirst-1,ilast+1
#  endif
         Zr(i,j)=  max(z_r(i,j,1)-z_w(i,j,0),z0sed(i,j)+1.E-4)
       enddo
      enddo

# ifdef MPI
      DO j=jfirst-1,jlast
        DO i=ifirst-2,ilast
# else   
      DO j=jfirst,jlast
        DO i=ifirst-1,ilast
#  endif
          raphbx(i,j)=ABS(u(i+1,j,1,nnew))/(ABS(u(i+1,j,1,nnew))+epsilon_MUSTANG)
#ifdef key_MUSTANG_tenfonUbar
          vit2=SQRT(0.0625_rsh*(vbar(i,j+1,nnew)+vbar(i+1,j+1,nnew)+  &
                  vbar(i,j,nnew)+vbar(i+1,j,nnew))**2         &
                  +ubar(i+1,j,nnew)**2) *ubar(i+1,j,nnew)
#ifdef MPI
          if (float(i +ii*Lm) .GE. 0) then
#else
          if (float(i       ) .GE. 0) then
#endif
          frofonx(i,j)=0.16_rsh*(LOG( ( z_w(i,j,N)-z_w(i,j,0))/   &
                  (z0sed(i,j)*2.718)))**(-2)*vit2
          else
          frofonx(i,j)=0.
          endif

#else
          vit2=SQRT(0.0625_rsh*(v(i,j+1,1,nnew)+v(i+1,j+1,1,nnew)+  &
                  v(i,j,1,nnew)+v(i+1,j,1,nnew))**2         &
                  +u(i+1,j,1,nnew)**2) *u(i+1,j,1,nnew)   
          frofonx(i,j)=0.16_rsh*(LOG(0.5*(Zr(i+1,j)+Zr(i,j))/   &
                  z0sed(i,j)))**(-2)*vit2

#endif
       enddo
      enddo
# ifdef MPI
      DO j=jfirst-2,jlast
        DO i=ifirst-1,ilast
# else   
      DO j=jfirst-1,jlast
        DO i=ifirst,ilast
#  endif  
          raphby(i,j)=ABS(v(i,j+1,1,nnew))/(ABS(v(i,j+1,1,nnew))+epsilon_MUSTANG)
#ifdef key_MUSTANG_tenfonUbar
          vit2=SQRT(0.0625_rsh*(ubar(i+1,j,nnew)+ubar(i+1,j+1,nnew)+  &
                  ubar(i,j+1,nnew)+ubar(i,j,nnew))**2         &
                  +vbar(i,j+1,nnew)**2) *vbar(i,j+1,nnew)
#ifdef MPI
          if (float(j +jj*Mm) .GE. 0) then
#else
          if (float(j       ) .GE. 0) then
#endif
          frofony(i,j)=0.16_rsh*(LOG( (z_w(i,j,N)-z_w(i,j,0))  /   &
                  (z0sed(i,j)*2.718)))**(-2)*vit2
          else
          frofony(i,j)=0.
          endif
#else
          vit2=SQRT(0.0625_rsh*(u(i+1,j,1,nnew)+u(i+1,j+1,1,nnew)+  &
                  u(i,j+1,1,nnew)+u(i,j,1,nnew))**2         &
                  +v(i,j+1,1,nnew)**2) *v(i,j+1,1,nnew)
          frofony(i,j)=0.16_rsh*(LOG(0.5*(Zr(i,j+1)+Zr(i,j))/   &
                  z0sed(i,j)))**(-2)*vit2
#endif
       enddo
      enddo

# ifdef MPI
      DO j=jfirst-1,jlast
        DO i=ifirst-1,ilast
# else   
      DO j=jfirst,jlast
        DO i=ifirst,ilast
#  endif  
# ifdef key_tenfon_upwind
            cff1=0.5*(1.0+SIGN(1.0,frofonx(i,j)))
            cff2=0.5*(1.0-SIGN(1.0,frofonx(i,j)))
            cff3=0.5*(1.0+SIGN(1.0,frofonx(i-1  ,j)))
            cff4=0.5*(1.0-SIGN(1.0,frofonx(i-1 ,j)))
            tenfonx(i,j)=cff3*(cff1*frofonx(i-1,j)+                     &
                       cff2*0.5*(frofonx(i-1,j)+frofonx(i,j)))+    &
                       cff4*(cff2*frofonx(i,j)+                   &
                       cff1*0.5*(frofonx(i-1,j)+frofonx(i,j)))

            cff1=0.5*(1.0+SIGN(1.0,frofony(i,j)))
            cff2=0.5*(1.0-SIGN(1.0,frofony(i,j)))
            cff3=0.5*(1.0+SIGN(1.0,frofony(i,j-1)))
            cff4=0.5*(1.0-SIGN(1.0,frofony(i,j-1)))
            tenfony(i,j)=cff3*(cff1*frofony(i,j-1)+                      &
                      cff2*0.5*(frofony(i,j-1)+frofony(i,j)))+      &
                       cff4*(cff2*frofony(i,j)+                     &
                       cff1*0.5*(frofony(i,j-1)+frofony(i,j)))

            tenfonc(i,j)=SQRT(tenfonx(i,j)**2+tenfony(i,j)**2)*(rho(i,j,1)+rho0)
# else
            tenfonc(i,j)=SQRT(((frofonx(i,j)*raphbx(i,j)+frofonx(i-1,j)        &
                     *raphbx(i-1,j))/(raphbx(i,j)                             &
                     +raphbx(i-1,j)+epsilon_MUSTANG))**2+((frofony(i,j)*raphby(i,j)   &
                     +frofony(i,j-1)*raphby(i,j-1))/(raphby(i,j)              &
                     +raphby(i,j-1)+epsilon_MUSTANG))**2)*(rho(i,j,1)+rho0)
# endif

# ifdef WAVE_OFFLINE
            fws2ij=fws2
            IF (l_fricwave .AND. Uwave(i,j)*Pwave(i,j) > 0.0_rsh .AND.  &
                     Pwave(i,j) > 0.001_rsh.AND.Uwave(i,j)>0.001_rsh) THEN
              fws2ij=0.5_rsh*1.39_rsh*(Uwave(i,j)*Pwave(i,j)/REAL(2.0_rlg*pi*z0sed(i,j),rsh))**(-0.52_rsh)
            ENDIF

            ! calculation of shear stress due to waves tenfonw 
            ! --------------------------
            tenfonw(i,j)=((rho(i,j,1)+rho0)*fws2ij*Uwave(i,j)**2)

            courant=SQRT(ubar(i,j,nnew)**2+vbar(i,j,nnew)**2)
            IF(tenfonc(i,j) > 0.0_rsh .AND. tenfonw(i,j) > 0.0_rsh .AND. courant> 0.0_rsh) THEN

            ! calculation of shear stress with the formula of Soulsby (1995)
            ! ======================================================
            ! calculation of  tenfonc (courrent) influenced by the waves
            ! ----------------------------------------------------
              tenfonc(i,j)=tenfonc(i,j)*(1+(1.2*(tenfonw(i,j)/(tenfonw(i,j)+tenfonc(i,j)))**3.2))

            ! calculating the difference in angle between the direction of the waves and the current
            ! ---------------------------------------------------------------------------
            ! calculating the direction of the current relative to the north
              alpha=ACOS(vbar(i,j,nnew)/courant)   ! en radians
            ! calculation of wave orientation relative to north
              beta=Dwave(i,j)   ! beta and Dwave in radians
            ! calculation of cos(alpha-beta) and sin(alpha-beta)
              cosamb=ABS(COS(alpha-beta))
              sinamb=ABS(SIN(alpha-beta))

            ! calculation of tenfon (waves + current)
            ! -----------------------------------
              tenfon(i,j)=SQRT((tenfonc(i,j)+tenfonw(i,j)*cosamb)**2 +(tenfonw(i,j)*sinamb)**2)

          ELSE
              tenfon(i,j)=tenfonw(i,j)+tenfonc(i,j)
          ENDIF


# else
          tenfon(i,j)=tenfonc(i,j)
# endif

# if defined WET_DRY && defined MASKING 
          tenfon(i,j)=tenfon(i,j)*rmask_wet(i,j)
#  endif
        enddo
      enddo

#  endif  /* end of #ifdef BBL */

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

#ifdef key_MARS
  SUBROUTINE bathy_actu_fromfile(h0)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE bathy_actu_fromfile  ***
   !&E
   !&E ** Purpose : actualization of bathymetry after a previous run with morphodynamic coupling 
   !&E              read filrepsed even if not l_repsed
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

   USE comvars2d,    ONLY : tsauv,t,saverestart_step,iscreenlog,l_restart_subs,iscreenlog
   USE parameters,   ONLY : riosh
   USE ionc4,        ONLY : ionc4_openr,ionc4_read_subxyt,ionc4_close

   !! * Arguments
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_h0),INTENT(INOUT)       :: h0

   !! * Local declarations
   INTEGER          :: i,j,k,iv,ivg
   REAL(KIND=riosh),DIMENSION(limin:limax,ljmin:ljmax)        ::  xind

   !!--------------------------------------------------------------------------   
   !! * Executable part


     CALL ionc4_openr(filrepsed,l_in_nc4par=.true.)

     CALL ionc4_read_subxyt(filrepsed,'H0_morpho',xind,limin,limax,ljmin,ljmax,1,0,0)
     DO j=ljmin,ljmax
     DO i=limin,limax
             h0(i,j)=xind(i,j)
     END DO
     END DO

     CALL ionc4_close(filrepsed)

  END SUBROUTINE bathy_actu_fromfile
    !!==============================================================================

  SUBROUTINE sed_obc_corflu
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_obc_corflu ***
   !&E
   !&E ** Purpose : treatment of horizontal flow corrections for the transport of sand in suspension
   !&E
   !&E ** Description : extrapolation at borders and MPI exchange between processors
   !&E
   !&E ** Called by : sed_MUSTANG_update
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used

   USE comvars2d,      ONLY : l_obc_cycl_x,l_obc_cycl_y
   USE parameters,     ONLY : liminm2,liminm1,limaxp1,limaxp2, &
                              ljminm2,ljminm1,ljmaxp1,ljmaxp2, &
                              imin,imax,jmin,jmax
   USE comsubstance,   ONLY : isand1,isand2
   USE_MPI toolmpi,    ONLY : ex_i_rsh,ex_j_rsh

   !! * Arguments

   !! * Local declarations
    INTEGER                              :: i,j,ivp


   !!--------------------------------------------------------------------------
   !! * Executable part

     ! hypothese gradient nul aux frontieres 
     IF_MPI (ljmax == jmax) THEN
#ifdef key_MARS 
        DO i=MAX0(limin,ig(jmax)+1),MIN0(limax,id(jmax)-1)
#else
        DO i=limin,limax
#endif
          DO ivp=isand1,isand2
             corflux(ivp,i,jmax)=corflux(ivp,i,jmax-1)
             corfluy(ivp,i,jmax)=corfluy(ivp,i,jmax-1)
          ENDDO
        ENDDO
     ENDIF_MPI
     IF_MPI (ljmin == jmin) THEN
#ifdef key_MARS 
        DO i=MAX0(limin,ig(jmin)+1),MIN0(limax,id(jmin)-1)
#else
        DO i=limin,limax
#endif
          DO ivp=isand1,isand2
             corflux(ivp,i,jmin+1)=corflux(ivp,i,jmin+2)
             corfluy(ivp,i,jmin+1)=corfluy(ivp,i,jmin+2)
          ENDDO
        ENDDO
     ENDIF_MPI
     IF_MPI (limin==imax) THEN
#ifdef key_MARS 
        DO j=MAX0(ljmin,jb(imax)+1),MIN0(ljmax,jh(imax)-1)
#else
        DO j=ljmin,ljmax
#endif
          DO ivp=isand1,isand2
             corflux(ivp,imax,j)=corflux(ivp,imax-1,j)
             corfluy(ivp,imax,j)=corfluy(ivp,imax-1,j)
          ENDDO
        ENDDO
     ENDIF_MPI
     IF_MPI (limin == imin) THEN
#ifdef key_MARS 
        DO j=MAX0(ljmin,jb(imin)+1),MIN0(ljmax,jh(imin)-1)
#else
        DO j=ljmin,ljmax
#endif
          DO ivp=isand1,isand2
             corflux(ivp,imin+1,j)=corflux(ivp,imin+2,j)
             corfluy(ivp,imin+1,j)=corfluy(ivp,imin+2,j)
          ENDDO
        ENDDO
     ENDIF_MPI

     ! les quatre coins
     IF_MPI (limax == imax .AND. ljmax == jmax) THEN
          DO ivp=isand1,isand2
             corflux(ivp,imax,jmax)=corflux(ivp,imax-1,jmax-1)
             corfluy(ivp,imax,jmax)=corfluy(ivp,imax-1,jmax-1)
          ENDDO
     ENDIF_MPI
     IF_MPI (limax==imax .AND. ljmin==jmin) THEN
          DO ivp=isand1,isand2
             corflux(ivp,imax,jmin+1)=corflux(ivp,imax-1,jmin+2)
             corfluy(ivp,imax,jmin+1)=corfluy(ivp,imax-1,jmin+2)
          ENDDO
     ENDIF_MPI
     IF_MPI (limin==imin .AND. ljmin==jmin) THEN
          DO ivp=isand1,isand2
             corflux(ivp,imin+1,jmin+1)=corflux(ivp,imin+2,jmin+2)
             corfluy(ivp,imin+1,jmin+1)=corfluy(ivp,imin+2,jmin+2)
          ENDDO
     ENDIF_MPI
     IF_MPI (limin==imin .AND. ljmax==jmax) THEN
          DO ivp=isand1,isand2
             corflux(ivp,imin+1,jmax)=corflux(ivp,imin+2,jmax-1)
             corfluy(ivp,imin+1,jmax)=corfluy(ivp,imin+2,jmax-1)
          ENDDO
     ENDIF_MPI

#ifdef key_toseelater
!  PROBLEME : corflux, corfluy se sont pas alloues de imin a imax et de jmin a jmax
!             il faut les echanger avec les processeurs extremes et les ranger dans m2 ou p2
!             voir comme dans tsobc3dapply par exemple
     ! condition de periodicite dans la direction x 
     IF(l_obc_cycl_x) THEN
        DO j=ljminp1,ljmax
          DO ivp=isand1,isand2
             corflux(ivp,imin+1,j)=corflux(ivp,imax-1,j)
             corfluy(ivp,imin+1,j)=corfluy(ivp,imax-1,j)
             corflux(ivp,imax,j)=corflux(ivp,imin+2,j)
             corfluy(ivp,imax,j)=corfluy(ivp,imin+2,j)
          ENDDO
        ENDDO
     ENDIF
     ! condition de periodicite dans la direction y
     IF(l_obc_cycl_y) THEN
        DO i=liminp1,limax
          DO ivp=isand1,isand2
             corflux(ivp,i,jmin+1)=corflux(ivp,i,jmax-1)
             corfluy(ivp,i,jmin+1)=corfluy(ivp,i,jmax-1)
             corflux(ivp,i,jmax)=corflux(ivp,i,jmin+2)
             corfluy(ivp,i,jmax)=corfluy(ivp,i,jmin+2)
          ENDDO
        ENDDO
     ENDIF
#endif


   PRINT_DBG*, 'FIN SED_OBC_corflu'

   
   END SUBROUTINE sed_obc_corflu


    !!==============================================================================

  SUBROUTINE sed_exchange_corflu_MARS
 
   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_corflu_MARS ***
   !&E
   !&E ** Purpose : treatment of horizontal flow corrections for the transport of sand in suspension
   !&E
   !&E ** Description : extrapolation at borders and MPI exchange between processors
   !&E
   !&E ** Called by : sed_MUSTANG_update
   !&E
   !&E ** External calls : 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE parameters,   ONLY : limin,limax,limaxp1,ljmin,ljmax,ljmaxp1
   USE_MPI toolmpi,    ONLY : ex_i_rsh,ex_j_rsh

   !! * Arguments

   !! * Local declarations


   !!--------------------------------------------------------------------------
   !! * Executable part
   

OMPMPI barrier
OMPMPI master
   CALL_MPI ex_i_rsh(0,1,nv_adv,limin,limaxp1,ljmin,ljmax  ,corflux(:,limin:limaxp1,ljmin:ljmax  ))
   CALL_MPI ex_j_rsh(0,1,nv_adv,limin,limax  ,ljmin,ljmaxp1,corfluy(:,limin:limax  ,ljmin:ljmaxp1))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(corflux,corfluy)


   PRINT_DBG*, 'FIN sed_exchange_corflu_MARS'

   
   END SUBROUTINE sed_exchange_corflu_MARS

  !!==============================================================================

  SUBROUTINE sed_exchange_w2s_MARS

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_w2s_MARS  ***
   !&E
   !&E ** Purpose :  MPI exchange of slip deposit flux between processors
   !&E
   !&E ** Description : used only if slopefac .NE. 0 (slip deposit if steep slope)
   !&E
   !&E ** Called by :  sed_MUSTANG_deposition
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limin,limax,limaxp1,ljminm1,ljmin,ljmax,ljmaxp1
 
   !! * Arguments


  !! * Local declarations

   !!---------------------------------------------------------------------------
   !! * Executable part


! slip depo
! echange MPI
    
OMPMPI barrier
OMPMPI master
    CALL_MPI ex_i_rsh(0,1,nvp ,liminm1,limaxp1,ljminm1,ljmaxp1,flx_w2s_corim1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_i_rsh(-1,0,nvp,liminm1,limaxp1,ljminm1,ljmaxp1,flx_w2s_corip1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_j_rsh(0,1,nvp,liminm1,limaxp1,ljminm1,ljmaxp1,flx_w2s_corjm1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_j_rsh(-1,0,nvp,liminm1,limaxp1,ljminm1,ljmaxp1,flx_w2s_corjp1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(flx_w2s_corim1,flx_w2s_corip1,flx_w2s_corjm1,flx_w2s_corjp1)


  END SUBROUTINE sed_exchange_w2s_MARS

  !!==============================================================================

  SUBROUTINE sed_exchange_s2w_MARS

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_MPI  ***
   !&E
   !&E ** Purpose :  MPI exchange of erosion lateral flux between processors
   !&E
   !&E ** Description : used only if coef_erolat .NE. 0. 
   !&E
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limin,limax,limaxp1,ljminm1,ljmin,ljmax,ljmaxp1
 
   !! * Arguments


  !! * Local declarations

   !!---------------------------------------------------------------------------
   !! * Executable part


! lateral erosion of dry cell
! echange MPI
    
OMPMPI barrier
OMPMPI master
    CALL_MPI ex_i_rsh(0,1,nv_adv+2 ,liminm1,limaxp1,ljminm1,ljmaxp1,flx_s2w_corim1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_i_rsh(-1,0,nv_adv+2,liminm1,limaxp1,ljminm1,ljmaxp1,flx_s2w_corip1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_j_rsh(0,1,nv_adv+2,liminm1,limaxp1,ljminm1,ljmaxp1,flx_s2w_corjm1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
    CALL_MPI ex_j_rsh(-1,0,nv_adv+2,liminm1,limaxp1,ljminm1,ljmaxp1,flx_s2w_corjp1(:,liminm1:limaxp1  ,ljminm1:ljmaxp1 ))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(flx_s2w_corim1,flx_s2w_corip1,flx_s2w_corjm1,flx_s2w_corjp1)


  END SUBROUTINE sed_exchange_s2w_MARS
    !!===========================================================================   
   subroutine sed_exchange_hxe_MARS(iwhat,xh0,xssh) 
     !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE sed_exchange_hxe_MARS  ***
  !&E
  !&E ** Purpose : MPI exchange  h0 and ssh for morpho betwwen processors
  !&E 
  !&E ** Description :
  !&E
  !&E ** Called by : morpho
  !&E
  !&E ** External calls : 
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
  USE parameters,   ONLY : liminm2,limaxp2,ljminm2,ljmaxp2
  USE comvars2d,      ONLY : hx,hy 

  !! * Arguments
  INTEGER,INTENT(IN)      :: iwhat
  REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_h0),INTENT(INOUT),OPTIONAL        :: xh0
  REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(INOUT),OPTIONAL :: xssh                         

  !! * Local declarations
  
  !!--------------------------------------------------------------------------
  !! * Executable part
  
   IF (iwhat==1) THEN
          ! exchange of h0 + ssh
OMPMPI barrier
OMPMPI master
       CALL_MPI ex_i_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,xh0(liminm2:limaxp2,ljminm2:ljmaxp2))
       CALL_MPI ex_j_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,xh0(liminm2:limaxp2,ljminm2:ljmaxp2))
       CALL_MPI ex_i_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,xssh(liminm2:limaxp2,ljminm2:ljmaxp2))
       CALL_MPI ex_j_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,xssh(liminm2:limaxp2,ljminm2:ljmaxp2))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(xssh,xh0)
          
   ELSE IF (iwhat==2) THEN
          ! exchange of hx,hy
OMPMPI barrier
OMPMPI master
     CALL_MPI ex_i_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,hx(liminm2:limaxp2,ljminm2:ljmaxp2))
     CALL_MPI ex_j_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,hx(liminm2:limaxp2,ljminm2:ljmaxp2))
     CALL_MPI ex_i_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,hy(liminm2:limaxp2,ljminm2:ljmaxp2))
     CALL_MPI ex_j_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,hy(liminm2:limaxp2,ljminm2:ljmaxp2))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(hx,hy)

   ENDIF

  PRINT_DBG*, 'END sed_exchange_hxe_MARS'
  END SUBROUTINE sed_exchange_hxe_MARS

  !!==============================================================================
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload && defined key_MPI_2D

  SUBROUTINE sed_exchange_maskbedload_MARS

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_maskbedload_MARS  ***
   !&E
   !&E ** Purpose :  exchange MPI mask
   !&E                 ATTENTION ON EST EN ZONE PARALLEL OMP
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limin,limax,limaxp1,ljminm1,ljmin,ljmax,ljmaxp1
 
   !! * Arguments


  !! * Local declarations

   !!---------------------------------------------------------------------------
   !! * Executable part


! bedload
! echange MPI
    
OMPMPI barrier
OMPMPI master
   CALL_MPI ex_i_rsh(-1,1,1,liminm1,limaxp1,ljminm1,ljmaxp1,sedimask_h0plusxe(liminm1:limaxp1,ljminm1:ljmaxp1))
   CALL_MPI ex_j_rsh(-1,1,1,liminm1,limaxp1,ljminm1,ljmaxp1,sedimask_h0plusxe(liminm1:limaxp1,ljminm1:ljmaxp1))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(sedimask_h0plusxe)


  END SUBROUTINE sed_exchange_maskbedload_MARS
  !!==============================================================================
  SUBROUTINE sed_exchange_flxbedload_MARS

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_exchange_flxbedload_MARS  ***
   !&E
   !&E ** Purpose :  exchange MPI bedload fluxes
   !&E                 ATTENTION ON EST EN ZONE PARALLEL OMP
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   USE_MPI toolmpi,  ONLY : ex_i_rsh,ex_j_rsh
   USE parameters,   ONLY : liminm1,limin,limax,limaxp1,ljminm1,ljmin,ljmax,ljmaxp1
 
   !! * Arguments


  !! * Local declarations

   !!---------------------------------------------------------------------------
   !! * Executable part


! bedload
! echange MPI
    
OMPMPI barrier
OMPMPI master
   CALL_MPI ex_i_rsh(-1,1,nvp,liminm1,limaxp1,ljmin,ljmax,flx_bx(:,liminm1:limaxp1,ljmin:ljmax))
   CALL_MPI ex_j_rsh(-1,1,nvp,limin,limax,ljminm1,ljmaxp1,flx_by(:,limin:limax,ljminm1:ljmaxp1))
OMPMPI end master
OMPMPI barrier
OMPMPI flush(flx_bx,flx_by)


  END SUBROUTINE sed_exchange_flxbedload_MARS
  
  ! end version V2+bedload+MPI2D
#endif
  !!==============================================================================
 
  SUBROUTINE sed_outsaverestart(h0)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE sed_outsaverestart  ***
   !&E
   !&E ** Purpose : saves sediment-related field for future runs
   !&E
   !&E ** Description : netcdf files
   !&E
   !&E ** Called by :  step 
   !&E
   !&E ** External calls : 
   !
   !objet    :  gestion de la reprise de calcul a partir d un etat
         !           enregistre dans un fichier (netcdf ou binaire acces direct)
         !           lors d un run anterieur
   !
   !     appelant         : stepu
         !                :
   !     appeles          : ionc4_createfile, ionc4_torigine, 
   !     appeles          : ionc4_createvar, ionc4_createvar
   !     appeles          : ionc4_createvar, ionc4_write_time
   !     appeles          : ionc4_write_xy, ionc4_write_xyt, ionc4_write_zxyt
         !                :
   !     canal d ecriture :
         !                :
   !     canal de lecture : 26, filrep
         !                :
   !     cle cpp          : key_notra
         !                :
   !     variables in     : t, u, uz, xu, v, vz, xv, sal, temp, tz
   !     variables in     : ect, boz, nz, kz, wz, vish_xe, vish_phi
   !     variables in     : filsauvsed, lon2d,lat2d, ntra
   !     variables in     : cdate, h0,sig
         !                :
   !     variables out    : 
         !                :
   !     variables inout  : tsauv
         !                :
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
   !!   ------------
   USE comvars2d,    ONLY :  lon2d_g,lat2d_g,lon2d,lat2d,l_saverestart_bydate,date_ref, &
                             tsauv,cdate,saverestart_step,name_out_h0,name_out_xe,t
   USE_MPI toolmpi,  ONLY : ADD_TAG
   USE output,       ONLY : l_out_nc4par
   USE ionc4,       ONLY : ionc4_openr,ionc4_read_subxyt,ionc4_read_subzxyt,ionc4_close,  &
                           ionc4_createfile,ionc4_gatt,ionc4_createvar,ionc4_write_sig,&
                           ionc4_vatt_fill,ionc4_torigine,ionc4_write_time,ionc4_write_xyt,&
                           ionc4_write_zxyt,ionc4_sync
                           
   !! * Arguments
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_h0),INTENT(IN)       :: h0                         

   !! * Local declarations
   !!   ------------------

   CHARACTER(LEN=8)      :: nomc
   CHARACTER(LEN=13)     :: prefsauv
   CHARACTER(LEN=lchain) :: chain,filsauvsed
   INTEGER               :: i,j,k,l,iv,ifois
   REAL(KIND=riosh),DIMENSION(limin:limax,ljmin:ljmax)          :: xind
   REAL(KIND=riosh),DIMENSION(ksdmin:ksdmax,limin:limax,ljmin:ljmax)   :: tlits
   REAL(KIND=riosh),DIMENSION(0:ksdmax+1)                       :: sigsed

   DATA ifois /0/
   
   !! * Executable part
   !!   ---------------

   PRINT_DBG*, 'enter sed_outsaverestart'

   ! initialisation of sigsed
   DO k=ksdmin,ksdmax
      sigsed(k)=k
   END DO
   sigsed(0)=0
   sigsed(ksdmax+1)=0
   
   ! creation
   ! ********

   IF (ifois.EQ.0 .AND. .NOT.(l_saverestart_bydate)) THEN
     ifois=1
     filsauvsed='./save_sedim.nc'
     IF_MPI(.NOT. l_out_nc4par) CALL ADD_TAG(filsauvsed,LEN(filsauvsed))
     IF ( l_out_nc4par ) THEN
        CALL ionc4_createfile(filsauvsed,lon2d_g,lat2d_g,imin,imax,1,          &
                              jmin,jmax,1,ksdmin,ksdmax,1,ksdmax,0,l_out_nc4par = .TRUE.)
     ELSE
        CALL ionc4_createfile(filsauvsed,lon2d,lat2d,limin,limax,1,          &
                              ljmin,ljmax,1,ksdmin,ksdmax,1,ksdmax,0)
     ENDIF

   ! save attribute of the spatial extension of the domain
   ! -----------------------------------------------------
     CALL ionc4_gatt(filsauvsed,'global_imin',imin)
     CALL ionc4_gatt(filsauvsed,'global_imax',imax)
     CALL ionc4_gatt(filsauvsed,'global_jmin',jmin)
     CALL ionc4_gatt(filsauvsed,'global_jmax',jmax)

   ! indispensable pour concatenation pas terrible a revoir
     CALL ionc4_createvar(filsauvsed,'SIG',' ','sigma variable',&
#if defined key_siggen || defined key_gencoord
                          standard_name= 'ocean_s_variable',&
#else
                          standard_name= 'ocean_sigma_variable',&
#endif
                          l_auxcoordinate=.false.,&
                          valid_min= -1.0_riosh,valid_max= 0.0_riosh,fill_value= rg_valmanq_io,&
                          dims= 'z', l_out_nc4par = l_out_nc4par)
     CALL ionc4_write_sig(filsauvsed,'SIG',sigsed,ksdmax)

     CALL ionc4_torigine(filsauvsed,date_ref)
     CALL ionc4_createvar(filsauvsed,'ksmi','s.o','ksmi',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min=0.0_riosh,valid_max=10.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'ksma','s.o','ksma',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min=0.0_riosh,valid_max=1000.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'dzsmax','s.o','ksma',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min=0.0_riosh,valid_max=1.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'dzs','m','dzs',fill_value=rg_valmanq_io, &
                          dims='zxyt',valid_min=0.0_riosh,valid_max=100.0_riosh,l_out_nc4par=l_out_nc4par)
     DO iv=-1,nv_adv     ! pour T et S
       IF (IV>0) WRITE(nomc,10) iv
       IF (iv==-1) nomc='concTEM'
       IF (iv==0) nomc='concSAL'
       CALL ionc4_createvar(filsauvsed,nomc,'u.s.i','concentration ',&
                            fill_value=rg_valmanq_io,dims='zxyt',l_out_nc4par=l_out_nc4par)
       IF (IV>0) WRITE(nomc,20) iv
       IF (iv==-1) nomc='fxw2sTEM'
       IF (iv==0) nomc='fxw2sSAL'
       CALL ionc4_createvar(filsauvsed,nomc,'u.s.i','flx_w2s ',&
                            fill_value=rg_valmanq_io,dims='xyt',l_out_nc4par=l_out_nc4par)
     END DO 
   ENDIF

   IF (l_saverestart_bydate.AND.t.GE.tsauv) THEN
     prefsauv='./save_sedim_'
     filsauvsed=prefsauv//cdate(7:10)//cdate(4:5)//  &
                    cdate(1:2)//cdate(12:13)//cdate(15:16)//cdate(18:19)//'.nc'
     IF_MPI(.NOT. l_out_nc4par) CALL ADD_TAG(filsauvsed,LEN(filsauvsed))
     IF ( l_out_nc4par ) THEN
        CALL ionc4_createfile(filsauvsed,lon2d_g,lat2d_g,imin,imax,1,          &
                              jmin,jmax,1,ksdmin,ksdmax,1,ksdmax,0,l_out_nc4par = .TRUE.)
     ELSE
        CALL ionc4_createfile(filsauvsed,lon2d,lat2d,limin,limax,1,          &
                              ljmin,ljmax,1,ksdmin,ksdmax,1,ksdmax,0)
     ENDIF

   ! save attribute of the spatial extension of the domain
   ! -----------------------------------------------------
     CALL ionc4_gatt(filsauvsed,'global_imin',imin)
     CALL ionc4_gatt(filsauvsed,'global_imax',imax)
     CALL ionc4_gatt(filsauvsed,'global_jmin',jmin)
     CALL ionc4_gatt(filsauvsed,'global_jmax',jmax)

   ! indispensable pour concatenation pas terrible a revoir
     CALL ionc4_createvar(filsauvsed,'SIG',' ','sigma variable',&
#if defined key_siggen || defined key_gencoord
                          standard_name= 'ocean_s_variable',&
#else
                          standard_name= 'ocean_sigma_variable',&
#endif
                          l_auxcoordinate=.false.,&
                          valid_min= -1.0_riosh,valid_max= 0.0_riosh,fill_value= rg_valmanq_io,&
                          dims= 'z', l_out_nc4par = l_out_nc4par)
     CALL ionc4_write_sig(filsauvsed,'SIG',sigsed,ksdmax)

     CALL ionc4_torigine(filsauvsed,date_ref)
     CALL ionc4_createvar(filsauvsed,'ksmi','s.o','ksmi',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min= 0.0_riosh,valid_max= 10.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'ksma','s.o','ksma',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min= 0.0_riosh,valid_max= 1000.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'dzsmax','s.o','ksma',fill_value=rg_valmanq_io, &
                          dims='xyt',valid_min=0.0_riosh,valid_max=1.0_riosh,l_out_nc4par=l_out_nc4par)
     CALL ionc4_createvar(filsauvsed,'dzs','m','dzs',fill_value=rg_valmanq_io, &
                          dims='zxyt',valid_min= 0.0_riosh,valid_max= 100.0_riosh,l_out_nc4par=l_out_nc4par)
     DO iv=-1,nv_adv     ! pour T et S
       IF (IV>0) WRITE(nomc,10) iv
       IF (iv==-1) nomc='concTEM'
       IF (iv==0) nomc='concSAL'
       CALL ionc4_createvar(filsauvsed,nomc,'u.s.i','concentration ',&
                            fill_value=rg_valmanq_io,dims='zxyt',l_out_nc4par=l_out_nc4par)
       IF (IV>0) WRITE(nomc,20) iv
       IF (iv==-1) nomc='fxw2sTEM'
       IF (iv==0) nomc='fxw2sSAL'
       CALL ionc4_createvar(filsauvsed,nomc,'u.s.i','flx_w2s ',&
                            fill_value=rg_valmanq_io,dims='xyt',l_out_nc4par=l_out_nc4par)
     END DO 
     IF (l_morphocoupl) THEN
       CALL ionc4_createvar(filsauvsed,'H0_morpho','u.s.i','H0_modify_morpho',&
                            fill_value=rg_valmanq_io,dims='xyt',l_out_nc4par=l_out_nc4par)       
     ENDIF
   ENDIF

   ! ecriture
   ! ********

   IF (t.GE.tsauv) THEN
     CALL ionc4_write_time(filsauvsed,1,t)
     DO j=ljmin,ljmax
       DO i=limin,limax
         xind(i,j)=REAL(ksmi(i,j))
       END DO
     END DO
     CALL ionc4_write_xyt(filsauvsed,'ksmi',xind,limin,limax,ljmin,ljmax,1,l_out_nc4par=l_out_nc4par)

     DO j=ljmin,ljmax
       DO i=limin,limax
         xind(i,j)=REAL(ksma(i,j))
       END DO
     END DO
     CALL ionc4_write_xyt(filsauvsed,'ksma',xind,limin,limax,ljmin,ljmax,1,l_out_nc4par=l_out_nc4par)
 
     DO j=ljmin,ljmax
       DO i=limin,limax
         xind(i,j)=REAL(dzsmax(i,j))
       END DO
     END DO
     CALL ionc4_write_xyt(filsauvsed,'dzsmax',xind,limin,limax,ljmin,ljmax,1,l_out_nc4par=l_out_nc4par)

     IF (l_morphocoupl) THEN
       DO j=ljmin,ljmax
         DO i=limin,limax
           xind(i,j)=REAL(h0(i,j))
         END DO
       END DO
       CALL ionc4_write_xyt(filsauvsed,'H0_morpho',xind,limin,limax,ljmin,ljmax,1,l_out_nc4par=l_out_nc4par)
     ENDIF
 
     DO j = ljmin,ljmax
       DO i = limin,limax
         DO k = 1,ksdmax
           tlits(k,i,j) = dzs(k,i,j)
         END DO
       END DO
     END DO
     CALL ionc4_write_zxyt(filsauvsed,'dzs',tlits,limin,limax,ljmin,ljmax,ksdmin,ksdmax,1,l_out_nc4par=l_out_nc4par)

     DO iv=-1,nv_adv     ! pour T et S
       IF (IV>0) WRITE(nomc,10) iv
       IF (iv==-1) nomc='concTEM'
       IF (iv==0) nomc='concSAL'
       DO j = ljmin,ljmax
         DO i = limin,limax
           DO k = 1,ksdmax
             tlits(k,i,j) = cv_sed(iv,k,i,j)
           ENDDO
         ENDDO
       ENDDO
       CALL ionc4_write_zxyt(filsauvsed,nomc,tlits,limin,limax,ljmin,ljmax,ksdmin,ksdmax,1,l_out_nc4par=l_out_nc4par)
       IF (IV>0) WRITE(nomc,20) iv
       IF (iv==-1) nomc='fxw2sTEM'
       IF (iv==0) nomc='fxw2sSAL'
       DO j = ljmin,ljmax
         DO i = limin,limax
           xind(i,j) = flx_w2s(iv,i,j)
         ENDDO
       ENDDO
       CALL ionc4_write_xyt(filsauvsed,nomc,xind,limin,limax,ljmin,ljmax,1,l_out_nc4par=l_out_nc4par)       
     ENDDO
     tsauv=t+saverestart_step

     IF (l_saverestart_bydate) THEN
       CALL ionc4_close(filsauvsed)
     ELSE
       ! To write the data on the disk and not loose data in case of run crash
       CALL ionc4_sync(filsauvsed)
     ENDIF

   ENDIF

   10 FORMAT('concen',i2.2)
   20 FORMAT('flxw2s',i2.2)

   PRINT_DBG*, 'exit sed_outsaverestart'

  END SUBROUTINE sed_outsaverestart
#endif

 !=========================================================================== 

#endif

END MODULE sed_MUSTANG_CROCO
