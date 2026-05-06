MODULE trajectools

   !!======================================================================
   !!                   ***  MODULE trajectools  ***
   !! Trajectories dynamics:  routines for trajectory estimates
   !! Several routines called by traject3d
   !!======================================================================
#include "cppdefs.h"
#include "toolcpp.h"

#if defined LAGRANGIAN || defined DEB_IBM

    !! * Modules used
    USE comtraj, ONLY       : imin,imax,jmin,jmax,kmax,rsh,rlg,riolg, &
                              htx,hty,hc_sig,dsigw,dsigu,valmanq,wz,&
                              dcusds,dcwsds,lonwest,latsouth,dlonr,dlatr,type_position
    USE module_lagrangian !, ONLY : sc_r,sc_w,h,zeta,kstp,theta_s,theta_b,We                        

    IMPLICIT NONE
    PRIVATE

    !! * Accessibility
    PUBLIC h0int,xeint,ksupkinf,loc_h0,kzprofile,splint,uint,vint,wint,dksdzint,interpvit
    PUBLIC update_htot,compute_dsig_dcuds,update_wz,set_htot_bc
    PUBLIC siggentoz,ztosiggen,hc_sigint,h_int_siggen,define_pos,CSF
    PUBLIC tool_ind2lat,tool_ind2lon,lonlat2ij,tool_latlon2i,tool_latlon2j

 CONTAINS


  !!======================================================================
  SUBROUTINE define_pos(pos)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE define_pos  ***
    !&E
    !&E ** Purpose : calcul indexes from global position for use in tab
    !&E
    !&E ** Called by :  LAGRANGIAN_update, avance, LAGRANGIAN_init, ibm_3d
    !&E                 get_Xdeb, ibm_nycth_mig, ibm_parameter_init,fish_move
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) For coupling LAGRANGIAN with CROCO
    !&E         
    !&E---------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(type_position),INTENT(inout) :: pos

    !compute indexes in local
    pos%idx = INT(pos%xp) - iminmpi + 1
    pos%idy = INT(pos%yp) - jminmpi + 1

    !indexes in float in local 
    pos%idx_r = pos%xp - REAL(iminmpi,kind=rsh) + 1.0_rsh
    pos%idy_r = pos%yp - REAL(jminmpi,kind=rsh) + 1.0_rsh
 
  END SUBROUTINE define_pos 



  !!======================================================================
  SUBROUTINE compute_dsig_dcuds
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE compute_dsig_dcuds  ***
    !&E
    !&E ** Purpose : Compute sigma parameters needed by traj 
    !&E
    !&E ** Called by :  LAGRANGIAN_init
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) For coupling LAGRANGIAN with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Local declaration
    INTEGER :: k 

    !dsigu,dsigw
    DO k=1,kmax
        dsigu(k) = sc_w(k) - sc_w(k-1)
    END DO
    DO k=1,kmax-1
        dsigw(k) = sc_r(k+1)-  sc_r(k)
    END DO

    ! sc_r(kmax+1)=-sc_r(kmax)
    dsigw(kmax) = -sc_r(kmax) - sc_r(kmax)

    !update hc_sig 
    hc_sig(:,:) = hc
    DO k=1,kmax
        dcusds(k) = (Cs_w(k) - Cs_w(k-1))/(sc_w(k) - sc_w(k-1))
    END DO
    
    DO k=1,kmax-1
        dcwsds(k) = (Cs_r(k+1) - Cs_r(k))/(sc_r(k+1) - sc_r(k))
    END DO
    dcwsds(kmax) = (0.0_rsh - Cs_r(kmax))/(0.0_rsh - sc_r(kmax)) 

  END SUBROUTINE compute_dsig_dcuds 
  


  !!======================================================================
  SUBROUTINE set_htot_bc(Istr,Iend,Jstr,Jend,IstrU,JstrV)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE set_htot_bc  ***
    !&E
    !&E ** Purpose : Update htot for Lagrangian
    !&E
    !&E ** Called by :  LAGRANGIAN_init, LAGRANGIAN_update
    !&E
    !&E ** External calls : exchange_u2d_tile,exchange_v2d_tile
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) For coupling LAGRANGIAN with CROCO
    !&E
    !&E---------------------------------------------------------------------

    INTEGER,     INTENT(in)  :: Istr,Iend,Jstr,Jend,IstrU,JstrV 
    INTEGER                  :: i,j
    INTEGER                  :: ISTRm,ISTRp,IENDp,JSTRm,JSTRp,JENDp

#ifndef EW_PERIODIC
    IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
            htx(Iend+1,j)=htx(Iend,j)
        END DO
        DO j=JstrV,Jend
            hty(Iend+1,j)=hty(Iend,j)
        END DO
    END IF

    IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
            htx(IstrU-1,j)=htx(IstrU,j)
        END DO
        DO j=JstrV,Jend
            hty(Istr-1,j)=hty(Istr,j)
        END DO
    END IF
#endif
#ifndef NS_PERIODIC
    IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
            htx(i,Jend+1) =htx(i,Jend)
        END DO
        DO i=Istr,Iend
            hty(i,Jend+1) =hty(i,Jend)
        END DO
    END IF
    IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
            htx(i,Jstr-1)=htx(i,Jstr)
        END DO
        DO i=Istr,Iend
            hty(i,JstrV-1)=hty(i,JstrV)
        END DO
    END IF
#endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
    ISTRm = Istr - 1
    ISTRp = Istr + 1
    IENDp = Iend + 1
    JSTRm = Jstr - 1
    JSTRp = Jstr + 1
    JENDp = Jend + 1

    IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        htx(Istr,JSTRm) = 0.5*(htx(ISTRp,JSTRm) + htx(Istr,Jstr))
        hty(ISTRm,Jstr) = 0.5*(hty(Istr,Jstr)   + hty(ISTRm,JSTRp))
    ENDIF
    IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        htx(IENDp,JSTRm) = 0.5*(htx(IENDp,Jstr)  + htx(Iend,JSTRm))
        hty(IENDp,Jstr)  = 0.5*(hty(IENDp,JSTRp) + hty(Iend,Jstr))
    ENDIF
    IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        htx(Istr,JENDp)  = 0.5*(htx(Istr,Jend)  + htx(ISTRp,JENDp))
        hty(ISTRm,JENDp) = 0.5*(hty(ISTRm,Jend) + hty(Istr,JENDp))
    ENDIF
    IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        htx(IENDp,JENDp) = 0.5*(htx(IENDp,Jend) + htx(Iend,JENDp))
        hty(IENDp,JENDp) = 0.5*(hty(IENDp,Jend) + hty(Iend,JENDp))
    ENDIF
#endif

! exchange variables
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
    call exchange_u2d_tile (Istr,Iend,Jstr,Jend,htx)
    call exchange_v2d_tile (Istr,Iend,Jstr,Jend,hty)
#endif
  END SUBROUTINE set_htot_bc


  !!======================================================================
  SUBROUTINE update_htot(limin,limax,ljmin,ljmax,liminu,ljminv)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE update_htot  ***
    !&E
    !&E ** Purpose : compute water level 
    !&E
    !&E ** Called by : LAGRANGIAN_init, LAGRANGIAN_update
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) For coupling LAGRANGIAN with CROCO
    !&E
    !&E---------------------------------------------------------------------

    !! * Arguments
    INTEGER,INTENT(in) :: limin,limax,ljmin,ljmax,liminu,ljminv

    !! * Local declarations
    INTEGER :: i,j

    DO j=ljmin,ljmax
        DO i=liminu,limax
            htx(i,j)= 0.5*(h(i-1,j)+z_w(i-1,j,N)+ &
                           h(i  ,j)+z_w(i  ,j,N))
        END DO
    END DO
    DO j=ljminv,ljmax
        DO i=limin,limax
            hty(i,j)=0.5*(h(i,j-1)+z_w(i,j-1,N)+  &
                          h(i,j  )+z_w(i,j  ,N))
        END DO
    END DO
  
  END SUBROUTINE update_htot



  !!======================================================================
  SUBROUTINE update_wz(limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE update_wz  ***
    !&E
    !&E ** Purpose : compute vertical current component 
    !&E
    !&E ** Description : We is a water flux through vertical interfaces (grid box) ; units is in m3/s
    !&E                  We is divided by dx*dy to get a velocity in m/s
    !&E                  Then We is divided by Hz_mars (ie Hz in mars framework == dz/dsig = Hz_croco/dsig) to have a Wz in sigma.s-1
    !&E                  (what is expected by the lagrangian module to compute vertical displacement)
    !&E
    !&E ** Called by   : LAGRANGIAN_init, LAGRANGIAN_update
    !&E
    !&E ** History :
    !&E       !  2024    (M. Caillaud) For coupling LAGRANGIAN with CROCO
    !&E
    !&E---------------------------------------------------------------------

    !! * Arguments
    INTEGER,   INTENT( in ) :: limin,limax,ljmin,ljmax

    !! * Local declarations
    INTEGER                 :: i,j,k
    
    DO k=1,kmax
        DO j=ljmin,ljmax
            DO i=limin,limax
                wz(i,j,k)=we(i,j,k)/(om_r(i,j)*on_r(i,j))
                wz(i,j,k)=wz(i,j,k)*dsigw(k)/Hz(i,j,k)
            END DO
        END DO
    END DO
  END SUBROUTINE update_wz



  !!======================================================================
  FUNCTION uint(xpos,ypos,spos,u,xe,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION uint  ***
    !&E
    !&E ** Purpose : Linear interpolation of zonal velocity (uint)
    !&E              at the location (xpos,ypos,spos)
    !&E
    !&E ** Called by      : avance.F90
    !&E
    !&E ** External calls : f_lag_sigz_uv,h_int_siggen,interpvit,loc_h0
    !&E
    !&E ** History :
    !&E       !  08-2002 (F. Dumas, P. Lazure) 
    !&E       !  01-2007 (V. Garnier) correction bug declaration px,py
    !&E       !  01-2011 (M. Huret) Adaptation to siggen (scale factor for vertical interpolation)
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                                     INTENT( in ) :: xpos,ypos,spos
    INTEGER,                                            INTENT( in ) :: limin,limax,ljmin,ljmax
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax),    INTENT( in ) :: u
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),         INTENT( in ) :: xe

    REAL(KIND=rsh)                                                   :: uint

    !! * Local declarations
    INTEGER                    :: igg,idd,jhh,jbb,k,kuv,kuvm
    INTEGER                    :: hlb,hlt,hrb,hrt  

    REAL(KIND=rsh)             :: uintm,uintp,px,py
    REAL(KIND=rsh)             :: f_lag_hzp,f_lag_hzm,f_lag_hz
    REAL(KIND=rsh)             :: xe_lag,hc_sig_lag,h0_lag

    !!----------------------------------------------------------------------
    !! * Executable part

    !---- reperage de la maille ou se situe la particule pour les U
    ! idd, igg relatifs a la grille U
    ! jbb ,jhh relatifs a la grille V
    ! xpos,ypos relatifs a la grille rho
    ! location of ssh : xpos=REAL(i), ypos=REAL(j)
    igg = nint(xpos)
    idd = igg+1
    jbb = int(ypos)
    jhh = jbb+1

    idd = min(max(idd,0),imax)
    igg = min(max(igg,0),imax)
    jbb = min(max(jbb,0),jmax)
    jhh = min(max(jhh,0),jmax)

    ! calcul des positions px et py dans la maille U
    px = xpos - real(igg) - 0.5_rsh
    py = ypos - real(jbb) 
 
    DO k = kmax,2,-1
        IF (spos <= sc_r(k).and.spos >= sc_r(k-1)) THEN
            kuv  = k
            kuvm = k-1
        END IF
    END DO
    IF (spos <= sc_r(1)) THEN
        kuv  = 1
        kuvm = 0
    END IF
    IF (spos >= sc_r(kmax)) THEN
        kuv  = kmax
        kuvm = kmax
    END IF

    !-- interpolation proprement dite sur la nappe sigma au-dessus
    uintp = interpvit(px,py,u(igg,jbb,kuv),u(idd,jbb,kuv), &
                          u(igg,jhh,kuv),u(idd,jhh,kuv), &
                          htx(igg,jbb),htx(idd,jbb),     &
                          htx(igg,jhh),htx(idd,jhh))

    !-- interpolation proprement dite sur la nappe sigma au-dessous
    IF (kuvm /= 0) THEN
        uintm = interpvit(px,py,u(igg,jbb,kuvm),u(idd,jbb,kuvm), &
                              u(igg,jhh,kuvm),u(idd,jhh,kuvm), &
                              htx(igg,jbb),htx(idd,jbb),       &
                              htx(igg,jhh),htx(idd,jhh))
    ELSE
        uintm = 0.0_rsh
    END IF

    !-- interpolation sur la verticale entre la valeur au-dessus (uintp)
    !-- et la valeur au-dessous (uintm)
    IF (spos >= sc_r(kmax)) THEN
        uint = uintm
    ELSE IF(kuv == 1) THEN
        uint = ((sc_r(kuv) - spos)*uintm + (spos + 1.0_rsh)*uintp)/(sc_r(kuv) + 1.0_rsh)
    ELSE
        f_lag_hzp = 1.0_rsh
        f_lag_hzm = 1.0_rsh
        f_lag_hz  = 1.0_rsh
        CALL loc_h0(xpos,ypos,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,limin,limax,ljmin,ljmax)
        CALL h_int_siggen(xe_lag,h0_lag,hc_sig_lag,xe,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,limin,limax,ljmin,ljmax)
        CALL f_lag_sigz_uv(spos,kuvm,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)

        uint = (f_lag_hzp*(sc_r(kuv)-spos)*uintm + f_lag_hzm*(spos - sc_r(kuvm))*uintp)/(f_lag_hz*dsigw(kuvm)) 
    END IF

  END FUNCTION uint



  !!======================================================================
  FUNCTION vint(xpos,ypos,spos,v,xe,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION vint  ***
    !&E
    !&E ** Purpose : Linear interpolation of meridional velocity (vint)
    !&E              at the location (xpos,ypos,spos)
    !&E
    !&E ** Called by      : avance
    !&E
    !&E ** External calls : f_lag_sigz_uv,h_int_siggen,interpvit,loc_h0
    !&E
    !&E ** History :
    !&E       !  08-2002 (F. Dumas, P. Lazure)
    !&E       !  01-2007 (V. Garnier) correction bug declaration px,py
    !&E       !  01-2011 (M. Huret) Adaptation to siggen (scale factor for vertical interpolation)
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                                  INTENT( in ) :: xpos,ypos,spos
    INTEGER,                                         INTENT( in ) :: limin,limax,ljmin,ljmax
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in ) :: v
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),      INTENT( in ) :: xe

    REAL(KIND=rsh)                                                :: vint   ! Result

    !! * Local declarations
    INTEGER                                                       :: idd,igg,jhh,jbb,k,kuv,kuvm
    INTEGER                                                       :: hlb,hlt,hrb,hrt
                                      
    REAL(KIND=rsh)                                                :: vintm,vintp,px,py
    REAL(KIND=rsh)                                                :: f_lag_hzp,f_lag_hzm,f_lag_hz
    REAL(KIND=rsh)                                                :: xe_lag,hc_sig_lag,h0_lag

    !!----------------------------------------------------------------------
    !! * Executable part

    !---- reperage de la maille ou se situe la particule pour les U
    ! idd, igg relatifs a la grille U
    ! jbb ,jhh relatifs a la grille V
    ! xpos,ypos relatifs a la grille rho
    ! location of ssh : xpos=REAL(i), ypos=REAL(j)
    igg = int(xpos)
    idd = igg+1
    jbb = nint(ypos)
    jhh = jbb+1

    idd = min(max(idd,0),imax)
    igg = min(max(igg,0),imax)
    jbb = min(max(jbb,0),jmax)
    jhh = min(max(jhh,0),jmax)

    ! calcul des positions px et py dans la maille V
    px = xpos - real(igg, kind=rsh)
    py = ypos - real(jbb, kind=rsh) - 0.5_rsh

    DO k=kmax,2,-1
        IF (spos <= sc_r(k).and.spos >= sc_r(k-1)) THEN
            kuv  = k
            kuvm = k-1
        END IF
    END DO
    IF (spos <= sc_r(1)) THEN
        kuv  = 1
        kuvm = 0
    END IF
    IF (spos >= sc_r(kmax)) THEN
        kuv  = kmax
        kuvm = kmax
    END IF

    !-- interpolation proprement dite sur la nappe sigma au-dessus
    vintp = interpvit(px,py,v(igg,jbb,kuv),v(idd,jbb,kuv), &
                            v(igg,jhh,kuv),v(idd,jhh,kuv), &
                            hty(igg,jbb),hty(idd,jbb),     &
                            hty(igg,jhh),hty(idd,jhh))

    !-- interpolation proprement dite sur la nappe sigma au-dessous
    IF (kuvm /= 0) THEN
        vintm = interpvit(px,py,v(igg,jbb,kuvm),v(idd,jbb,kuvm), &
                                v(igg,jhh,kuvm),v(idd,jhh,kuvm), &
                                hty(igg,jbb),hty(idd,jbb),       &
                                hty(igg,jhh),hty(idd,jhh))
    ELSE
        vintm = 0.0_rsh
    END IF

    !-- interpolation sur la verticale entre la valeur au-dessus (vintp)
    !-- et la valeur au-dessous (vintm).
    IF (spos > sc_r(kmax)) THEN
        vint = vintm
    ELSE IF(kuv == 1) THEN
        vint = ((sc_r(kuv) - spos)*vintm + (spos + 1.0_rsh)*vintp)/(sc_r(kuv) + 1.0_rsh)
    ELSE
        f_lag_hzp = 1.0_rsh
        f_lag_hzm = 1.0_rsh
        f_lag_hz  = 1.0_rsh
        CALL loc_h0(xpos,ypos,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,limin,limax,ljmin,jmax)
        CALL h_int_siggen(xe_lag,h0_lag,hc_sig_lag,xe,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,limin,limax,ljmin,ljmax)
        CALL f_lag_sigz_uv(spos,kuvm,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)

        vint = (f_lag_hzp*(sc_r(kuv) - spos)*vintm + f_lag_hzm*(spos - sc_r(kuvm))*vintp)/(f_lag_hz*dsigw(kuvm))
    END IF
 
  END FUNCTION vint



  !!======================================================================
  FUNCTION interpvit(x,y,v00,v10,v01,v11,rh00,rh10,rh01,rh11)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION interpvit  ***
    !&E
    !&E ** Purpose : interpolation horizontale d une composante de vitesse
    !&E              tenant compte des profondeurs aux points encadrant la
    !&E              position ou l on doit interpoler.
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : uint,vint
    !&E
    !&E ** External calls : 
    !&E
    !&E ** History :
    !&E       !  08-2002 (F. Dumas, P. Lazure)
    !&E       !  2024    (M. Caillaud) Coupled with CROCO  
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh), INTENT( in )      :: x,y,v00,v10,v01,v11,rh00,rh10,rh01,rh11
    REAL(KIND=rsh)                    :: interpvit

    !! * Local declarations
    INTEGER                           :: icas
    REAL(KIND=rsh)                    :: xref,yref

    !!----------------------------------------------------------------------
    !! * Executable part

    IF ( rh00 > 0.0_rsh .and. rh10 > 0.0_rsh .and. rh01 > 0.0_rsh .and. rh11 > 0.0_rsh ) THEN
        interpvit = x*y*v11 + x*(1.0_rsh - y)*v10 + (1.0_rsh - x)*y*v01 + (1.0_rsh - x)*(1.0_rsh - y)*v00
        icas      = 1
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = 0.0_rsh
        xref      = 0.5_rsh
        yref      = 0.5_rsh
        icas      = 0
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = v11
        xref      = 1.0_rsh
        yref      = 1.0_rsh
        icas      = 2
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = v10
        xref      = 1.0_rsh
        yref      = 0.0_rsh
        icas      = 4
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = (1.0_rsh - y)*v10 + y*v11
        xref      = 1.0_rsh
        yref      = 0.5_rsh
        icas      = 3
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = v01
        xref      = 0.0_rsh
        yref      = 1.0_rsh
        icas      = 5
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = ((v10 - v01)/2.0_rsh)*(x - y) + (v10 + v01)/2.0_rsh
        xref      = 0.5_rsh
        yref      = 0.5_rsh
        icas      = 6
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit =(1.0_rsh - x)*v01 + x*v11
        xref      = 0.5_rsh
        yref      = 1.0_rsh
        icas      = 7
    ELSE IF (rh00 <= 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = (v11 - v01)*x + (v11 - v10)*y + v10 - (v11 - v01)
        xref     = 1.0_rsh
        yref     = 1.0_rsh
        icas     = 8
    ELSE IF (rh00 > 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = v00
        xref      = 0.0_rsh
        yref      = 0.0_rsh
        icas      = 9
    ELSE IF (rh00 > 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = ((v11 - v00)/2.0_rsh)*(x + y) + v00
        xref      = 0.5_rsh
        yref      = 0.5_rsh
        icas      = 10
    ELSE IF (rh00 > 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = (1.0_rsh - y)*v00 + y*v01
        xref      = 0.0_rsh
        yref      = 0.5_rsh
        icas      = 11
    ELSE IF (rh00 > 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 <= 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = (v11 - v01)*x + (v01 - v00)*y + v00
        xref      = 0.0_rsh
        yref      = 1.0_rsh
        icas      = 12
    ELSE IF (rh00 > 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = (1.0_rsh - x)*v00 + x*v10
        xref      = 0.5_rsh
        yref      = 0.0_rsh
        icas      = 13
    ELSE IF (rh00 > 0.0_rsh .and. rh01 > 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 <= 0.0_rsh) THEN
        interpvit = (v10 - v00)*x + (v01 - v00)*y + v00
        xref      = 0.0_rsh
        yref      = 0.0_rsh
        icas      = 14
    ELSE IF (rh00 > 0.0_rsh .and. rh01 <= 0.0_rsh .and. rh10 > 0.0_rsh .and. rh11 > 0.0_rsh) THEN
        interpvit = (v10 - v00)*x + (v11 - v10)*y + v00
        xref      = 1.0_rsh
        yref      = 0.0_rsh
        icas      = 15
    END IF
  
  END FUNCTION interpvit



  !!======================================================================
  SUBROUTINE ksupkinf(pos,refpos,sizerefpose,kw,kwm,appel)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE ksupkinf ***
    !&E
    !&E ** Purpose : Seek neighboor index in the vertical for a particle position
    !&E
    !&E ** Called by   : LAGRANGIAN_update, ibm_profmean, ibm_profmean,
    !&E                  wint, dksdzint, splint, ztosiggen
    !&E
    !&E ** History :
    !&E       !  ( ? )
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    
    REAL(KIND=rsh),                            INTENT(in)    :: pos
    INTEGER,                                   INTENT(in)    :: sizerefpose,appel
    REAL(KIND=rsh),DIMENSION(0:sizerefpose-1), INTENT(in)    :: refpos 
    INTEGER,                                   INTENT(inout) :: kw,kwm

    !! * Local declarations
    INTEGER                                                  :: k
    INTEGER                                                  :: IERR_MPI

    !!----------------------------------------------------------------------
    !! * Executable part

    DO WHILE (kw-kwm > 1) 
        k=int((kw+kwm)/2)  
        IF (refpos(k) > pos) THEN
            kw=k
        ELSE
            kwm=k
        ENDIF
    ENDDO

    IF ((refpos(kw) < pos) .OR. (refpos(kwm) > (pos+0.00001_rsh))) THEN
        print*,'erreur subroutine ksupkinf ',appel
        print*,'kw',kw
        print*,'kwm',kwm
        print*,'pos',pos
        print*,'refpos(1)',refpos(1)
        print*,'refpos(kw)',refpos(kw)
        print*,'refpos(kwm)',refpos(kwm)
        print*,'refpos',refpos
        print*,'refpos',appel
#ifdef MPI
        CALL MPI_FINALIZE(ierr_mpi)
#endif
        STOP
    END IF

 END SUBROUTINE ksupkinf



  !!======================================================================
  SUBROUTINE loc_h0(xpos,ypos,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  SUBROUTINE location  ***
    !&E
    !&E ** Purpose : search index of the cells around particle
    !&E              and location of the particle inside the cell 
    !&E              for later interpolation 
    !&E              Valid relative to h0,xe,wz,kz... 
    !&E
    !&E ** Called by : LAGRANGIAN_update, LAGRANGIAN_init, ibm_3d,
    !&E                uint, vint
    !&E
    !&E ** History :
    !&E       !  2012-07 (M. Huret)
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh), INTENT( in )        :: xpos,ypos
    INTEGER,        INTENT( in )        :: limin,limax,ljmin,ljmax
    INTEGER,        INTENT( out )       :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    REAL(KIND=rsh), INTENT( out )       :: px,py 

    !!----------------------------------------------------------------------
    !! * Executable part

    ! indexes of the cell where the particle is
    ! reperage de la maille XE ou se situe la particule
    ! location of ssh : xpos=REAL(i), ypos=REAL(j)
    i1 = INT(xpos)
    i2 = INT(xpos)+1
    j1 = INT(ypos)
    j2 = INT(ypos)+1

    i1 = MIN(MAX(i1,limin-1),limax+1)
    i2 = MIN(MAX(i2,limin-1),limax+1)
    j1 = MIN(MAX(j1,ljmin-1),ljmax+1)
    j2 = MIN(MAX(j2,ljmin-1),ljmax+1)

    ! location of the particle inside the cell
    ! calcul des positions px et py dans la maille XE
    px = xpos - INT(xpos)
    py = ypos - INT(ypos)

    ! check for next cells in water
    rh11 = 1
    rh12 = 1
    rh21 = 1
    rh22 = 1

    !IF (h(i1,j1)==-valmanq) THEN
    !    rh11=0
    !END IF
    !IF (h(i1,j2)==-valmanq) THEN
    !    rh12=0
    !END IF
    !IF (h(i2,j1)==-valmanq) THEN
    !    rh21=0
    !END IF
    !IF (h(i2,j2)==-valmanq) THEN
    !    rh22=0
    !END IF
#ifdef MASKING
    IF (rmask(i1,j1)==0)  rh11 = 0
    IF (rmask(i1,j2)==0)  rh12 = 0
    IF (rmask(i2,j1)==0)  rh21 = 0
    IF (rmask(i2,j2)==0)  rh22 = 0
#endif

  END SUBROUTINE loc_h0



  !!======================================================================
  FUNCTION h0int(px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION h0int  ***
    !&E
    !&E ** Purpose : calcul par interpolation lineaire de la bathymetry
    !&E
    !&E ** Called by : LAGRANGIAN_update, LAGRANGIAN_init, ibm_3d
    !&E
    !&E ** History :
    !&E       !          (P. Lazure, F. Dumas)  Original code
    !&E       !  08-2008 (F. Dumas ?)
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                                 INTENT( in )  :: px,py
    INTEGER,                                        INTENT( in )  :: i1,i2,j1,j2,rh11,rh12,rh21,rh22

    REAL(KIND=rsh)                                                :: h0int
    
    !!----------------------------------------------------------------------
    !! * Executable part

    ! spatial interpolation
    h0int = px*(1.0_rsh-py)*h(i2,j1)*real(rh21)+ px*py*h(i2,j2)*real(rh22) &
            +(1.0_rsh-px)*py*h(i1,j2)*real(rh12)+(1.0_rsh-px)*(1.0_rsh-py)*h(i1,j1)*real(rh11)

  END FUNCTION h0int



  !!======================================================================
  FUNCTION xeint(xe,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION xe_int  ***
    !&E
    !&E ** Purpose : calcul par interpolation lineaire de l elevation hauteur d eau
    !&E
    !&E ** Called by : LAGRANGIAN_update, LAGRANGIAN_init, ibm_3d
    !&E
    !&E ** History :
    !&E       !          (P. Lazure, F. Dumas)  Original code
    !&E       !  08-2008 (F. Dumas ?)
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                                  INTENT( in ) :: px,py
    INTEGER,                                         INTENT( in ) :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),      INTENT( in ) :: xe
    INTEGER,                                         INTENT( in ) :: limin,limax,ljmin,ljmax

    REAL(KIND=rsh)                                                :: xeint

    !!----------------------------------------------------------------------
    !! * Executable part

    ! spatial interpolation
    xeint = px*(1.0_rsh-py)*xe(i2,j1)*REAL(rh21,rsh)+ px*py*xe(i2,j2)*REAL(rh22,rsh) &
            +(1.0_rsh-px)*py*xe(i1,j2)*REAL(rh12,rsh)+(1.0_rsh-px)*(1.0_rsh-py)*xe(i1,j1)*REAL(rh11,rsh)

  END FUNCTION xeint



  !!======================================================================
  FUNCTION interp(varint,k,i1,i2,j1,j2,px,py,rh11,rh21,rh12,rh22,&
                  limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION interp  ***
    !&E
    !&E ** Purpose : Interpolate at H0 or tracer location
    !&E              with check of terrestrial boundary 
    !&E
    !&E ** Called by      : wint, kzprofile
    !&E
    !&E ** External calls : varint
    !&E
    !&E ** History :
    !&E       !  01-2008 (M. Huret)
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments

    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,0:kmax),INTENT( in )  :: varint  
    REAL(KIND=rsh),                                   INTENT( in )  :: px,py
    INTEGER,                                          INTENT( in )  :: k
    INTEGER,                                          INTENT( in )  :: i1,j1,i2,j2
    INTEGER,                                          INTENT( in )  :: rh11,rh21,rh12,rh22
    INTEGER,                                          INTENT( in )  :: limin,limax,ljmin,ljmax

    REAL(KIND=rsh)                                                  :: interp

    !!----------------------------------------------------------------------
    !! * Executable part

    ! check if boundary condition
    IF (rh11 > 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        interp = px*(1.0_rsh - py)*varint(i2,j1,k) + px*py*varint(i2,j2,k) + &
                 py*(1.0_rsh - px)*varint(i1,j2,k) + (1.0_rsh - px)*(1.0_rsh - py)*varint(i1,j1,k)

    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        interp = 0

    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        interp = varint(i2,j2,k)

    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        interp = varint(i1,j2,k)

    ELSE IF (rh11 <= 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        interp = px*varint(i2,j2,k) + (1.0_rsh-px)*varint(i1,j2,k)

    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        interp = varint(i2,j1,k)

    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        IF (px > 0.5_rsh) THEN
            interp = varint(i2,j1,k) 
        ELSE
            interp = varint(i1,j2,k)
        END IF

    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
      interp = (1.0_rsh - py)*varint(i2,j1,k) + py*varint(i2,j2,k)

    ELSE IF (rh11 <= 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        interp = px*(varint(i2,j2,k) - varint(i1,j2,k)) + py*(varint(i2,j2,k) - varint(i2,j1,k)) & 
                 - varint(i2,j2,k) + varint(i1,j2,k) + varint(i2,j1,k)
    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
        interp = varint(i1,j1,k)

    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        IF (px > 0.5_rsh) THEN
            interp = varint(i2,j2,k)
        ELSE
            interp = varint(i1,j1,k)
        END IF

    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 <= 0) THEN
      interp = px*varint(i2,j1,k) + (1.0_rsh - px)*varint(i1,j1,k)

    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 <= 0 .and. rh22 > 0) THEN
        interp = px*(varint(i2,j1,k) - varint(i1,j1,k)) + py*(varint(i2,j2,k) - varint(i2,j1,k)) + varint(i1,j1,k)

    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        interp = py*varint(i1,j2,k) + (1.0_rsh-py)*varint(i1,j1,k)

    ELSE IF (rh11 > 0 .and. rh21 > 0 .and. rh12 > 0 .and. rh22 <= 0) THEN
        interp = px*(varint(i2,j1,k) - varint(i1,j1,k)) + py*(varint(i1,j2,k) - varint(i1,j1,k)) + varint(i1,j1,k)

    ELSE IF (rh11 > 0 .and. rh21 <= 0 .and. rh12 > 0 .and. rh22 > 0) THEN
        interp = px*(varint(i2,j2,k) - varint(i1,j2,k)) + py*(varint(i1,j2,k) - varint(i1,j1,k)) + varint(i1,j1,k)
    END IF

  END FUNCTION interp



  !!======================================================================
  FUNCTION wint(spos_s,xelag,h0lag,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION wint  ***
    !&E
    !&E ** Purpose : Linear interpolation of vertical velocity (wint)
    !&E              at the location (spos_s) between indices i1,i2,j1,j2
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by      : LAGRANGIAN_update
    !&E
    !&E ** External calls : ksupkinf,interp,f_lag_sigz_w, hc_sigint
    !&E
    !&E ** History :
    !&E       !  08-2002 (F. Dumas, P. Lazure) 
    !&E       !  01-2011 (M. Huret) Adaptation to siggen 
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),         INTENT( in ) :: spos_s,xelag,h0lag,px,py
    INTEGER,                INTENT( in ) :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    INTEGER,                INTENT( in ) :: limin,limax,ljmin,ljmax

    REAL(KIND=rsh)                       :: wint

    !! * Local declarations
    INTEGER                              :: IERR_MPI
    INTEGER                              :: kw,kwm
    REAL(KIND=rsh)                       :: wintm,wintp,f_lag_hzp,f_lag_hzm,f_lag_hz,hc_siglag
    
    !!----------------------------------------------------------------------
    !! * Executable part

    kwm = 0 
    kw  = kmax
    CALL ksupkinf(spos_s,sc_w,kmax+1,kw,kwm,5)
    
    IF ((kw > kmax) .or. (kwm > kmax) .or. (kw < 0) .or. (kwm < 0)) THEN
        WRITE(*,*) 'mauvais sigma fonction wint'
#ifdef MPI
        CALL MPI_FINALIZE(ierr_mpi)
#endif
        STOP
    END IF
    
    !-- interpolation proprement dite sur la nappe sigma au-dessus
    wintp = interp(wz,kw,i1,i2,j1,j2,px,py,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)

    !-- interpolation proprement dite sur la nappe sigma au-dessous
    IF (kwm /= 0) THEN
        wintm = interp(wz,kwm,i1,i2,j1,j2,px,py,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)
    ELSE
        wintm = 0.0_rsh
    END IF

    !-- interpolation sur la verticale entre la valeur au-dessus (wintp)
    !-- et la valeur au-dessous (wintm).
    f_lag_hzp = 1.0_rsh
    f_lag_hzm = 1.0_rsh
    f_lag_hz  = 1.0_rsh
    hc_siglag = hc_sigint(px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22)
    CALL f_lag_sigz_w(spos_s,kw,xelag,h0lag,hc_siglag,f_lag_hzp,f_lag_hzm,f_lag_hz)
    wint=(f_lag_hzp*(sc_w(kw) - spos_s)*wintm + f_lag_hzm*(spos_s - sc_w(kwm))*wintp) &
         /(f_lag_hz*dsigu(kw))

  END FUNCTION wint



  !!======================================================================
  FUNCTION dksdzint(zpos_s,zpos,kzsmth,dkz2)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE dksdzint  ***
    !&E
    !&E ** Purpose : calcul de la derive de la composante verticale
    !&E              du melange turbulent (notee dksdzint) pour le profil
    !&E              kz (kzsmth) lisse calcule dans kzprofile dont derive
    !&E              seconde est dkz2
    !&E              adapte de splint
    !&E
    !&E ** Called by : LAGRANGIAN_update
    !&E
    !&E ** History :
    !&E       !  2007-12  (M. Huret)
    !&E       !  2010-12  (M. Huret) Switch to z
    !&E       !  2024     (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                     INTENT( in )  :: zpos_s
    REAL(KIND=rsh), DIMENSION(0:kmax),  INTENT( in )  :: kzsmth,dkz2,zpos

    REAL(KIND=rsh)                                    :: dksdzint

    !! * Local declarations
    REAL(KIND=rsh)                                    :: a,b,dzw_local
    INTEGER                                           :: k,khi,klo
    INTEGER                                           :: IERR_MPI

    !!----------------------------------------------------------------------
    !! * Executable part
    klo = 0
    khi = kmax
    CALL ksupkinf(zpos_s,zpos,kmax+1,khi,klo,3)

    IF ((khi > kmax) .OR. (klo > kmax) .OR. (khi < 0) .OR. (klo < 0)) THEN
        WRITE(*,*) 'mauvais sigma fonction dksdzint'
#ifdef MPI
        CALL MPI_FINALIZE(ierr_mpi)
#endif
        STOP
    END IF
    
    dzw_local = zpos(khi) - zpos(klo)
    IF (dzw_local <= 0.0_rsh) THEN
        WRITE(*,*) khi,klo,zpos_s,'bad sigw input in splint'
        READ(*,*) 
    END IF
    
    a = (zpos(khi) - zpos_s)/dzw_local
    b = (zpos_s - zpos(klo))/dzw_local

    dksdzint = (-kzsmth(klo) + kzsmth(khi))/dzw_local + ((-3.0_rsh*a**2 + 1)*dkz2(klo) + &
               (3.0_rsh*b**2 - 1)*dkz2(khi))*dzw_local/6.0_rsh
     
  END FUNCTION dksdzint



  !!======================================================================
  SUBROUTINE kzprofile(kzsmth,dkz2,zpos,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,limin,limax,ljmin,ljmax)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE kzprofile  ***
    !&E
    !&E ** Purpose : Calcul les derivees secondes du coeff de diffusivite 
    !&E              vertical en sigw a l aide de spline apres avoir filtre
    !&E              le profil kz
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : LAGRANGIAN_update
    !&E
    !&E ** External calls : interp, spline
    !&E
    !&E ** History :
    !&E       !  2012-07 (M. Huret)
    !&E       !  2024    (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(0:kmax),       INTENT( out ) :: dkz2,kzsmth
    REAL(KIND=rsh), DIMENSION(0:kmax),       INTENT( in )  :: zpos
    INTEGER,                                 INTENT( in )  :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    REAL(KIND=rsh),                          INTENT( in )  :: px,py
    INTEGER,                                 INTENT( in )  :: limin,limax,ljmin,ljmax

    !! * Local declarations
    REAL(KIND=rsh), DIMENSION(0:kmax)    :: kzint
    REAL(kind=rsh)                       :: big

    INTEGER                              :: k
    INTEGER                              :: IERR_MPI

    !!----------------------------------------------------------------------
    !! * Executable part
    big = 1.0E30

    ! horizontal interpolation at each sigma level
    kzint(0)    = 0.0_rsh
    kzint(kmax) = 0.0_rsh

   DO k=1,kmax-1
        kzint(k) = interp(Akt(:,:,0:kmax,nstp),k,i1,i2,j1,j2,px,py,rh11,rh21,rh12,rh22, &
                          limin,limax,ljmin,ljmax)
   END DO

   ! smooth of the kz profile with weighting by the distance
   IF (kmax <= 3) THEN
        PRINT*,'float interp code cant handle'
        PRINT*,'kmax.le.3'
#ifdef MPI
        CALL MPI_FINALIZE(ierr_mpi)
#endif
        STOP
   ENDIF

   kzsmth(0)    = kzint(0) 
   kzsmth(kmax) = kzint(kmax)

   DO k = 1,kmax-1
        kzsmth(k) = (4.0_rsh/6.0_rsh)*kzint(k) + (2.0_rsh/6.0_rsh)*(kzint(k+1)*(zpos(k) - zpos(k-1)) + &
                     kzint(k-1)*(zpos(k+1) - zpos(k)))/(zpos(k+1) - zpos(k-1))
   ENDDO

   CALL spline(zpos(0:kmax),kzsmth,kmax+1,big,big,dkz2) 

   END SUBROUTINE kzprofile

   
  !!======================================================================
  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE spline  ***
    !&E
    !&E     From numerical recipies vol 2, but modified so that nmax=50
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    real(kind=rsh),                 INTENT( in )  :: yp1,ypn
    real(kind=rsh),dimension(n),    INTENT( in )  :: x,y
    integer,                        INTENT( in )  :: n
    real(kind=rsh),dimension(n),    INTENT( out ) :: y2

    !! * Local declarations
    integer                                       :: i,k
    real(kind=rsh),dimension(n)                   :: u
    real(kind=rsh)                                :: p,qn,sig,un

    ! executable statement
    !---------------------------------------------------------
    IF (yp1 > .99e30) THEN
        y2(1)=0.0
        u(1)=0.0
    ELSE
        y2(1)=-0.5
        u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF

    DO i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0
        y2(i)=(sig-1.0)/p
        u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    END DO

    IF (ypn > 0.99e30) THEN
        qn=0.0
        un=0.0
    ELSE
        qn=0.5
        un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)

    DO  k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    END DO

  END SUBROUTINE spline



  !!======================================================================
  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE splint  ***
    !&E
    !&E     From numerical recipies vol 2
    !&E---------------------------------------------------------------------

    !! * Arguments
    integer,                         INTENT( IN )  :: n
    real(kind=rsh),                  INTENT( IN )  :: x
    real(kind=rsh),                  INTENT( OUT ) :: y
    real(kind=rsh),dimension(0:n-1), INTENT( IN )  :: xa,y2a,ya
    
    !! * Local declarations
    integer                     :: k,khi,klo
    real(kind=rsh)              :: a,b,h
    INTEGER                     :: IERR_MPI
    
    !!----------------------------------------------------------------------
    !! * Executable statment
    
    klo = 0
    khi = n-1
    CALL ksupkinf(x,xa,n,khi,klo,4)

    IF ((khi > kmax) .or. (klo > kmax) .or. (khi < 0) .or. (klo < 0)) THEN
        WRITE(*,*) 'mauvais sigma fonction splint'
        WRITE(*,*) khi,klo,kmax
        WRITE(*,*) x
        WRITE(*,*) xa
        CALL_MPI MPI_FINALIZE(ierr_mpi)
        STOP
    END IF
    
    h=xa(khi)-xa(klo)
    IF (h==0.0) THEN
        WRITE(*,*) 'bad xa input in splint'
        READ(*,*)
    END IF
    a = (xa(khi) - x)/h
    b = (x - xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3  -b)*y2a(khi))*(h**2)/6.0
    
  END SUBROUTINE splint



  !!======================================================================
  SUBROUTINE f_lag_sigz_uv(spos,k,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE f_lag  ***
    !&E
    !&E ** Purpose : Scale factor from generalized sigma to z (relative to u and v location
    !&E              for vertical interpolation at particle location
    !&E
    !&E ** Description :
    !&E
    !&E ** Called by : None
    !&E
    !&E ** History : 
    !&E       !  2010-12  (M. Huret)
    !&E       !  2024     (M. Caillaud) Coupled with CROCO 
    !&E---------------------------------------------------------------------
    !! * Modules used
    !USE comsiggen, ONLY : b_sig,theta_sig,Csu_sig,dcwsds

    !! * Arguments
    REAL(KIND=rsh), INTENT( in ) :: spos,xe_lag,h0_lag,hc_sig_lag
    REAL(KIND=rsh), INTENT(out ) :: f_lag_hzp,f_lag_hzm,f_lag_hz
    INTEGER,        INTENT(IN)   :: k

    !! * Local declarations
    REAL(KIND=rsh)      :: Csusig, dcwsdsp,dcwsdsm,hinv
    
    !!----------------------------------------------------------------------
    !! * Executable part

    Csusig=CSF(spos,theta_s,theta_b)
    hinv=1./(h0_lag+hc_sig_lag)
    dcwsdsp=0.0_rsh
    dcwsdsm=0.0_rsh

    IF (spos /= sc_r(k+1)) THEN
        dcwsdsp = (Cs_r(k+1)-Csusig)/(sc_r(k+1)-spos)
    ENDIF
    IF (spos /= sc_r(k)) THEN
        dcwsdsm = (Csusig-Cs_r(k))/(spos-sc_r(k))
    ENDIF

    f_lag_hzp = hinv*(hc_sig_lag + h0_lag*dcwsdsp)*(h0_lag + xe_lag)
    f_lag_hzm = hinv*(hc_sig_lag + h0_lag*dcwsdsm)*(h0_lag + xe_lag)
    f_lag_hz  = hinv*(hc_sig_lag + h0_lag*dcwsds(k))*(h0_lag + xe_lag)

  END SUBROUTINE f_lag_sigz_uv


  !!======================================================================
  SUBROUTINE h_int_siggen(xeint,h0int,hc_sigint,xe,px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22,&
                          limin,limax,ljmin,ljmax)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE h_int  ***
    !&E
    !&E ** Purpose : calcul par interpolation lineaire de la hauteur d eau
    !&E              a la position (xpos_s,ypos_s,spos_s)
    !&E              (also return xe and h0 interpolated)
    !&E
    !&E ** Called by : uint, vint, ibm_traint
    !&E
    !&E ** History :
    !&E       !          (P. Lazure, F. Dumas)  Original code
    !&E       !  08-2008 (F. Dumas ?)
    !&E       !  12/2010 (M. huret)
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(liminm2:limaxp2,ljminm2:ljmaxp2)
    !&E       !  2024     (M. Caillaud) Coupled with CROCO 
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                                 INTENT( in )  :: px,py
    INTEGER,                                        INTENT( in )  :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    INTEGER,                                        INTENT( in )  :: limin,limax,ljmin,ljmax
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),     INTENT( in )  :: xe
    REAL(KIND=rsh),                                 INTENT( out ) :: xeint,h0int,hc_sigint

    !! * Local declarations

    !!----------------------------------------------------------------------
    !! * Executable part

    ! spatial interpolation
    xeint = px*(1.0-py)*xe(i2,j1)*rh21 + px*py*xe(i2,j2)*rh22 +                   &
            (1.0-px)*py*xe(i1,j2)*rh12 + (1.0-px)*(1.0-py)*xe(i1,j1)*rh11
    h0int = px*(1.0-py)*h(i2,j1)*rh21 + px*py*h(i2,j2)*rh22 +                     &
            (1.0-px)*py*h(i1,j2)*rh12 + (1.0-px)*(1.0-py)*h(i1,j1)*rh11
    hc_sigint = px*(1.0-py)*hc_sig(i2,j1)*rh21 + px*py*hc_sig(i2,j2)*rh22 +       &
            (1.0-px)*py*hc_sig(i1,j2)*rh12 + (1.0-px)*(1.0-py)*hc_sig(i1,j1)*rh11

  END SUBROUTINE h_int_siggen



  !!======================================================================
  FUNCTION hc_sigint(px,py,i1,i2,j1,j2,rh11,rh21,rh12,rh22)

    !&E                 ***  ROUTINE hc_sigint  ***
    !&E---------------------------------------------------------------------
    !&E
    !&E ** Purpose : Linear interpolation of hc_sig at the particle location 
    !&E
    !&E ** Called by : LAGRANGIAN_update, LAGRANGIAN_init, ibm_3d, wint
    !&E
    !&E
    !&E ** History :
    !&E       !  2010-12  (M. Huret)
    !&E       !  2024     (M. Caillaud) Coupled with CROCO 
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),                   INTENT( in )  :: px,py
    INTEGER,                          INTENT( in )  :: i1,i2,j1,j2,rh11,rh12,rh21,rh22
    REAL(KIND=rsh)                                  :: hc_sigint

    !! * Local declarations

    !!----------------------------------------------------------------------
    !! * Executable part

    ! spatial interpolation
    hc_sigint = px*(1.0-py)*hc_sig(i2,j1)*rh21 + px*py*hc_sig(i2,j2)*rh22  + &
                (1.0-px)*py*hc_sig(i1,j2)*rh12 + (1.0-px)*(1.0-py)*hc_sig(i1,j1)*rh11

   END FUNCTION hc_sigint



  !!======================================================================
  SUBROUTINE f_lag_sigz_w(spos,k,xe_lag,h0_lag,hc_sig_lag,f_lag_hzp,f_lag_hzm,f_lag_hz)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE f_lag  ***
    !&E
    !&E ** Purpose : Scale factor from generalized sigma to z
    !&E              for vertical interpolation at particle location
    !&E
    !&E ** Called by : None
    !&E
    !&E ** External calls : CSF
    !&E
    !&E ** History :
    !&E       !  2010-12  (M. Huret)
    !&E       !  2024     (M. Caillaud) Coupled with CROCO 
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    !USE comsiggen, ONLY : b_sig,theta_sig,dcusds,Csw_sig
    !USE comvars3d,   ONLY : sigw

    !! * Arguments
    REAL(KIND=rsh),              INTENT( in ) :: spos,xe_lag,h0_lag,hc_sig_lag
    REAL(KIND=rsh),              INTENT(out ) :: f_lag_hzp,f_lag_hzm,f_lag_hz
    INTEGER,                     INTENT( in)  :: k

    !! * Local declarations
    REAL(KIND=rsh)      :: Cswsig, dcusdsp,dcusdsm,hinv
    
    !!----------------------------------------------------------------------
    !! * Executable part

    Cswsig  = CSF(spos,theta_s,theta_b)
    hinv    = 1/(h0_lag+hc_sig_lag)
    dcusdsp = 0.0_rsh
    dcusdsm = 0.0_rsh

    IF (spos /= sc_w(k)) THEN
        dcusdsp = (Cs_w(k) - Cswsig)/(sc_w(k) - spos)
    ENDIF
    IF (spos /= sc_w(k-1)) THEN
        dcusdsm = (Cswsig - Cs_w(k-1))/(spos - sc_w(k-1))
    ENDIF

    f_lag_hzp = hinv*(hc_sig_lag + h0_lag*dcusdsp)*(h0_lag + xe_lag)
    f_lag_hzm = hinv*(hc_sig_lag + h0_lag*dcusdsm)*(h0_lag + xe_lag)
    f_lag_hz  = hinv*(hc_sig_lag + h0_lag*dcusds(k))*(h0_lag + xe_lag)

  END SUBROUTINE f_lag_sigz_w



  !!======================================================================
  SUBROUTINE siggentoz(zpos,spos,xe_lag,h0_lag,hc_sig_lag)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE siggen_int  ***
    !&E
    !&E ** Purpose : Calculate z from sigma (for siggen case)
    !&E
    !&E ** Called by : LAGRANGIAN_update, ztosiggen, ibm_profmean
    !&E
    !&E ** History :
    !&E       !  2024     (M. Caillaud) Changed for CROCO way to calculate sigma
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),             INTENT( in )  :: spos,xe_lag,hc_sig_lag,h0_lag
    REAL(KIND=rsh),             INTENT( out ) :: zpos

    !! * Local declarations
    REAL(KIND=rsh)                            :: Cssig,hinv,z_r0

    !!----------------------------------------------------------------------
    !! * Executable part

    Cssig = CSF(spos,theta_s,theta_b) 
    hinv  = 1./(h0_lag + hc_sig_lag)     
    z_r0  = hc_sig_lag*spos + Cssig*h0_lag
    zpos  = z_r0*h0_lag*hinv + xe_lag*(1. + z_r0*hinv)

  END SUBROUTINE siggentoz



  !!======================================================================
  SUBROUTINE ztosiggen(zpos_s,spos,xelag,h0lag,hc_siglag)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE ztosiggen  ***
    !&E
    !&E ** Purpose : Calcul du sigma a partir du z (en sigma generalises)
    !&E
    !&E ** Called by : LAGRANGIAN_update, LAGRANGIAN_init, ibm_3d, ibm_parameter_init
    !&E
    !&E ** External calls : siggentoz, ksupkinf
    !&E
    !&E ** History : 
    !&E       ! 2010-12   (M. Huret) (to be improved, only linear interpolation
    !&E                   should use a correction factor considering Ds/Dz
    !&E       !  2024     (M. Caillaud) Changed for CROCO way to calculate sigma
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(KIND=rsh),               INTENT( in )  :: zpos_s,xelag,h0lag,hc_siglag
    REAL(KIND=rsh),               INTENT( out ) :: spos

    !! * Local declarations
    INTEGER                                     :: kw,kwm,l
    REAL(KIND=rsh), DIMENSION(0:kmax)           :: zpos
    
    !!----------------------------------------------------------------------
    !! * Executable part

    DO l=0,kmax
        CALL siggentoz(zpos(l),sc_w(l),xelag,h0lag,hc_siglag)
    END DO

    kw=kmax
    kwm=0
    CALL ksupkinf(zpos_s,zpos,kmax+1,kw,kwm,1)

    spos = ((zpos_s - zpos(kwm))*sc_w(kw) + (zpos(kw) - zpos_s)*sc_w(kwm))/(zpos(kw) - zpos(kwm))  !!Approximation

  END SUBROUTINE ztosiggen



  !!======================================================================
  FUNCTION CSF(sc,theta_s,theta_b)  
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION CSF   ***
    !&E
    !&E ** Purpose : compute of Cs sigma parameter given sigma position and
    !&E              stretching parameters 
    !&E
    !&E ** Called by : siggentoz, f_lag_sigz_w, f_lag_sigz_uv
    !&E
    !&E ** History : 
    !&E       !  (??)
    !&E       !  2024     (M. Caillaud) Changed for CROCO way to calculate sigma
    !&E              
    !&E---------------------------------------------------------------------
    !! * Arguments
    IMPLICIT NONE  
    REAL(kind=rsh),INTENT(in)   :: sc, theta_s,theta_b 
    REAL(kind=rsh)              :: CSF,csrf 

    IF (theta_s.gt.0.D0) THEN
        csrf = (1.D0 - cosh(theta_s*sc))/(cosh(theta_s) - 1.D0)
    ELSE
        csrf = -sc**2
    ENDIF
    IF (theta_b.gt.0.D0) THEN
        CSF = (exp(theta_b*csrf) - 1.D0)/(1.D0 - exp(-theta_b))
    ELSE
        CSF = csrf
    ENDIF
  END FUNCTION CSF


  !!======================================================================
  FUNCTION tool_ind2lon(x,y)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_ind2lon   ***
    !&E
    !&E ** Purpose : compute grid indexes into a longitude
    !&E
    !&E ** Called by : traj_save3d, ibm_save, eggs_grid
    !&E
    !&E ** History : 
    !&E       !  (??)
    !&E       !  2024     (M. Caillaud)
    !&E              
    !&E---------------------------------------------------------------------
    !! * Function declaration
    REAL(kind=rlg)             :: tool_ind2lon

    !! * Arguments
    REAL(kind=rsh),INTENT(in)  :: x,y

    tool_ind2lon = lonwest + (REAL(x,riolg) - 1.0_riolg)*dlonr

  END FUNCTION tool_ind2lon



  
  !!======================================================================
  FUNCTION tool_ind2lat(x,y)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_ind2lon   ***
    !&E
    !&E ** Purpose : compute grid indexes into a latitude
    !&E
    !&E ** Called by : traj_save3d, ibm_save, eggs_grid
    !&E
    !&E ** History : 
    !&E       !  2024     (M. Caillaud) 
    !&E              
    !&E---------------------------------------------------------------------
    !! * Function declaration
    REAL(kind=rlg)             :: tool_ind2lat

    !! * Arguments
    REAL(kind=rsh),INTENT(in)  :: x,y

    tool_ind2lat=latsouth+(REAL(y,riolg)-1.0_riolg)*dlatr

  END FUNCTION tool_ind2lat



  !!======================================================================
  SUBROUTINE lonlat2ij(lon_array,lat_array,target_lon,target_lat,indexi,indexj)
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_ind2lon   ***
    !&E
    !&E ** Purpose : find indexes in longitude and latitude arrays from given position
    !&E
    !&E ** Called by : None
    !&E
    !&E ** History : 
    !&E       !  2024     (M. Caillaud) 
    !&E              
    !&E---------------------------------------------------------------------

    !!* Arguments
    real, dimension(:,:),intent(in)     :: lon_array, lat_array
    real ,intent(in)                    :: target_lon, target_lat
    integer,intent(inout)               :: indexi,indexj
    !!* local declarations
    real,dimension(:,:),allocatable     :: lon_diff, lat_diff
    integer,dimension(2)                :: shp,tmp

    shp = shape(lon_array)
    ALLOCATE(lon_diff(shp(1),shp(2)))
    ALLOCATE(lat_diff(shp(1),shp(2)))

    lon_diff = abs(lon_array - target_lon)
    lat_diff = abs(lat_array - target_lat)

    ! Find the indices of the minimum differences
    tmp     = minloc(lon_diff)
    indexi  = tmp(1)
    tmp     = minloc(lat_diff)
    indexj  = tmp(2)

    DEALLOCATE(lon_diff,lat_diff)

  END SUBROUTINE lonlat2ij



  !!======================================================================
  FUNCTION tool_latlon2i(longitude,latitude)
 
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_latlon2i  ***
    !&E
    !&E ** Purpose : convert the geographic coordinates into an i grid index
    !&E
    !&E ** Called by : ibm_3d, LAGRANGIAN_init
    !&E
    !&E ** Modified variables : tool_latlon2i
    !&E
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E         !  2004-10-28
    !&E         !  2024     (M. Caillaud)  Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    IMPLICIT NONE
 
    !! * Function declaration
    REAL(kind=rsh)             :: tool_latlon2i
 
    !! * Arguments
    REAL(kind=rlg),INTENT(in)  :: longitude,latitude
     
    !! * Local declarations
    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! aie et comment fait on pour grilles non regulieres ????
    tool_latlon2i = (REAL(longitude,rlg) - lonwest)/dlonr + 1.0_rlg
 
   END FUNCTION tool_latlon2i


  !!======================================================================  
  FUNCTION tool_latlon2j(longitude,latitude)
 
    !&E---------------------------------------------------------------------
    !&E                 ***  FUNCTION tool_latlon2j  ***
    !&E
    !&E ** Purpose : convert the geographic coordinates into a j grid index
    !&E
    !&E ** Called by : ibm_3d, LAGRANGIAN_init
    !&E
    !&E ** Modified variables : tool_latlon2j
    !&E
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E         !  2004-10-28
    !&E         !  2024     (M. Caillaud)  Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    IMPLICIT NONE
 
    !! * Function declaration
    REAL(kind=rsh)             :: tool_latlon2j
 
    !! * Arguments
    REAL(kind=rlg),INTENT(in)  :: longitude,latitude
 
    !!----------------------------------------------------------------------
    !! * Executable part
 
    ! aie et comment fait on pour grilles non regulieres ????
    tool_latlon2j = (REAL(latitude,rlg) - latsouth)/dlatr + 1.0_rlg
 
   END FUNCTION tool_latlon2j

#endif

 END MODULE
