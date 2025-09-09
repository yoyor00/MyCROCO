MODULE traject3d
    !!======================================================================
    !!                   ***  MODULE traject3d  ***
    !! Trajectories dynamics:  time evolution of a patch in 3D
    !&E
    !&E Description :
    !&E
    !&E       ! At the location of ssh, the middle of the cell, the position of the particule
    !&E       ! (xpos ypos variables) is equal to the indexes i j.
    !&E       ! location of ssh : xpos=REAL(i)
    !&E       ! in a cell, REAL(i)-0.5 <= xpos <= REAL(i) + 0.5 - epsilon
    !&E       ! in a cell, REAL(j)-0.5 <= ypos <= REAL(j) + 0.5 - epsilon
    !&E
    !&E ** History :
    !&E       !  2010-05 (M. Sourisseau) Introduction of derivate type
    !&E       !  2011-01 (M. Huret) Review and introduction of generalized sigma
    !&E       !  2011-01 (R. Ramel Alyotech) Ionc library in fortran 90/95, netcdf4 and MPI
    !&E       !  2011-01 (M. Huret) call from step (if key_trajec3d) or from
    !&E       !           module ibm (if key_ibm in case of particle tracking with biology)
    !&E       !  2011-09 (T. Odaka, V. Garnier) MPI
    !&E       !  2011-10 (T. Odaka, M. Huret, V. Garnier) Introduction of MPI
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&
    !!======================================================================

#include "cppdefs.h"
#include "toolcpp.h"
#ifdef MPI
   USE mpi
#endif

#if defined LAGRANGIAN || defined DEB_IBM 
   !! * Modules used
   USE comtraj,  ONLY : imin,imax,jmin,jmax,kmax,rsh,rlg,riosh,lchain,valmanq

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC :: LAGRANGIAN_update    ! routine called by step3dxy.F90 step3dyx.F90, ibm_3d

   !! * Shared module variables

   !! * Private variables
   REAL(KIND=rsh), parameter :: r = 1.0_rsh/3.0_rsh

 CONTAINS

!!======================================================================

 SUBROUTINE LAGRANGIAN_update(xe,uz,vz,Istr,Iend,Jstr,Jend)
 
    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE traj_3d  ***
    !&E
    !&E ** Purpose : calcul de l advection dispersion d un nuage de
    !&E              particules.
    !&E
    !&E ** Description : Methode de random walk version non naive
    !&E
    !&E ** Called by : LAGRANGIAN_update_main, ibm_3d
    !&E
    !&E ** External calls : traj_save3d, h0int,xeint,ksupkinf,loc_h0,wint,dksdzint
    !&E                     splint,kzprofile,update_htot,update_wz,siggentoz, ztosiggen
    !&E                     hc_sigint,set_htot_bc,define_pos,ex_traj
    !&E
    !&E ** Reference : Visser A. W., Marine ecology progress series,
    !&E                vol 158, pp275-281, annee 1997 : Using random
    !&E                walk models to simulate the vertical distribution of
    !&E                particles in a turbulent water column
    !&E
    !&E ** History :
    !&E       !          (P. Lazure, F. Dumas)  Original code
    !&E       !  2008-08 (F. Dumas ?)
    !&E       !  2007-12 (M. Huret) ajout modif du profil kz (smooth + spline)
    !&E       !                     see (Ross and Sharples, 2004, Limnol.Oceanogr.) and 
    !&E       !                     North et al. 2006. to garantee the Well Mixed Condition test
    !&E       !                     where sharp kz gradients occur
    !&E       !  2010-01 (F Batifoulier, M Sourisseau) 
    !&E       !  2010-05 (M. Sourisseau)  Introduction of derivate type
    !&E       !  2011-01 (M. Huret) Move of biology and behaviour in module ibm
    !&E       !  2011-01 (M. Huret) Adaptation to siggen, RW in z
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(Istrm2:Iendp2,Jstrm2:Jendp2)
    !&E       !  2014-12 (M. Honnorat) Adapt for IBM upgrade
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE module_lagrangian ! sc_w,h,dt,time,om_r(dx),on_r(dy)
    USE trajinitsave,   ONLY : traj_save3d
    USE trajectools,    ONLY : h0int,xeint,ksupkinf,loc_h0,wint,dksdzint, &
                               splint,kzprofile,update_htot,update_wz,    &
                               siggentoz, ztosiggen, hc_sigint,set_htot_bc,define_pos
#ifdef MPI
    USE toolmpi,        ONLY : ex_traj
    USE comtraj,        ONLY : down_give, up_give, right_give, left_give
#endif
    USE comtraj,        ONLY : patches, type_patch, type_particle, type_position,ndtz,htx,hty,wz

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,4),      INTENT( in )    :: xe
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax,3), INTENT( in )    :: uz,vz
    INTEGER,                                           INTENT( in )    ::  Istr,Iend,Jstr,Jend

    !! * Local declarations
    ! Indexes on number of patch and particles
    INTEGER                                     :: npa,npart,ndt,nb_patch

    ! Old and intermediate position of particles, as bathymetry at old and intermediate position of particles
    REAL(KIND=rsh)                              :: x_old,y_old,s_old,z_old,z_int,s_int,slag,     &
                                                   x_mid,y_mid,z_mid,d3_mid,xe_mid,h0_mid,       &
                                                   tir,dksdz,kzz,ds_adv,ds_dif,d3,d3_final,      &
                                                   xe_final,h0_final,d3avt,xeavt,hcavt,h0avt,    &
                                                   hc_sig_mid,hc_sig_final,zposf,zlag

    ! For Random walk
    REAL(KIND=rsh)                              :: lb,lt,lbs,lts,ds_bio,wbio
    REAL(KIND=rsh), DIMENSION(0:kmax)           :: dkz2,kzsmth                  ! For spline
    REAL(KIND=rsh), DIMENSION(0:kmax)           :: zpos
    INTEGER                                     :: kw, kwm

    ! To save a local and global position of particle for MPI and Sequential compatibility
    TYPE(type_position)                         :: pos_old,pos_mid,pos_temp
    
    ! Integer for loops
    INTEGER                                     :: l

    ! Model time step
    REAL(KIND=rlg)                              :: dtm,dtz
    INTEGER                                     :: time_step

    TYPE(type_patch),                  POINTER  :: patch    => NULL()
    TYPE(type_particle),               POINTER  :: particle => NULL()
    

    ! Parameters for interpolation at particle location
    REAL(KIND=rsh)       :: px,py
    INTEGER              :: igg,idd,jbb,jhh,hlb,hlt,hrb,hrt,ipos,jpos

    ! Integer for MPI errors
    INTEGER              :: ierr_mpi

# include "compute_auxiliary_bounds.h"
    !!----------------------------------------------------------------------
    !! * Executable part

    ! Initialization
    !------------
    dtm = dt
    nb_patch = patches % nb
    time_step=nrhs

    !Update Htot
    !------------
    call update_htot(Istr,Iend,Jstr,Jend,IstrU,JstrV)
    call set_htot_bc(Istr,Iend,Jstr,Jend,IstrU,JstrV)

    !Update wz
    !------------------
    CALL update_wz(Istr,Iend,Jstr,Jend)
#ifdef MPI
    call exchange_w3d_tile (Istr,Iend,Jstr,Jend,wz(START_2D_ARRAY,0))
#endif

    ! Implement intermediate time step, for Random Walk
    IF ( nb_patch /= 0 ) dtz = dtm/real(ndtz,kind=rlg)

    patch => patches%first
    DO npa = 1,nb_patch
   
        IF ( (time < patch % t_beg) .OR. (time > patch % t_end) ) THEN
            patch => patch % next
            CYCLE
        END IF

        ! We modify spos from zpos because it was first calculated at simulation intial time step with another xe...
        ! to be removed when reading of positions at debtraj time step and not beginning of run
        IF ( (time-patch%t_beg) <= dtm ) THEN ! first pass for this patch
            
            ! We modify spos from zpos because it was first calculated at simulation intial time step with another xe...
            DO npart = 1, patch % nb_part_alloc
                particle => patch % particles(npart)

                IF ( .NOT. particle % active ) CYCLE

                ! Get local position if MPI, doesn't change anything in sequential
                pos_temp%xp = particle%xpos; pos_temp%yp = particle%ypos
                call define_pos(pos_temp)

                ! total depth at particle s location
                CALL loc_h0(pos_temp%idx_r,pos_temp%idy_r,px,py,igg,idd,jbb,jhh, &
                            hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
                particle%xe = xeint(xe(:,:,time_step),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                                    Istr,Iend,Jstr,Jend)
                particle%d3 = particle%h0 + particle%xe
                IF ( particle%d3 > particle%zpos ) THEN                
                    zlag = -(particle%zpos) + particle%xe ! switch from immersion to real z
                    CALL ztosiggen(zlag, slag, particle%xe, particle%h0, particle%hc)
                    particle  % spos = slag
                ELSE
                    particle%zpos = particle%h0
                    particle%spos = -1.0_rsh
                END IF
            END DO
        END IF

        patch => patch % next
            
    END DO

!**********************************
! Ecriture des sorties (before any advection to get the exact release location in the output file
!**********************************
    CALL traj_save3d

#ifdef MPI
    down_give  = 0
    up_give    = 0
    right_give = 0
    left_give  = 0
#endif

   !***************************
   ! Start loop on patches 
   !***************************
    patch => patches % first
    DO npa = 1,nb_patch

        IF ( (time < patch % t_beg) .OR. (time > patch % t_end) ) THEN
           patch => patch % next
           CYCLE
        END IF
    
        DO npart = 1, patch % nb_part_alloc
            particle => patch % particles(npart)
            ! Skip if particle is inactive.
            IF ( .NOT. particle % active ) CYCLE
#ifdef IBM_SPECIES
            ! Skip if species stage not appropriate
            IF ( particle%stage >= 5 .OR. particle%super <= 0.0_rsh )  CYCLE
#endif
            ! Skip if flag is missing...
            IF ( particle % flag == -valmanq )  CYCLE

            ! Get local position if MPI, doesn't change anything in sequential
            pos_temp%xp = particle%xpos; pos_temp%yp = particle%ypos
            call define_pos(pos_temp)

            ! depth of cell in which particle is located
            d3 = h(pos_temp%idx,pos_temp%idy) + xe(pos_temp%idx,pos_temp%idy,time_step)
            IF ( d3 <= 0.0_rsh ) CYCLE

            !***************************
            ! horizontal advection 
            !***************************

            ! Save former position
            pos_old = pos_temp 
            s_old =   particle % spos
            d3avt =   particle % d3
            xeavt =   particle % xe
            h0avt =   particle % h0
            hcavt =   particle % hc

            ! along-sigma advection
            CALL avance(uz(:,:,:,time_step),vz(:,:,:,time_step),xe(:,:,time_step),&
                        dtm,pos_temp,particle%spos,particle%flag,&
                        Istr,Iend,Jstr,Jend) 
            CALL loc_h0(pos_temp%idx_r,pos_temp%idy_r,px,py,igg,idd,jbb,jhh, &
                        hlb,hrb,hlt,hrt,Istr,Iend,Jstr,Jend)
            xe_final = xeint(xe(:,:,time_step),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                             Istr,Iend,Jstr,Jend)
            h0_final     = h0int(    px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
            d3_final     = h0_final + xe_final
            hc_sig_final = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

            ! Mid-position between old(s_old) and new (particle%xpos)
            pos_mid%xp = 0.5_rsh*(pos_temp%xp + pos_old%xp)
            pos_mid%yp = 0.5_rsh*(pos_temp%yp + pos_old%yp)
            call define_pos(pos_mid)

            CALL loc_h0(pos_mid%idx_r,pos_mid%idy_r,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt, & 
                        Istr,Iend,Jstr,Jend)
            xe_mid = xeint(xe(:,:,time_step),px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                           Istr,Iend,Jstr,Jend)
            h0_mid     = h0int(    px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)
            d3_mid     = h0_mid + xe_mid
            hc_sig_mid = hc_sigint(px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt)

            particle%xpos = pos_temp%xp; particle%ypos = pos_temp%yp
            IF ( (d3_final > 0.0_rsh).AND.(d3_mid > 0.0_rsh) ) THEN
                particle % d3 = d3_final
                particle % xe = xe_final
                particle % h0 = h0_final
                particle % hc = hc_sig_final
            ELSE
                write(*,*) 'WARNING: Particule a la cote'
                write(*,*) d3_final,h0_final,d3_mid,h0_mid,npa,npart
                particle % xpos = pos_old%xp ! on remet dans l eau
                particle % ypos = pos_old%yp 
                stop
            endif
#ifdef MPI
            ! check if we need to exchange or not
            ! -----------------------------------
            ipos = NINT(particle%xpos)
            jpos = NINT(particle%ypos)

            ! Here, we verify if the particle stays in your zone
            ! although it gets out of its zone; it will only move max to i+-1 and j+-1 
            ! thus we can compute the particle s new value in this time step without pb
            !
            ! If the particle went out the zone, we check on which side it went out and
            ! count it in the proper variable 'down_give', 'up_give', 'right_give' or 'left_give'
            ! (used to call ex_traj)
            !
            ! The neighbor domain where the particle went is coded in variable 'limitbye' according
            ! to this scheme :
            !         7  |  6  |  5
            !       -----+-----+-----
            !         8  |  0  |  4
            !       -----+-----+-----
            !         1  |  2  |  3

            IF ( jpos < jminmpi ) THEN
                particle % limitbye = 2
                down_give = down_give+1
                IF ( ipos < iminmpi ) particle % limitbye = 1
                IF ( ipos > imaxmpi ) particle % limitbye = 3
            ELSE IF ( jpos > jmaxmpi ) THEN
                particle % limitbye = 6
                up_give = up_give+1
                IF ( ipos < iminmpi ) particle % limitbye = 7
                IF ( ipos > imaxmpi ) particle % limitbye = 5
            ELSE IF ( ipos < iminmpi ) THEN
                particle % limitbye = 8
                left_give = left_give+1
            ELSE IF ( ipos > imaxmpi ) THEN
                particle % limitbye = 4
                right_give = right_give+1
            END IF
#endif

            !**********************************
            ! Vertical advection
            ! is estimated from the environment at the averaged location x_mid,y_mid
            !**********************************   
            IF ( (d3_final <= 0.0_rsh).OR.(d3_mid <= 0.0_rsh) )  CYCLE

            IF ( particle%itypevert == 0 ) THEN

                !**********************************
                ! No vertical advection (float trajectory at constant immersion)
                !**********************************    
                particle % zpos = MIN(particle%zpos,d3_final)   ! when particle wants to go in shallower areas than immersion 
                zlag = -particle%zpos + xe_final                ! switch from immersion to real z
                CALL ztosiggen(zlag,slag,xe_final,h0_final,hc_sig_final)
                particle%spos = MAX(slag,-1.0_rsh)              ! when particle wants to go in shallower areas than immersion  

            ELSE   
                !**************************************
                ! Vertical advection + random walk
                !**************************************

                ! vertical displacement induced by advection : ds_adv
                ds_adv = dtm*wint(s_old,xe_mid,h0_mid,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                                  Istr,Iend,Jstr,Jend)
            
                IF ( particle%itypevert > 0 ) THEN ! Random walk (turbulent diffusion)
                
                    !*********************************************************
                    ! + Random walk (turbulent diffusion)
                    !*********************************************************
                    ! Switch to z coordinate for RW
                    ! allow gain in computing time when in siggen
                    ! because then no switch of coordinate in the RW loop
                
                    DO l=0,kmax
                       CALL siggentoz(zpos(l),sc_w(l),xe_mid,h0_mid,hc_sig_mid)
                    END DO
                
                    ! smooth and spline the kz profile at x_mid,y_mid location
                    CALL kzprofile(kzsmth,dkz2,zpos,px,py,igg,idd,jbb,jhh,hlb,hrb,hlt,hrt,&
                                   Istr,Iend,Jstr,Jend)
                
                    ! vertical displacement induced by diffusion : ds_dif
                    ! displacement linked to kz intensity and structure
                    ! ( kz*(dkz/dz) see Visser)
                    ! loop on t to use a smaller time step
                    ! sous-boucle de calcul du deplacement vertical
                    ! lie a kz + structure de  kz (dkz/dz)
                    ! deplacement net du centre de masse (voir Visser) : ds_dif
                    ! s_int positions intermediaires successives
                    ! s_mid position mediane entre deplacement du centre de masse et s_int

                    CALL siggentoz(z_old,s_old,xe_mid,h0_mid,hc_sig_mid)
                    z_int = z_old
                
                    IF ( d3_mid > 5.0_rsh ) THEN   ! attention : lb+lt for boundary layers should
                                                   ! be lower than criteria here
                        DO ndt=1,ndtz
                            ! get the first derivative of Kz at z_int location
                            dksdz  = dksdzint(z_int, zpos, kzsmth, dkz2)                                       
                            ds_dif = dksdz*dtz 
                
                            z_mid = MIN(MAX(z_int + ds_dif*0.5_rsh, zpos(0)), zpos(kmax))
                            ! get the kz value at s_mid location entre son centre de masse et sa position initiale
                            ! ancien calcul         kzz=kzint(x_mid,y_mid,s_mid)    
                            CALL splint(zpos, kzsmth, dkz2, kmax+1, z_mid, kzz)
                
                            IF ( kzz < 0.0 ) kzz = 0.0                  ! kzz should be > 0 but...
                            
                            ! random walk component of vertical diffusion (between 0 and 1)
                            CALL random_number(tir)
                
                            ! s_int=s_int+ds_dif+(2.0_rsh*tir-1.0_rsh)*sqrt(2.0_rsh*dtz*kzz/(r*d3_mid*d3_mid))  en sigma
                            z_int = z_int + (ds_dif + (2.0_rsh*tir - 1.0_rsh)*sqrt(2.0_rsh*dtz*kzz/r))
                            
                            ! handle boundary conditions
                            ! s_int=MIN(MAX(s_int,-1.0_rsh),0.0_rsh)         ! particle set on boundary
                            ! s_int=MIN(MAX(s_int,-(2_rsh+s_int)),-s_int)  ! reflective boundary condition
                            ! random mixed layer to avoid accumulation (see Ross and Sharples, 2004)
                            ! lb=0.5_rsh        ! distance in meters for the random boundary layer
                            ! lt=0.5_rsh
                            ! lts=lt/d3_mid
                            ! lbs=lb/d3_mid
                            ! CALL random_number(tir)
                            ! IF (s_int > -lts) s_int=-tir*lts
                            ! IF (s_int < -1_rsh+lbs) s_int=tir*lbs-1_rsh
                
                            ! handle boundary conditions
                            ! z_int=MIN(MAX(z_int,-h0_mid),xe_mid)         ! particle set on boundary
                            ! z_int=MIN(MAX(z_int,-(2_rsh*h0_mid+z_int)),2_rsh*xe_mid-z_int)  ! reflective boundary condition
                            ! random mixed layer to avoid accumulation (see Ross and Sharples, 2004)
                            lb = 2.0_rsh        ! distance in meters for the random boundary layer
                            lt = 2.0_rsh        ! twice the boundary layer should be ok
                            CALL random_number(tir)
                            IF (z_int > zpos(kmax)-lt) z_int = zpos(kmax) - tir*lt
                            IF (z_int < zpos(0)+lb)    z_int = zpos(0)    + tir*lb
                
                        END DO   ! end loop on RW time
                    ENDIF  ! if d3_mid > 5.0

                    ! back to sigma position
                    kw  = kmax
                    kwm = 0
                    CALL ksupkinf(z_int, zpos, kmax+1, kw, kwm, 2)

                    particle%spos = ((z_int-zpos(kwm))*sc_w(kw)+ &
                                     (zpos(kw)-z_int)*sc_w(kwm))/(zpos(kw)-zpos(kwm))  !! Interpolation lineaire
                
                ENDIF ! end of RW (itypevert > 0)

                !**************************************
                ! Vertical behavior (larval behavior)
                !**************************************

                ! Behaviour: speed of vertical migration
                wbio = 0.0_rsh ! may be done here (if simple) or in module IBM

!#ifdef key_ibm_behav
!                  CALL w_behav_part(n,m,sal,temp)
!                  wbio = particle%w
!#endif
                ds_bio = wbio*dtm/d3_mid  ! approximation
                ! End of vertical migration
                !**************************************

                ! Estimate of the new vertical position
                particle%spos = MIN(MAX(particle%spos + ds_adv + ds_bio, -1.0_rsh),0.0_rsh) 
                !particle%spos = MIN(MAX(particle%spos+ds_bio,-1.0_rsh),0.0_rsh) 

                CALL siggentoz(zposf,particle%spos, xe_final, h0_final, hc_sig_final)
                particle%zpos = -zposf + xe_final           ! immersion
            END IF     ! ends test on vertical movement (itypevert)

        END DO   ! ends loop on number of particles

        patch => patch % next

    END DO   ! Ends loop on number of patches

    CALL_MPI ex_traj(down_give, up_give, right_give, left_give)

 END SUBROUTINE LAGRANGIAN_update 





 !!====================================================================== 
 SUBROUTINE avance(uz,vz,xe,deltat,pos,sig0,statp,Istr,Iend,Jstr,Jend)

    !&E---------------------------------------------------------------------
    !&E                 ***  ROUTINE avance  ***
    !&E
    !&E ** Purpose : calcul du deplacement cinematique d une particule
    !&E              sous l action d un courant et d une composante
    !&E              aleatoire de deplacement (random walk), dans
    !&E              le plan horizontal.
    !&E              Le deplacement sous l'action du courant est calcule avec 
    !&E              un schema de Runge-Kutta 2
    !&E
    !&E
    !&E ** Called by : LAGRANGIAN_update
    !&E
    !&E ** External calls : uint,vint,define_pos,MPI_glob2loc
    !&E
    !&E ** Reference :
    !&E
    !&E ** History :
    !&E       !  08-2002 (F. Dumas, P. Lazure)
    !&E       !  2011-11 (V. Garnier) Spatial extension of ssh(Istrm2:Iendp2,Jstrm2:Jendp2)
    !&E       !  2024    (M. Caillaud, D. Gourves) Coupled with CROCO
    !&E
    !&E---------------------------------------------------------------------
    !! * Modules used
    USE comtraj,            ONLY : valmanq, type_position
    USE module_lagrangian ! on_r,om_r
    USE trajectools,        ONLY: uint,vint,define_pos
#ifdef MPI
    USE toolmpi,            ONLY : MPI_glob2loc
#endif

    !! * Arguments
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY,kmax), INTENT( in )               :: uz,vz
    REAL(KIND=rsh), DIMENSION(GLOBAL_2D_ARRAY),      INTENT( in )               :: xe
    TYPE(type_position),                             INTENT( inout )            :: pos
    REAL(KIND=rlg),                                  INTENT( in )               :: deltat
    REAL(KIND=rsh),                                  INTENT( in )               :: sig0
    INTEGER,                                         INTENT( in )               :: Istr,Iend,Jstr,Jend
    REAL(KIND=rsh),                                  INTENT( inout ),OPTIONAL   :: statp

    !! * Local declarations
    INTEGER                                  :: j0,jst,i0,ist,i1,j1
    REAL(KIND=rsh)                           :: uxp,uyp,ux0,uy0,xst,yst

    ! For randow movement
    REAL(KIND=rsh)                           :: a,b,tir1,tir2

    REAL(KIND=rsh)                           :: dkx,disp
    REAL(KIND=rsh)                           :: dmaill,tetha

    TYPE(type_position)                      :: pos0,pos1,posint

    !!----------------------------------------------------------------------
    !! * Executable part
    dmaill = 0.01_rsh ! plus eleve plus jump en z, car advecte le long des sigmas, martin
    disp   = 0.0_rsh
    pos0   = pos

    !interpolate current components on particle position
    !----------------------------------------------------
    ux0=uint(pos0%idx_r,pos0%idy_r,sig0,uz,xe,Istr,Iend,Jstr,Jend)
    uy0=vint(pos0%idx_r,pos0%idy_r,sig0,vz,xe,Istr,Iend,Jstr,Jend)

    posint%xp = pos0%xp + ux0*deltat/om_r(nint(pos0%idx_r),nint(pos0%idy_r))
    posint%yp = pos0%yp + uy0*deltat/on_r(nint(pos0%idx_r),nint(pos0%idy_r)) 
    CALL define_pos(posint) 
    uxp = uint(posint%idx_r,posint%idy_r,sig0,uz,xe,Istr,Iend,Jstr,Jend)
    uyp = vint(posint%idx_r,posint%idy_r,sig0,vz,xe,Istr,Iend,Jstr,Jend)

    pos1%xp = pos0%xp + 0.5_rsh*deltat*(ux0/om_r(nint(pos0%idx_r),nint(pos0%idy_r))+uxp/om_r(nint(posint%idx_r),nint(posint%idy_r)))
    pos1%yp = pos0%yp + 0.5_rsh*deltat*(uy0/on_r(nint(pos0%idx_r),nint(pos0%idy_r))+uyp/on_r(nint(posint%idx_r),nint(posint%idy_r)))
    CALL define_pos(pos1)

    ! add random dispersion 
    IF (disp /= 0.0_rsh) then
        call random_number(tir1)
        tetha = 2.0_rsh*pi*tir1
        call random_number(tir2)
        dkx  = sqrt(2.0_rsh*disp*deltat*3.0_rsh)*tir2
        IF (0.5_rsh*(ux0+uxp) /= 0.0_rsh .and. 0.5_rsh*(uy0+uyp) /= 0.0_rsh) THEN
            pos1%xp = pos1%xp + dkx*cos(tetha)/om_r(nint(pos1%idx_r),nint(pos1%idy_r))
            pos1%yp = pos1%yp + dkx*sin(tetha)/on_r(nint(pos1%idx_r),nint(pos1%idy_r))
        END IF
    ENDIF

    !     detection d un franchissement de mur hy
    !--------------------------------------------------
    j0 = nint(pos0%yp)
    j1 = nint(pos1%yp)
    IF (j0 /= j1) THEN
        jst = min (j0,j1)
        yst = real(jst,kind=rsh) + 0.5_rsh
        xst = pos0%xp
        IF (pos1%xp /= pos0%xp) THEN
            a = (pos1%yp - pos0%yp)/(pos1%xp - pos0%xp)
            b =  pos0%yp - a*pos0%xp
            xst = (yst-b)/a
        END IF
        ist = nint(xst)

#ifdef MPI
        call MPI_glob2loc(ist,jst)
#endif
        ! To catch the good v current in the croco grid.
        ! Otherwise it takes a v current not at the wall but the cell before on eta_rho
        IF (vz(ist,jst+1,1) == 0.0_rsh) THEN

            !------------------------------------------------------------
            !    real(j1-j0) est juste la pour donner le signe
            !    et savoir si l on s arrete au-dessus ou au-dessous
            !    du mur.
            !    Pourquoi xst n est-il pas modifie ?
            !    Reflection a la normale de la position finale sur terre
            !------------------------------------------------------------

            pos1%yp = yst - dmaill*real((j1-j0),kind=rsh)
            !update index for pos1
            call define_pos(pos1) 
            statp = statp + 1
        END IF
    END IF

   !     detection d un franchissement de mur hx
   !--------------------------------------------------

    i0 = nint(pos0%xp)
    i1 = nint(pos1%xp)
    IF (i0 /= i1) THEN
        ist = min (i0,i1)
        xst = real(ist,kind=rsh) + 0.5_rsh
        yst = pos0%yp
        IF (pos1%yp /= pos0%yp) THEN
            a = (pos1%xp - pos0%xp)/(pos1%yp - pos0%yp)
            b =  pos0%xp - a*pos0%yp
            yst = (xst-b)/a
        END IF
        jst = nint(yst)
#ifdef MPI
        call MPI_glob2loc(ist,jst)
#endif
	! To catch the good u current in the croco grid.
	! Otherwise it takes a u current not at the wall but the cell before on xi_rho
        IF (uz(ist+1,jst,1) == 0.0_rsh) THEN
            !------------------------------------------------------------
            !    real(i1-i0) est juste la pour donner le signe
            !    et savoir si l on s arrete au-dessus ou au-dessous
            !    du mur.
            !    Pourquoi yst n est-il pas modifie ?
            !------------------------------------------------------------

            pos1%xp = xst - dmaill*real((i1-i0),kind=rsh)
            call define_pos(pos1)
            statp = statp + 1
        END IF
    END IF

    ! detection d un franchissement limite de domaine
    IF (pos1%yp <= 1.0_rsh .or. pos1%yp >= real(jmax,kind=rsh)) THEN
        statp = -valmanq ! on prefere flagger la particule et la laisser ou elle est (ne bougera plus) 
    END IF
    IF (pos1%xp <= 1.0_rsh .or. pos1%xp >= real(imax,kind=rsh)) THEN
        statp = -valmanq ! on prefere flagger la particule et la laisser ou elle est (ne bougera plus) 
    END IF

    pos = pos1

    IF (h(pos%idx,pos%idy) == -valmanq) THEN ! pour cas special ou passage exact sur coin de cellule, avec un mur unique partant de la section, tests precedents ne marche pas forcement
        pos   = pos0
        statp = statp + 1
    ENDIF

   END SUBROUTINE avance

#endif

END MODULE
