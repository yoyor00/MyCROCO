MODULE OBSTRUCTIONS

#include "cppdefs.h"

#ifdef OBSTRUCTIONS
   !&E==========================================================================
   !&E                   ***  MODULE  OBSTRUCTIONS  ***
   !&E
   !&E
   !&E ** Purpose : concerns all subroutines related to obstructions interactions
   !&E    with hydrodynamics
   !&E
   !&E ** Description :
   !&E     
   !&E     subroutine OBSTRUCTIONS_update             ! Main subroutine controlling updates of
   !&E                                                ! obstructions characteristics at each half dt
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_uzvz          ! Computes 3D or pseudo-3D velocities
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_abdelposture  ! Computes flexible obstructions posture
   !&E                                                ! based on Abdelhrman 2007
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_height        ! Computes the height of obstructions 
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_theta         ! Computes the bending angle of flexible 
   !&E                                                ! obstructions from vertical
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_fracz         ! Computes the vertical fraction of layer 
   !&E                                                ! effectively occupied by obstructions
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_diam          ! Computes width and thickness of obstructions 
   !&E                                                ! resulting from bending
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_distrib       ! Computes the vertical distribution of 
   !&E                                                ! obstructions density
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_fracxy        ! Computes correction term for obstruction
   !&E                                                ! coverage (i.e. fragmentation) within one cell
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_projarea      ! Computes obstruction vertical area facing flow (called s)
   !&E                                                ! and horizontal area (called a)
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_obstroughness ! Computes obstruction induced roughness length
   !&E                                                ! used to correct velocity profile in 3D-small-depth (with
   !&E                                                ! turbulence model), or for fully simplified formulation
   !&E
   !&E     subroutine OBSTRUCTIONS_comp_hydroparam    ! Computes source terms used within 
   !&E                                                ! momentum equation and tubulence closure
   !&E
   !&E
   !&E ** History :
   !&E     ! 2011-03 (F. Ganthy) Original code
   !&E     ! 2012-04 (F. Ganthy) Adaptation for different kinds of obstructions 
   !&E                           and simplification for an easier use
   !&E     ! 2014-02 (F. Ganthy) Add subroutine obst_write_summary
   !&E     ! 2014-02 (F. Ganthy) Update for Ionc library netcdf 4
   !&E     ! 2014-02 (F. Ganthy) Add the computation of bottom shear stress
   !&E     ! 2014-08 (F. Ganthy) Some modifications, cleaning, debugging and armoring
   !&E     ! 2014-10 (F. Ganthy) More modifications + computation of obstructions
   !&E                           posture following Abdelrhman 2007
   !&E     ! 2015-09 (F. Ganthy) Some fixes
   !&E     ! 2016-03 (F. Ganthy) Add more outputs
   !&E     ! 2016-03 (F. Ganthy) Some fixes and reorganizations:
   !&E                             - Limit obst_kmin and obst_kmax to (1,kmax) within obst_comp_abdelposture
   !&E                             - Remove computation of obst_theta in obst_comp_diam (already done within
   !&E                               obst_comp_height or obst_comp_abdelposture
   !&E                             - Change method for computation of obst_thick2d and obst_with2d in 3D (was geometrical
   !&E                               averaged, replaced by direct computation as for 2D model)
   !&E                             - Change computation method (and associate logical name) for obstructions drag coefficient
   !&E                             - Add logical to choose turbulence coefficient values (default or user-defined)
   !&E                             - Modification to better taken into account for obstructions parameters in 3D when
   !&E                               a sigma layer is not completely filled (in z direction) by obstructions
   !&E     ! 2016-03 (F. Ganthy) Remove useless option l_obst_2donly
   !&E     ! 2016-09 (F. Ganthy) Some formatting in subroutines history, PRINT_DEBUG, error and warning message,
   !&E                           as well as messages in simulog.
   !&E     ! 2016-09 (F. Ganthy) New and optimized computation of flexible obstruction height using Abdelrhman's 2007
   !&E                           method
   !&E     ! 2017-02 (F. Ganthy) Some modifications:
   !&E                             - Allowing multiple obstructions type in a single grid cell
   !&E                               --> allowed multispecific computation.
   !&E                               This imply that some tables must be allocated depending on (iv,k,i,j) or (iv,i,j)
   !&E                               These allocations are done within obst_alloc_xyz (because tables allocated within
   !&E                               obst_alloc_nbvar are only those read within namelist).
   !&E                             - Changes on instantaneous obstruction state variables for future coupling with Zostera growth module
   !&E                             - Differenciation of cylindric / parallelepipedic structures
   !&E                             - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E                             - Cleaning (removing useless parameters and tests)
   !&E     ! 2017-04 (F. Ganthy) Change initialization of turbulence coefficients cmu and c2turb
   !&E     ! 2017-04 (F. Ganthy) Add subroutine obst_comp_obstroughness
   !&E     ! 2017-04 (F. Ganthy) Add subroutine obst_comp_projarea
   !&E     ! 2017-04 (F. Ganthy) Remove previous allocation related to vertical growth of zostera root system (not purpose of obstruction module)
   !&E     ! 2017-04 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E     ! 2017-04 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E                           --> Add subroutine obst_comp_uzvz to compute uz and vz which will be used
   !&E                               during obstructions procedure (in 3D, 3D-small-depth and 2D)
   !&E     ! 2017-10 (F. Ganthy) Re-organization and cleaning related to 3D obstructions, multiple obstructions within
   !&E                           one grid cell, and add simplified obstructions kind (based on roughness length)
   !&E     ! 2017-11 (F. Ganthy) Remove useless variables related to turbulence coefficients
   !&E     ! 2018-01 (F. Ganthy) Some changes related to coupling with sediment module
   !&E     ! 2018-02 (F. Ganthy) Bug corrected within comp_distrib subroutine
   !&E     ! 2018-02 (F. Ganthy) Changes on parameterization for sedimentological roughness length within obstructions
   !&E     ! 2018-04 (F. Ganthy) Moved variable declaration and allocation within comobstructions.F90
   !&E     ! 2018-04 (F. Ganthy) Parameterization of obstructions roughness length (no turb) using Abdelrhman 2003 formulations
   !&E     ! 2018-04 (F. Ganthy) Theoretical velocity profiles (2D and 3D small-depth) using Abdelrhman 2003 formulations
   !&E     ! 2021-10 (F. Ganthy) Initialization routines extracted to initOBSTRUCTIONS.F90
   !&E==========================================================================

   !! * Modules used
   USE comOBSTRUCTIONS


   IMPLICIT NONE

   !! * Accessibility

   ! function & routines of this module, called outside :
   ! PUBLIC functions
   PUBLIC OBSTRUCTIONS_update
   PUBLIC OBSTRUCTIONS_comp_uzvz
   PUBLIC OBSTRUCTIONS_comp_abdelposture
   PUBLIC OBSTRUCTIONS_comp_height
   PUBLIC OBSTRUCTIONS_comp_theta
   PUBLIC OBSTRUCTIONS_comp_fracz
   PUBLIC OBSTRUCTIONS_comp_diam
   PUBLIC OBSTRUCTIONS_comp_distrib
   PUBLIC OBSTRUCTIONS_comp_fracxy
   PUBLIC OBSTRUCTIONS_comp_projarea
   PUBLIC OBSTRUCTIONS_comp_obstroughness
   PUBLIC OBSTRUCTIONS_comp_hydroparam
   PUBLIC OBSTRUCTIONS_comp_bedroughness

   PRIVATE

   !! * Shared or public module variables

   !! * Private variables
    
  CONTAINS


   !!==========================================================================================================
 
   SUBROUTINE OBSTRUCTIONS_update(limin, limax, ljmin, ljmax, cm0, h, &
                        z0b, ssh, ubar, vbar, u, v)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_update  ***
   !&E
   !&E ** Purpose : Update of obstruction parameters at each time step
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08    (F. Ganthy) Simplification and re-organization
   !&E       ! 2014-10    (F. Ganthy) Some modifications + computation of obstructions
   !&E                                posture following Abdelrhman 2007
   !&E       ! 2015-09    (F. Ganthy) Minor correction (obst_cmu)
   !&E       ! 2016-03    (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2017-02    (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       ! 2017-04    (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04    (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10    (F. Ganthy) Optimization for full 3D obstructions
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
   USE initOBSTRUCTIONS, ONLY : OBSTRUCTIONS_readfile_char
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax
   INTEGER(KIND=rsh), INTENT(IN) :: cm0
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: h 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: z0b 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: ssh 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: ubar
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: vbar 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax, kmax),INTENT(IN) :: u
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax, kmax),INTENT(IN) :: v

   !! * Local declaration

   !!----------------------------------------------------------------------
   !! * Executable part
   ! ******************************
   ! * READING CHARACTERISTICS FILE
   ! ******************************
   CALL OBSTRUCTIONS_readfile_char(limin,limax,ljmin,ljmax)
   ! ************************
   ! * COMPUTES 3D VELOCITIES
   ! ************************
   CALL OBSTRUCTIONS_comp_uzvz(limin,limax,ljmin,ljmax)
   ! *********************************
   ! * UPDATING OBSTRUCTION PARAMETERS
   ! *********************************
   !---------------------
   ! * Obstruction height
   !---------------------
   CALL OBSTRUCTIONS_comp_abdelposture(limin,limax,ljmin,ljmax)
   CALL OBSTRUCTIONS_comp_height(limin,limax,ljmin,ljmax)
   !----------------------------
   ! * Obstruction bending angle
   !----------------------------
   CALL OBSTRUCTIONS_comp_theta(limin,limax,ljmin,ljmax)
   !--------------------------------------
   ! * Obstruction fraction of sigma layer
   !--------------------------------------
   CALL OBSTRUCTIONS_comp_fracz(limin,limax,ljmin,ljmax)
   !----------------------
   ! * Obstruction density
   !----------------------
   CALL OBSTRUCTIONS_comp_distrib(limin,limax,ljmin,ljmax)
   !------------------------------
   ! * Obstruction width/thickness
   !------------------------------
   CALL OBSTRUCTIONS_comp_diam(limin,limax,ljmin,ljmax)
   !---------------------------------------------------
   ! * Obstruction density correction for fragmentation
   !---------------------------------------------------
   CALL OBSTRUCTIONS_comp_fracxy(limin,limax,ljmin,ljmax)
   !---------------------------------------------
   ! * Obstruction projected area (obst_a,obst_s)
   !---------------------------------------------
   CALL OBSTRUCTIONS_comp_projarea(limin,limax,ljmin,ljmax)
   !---------------------------------------
   ! * Obstruction induced roughness length
   !---------------------------------------
   CALL OBSTRUCTIONS_comp_obstroughness(limin,limax,ljmin,ljmax)
   ! ***********************************
   ! * UPDATING OBSTRUCTION SOURCE TERMS
   ! ***********************************
   CALL OBSTRUCTIONS_comp_hydroparam(limin,limax,ljmin,ljmax)
   ! **********************************
   ! * UPDATING OBSTRUCTION OTHER TERMS
   ! **********************************
   !-------------------------------------------
   ! * ROUGHNESS LENGTH FOR BOTTOM SHEAR STRESS
   !-------------------------------------------
! * FG-OBST : Improve coupling with sediment
#ifndef MUSTANG
   ! Only used for running the model without sediment
   ! For sediment run, this the function is within tools_obstructions
   ! and median grain size is used, the call to the function is performed
   ! within sed_MUSTANG_MARS.F90/sed_skinstress
   IF(l_obst_z0bstress_tot) CALL OBSTRUCTIONS_comp_bedroughness(limin,limax,ljmin,ljmax)
#endif
! * FG-OBST-END
   ! *********************************
   END SUBROUTINE OBSTRUCTIONS_update 
 
   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_uzvz(limin,limax,ljmin,ljmax,h0,ssh,uz,vz)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_uzvz  ***
   !&E
   !&E ** Purpose : Computes 3D (or pseudo-3D) velocities used over the whole
   !&E              obstructions procedure
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : ssh, u, v, h0, uz, vz
   !&E
   !&E ** Modified variables : obst_zc,obst_dz,obst_uz,obst_vz
   !&E
   !&E ** Reference : 
   !&E
   !&E ** History :
   !&E       ! 2017-04-14 (F. Ganthy) Original code
   !&E       ! 2018-04-19 (F. Ganthy) Modification to provide theoretical velocity
   !&E                                profiles for obstructions based on Abdelrhman's
   !&E                                 (2003, MEPS) method
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE
   !! * Arguments
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: h0 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax),INTENT(IN) :: ssh 
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax, kmax),INTENT(IN) :: uz
   REAL(KIND=rsh),DIMENSION(imin:imax, jmin:jmax, kmax),INTENT(IN) :: vz
   !! * Local declaration
   INTEGER                   :: i,j,k
   REAL(KIND=rsh)            :: hwat
   !!----------------------------------------------------------------------
   !! * Executable part
   !---------------------------------
   obst_zc(:,:,:) = 0.0_rsh
   obst_dz(:,:,:) = 0.0_rsh
   obst_uz(:,:,:) = 0.0_rsh
   obst_vz(:,:,:) = 0.0_rsh
   DO j=ljmin,ljmax
     DO i=limin,limax
       hwat = h0(i,j) + ssh(i,j)  ! z_w(i,j,N)+h(i,j)
       IF(hwat.GT.h0fond) THEN
         !--------------------------
         ! * Height above bed, dz, uz, vz
         !--------------------------
         DO k=1,obst_kmax
            obst_zc(k,i,j) = (1.0_rsh+obst_sig(k))*hwat !TODO : compute sig and dsig with generalized sigmas of croco
            obst_dz(k,i,j) = obst_dsig(k)*hwat
            obst_uz(k,i,j) = (uz(k,i,j) + uz(k,i+1,j)) / 2.
            obst_vz(k,i,j) = (vz(k,i,j) + uz(k,i,j+1)) / 2.
         ENDDO
       ENDIF ! * END TEST ON HWAT
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !---------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_uzvz

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_abdelposture(limin,limax,ljmin,ljmax,ssh,h0,h0fond)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_abdelposture  ***
   !&E
   !&E ** Purpose : Computes obstruction posture (height, diameters and bending angle)
   !&E
   !&E ** Description : Based on the balance between forces of drag, lift, friction,
   !&E                  weight and buoyancy of single obstruction element discretized
   !&E                  along the vertical, the Abdelrhman's model defined the bending
   !&E                  of obstruction elements
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : ssh, u, v, h0, obst_uz, obst_vz
   !&E
   !&E ** Modified variables : obst_width3d, obst_thick3d,obst_theta3d,obst_height,
   !&E                         obst_kmin,obst_kmax
   !&E
   !&E ** Reference : Abdelrhman M.A. (2007) Modeling coupling between eelgrass
   !&E                Zostera marina and water flow. Mar. Ecol. Prog. Ser. 338:81-96
   !&E
   !&E ** History :
   !&E       ! 2014-09-12 (F. Ganthy) Original code
   !&E       ! 2016-03-11 (F. Ganthy) Small fixe
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2016-09-20 (F. Ganthy) Completely re-coded with lot of optimizations and stability criteria
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       ! 2017-04-05 (F. Ganthy) Modification for case of 2D model
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-12-20 (F. Ganthy) Add test for small depth (when water depth LE obst_height)
   !&E       ! 2019-03-21 (F. Ganthy remove error on theta research (remplaced by warning) + small bug correction
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE
   !! * Arguments
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

   !! * Local declaration
   LOGICAL                                     :: IERR_MPI,l_theta
   INTEGER                                     :: i,j,k,iv,niter,niter_eff,s
   INTEGER                                     :: k0,k1
   REAL(KIND=rsh)                              :: z0,z1,zc,zt0,zt1,dz,udz,uv,htmp
   REAL(KIND=rsh)                              :: hwat,lseg,cd,cl,cf,dtheti,phi,re,cshelt
   REAL(KIND=rsh)                              :: th0,th0p,th0m,th1,th1p,th1m,sm0,sm0p,sm0m,sm1,sm1p,sm1m
   REAL(KIND=rsh)                              :: drag_xz_0,drag_xz_1m,drag_xz_1p
   REAL(KIND=rsh)                              :: lift_xz_0,lift_xz_1m,lift_xz_1p
   REAL(KIND=rsh)                              :: buoy_xz_0,buoy_xz_1m,buoy_xz_1p
   REAL(KIND=rsh)                              :: drag_x,fric_x,lift_z,buoy_z,fric_z
   REAL(KIND=rsh),DIMENSION(0:obst_kmax+1)     :: zz,uvz
   REAL(KIND=rsh),PARAMETER                    :: g=9.81_rsh,rw=1025.0_rsh,nu=1.4E-6_rsh,     &
                                                  gamma=1.0e-6_rsh
   INTEGER,PARAMETER                           :: niter_max   = 25                    ! Maximum number of iterations
   REAL(KIND=rsh),PARAMETER                    :: dtheta_max  = 0.25_rsh              ! Angle variation (degree) to reach stability
   REAL(KIND=rsh),PARAMETER                    :: dtheta_iter = 10.0_rsh*pi/180.0_rsh ! Searching step 
   REAL(KIND=rsh),PARAMETER                    :: phi_max     = 10.0_rsh*pi/180.0_rsh ! Max angle for lift
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   ! FG_TO_DO : will see if possible to use this procedure also for downward-flexible obstructions
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         IF((l_obst_flexible(iv)).AND.(l_obst_abdelposture(iv)))THEN
           IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
             hwat = h0(i,j) + ssh(i,j)
             IF(hwat.GT.h0fond)THEN
               IF(hwat.LE.obst_height_inst(iv,i,j))THEN
                 obst_height(iv,i,j)    = hwat
                 obst_theta3d(iv,:,i,j) = ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j))
               ELSE
                 !--------------------------
                 ! **** Initializations ****
                 !--------------------------
                 lseg  = obst_height_inst(iv,i,j)/REAL(obst_c_abdel_nmax(iv),rsh)
                 obst_abdel_zcent(iv,:)   = 0.0_rsh
                 obst_abdel_t0cent(iv,:)  = 0.0_rsh
                 obst_abdel_t1cent(iv,:)  = 0.0_rsh
                 obst_abdel_dtheta(iv,:)  = 0.0_rsh
                 obst_abdel_uvcent(iv,:)  = 0.0_rsh
                 obst_abdel_zn(iv,:)      = 0.0_rsh
                 obst_abdel_tn(iv,:)      = 0.0_rsh
                 !---------------------------------------
                 ! **** Computes velocity k=0:kmax+1 ****
                 !---------------------------------------
                 ! For k=0
                 k=0
                 zz(k)  = 0.0_rsh
                 uvz(k) = 0.0_rsh
                 ! For k=1,obst_kmax
                 DO k=1,obst_kmax
                   uvz(:) = SQRT(obst_uz(k,i,j)**2.0_rsh + obst_vz(k,i,j)**2.0_rsh)
                   zz(k)  = (1.0_rsh+obst_sig(k))*hwat + obst_dz(k,i,j)*0.5_rsh
                 ENDDO
                 ! For k=kmax+1
                 k=kmax+1
                 zz(k)  = hwat
                 uvz(k) = uvz(k-1)
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! **** STARTING ITERATIVE PROCEDURE OVER ALL SEGMENTS **** !
                 ! ****      TO REACH A STABLE OBSTRUCTION STATE       **** !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 niter = 0
                 DO WHILE (niter.LT.niter_max)
                   !-------------------------
                   ! Some re-initializations
                   !------------------------
                   obst_abdel_fx(iv,:)      = 0.0_rsh
                   obst_abdel_fz(iv,:)      = 0.0_rsh
                   !---------------------------------------------------------------------------!
                   ! *** COMPUTE HEIGHT AT EACH HALF-SEGMENT THROUGOUT AN UPWARD PROCEDURE *** !
                   !---------------------------------------------------------------------------!
                   DO s=1,obst_c_abdel_nmax(iv)
                     ! Height at half segment
                     IF (s.EQ.1) THEN ! First segment
                       obst_abdel_zcent(iv,s) = (lseg/2.0_rsh)*COS(obst_abdel_t0cent(iv,s))
                     ELSE             ! Other segments
                       obst_abdel_zcent(iv,s) = (obst_abdel_zcent(iv,s-1) + (lseg/2.0)*COS(obst_abdel_t0cent(iv,s-1))) + &
                                                (lseg/2.0_rsh)*COS(obst_abdel_t0cent(iv,s))
                     ENDIF
                     ! Velocity at height segment
                     k = 0
                     DO WHILE ((k.LT.obst_kmax) .AND. (zz(k).LE.obst_abdel_zcent(iv,s)))
                       k = k+1
                     ENDDO
                     k0 = MAX(MIN(k-1,obst_kmax-1),0)
                     k1 = MIN(MAX(k,1),obst_kmax)
                     obst_abdel_uvcent(iv,s) = uvz(k0) + (obst_abdel_zcent(iv,s) - zz(k0)) * &
                                                         (uvz(k1)-uvz(k0))/(zz(k1)-zz(k0))
                   ENDDO
                   !----------------------------------------------------------!
                   ! *** GENERAL DOWNWARD PROCEDURE FOR EACH LEAF SEGMENT *** !
                   !----------------------------------------------------------!
                   DO s=obst_c_abdel_nmax(iv),1,-1 ! Downward loop on each leaf segment
                     !-------------------------------------------------------------!
                     ! * Starting iterative procedure to search theta by solving * !
                     ! *  iteratively Eq. 7 or Eq. 10 of Abdelrhman 2007 with    * !
                     ! * with values starting at theta=0 and ending at pi/2 with * !
                     ! * its symetry at -pi/2, respectively for th0m,th0p,th1m,  * !
                     ! *              th1p and sm0p,sm0m,sm1p,sm1m               * !
                     !-------------------------------------------------------------!
                     th0     = 0.0_rsh ! Starting angle = 0 (vertical segment)
                     l_theta = .TRUE.  ! Condition to stop iterative research of theta value
                     !------------------------------------------------------------!
                     ! Computing xz contributions and summing all moments for th0 !
                     !------------------------------------------------------------!
                     phi = ABS((pi/2.0_rsh) - th0)
                     cshelt = MIN(MAX(1.0_rsh,obst_c_shelter(iv)*obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_thick_inst(iv,i,j)*lseg/(lseg*COS(th0)) - 1.0_rsh),4.0_rsh)
                     ! For drag
                     cd        = (phi*obst_c_drag(iv)/(pi/2.0_rsh)) / cshelt
                     drag_xz_0 = 0.5_rsh * cd * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                 obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (COS(th0)**2.0_rsh)
                     ! For lift
                     IF (phi.LE.phi_max) THEN
                       cl      = (obst_c_lift(iv)*(phi/phi_max)) / cshelt
                     ELSE
                       cl      = 0.0_rsh
                     ENDIF
                     lift_xz_0 = 0.5_rsh * cl * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                 obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (SIN(th0)**2.0_rsh)
                     ! For buoyancy
                     buoy_xz_0 = (rw - obst_c_rho(iv)) * g * obst_width_inst(iv,i,j) *        &
                                 obst_thick_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * SIN(th0)
                     ! Now summing all moments
                     IF (s.EQ.obst_c_abdel_nmax(iv)) THEN
                       sm0p = drag_xz_0 + lift_xz_0 - buoy_xz_0
                     ELSE
                       sm0p = drag_xz_0 + lift_xz_0 - buoy_xz_0 +                         &
                              obst_abdel_fx(iv,s+1)*lseg*COS(th0) -                       &
                              obst_abdel_fz(iv,s+1)*lseg*SIN(th0)
                     ENDIF
                     sm0m = sm0p
                     th0m = th0
                     th0p = th0
                     IF (obst_abdel_uvcent(iv,s).LE.0.001) THEN
                       !------------------------------!
                       ! Not enough velocity, theta=0 !
                       !------------------------------!
                       obst_abdel_t1cent(iv,s) = 0.0_rsh
                       obst_abdel_fx(iv,s)     = 0.0_rsh
                       obst_abdel_fz(iv,s)     = 0.0_rsh
                     ELSE
                       !----------------------------------------!
                       ! Effective start of iterative procedure !
                       !----------------------------------------!
                       DO WHILE (l_theta)
                         ! Incrementation of new theta values
                         th1m = th0m - dtheta_iter
                         th1p = th0p + dtheta_iter
                         !-------------------------------------------------------------!
                         ! Computing xz contributions and summing all moments for th1m !
                         !-------------------------------------------------------------!
                         phi = ABS((pi/2.0_rsh) - th1m)
                         cshelt = MIN(MAX(1.0_rsh,obst_c_shelter(iv)*obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_thick_inst(iv,i,j)*lseg/(lseg*COS(th1m)) - 1.0_rsh),4.0_rsh)
                         ! For drag
                         cd         = (phi*obst_c_drag(iv)/(pi/2.0_rsh)) / cshelt
                         drag_xz_1m = 0.5_rsh * cd * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                      obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (COS(th1m)**2.0_rsh)
                         ! For lift
                         IF (phi.LE.phi_max) THEN
                           cl       = (obst_c_lift(iv)*(phi/phi_max)) / cshelt
                         ELSE
                           cl       = 0.0_rsh
                         ENDIF
                         lift_xz_1m = 0.5_rsh * cl * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                      obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (SIN(th1m)**2.0_rsh)
                         ! For buoyancy
                         buoy_xz_1m = (rw - obst_c_rho(iv)) * g * obst_width_inst(iv,i,j) *    &
                                      obst_thick_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * SIN(th1m)
                         ! Now summing all moments
                         IF (s.EQ.obst_c_abdel_nmax(iv)) THEN
                           sm1m = drag_xz_1m + lift_xz_1m - buoy_xz_1m
                         ELSE
                           sm1m = drag_xz_1m + lift_xz_1m - buoy_xz_1m +                       &
                                  obst_abdel_fx(iv,s+1)*lseg*COS(th1m) -                       &
                                  obst_abdel_fz(iv,s+1)*lseg*SIN(th1m)
                         ENDIF
                         !-------------------------------------------------------------!
                         ! Computing xz contributions and summing all moments for th1p !
                         !-------------------------------------------------------------!
                         phi = ABS((pi/2.0_rsh) - th1p)
                         cshelt = MIN(MAX(1.0_rsh,obst_c_shelter(iv)*obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_thick_inst(iv,i,j)*lseg/(lseg*COS(th1p)) - 1.0_rsh),4.0_rsh)
                         ! For drag
                         cd         = (phi*obst_c_drag(iv)/(pi/2.0_rsh)) / cshelt
                         drag_xz_1p = 0.5_rsh * cd * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                      obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (COS(th1p)**2.0_rsh)
                         ! For lift
                         IF (phi.LE.phi_max) THEN
                           cl       = (obst_c_lift(iv)*(phi/phi_max)) / cshelt
                         ELSE
                           cl       = 0.0_rsh
                         ENDIF
                         lift_xz_1p = 0.5_rsh * cl * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                      obst_width_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * (SIN(th1p)**2.0_rsh)
                         ! For buoyancy
                         buoy_xz_1p = (rw - obst_c_rho(iv)) * g * obst_width_inst(iv,i,j) *    &
                                      obst_thick_inst(iv,i,j) * ((lseg**2.0_rsh)/2.0_rsh) * SIN(th1p)
                         ! Now summing all moments
                         IF (s.EQ.obst_c_abdel_nmax(iv)) THEN
                           sm1p = drag_xz_1p + lift_xz_1p - buoy_xz_1p
                         ELSE
                           sm1p = drag_xz_1p + lift_xz_1p - buoy_xz_1p +                         &
                                  obst_abdel_fx(iv,s+1)*lseg*COS(th1p) -                         &
                                  obst_abdel_fz(iv,s+1)*lseg*SIN(th1p)
                         ENDIF
                         !------------------------------------------!
                         ! Now testing if there sign change between !
                         ! sm0p and sm1p or between sm0m and sm1m   !
                         !------------------------------------------!
                         IF (((sm0m.GE.0.0_rsh) .AND. (sm1m.LE.0.0_rsh)) .OR. ((sm0m.LE.0.0_rsh) .AND. (sm1m.GE.0.0_rsh))) THEN
                           l_theta = .FALSE.
                           th0     = ABS(th0m)
                           th1     = ABS(th1m)
                           sm0     = sm0m
                           sm1     = sm1m
                         ELSEIF (((sm0p.GE.0.0_rsh) .AND. (sm1p.LE.0.0_rsh)) .OR. ((sm0p.LE.0.0_rsh) .AND. (sm1p.GE.0.0_rsh))) THEN
                           l_theta = .FALSE.
                           th0     = ABS(th0p)
                           th1     = ABS(th1p)
                           sm0     = sm0p
                           sm1     = sm1p
                         ELSEIF (th1p.GT.(pi/2.0_rsh)+dtheta_iter) THEN
                           IF(obst_abdel_uvcent(iv,s).LE.0.1_rsh) THEN
                               th0 = -(pi/4.0_rsh)
                               th1 =  (pi/4.0_rsh)
                               sm0 = -1.0_rsh
                               sm1 =  1.0_rsh
                           ELSE
                               th0 = (pi/2.0_rsh)-(pi/4.0_rsh)
                               th1 = (pi/2.0_rsh)+(pi/4.0_rsh)
                               sm0 = -1.0_rsh
                               sm1 =  1.0_rsh
                           ENDIF
                           l_theta = .FALSE.
                           WRITE(iwarnlog,*) ' '
                           WRITE(iwarnlog,*) ' '
                           WRITE(iwarnlog,*) '**************************************************************************'
                           WRITE(iwarnlog,*) '***** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_COMP_ABDELPOSTURE *****'
                           WRITE(iwarnlog,*) ' WARNING : no solution was found for the sum of moment'
                           WRITE(iwarnlog,*) '           At iv',iv,'i',i,'j',j,'s',s
                           WRITE(iwarnlog,*) '           hwat',hwat,'obst_height_inst',obst_height_inst(iv,i,j)
                           WRITE(iwarnlog,*) '           uzvz',obst_abdel_uvcent(iv,s)
                           WRITE(iwarnlog,*) ' --> Th=0 or Th=pi/2 applied depending on local velocity !!! '
                           WRITE(iwarnlog,*) '**************************************************************************'
                         ELSE
                           th0m = th1m
                           th0p = th1p
                           sm0m = sm1m
                           sm0p = sm1p
                         ENDIF
                       ENDDO ! End of iterative loop for theta value
                       !-------------------------------------------!
                       ! Linear interpolation of final theta value !
                       !-------------------------------------------!
                       obst_abdel_t1cent(iv,s) = th0 - sm0 * (th1-th0)/(sm1-sm0)
                       !-----------------------------------------------------!
                       ! Apply this good theta value to compute forces which !
                       !     will acts to the next obstruction segment       !
                       !-----------------------------------------------------!
                       IF (obst_c_abdel_nmax(iv).GT.1) THEN
                         phi = ABS((pi/2.0_rsh) - obst_abdel_t1cent(iv,s))
                         cshelt = MIN(MAX(1.0_rsh,obst_c_shelter(iv)*obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_thick_inst(iv,i,j)*lseg/(lseg*COS(obst_abdel_t1cent(iv,s))) - 1.0_rsh),4.0_rsh)
                         ! Computing drag
                         cd     = (phi*obst_c_drag(iv)/(pi/2.0_rsh)) / cshelt
                         drag_x = 0.5_rsh * cd * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                  obst_width_inst(iv,i,j) * lseg * COS(obst_abdel_t1cent(iv,s))
                         ! Computing friction
                         re     = (obst_abdel_uvcent(iv,s)*obst_width_inst(iv,i,j))/nu ! reynolds number
                         cf     = (0.074_rsh*(re**(-1.0_rsh/5.0_rsh))) / cshelt
                         fric_x = cf * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) *           &
                                  obst_width_inst(iv,i,j) * lseg *                         &
                                  (SIN(obst_abdel_t1cent(iv,s))**3.0_rsh)
                         fric_z = cf * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) *           &
                                  obst_width_inst(iv,i,j) * lseg *                         &
                                  (SIN(obst_abdel_t1cent(iv,s))**2.0_rsh) * COS(obst_abdel_t1cent(iv,s))
                         ! Computing lift
                         IF (phi.LE.phi_max) THEN
                           cl   = (obst_c_lift(iv)*(phi/phi_max)) / cshelt
                         ELSE
                           cl   = 0.0_rsh
                         ENDIF
                         lift_z = 0.5 * cl * rw * (obst_abdel_uvcent(iv,s)**2.0_rsh) * &
                                  obst_width_inst(iv,i,j) * lseg * SIN(obst_abdel_t1cent(iv,s))
                         ! Computing buoyancy
                         buoy_z = (rw - obst_c_rho(iv)) * g * obst_width_inst(iv,i,j) * &
                                   obst_thick_inst(iv,i,j) * lseg
                         ! Now summing all forces in x and z directions
                         IF (s.EQ.obst_c_abdel_nmax(iv)) THEN
                           obst_abdel_fx(iv,s) = -(  drag_x + fric_x) 
                           obst_abdel_fz(iv,s) = -( -lift_z + fric_z + buoy_z)
                         ELSE
                           obst_abdel_fx(iv,s) = -(  drag_x + fric_x + obst_abdel_fx(iv,s+1))
                           obst_abdel_fz(iv,s) = -( -lift_z + fric_z + buoy_z + obst_abdel_fz(iv,s+1))
                         ENDIF
                       ENDIF
                     ENDIF ! End test on enough velocity
                   ENDDO ! End of downward loop on each leaf segment
                   !---------------------------------------------------!
                   ! *** TEST FOR CONVERGENCE OF OBSTRUCTION STATE *** !
                   !---------------------------------------------------!
                   ! Computing difference in bending angle between current iteration and previous one
                   DO s=1,obst_c_abdel_nmax(iv)
                     obst_abdel_dtheta(iv,s) = ABS(obst_abdel_t0cent(iv,s)-obst_abdel_t1cent(iv,s))*180.0_rsh/pi
                   ENDDO
                   ! Stop iterative procedure if maximum bendig difference is less than dtheta_max
                   IF (MAXVAL(obst_abdel_dtheta(iv,:)).LE.dtheta_max) THEN
                     niter_eff = niter
                     niter     = niter_max
                   ELSE
                     niter_eff = niter
                     niter     = niter + 1
                   ENDIF
                   !-------------------------------------------------------------------------!
                   ! *** APPLY NEW THETA VALUE TO OLD ONE FOR NEXT ITERATION OR AVERAGED *** !
                   ! ***     THETA BETWEEN LAST ITERATION AND PREVIOUS ONE IF STABLE     *** !
                   ! ***   IF STABLE OBSTRUCTION STATE WAS NOT ACHIEVED (NO CONVERGENCE) *** !
                   !-------------------------------------------------------------------------!
                   IF (niter_eff.EQ.niter_max-1) THEN
                     obst_abdel_t0cent(iv,:) = 0.5_rsh*(obst_abdel_t0cent(iv,:)+obst_abdel_t1cent(iv,:))
                   ELSE
                     obst_abdel_t0cent(iv,:) = obst_abdel_t1cent(iv,:)
                   ENDIF
                 ENDDO ! End of iterative loop for osbtruction state
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! ****                   FINALIZATION                 **** !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !------------------------------------------------------------------!
                 ! *** WRITING INTO WARNING LOG IF CONVERGENCE WAS NOT ACHIEVED *** !
                 !------------------------------------------------------------------!
                 IF (niter_eff.EQ.niter_max-1) THEN
                   WRITE(iwarnlog,*) ' '
                   WRITE(iwarnlog,*) ' '
                   WRITE(iwarnlog,*) '************************************************************************'
                   WRITE(iwarnlog,*) '**** module OBSTRUCTIONS, subroutine OBSTRUCTIONS_COMP_ABDELPOSTURE ****'
                   WRITE(iwarnlog,*) ' WARNING : convergence for flow-obstruction coupling was not reached'
                   WRITE(iwarnlog,*) '           At iv',iv,'i',i,'j',j
                   WRITE(iwarnlog,*) '           Niter',niter_eff,'Niter_max',niter_max
                   WRITE(iwarnlog,*) '           Dtheta',MAXVAL(obst_abdel_dtheta(iv,:)),'Dtheta_max',dtheta_max
                   WRITE(iwarnlog,*) '           This may significantly alter hydrodynamic results and'
                   WRITE(iwarnlog,*) '           produce canopy oscillations...'
                   WRITE(iwarnlog,*) '           To prevent this, bending angles obtain at the last'
                   WRITE(iwarnlog,*) '           iteration and the preivous one are averaged'
                   WRITE(iwarnlog,*) '************************************************************************'
                 ENDIF
                 !--------------------------------------!
                 ! *** Updating obstructions height *** !
                 !--------------------------------------!
                 obst_height(iv,i,j) = 0.0_rsh
                 DO s=1,obst_c_abdel_nmax(iv)
                   obst_height(iv,i,j) = obst_height(iv,i,j) + lseg*COS(obst_abdel_t0cent(iv,s))
                 ENDDO
                 obst_height(iv,i,j) = MIN(MAX(obst_height(iv,i,j),obst_p_hmin*obst_height_inst(iv,i,j)),0.99_rsh*obst_height_inst(iv,i,j))
                 !--------------------------------------------------------------!
                 ! **** INTERPOLATION OF REAL BENDING ANGLES ON SIGMA GRID **** !
                 !--------------------------------------------------------------!
                 obst_theta3d(iv,:,i,j) = 0.0_rsh
                 s=0
                 obst_abdel_zn(iv,s) = 0.0_rsh
                 obst_abdel_tn(iv,s) = 0.0_rsh
                 DO s=1,obst_c_abdel_nmax(iv)
                   obst_abdel_zn(iv,s) = obst_abdel_zcent(iv,s)
                   obst_abdel_tn(iv,s) = obst_abdel_t0cent(iv,s)
                 ENDDO
                 s=obst_c_abdel_nmax(iv)+1
                 obst_abdel_zn(iv,s) = obst_height(iv,i,j)
                 obst_abdel_tn(iv,s) = obst_abdel_t0cent(iv,s-1)
                 DO k = 1,obst_kmax
                   z0  = obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh
                   z1  = obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh
                   IF(z1.LE.obst_height(iv,i,j))THEN
                     DO s=0,obst_c_abdel_nmax(iv)-1
                       zt0 = obst_abdel_zn(iv,s)
                       zt1 = obst_abdel_zn(iv,s+1)
                       th0 = obst_abdel_tn(iv,s)
                       th1 = obst_abdel_tn(iv,s+1)
                       IF((obst_zc(k,i,j).GE.zt0).AND.(obst_zc(k,i,j).LE.zt1))THEN
                         obst_theta3d(iv,k,i,j) = th0 + (obst_zc(k,i,j)-zt0)*(th1-th0)/(zt1-zt0)
                       ENDIF
                     ENDDO ! * END LOO ON S
                   ELSEIF((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                      obst_theta3d(iv,k,i,j) = obst_abdel_t0cent(iv,obst_c_abdel_nmax(iv))
                   ENDIF ! * END TEST ON Z0 AND Z1
                 ENDDO ! * END LOOP ON K
                 ! Additional check
                 DO k = 2,obst_kmax-1
                   IF(obst_theta3d(iv,k,i,j).EQ.0.0_rsh)THEN
                     IF((obst_theta3d(iv,k-1,i,j).GT.0.0_rsh).AND.(obst_theta3d(iv,k+1,i,j).GT.0.0_rsh))THEN
                       obst_theta3d(iv,k,i,j) = 0.5*(obst_theta3d(iv,k-1,i,j)+obst_theta3d(iv,k+1,i,j))
                     ENDIF
                   ENDIF
                 ENDDO
               ENDIF ! * END TEST ON HWAT AND OBST_HEIGHT
             ELSE ! * Not enough water
               obst_height(iv,i,j)    = 0.99_rsh*obst_height_inst(iv,i,j)
               obst_theta3d(iv,:,i,j) = 0.0_rsh
             ENDIF ! * END TEST ON Hwat>H0fond
           ELSE ! * No obstruction within this cell
             obst_height(iv,i,j)    = 0.0_rsh
             obst_theta3d(iv,:,i,j) = 0.0_rsh
           ENDIF ! END TEST ON PRESENCE OF OBSTRUCTIONS
         ENDIF ! END TEST ON RIGID/FLEXIBLE OBSTRUCTIONS AND ON ABDELRHMAN PROCEDURE
       ENDDO ! END LOOP ON iv
     ENDDO ! END LOOP ON i
   ENDDO ! END LOOP ON j
   !--------------------------------------- 
   END SUBROUTINE OBSTRUCTIONS_comp_abdelposture

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_height(limin, limax, ljmin, ljmax, ssh, h0, h0fond)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_height  ***
   !&E
   !&E ** Purpose : Computes obstructions height depending on choosen 
   !&E              parameterization
   !&E
   !&E ** Description : For rigid obstructions, obst_height is not changed
   !&E                  only obst_kmin and obst_kmax are computed, if 
   !&E                  obst_height < dsigu(k) obst_kmin=obst_kmax
   !&E                  For flexible obstructions, partial depth averaged velocity
   !&E                  (from bottom or surface to 2.5*previous obst_height)
   !&E                  corresponding to unconfined canopy is computed and used
   !&E                  to compute obst_height. Then obst_kmin and obst_kmax are
   !&E                  computed similarily than for rigid obstructions.
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_height
   !&E                     ssh, u, v, h0, uz, vz        
   !&E
   !&E ** Modified variables : obst_height
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08-26 (F. Ganthy) Computation for rigid obstruction moved here,
   !&E       !                        Re-organisation and optimization of the routine
   !&E       !                        Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       ! 2017-04-05 (F. Ganthy) Modification for case of 2D model
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. GAnthy) Optimization for full 3D obstructions
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE
   !! * Arguments
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

   !! * Local declaration
   INTEGER                               :: i,j,k,iv
   REAL(KIND=rsh)                        :: htmp,dz,udz,us
   REAL(KIND=rsh)                        :: uv,hwat
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBST_COMP_HEIGHT'
   !-----------------------------------------
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         hwat = h0(i,j) + ssh(i,j)
         !***************************
         ! * TEST ON OBSTRUCTION KIND
         !***************************
         IF(l_obst_flexible(iv)) THEN ! * FLEXIBLE OBSTRUCTIONS
           IF(.NOT.l_obst_abdelposture(iv)) THEN
             IF(l_obst_param_height(iv)) THEN ! * USE EXPONENTIAL PARAMETERIZATION
               IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
                 IF(hwat.GT.h0fond)THEN
                   !-------------------------------------------------
                   ! * COMPUTATION OF PARTIAL DEPTH-AVERAGED VELOCITY
                   !-------------------------------------------------
                   IF(obst_height(iv,i,j).EQ.0.0_rsh)THEN
                     htmp = obst_c_paramhuv * obst_height_inst(iv,i,j) ! At nbouc=0
                   ELSE
                     htmp = obst_c_paramhuv * obst_height(iv,i,j)
                   ENDIF
                   dz  = 0.0_rsh
                   udz = 0.0_rsh
                   IF(l_obst_downward(iv))THEN ! * DOWNWARD
                     k=obst_kmax
                     DO WHILE((k.GE.1).AND.(hwat-(obst_zc(k,i,j)-obst_dz(k,i,j)/2.0_rsh).LT.htmp))
                       dz  =  dz + obst_dz(k,i,j)
                       udz = udz + SQRT(obst_uz(k,i,j)**2.0_rsh + obst_vz(k,i,j)**2.0_rsh)*obst_dz(k,i,j)
                       k=k-1
                     ENDDO
                   ELSE ! * UPWARD
                     dz  = 0.0_rsh
                     udz = 0.0_rsh
                     k = 1
                     DO WHILE((k.LE.obst_kmax).AND.(obst_zc(k,i,j)+obst_dz(k,i,j)/2.0_rsh.LT.htmp))
                       dz  =  dz + obst_dz(k,i,j)
                       udz = udz + SQRT(obst_uz(k,i,j)**2.0_rsh + obst_vz(k,i,j)**2.0_rsh)*obst_dz(k,i,j)
                       k=k+1
                     ENDDO
                   ENDIF ! END TEST ON UPWARD/DOWNWARD
                   uv = udz/dz
                   !-----------------------------------------
                   ! * COMPUTATION OF BENT OBSTRUCTION HEIGHT
                   !-----------------------------------------
                   htmp = obst_c_height_x0(iv) * obst_height_inst(iv,i,j)* EXP(obst_c_height_x1(iv)*uv)
                   obst_height(iv,i,j) = MIN(MAX(htmp,obst_p_hmin*obst_height_inst(iv,i,j)),0.99_rsh*obst_height_inst(iv,i,j))
                 ELSE ! * Not enough water depth
                   obst_height(iv,i,j) = 0.99_rsh*obst_height_inst(iv,i,j)
                 ENDIF! * END TEST ON WATER DEPTH
               ELSE ! * No obstruction within this cell
                 obst_height(iv,i,j) = 0.0_rsh
               ENDIF ! * END TEST ON OBSTRUCTION POSITION
             ELSE ! * USE ONLY FIRST COEFFICIENT
               IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
                 IF(hwat.GT.h0fond)THEN
                   !-------------------------------------------
                   ! * COMPUTATION OF UNBENT OBSTRUCTION HEIGHT
                   !-------------------------------------------
                   obst_height(iv,i,j) = obst_height_inst(iv,i,j)*obst_c_height_x0(iv)
                 ELSE ! * Not enough water depth
                   obst_height(iv,i,j) = 0.99_rsh* obst_height_inst(iv,i,j)
                 ENDIF! * END TEST ON WATER DEPTH
               ELSE ! * No obstruction within this cell
                 obst_height(iv,i,j) = 0.0_rsh
               ENDIF ! * END TEST ON OBSTRUCTION POSITION
             ENDIF ! * END TEST ON PARAMETERIZATION
           ENDIF ! * END TEST ON ABDELRHMAN PROCEDURE
         ELSE ! * RIDIG OBSTRUCTIONS
           !------------------------------------------
           ! * COMPUTATION OF RIGID OBSTRUCTION HEIGHT
           !------------------------------------------
           IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
             IF(hwat.GT.h0fond)THEN
               obst_height(iv,i,j) = obst_height_inst(iv,i,j)*obst_c_height_x0(iv)
             ELSE ! * Not enough water depth
               obst_height(iv,i,j) = 0.99_rsh*obst_height_inst(iv,i,j)
             ENDIF! * END TEST ON WATER DEPTH
           ELSE ! * No obstruction within this cell
             obst_height(iv,i,j) = 0.0_rsh
           ENDIF ! * END TEST ON OBSTRUCTION POSITION
         ENDIF ! * END TEST ON FLEXIBLE/RIGID OBSTRUCTION
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-----------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_height

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_theta(limin, limax, ljmin, ljmax, ssh, h0, h0fond)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_theta  ***
   !&E
   !&E ** Purpose : Computes obstructions height depending on choosen 
   !&E              parameterization
   !&E
   !&E ** Description : For rigid obstructions, obst_height is not changed
   !&E                  only obst_kmin and obst_kmax are computed, if 
   !&E                  obst_height < dsigu(k) obst_kmin=obst_kmax
   !&E                  For flexible obstructions, partial depth averaged velocity
   !&E                  (from bottom or surface to 2.5*previous obst_height)
   !&E                  corresponding to unconfined canopy is computed and used
   !&E                  to compute obst_height. Then obst_kmin and obst_kmax are
   !&E                  computed similarily than for rigid obstructions.
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_height
   !&E                     ssh,h0      
   !&E
   !&E ** Modified variables :  obst_theta3d
   !&E
   !&E ** Reference :
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08-26 (F. Ganthy) Computation for rigid obstruction moved here,
   !&E       !                        Re-organisation and optimization of the routine
   !&E       !                        Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. Ganthy) Optimization for full 3D obstructions
   !&E       ! 2020-04-23 (F. Ganthy) Small bug correstion
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE
   !! * Arguments
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

   !! * Local declaration
   LOGICAL                          :: IERR_MPI
   INTEGER                          :: i,j,k,iv
   REAL(KIND=rsh)                   :: hwat,z0,z1
   !!----------------------------------------------------------------------
   !! * Executable part
   PRINT_DBG*, 'ENTER OBST_COMP_THETA'
   !-----------------------------------------
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         IF(l_obst_flexible(iv))THEN
           IF(.NOT.l_obst_abdelposture(iv))THEN
             IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
               !*********************************
               ! *** COMPUTATION OF BENDING ANGLE
               !*********************************
               IF(l_obst_downward(iv))THEN ! * DOWNWARD
                 k=obst_kmax
                 DO WHILE ((k.GE.1).AND.(hwat-(obst_zc(k,i,j)+obst_dz(k,i,j)/2.0_rsh).LT.obst_height(iv,i,j)))
                   obst_theta3d(iv,k,i,j) = ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j))
                   k=k-1
                 ENDDO
               ELSE ! * UPWARD
                 k=1
                 !DO WHILE ((k.LE.obst_kmax).AND.(obst_zc(k,i,j)+obst_dz(k,i,j)/2.0_rsh.LT.obst_height(iv,i,j)))
                 DO WHILE ((k.LE.obst_kmax).AND.(obst_zc(k,i,j)-obst_dz(k,i,j)/2.0_rsh.GE.obst_height(iv,i,j)))
                   obst_theta3d(iv,k,i,j) = ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j))
                   k=k+1
                 ENDDO
               ENDIF ! * END TEST ON UPWARD/DOWNWARD
             ELSE ! * NO OBSTRUCTION WITHIN THIS CELL
               obst_theta3d(iv,:,i,j) = 0.0_rsh
             ENDIF ! * END TEST ON OBSTRUCTION POSITION
           ENDIF ! * END TEST ON ABDELHRMAN PROCEDURE
         ELSE ! * FOR RIGID OBSTRUCTION
           obst_theta3d(iv,:,i,j) = 0.0_rsh
         ENDIF ! * END TEST ON FLEXIBLE/RIGID OBSTRUCTION
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-----------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_theta

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_fracz(limin, limax, ljmin, ljmax, ssh, h0)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_fracz ***
   !&E
   !&E ** Purpose : Computes obstructions height depending on choosen 
   !&E              parameterization
   !&E
   !&E ** Description : For rigid obstructions, obst_height is not changed
   !&E                  only obst_kmin and obst_kmax are computed, if 
   !&E                  obst_height < dsigu(k) obst_kmin=obst_kmax
   !&E                  For flexible obstructions, partial depth averaged velocity
   !&E                  (from bottom or surface to 2.5*previous obst_height)
   !&E                  corresponding to unconfined canopy is computed and used
   !&E                  to compute obst_height. Then obst_kmin and obst_kmax are
   !&E                  computed similarily than for rigid obstructions.
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_height
   !&E                     ssh,h0        
   !&E
   !&E ** Modified variables : obst_fracz
   !&E
   !&E ** Reference :
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-08-26 (F. Ganthy) Computation for rigid obstruction moved here,
   !&E       !                        Re-organisation and optimization of the routine
   !&E       !                        Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. Ganthy) Optimization for full 3D obstructions
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE
   !! * Arguments
   INTEGER, INTENT(IN) :: limin, limax, ljmin, ljmax

   !! * Local declaration
   LOGICAL                          :: IERR_MPI
   INTEGER                          :: i,j,k,iv
   REAL(KIND=rsh)                   :: hwat,z0,z1
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   obst_fracz3d(:,:,:,:) = 0.0_rsh
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         hwat = h0(i,j) + ssh(i,j)
         IF((hwat.GT.h0fond).AND.(obst_position(iv,i,j).GT.0.0_rsh))THEN
           IF(l_obst_downward(iv))THEN ! * DOWNWARD
             DO k=obst_kmax,1,-1
               z0  = hwat - (obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh)
               z1  = hwat - (obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh)
               IF(z1.LE.obst_height(iv,i,j)) THEN
                 obst_fracz3d(iv,k,i,j) = 1.0_rsh
               ELSEIF ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                 obst_fracz3d(iv,k,i,j) = (obst_height(iv,i,j)-z0)/(z1-z0)
               ENDIF
             ENDDO
           ELSE ! * UPWARD
             IF(l_obst_3dobst(iv))THEN
               DO k=1,obst_kmax
                 z0   = obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh
                 z1   = obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh
                 IF(z1.LE.obst_height(iv,i,j)) THEN
                   IF(obst_dens3d(iv,k,i,j).EQ.0.0_rsh)THEN
                     obst_fracz3d(iv,k,i,j) = 0.0_rsh
                   ELSE
                     obst_fracz3d(iv,k,i,j) = 1.0_rsh
                   ENDIF
                 ELSEIF ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                   obst_fracz3d(iv,k,i,j) = (obst_height(iv,i,j)-z0)/(z1-z0)
                 ENDIF
               ENDDO
             ELSE               
               DO k=1,obst_kmax
                 z0   = obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh
                 z1   = obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh
                 IF(z1.LE.obst_height(iv,i,j)) THEN
                   obst_fracz3d(iv,k,i,j) = 1.0_rsh
                 ELSEIF ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                   obst_fracz3d(iv,k,i,j) = (obst_height(iv,i,j)-z0)/(z1-z0)
                 ENDIF
               ENDDO
             ENDIF ! * END TEST ON 3D OBSTRUCTION
           ENDIF ! * END TEST ON UPWARD/DOWNWARD
         ENDIF ! * END TEST ON HWAT AND PRESENCE OF OBSTRUCTIONS
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-----------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_fracz

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_diam(limin, limax, ljmin, ljmax, ssh, h0)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_diam  ***
   !&E
   !&E ** Purpose : Computes obstruction diameters (width and thickness)
   !&E
   !&E ** Description : Computes obstruction width and thickness according to
   !&E                  bending angle (for flexible obstructions)
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_kmin, obst_kmax, obst_width3d, obst_thick3d, obst_thick3dd
   !&E                     ssh, u, v, h0, uz, vz
   !&E
   !&E ** Modified variables : obst_width3d, obst_thick3d, obst_thick3dd
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-10-06 (F. Ganthy) Computation for rigid obstruction moved here,
   !&E       !                        Re-organisation and optimization of the routine
   !&E       !                        Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2016-03-11 (F. Ganthy) Change computation of obst_thick2d and obst_width2d
   !&E       ! 2016-03-11 (F. Ganthy) Minor bug correction
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       !                        and removing useless parameters and tests
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. Ganthy) Optimization for full 3D obstructions
   !&E       ! 2019-03-21 (F. Ganthy) Small modification of thickness computation to be more realistic
   !&E       ! 2020-04-23 (F. Ganthy) More modification of thickness computation to be more realistic
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Arguments
   !! * Local declaration
   LOGICAL                          :: IERR_MPI
   INTEGER                          :: i,j,k,iv
   REAL(KIND=rsh)                   :: hwat,lsegl,ltotl
 
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   obst_width3d(:,:,:,:) = 0.0_rsh
   obst_thick3d(:,:,:,:) = 0.0_rsh
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         hwat = h0(i,j) + ssh(i,j)
         IF((hwat.GT.h0fond).AND.(obst_position(iv,i,j).GT.0.0_rsh))THEN
           !***************************************
           ! *** COMPUTATION OF WIDTH AND THICKNESS
           !***************************************
           obst_width3d(iv,:,i,j) = 0.0_rsh
           obst_thick3d(iv,:,i,j) = 0.0_rsh
           IF(l_obst_downward(iv))THEN ! * DOWNWARD
             ltotl = 0.0_rsh
             DO k=obst_kmax,1,-1
               IF(obst_dens3d(iv,k,i,j).GT.0.0_rsh)THEN
                 obst_width3d(iv,k,i,j) = obst_width_inst(iv,i,j)
                 IF(l_obst_cylindre(iv))THEN ! * CYLINDER
                     IF(obst_theta3d(iv,k,i,j).EQ.0.0_rsh)THEN
                       obst_thick3d(iv,k,i,j) = obst_width_inst(iv,i,j)
                     ELSE
                       lsegl = MIN(MAX(obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)/COS(obst_theta3d(iv,k,i,j)),obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)),obst_height_inst(iv,i,j)-ltotl)
                       obst_thick3d(iv,k,i,j) = MIN(lsegl,MAX(obst_width_inst(iv,i,j),lsegl*SIN(obst_theta3d(iv,k,i,j))))
                       ltotl = ltotl+lsegl
                     ENDIF
                 ELSE ! * PARALELEPIPED
                     IF(obst_theta3d(iv,k,i,j).EQ.0.0_rsh)THEN
                       obst_thick3d(iv,k,i,j) = obst_thick_inst(iv,i,j)
                     ELSE
                       lsegl = MIN(MAX(obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)/COS(obst_theta3d(iv,k,i,j)),obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)),obst_height_inst(iv,i,j)-ltotl)
                       obst_thick3d(iv,k,i,j) = MIN(lsegl,MAX(obst_thick_inst(iv,i,j),lsegl*SIN(obst_theta3d(iv,k,i,j))))
                       ltotl = ltotl+lsegl
                     ENDIF
                 ENDIF ! * END TEST CYLINDER
               ENDIF
             ENDDO
           ELSE ! * UPWARD
             ltotl = 0.0_rsh
             DO k=1,obst_kmax
               IF(obst_dens3d(iv,k,i,j).GT.0.0_rsh)THEN
                 obst_width3d(iv,k,i,j) = obst_width_inst(iv,i,j)
                 IF(l_obst_cylindre(iv))THEN ! * CYLINDER
                   IF(obst_theta3d(iv,k,i,j).EQ.0.0_rsh)THEN
                     obst_thick3d(iv,k,i,j) = obst_width_inst(iv,i,j)
                   ELSE
                     lsegl = MIN(MAX(obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)/COS(obst_theta3d(iv,k,i,j)),obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)),obst_height_inst(iv,i,j)-ltotl)
                     obst_thick3d(iv,k,i,j) = MIN(lsegl,MAX(obst_width_inst(iv,i,j),lsegl*SIN(obst_theta3d(iv,k,i,j))))
                     ltotl = ltotl+lsegl
                   ENDIF
                 ELSE ! * PARALELEPIPED
                   IF(obst_theta3d(iv,k,i,j).EQ.0.0_rsh)THEN
                     obst_thick3d(iv,k,i,j) = obst_thick_inst(iv,i,j)
                   ELSE
                     lsegl = MIN(MAX(obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)/COS(obst_theta3d(iv,k,i,j)),obst_fracz3d(iv,k,i,j)*obst_dz(k,i,j)),obst_height_inst(iv,i,j)-ltotl)
                     obst_thick3d(iv,k,i,j) = MIN(lsegl,MAX(obst_thick_inst(iv,i,j),lsegl*SIN(obst_theta3d(iv,k,i,j))))
                     ltotl = ltotl+lsegl
                   ENDIF
                 ENDIF ! * END TEST CYLINDER
               ENDIF
             ENDDO
           ENDIF ! * END TEST ON UPWARD/DOWNWARD
         ENDIF ! * END TEST ON HWAT AND PRESENCE OF OBSTRUCTIONS
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-----------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_diam

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_distrib(limin,limax,ljmin,ljmax,h0,ssh)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_distrib  ***
   !&E
   !&E ** Purpose : Computes vertical distribution of obstructions densities
   !&E
   !&E ** Description : In all case, obst_dens3d(obst_kmin) = obst_dens_inst (100%)
   !&E                  if obstructions are flexible
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_position, obst_dens3d, obst_kmin,
   !&E                     obst_kmax, ssh, h0
   !&E
   !&E ** Modified variables : obst_dens3d
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-10-06 (F. Ganthy) Computation for rigid obstruction moved here,
   !&E       !                        Re-organisation and optimization of the routine
   !&E       !                        Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2015-09-10 (F. Ganthy) Minor corrections
   !&E       ! 2016-03-11 (F. Ganthy) Minor bug correction
   !&E       ! 2017-02-16 (F. Ganthy) Modification for multispecific obstructions within single grid cell
   !&E       !                        and Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. Ganthy) Optimization for full 3D obstructions
   !&E       ! 2018-03-09 (F. Ganthy) Small bug correction for variable distribution when obst_height is very small
   !&E       !                        compared with layer thickness
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Arguments
   !! * Local declaration
   LOGICAL                                  :: IERR_MPI
   INTEGER                                  :: iv,i,j,k,kk
   REAL(KIND=rsh)                           :: hwat,dhnn,dnn,shtot
   REAL(KIND=rsh)                           :: z0,z1,zc,zn0,zn1,dn0,dn1
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   obst_dens3d(:,:,:,:) = 0.0_rsh
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         hwat = h0(i,j) + ssh(i,j)
         IF((hwat.GT.h0fond).AND.(obst_position(iv,i,j).GT.0.0_rsh))THEN
           !****************************************
           ! *** COMPUTATION OF DENSITY DISTRIBUTION
           !****************************************
           IF(l_obst_filedistri(iv)) THEN ! * VARIABLE DISTRIBUTION
             !----------------------------------------
             ! *** COMPUTATION OF DENSITY DISTRIBUTION
             !----------------------------------------
             obst_dens3d(iv,:,i,j) = 0.0_rsh
             IF(l_obst_downward(iv))THEN ! * DOWNWARD
               DO k=obst_kmax,1,-1
                 z0  = MAX(hwat - (obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh) , 0.0_rsh)
                 z1  = MIN(hwat - (obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh) , hwat)
                 IF((z1.LE.obst_height(iv,i,j)) .OR. ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j))))THEN
                   IF(z1.LE.obst_height(iv,i,j)) THEN
                     zc = hwat - obst_zc(k,i,j)
                   ELSEIF((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                     zc = z0 + 0.5_rsh*(obst_height(iv,i,j)-z0)
                   ENDIF
                   DO kk=1,obst_nbhnorm(iv)-1
                     !--------------------------------------------------
                     ! *** CONVERSION OF NORMALIZED HEIGHT AND DENSITY
                     !--------------------------------------------------
                     zn0 = MIN(MAX(obst_height_norm(iv,kk)   * obst_height(iv,i,j) / 100.0_rsh,0.0_rsh),obst_height(iv,i,j))
                     zn1 = MIN(MAX(obst_height_norm(iv,kk+1) * obst_height(iv,i,j) / 100.0_rsh,0.0_rsh),obst_height(iv,i,j))
                     dn0 = MIN(MAX(obst_dens_norm(iv,kk)   * obst_dens_inst(iv,i,j) / 100.0_rsh,0.0_rsh),obst_dens_inst(iv,i,j))
                     dn1 = MIN(MAX(obst_dens_norm(iv,kk+1) * obst_dens_inst(iv,i,j) / 100.0_rsh,0.0_rsh),obst_dens_inst(iv,i,j))
                     IF((zc.GE.zn0).AND.(zc.LE.zn1))THEN
                       obst_dens3d(iv,k,i,j) = dn0 + (zc-zn0)*(dn1-dn0)/(zn1-zn0)
                     ENDIF
                   ENDDO
                 ENDIF
               ENDDO
             ELSE ! * UPWARD
               DO k=1,obst_kmax
                 z0   = MAX(obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh , 0.0_rsh)
                 z1   = MIN(obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh , hwat)
                 IF((z1.LE.obst_height(iv,i,j)) .OR. ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j))))THEN
                   IF(z1.LE.obst_height(iv,i,j)) THEN
                     zc = obst_zc(k,i,j)
                   ELSEIF((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j)))THEN
                     zc = z0 + 0.5_rsh*(obst_height(iv,i,j)-z0)
                   ENDIF
                   DO kk=1,obst_nbhnorm(iv)-1
                     zn0 = MIN(MAX(obst_height_norm(iv,kk)   * obst_height(iv,i,j) / 100.0_rsh,0.0_rsh),obst_height(iv,i,j))
                     zn1 = MIN(MAX(obst_height_norm(iv,kk+1) * obst_height(iv,i,j) / 100.0_rsh,0.0_rsh),obst_height(iv,i,j))
                     dn0 = MIN(MAX(obst_dens_norm(iv,kk)   * obst_dens_inst(iv,i,j) / 100.0_rsh,0.0_rsh),obst_dens_inst(iv,i,j))
                     dn1 = MIN(MAX(obst_dens_norm(iv,kk+1) * obst_dens_inst(iv,i,j) / 100.0_rsh,0.0_rsh),obst_dens_inst(iv,i,j))
                     IF((zc.GE.zn0).AND.(zc.LE.zn1))THEN
                       obst_dens3d(iv,k,i,j) = dn0 + (zc-zn0)*(dn1-dn0)/(zn1-zn0)
                     ENDIF
                   ENDDO
                 ENDIF
               ENDDO
             ENDIF ! * END TEST ON UPWARD/DOWNWARD
           ELSE ! * CONSTANT DISTRIBUTION
             IF(l_obst_downward(iv))THEN ! * DOWNWARD
               DO k=obst_kmax,1,-1
                 z0  = MAX(hwat - (obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh) , 0.0_rsh)
                 z1  = MIN(hwat - (obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh) , hwat)
                 IF((z1.LE.obst_height(iv,i,j)) .OR. ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j))))THEN
                   obst_dens3d(iv,k,i,j) = obst_dens_inst(iv,i,j)
                 ENDIF
               ENDDO
             ELSE ! * UPWARD
               DO k=1,obst_kmax
                 z0   = MAX(obst_zc(k,i,j) - obst_dz(k,i,j)/2.0_rsh , 0.0_rsh)
                 z1   = MIN(obst_zc(k,i,j) + obst_dz(k,i,j)/2.0_rsh , hwat)
                 IF((z1.LE.obst_height(iv,i,j)) .OR. ((z0.LE.obst_height(iv,i,j)).AND.(z1.GT.obst_height(iv,i,j))))THEN
                   obst_dens3d(iv,k,i,j) = obst_dens_inst(iv,i,j)
                 ENDIF
               ENDDO
             ENDIF ! * END TEST ON UPWARD/DOWNWARD
           ENDIF ! * END TEST ON CONSTANT DISTRIBUTION
         ENDIF ! * END TEST ON HWAT AND PRESENCE OF OBSTRUCTIONS
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !----------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_distrib

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_fracxy(limin,limax,ljmin,ljmax,h0,ssh)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_fracxy  ***
   !&E
   !&E ** Purpose : Computes obstructions correction term for coverage (fragmentation)
   !&E              within one single grid cell
   !&E ** Description : 
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_height
   !&E                     ssh, u, v, h0, uz, vz        
   !&E
   !&E ** Modified variables : 
   !&E ** Reference : Laporte-Fauret and Ganthy, ISOBAY 2016
   !&E
   !&E ** History :
   !&E       ! 2012-03-20 (F. Ganthy) Original code
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-10 (F. GAnthy) Optimization for full 3D obstructions
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Arguments

   !! * Local declaration
   INTEGER            :: i,j,iv
   REAL(KIND=rsh)     :: hwat,kv
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   obst_fracxy(:,:,:) = 0.0_rsh
   DO j=ljmin,ljmax
     DO i=limin,limax
       DO iv=1,obst_nbvar
         hwat = h0(i,j) + ssh(i,j)
         IF((hwat.GT.h0fond).AND.(obst_position(iv,i,j).GT.0.0_rsh))THEN
           !--------------------------------------------------
           ! *** COMPUTATION OF CORRECTION FACTOR FOR COVERAGE
           !--------------------------------------------------
           IF (l_obst_fracxy(iv))THEN
             IF (obst_fracxy_type(iv).EQ.0)THEN
               !------------------
               ! Linear correction
               !------------------
               obst_fracxy(iv,i,j) = obst_position(iv,i,j)
             ELSEIF (obst_fracxy_type(iv).EQ.1)THEN
               !------------------------------
               ! Simple exponential correction
               !------------------------------
               obst_fracxy(iv,i,j) = (EXP(obst_position(iv,i,j)*obst_c_fracxy_k0(iv))-1.0_rsh) / &
                                     (EXP(obst_c_fracxy_k0(iv))-1.0_rsh)
             ELSEIF (obst_fracxy_type(iv).EQ.2)THEN
               !-------------------------------
               ! Complex exponential correction
               !-------------------------------
               kv = obst_c_fracxy_k0(iv) + obst_c_fracxy_k1(iv)*EXP(-(1.0_rsh-obst_position(iv,i,j))*obst_c_fracxy_l(iv))
               obst_fracxy(iv,i,j) = (EXP(obst_position(iv,i,j)*kv)-1.0_rsh) / &
                                     (EXP(kv)-1.0_rsh)
             ELSEIF (obst_fracxy_type(iv).EQ.3)THEN
               !---------------------------------------------------------
               ! Constant correction factor (for parameteriation purpose)
               !---------------------------------------------------------
               obst_fracxy(iv,i,j) = obst_c_fracxy_k0(iv) * obst_position(iv,i,j)
             ENDIF ! * END TEST ON CORRECTION TYPE
           ELSE
               !--------------
               ! No correction
               !--------------
               obst_fracxy(iv,i,j) = 1.0_rsh
           ENDIF ! * END TEST ON CORRECTION TYPE
         ENDIF ! END TEST ON HWAT AND PRESENCE OF OBSTRUCTIONS
       ENDDO ! * END LOOP ON iv
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !----------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_fracxy

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_projarea(limin,limax,ljmin,ljmax,h0,ssh)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_projarea  ***
   !&E
   !&E ** Purpose : Computes obstuctions projected horizontal and vertical area
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2017-04-12 (F. Ganthy) Original code
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-04 (F. Ganthy) Correction of obst_s3d computation and computation of obst_a2d and obst_s2d
   !&E                                for the three additional types (NoTurb, Turb and Tot) with respect to physics
   !&E                                behind
   !&E       ! 2019-03-21 (F. Ganthy) Small correction of computation of obst_a3d and obst_s3d to be more realistic
   !&E                                + change method for averaging (average over the obstruction height instead of
   !&E                                average over the water depth : more meaningfull for further application)
   !&E
   !&E---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Arguments
   !! * Local declaration
   INTEGER                  :: iv,i,j,k
   REAL(KIND=rsh)           :: hwat,dzas,spos,wtmp,dtmp,sfz
   REAL(KIND=rsh),PARAMETER :: gamma=1.0_rsh-1.0e-6_rsh
   !!----------------------------------------------------------------------
   !! * Executable part
   ! ********************************
   ! * RE-INITIALIZATION OF VARIABLES
   ! ********************************
   obst_a2d(:,:,:)       = 0.0_rsh
   obst_s2d(:,:,:)       = 0.0_rsh
   obst_a3d(:,:,:,:)     = 0.0_rsh
   obst_s3d(:,:,:,:)     = 0.0_rsh
   obst_dens_mean(:,:)   = 0.0_rsh
   obst_width_mean(:,:)  = 0.0_rsh
   obst_height_mean(:,:) = 0.0_rsh
   !*************************
   DO j=ljmin,ljmax
     DO i=limin,limax
       hwat = h0(i,j) + ssh(i,j)
       IF(hwat.GT.h0fond) THEN
         ! ******************** !
         ! *** 3D VARIABLES *** !
         ! ******************** !
         DO k=1,obst_kmax
           !----------------------
           ! * LOOP ON OBST_NB_VAR
           !----------------------
           DO iv=1,obst_nbvar
             IF((obst_position(iv,i,j).GT.0.0_rsh).AND.(obst_dens3d(iv,k,i,j).GT.0.0_rsh))THEN
               !------------------------------------------
               ! * Horizontal surface area of obstructions
               !------------------------------------------
               IF(l_obst_cylindre(iv))THEN ! Cylindric/Ellipse obstruction
                 obst_a3d(iv,k,i,j) = MIN(pi*(0.5_rsh*obst_width3d(iv,k,i,j)*0.5_rsh*obst_thick3d(iv,k,i,j))* &
                                          obst_dens3d(iv,k,i,j)*obst_fracxy(iv,i,j),gamma)
               ELSE
                 obst_a3d(iv,k,i,j) = MIN(obst_width3d(iv,k,i,j)*obst_thick3d(iv,k,i,j)*                      &
                                          obst_dens3d(iv,k,i,j)*obst_fracxy(iv,i,j),gamma)
               ENDIF
               !----------------------------------------
               ! * Vertical surface area of obstructions
               !----------------------------------------
               obst_s3d(iv,k,i,j) = MIN(obst_fracz3d(iv,k,i,j)*obst_width3d(iv,k,i,j)*        &
                                        SQRT(obst_dens3d(iv,k,i,j))*obst_fracxy(iv,i,j),gamma)
             ENDIF ! * END TEST ON POSITION AND DENSITY
           ENDDO ! * END LOOP ON OBST_NBVAR
           !------------------------------------------------------------------------
           ! * Computes obst_a3d and obst_s3d for NoTurb/Turb and Total obstructions
           !------------------------------------------------------------------------
           DO iv=1,obst_nbvar
             IF(l_obst_noturb(iv)) THEN
               !-------------------------------------------
               ! First for NoTurb (simplified obstructions)
               !-------------------------------------------
               obst_a3d(obst_nbvar+1,k,i,j) = MIN(obst_a3d(obst_nbvar+1,k,i,j) + obst_a3d(iv,k,i,j),gamma)
               obst_s3d(obst_nbvar+1,k,i,j) = MIN(obst_s3d(obst_nbvar+1,k,i,j) + obst_s3d(iv,k,i,j),gamma)
             ELSE
               !---------------------------------------
               ! Then for Turb (turbulent obstructions)
               !---------------------------------------
               obst_a3d(obst_nbvar+2,k,i,j) = MIN(obst_a3d(obst_nbvar+2,k,i,j) + obst_a3d(iv,k,i,j),gamma)
               obst_s3d(obst_nbvar+2,k,i,j) = MIN(obst_s3d(obst_nbvar+2,k,i,j) + obst_s3d(iv,k,i,j),gamma)
             ENDIF
             !--------------------------------
             ! Then for Tot (all obstructions)
             !--------------------------------
             obst_a3d(obst_nbvar+3,k,i,j) = MIN(obst_a3d(obst_nbvar+3,k,i,j) + obst_a3d(iv,k,i,j),gamma)
             obst_s3d(obst_nbvar+3,k,i,j) = MIN(obst_s3d(obst_nbvar+3,k,i,j) + obst_s3d(iv,k,i,j),gamma)
           ENDDO ! * END LOOP ON OBST_NBVAR
         ENDDO ! * END LOOP ON K
         !----------------------------------
         ! Now computes depth-averaged value
         !----------------------------------
         DO iv=1,obst_nbvar
           IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
             dzas = 0.0_rsh
             DO k=1,obst_kmax
                 !-----------
                 ! For obst_a
                 !-----------
                 obst_a2d(iv,i,j) = obst_a2d(iv,i,j) + (obst_a3d(iv,k,i,j)*obst_dz(k,i,j))
                 !-----------
                 ! For obst_s
                 !-----------
                 obst_s2d(iv,i,j) = obst_s2d(iv,i,j) + (obst_s3d(iv,k,i,j)*obst_dz(k,i,j))
                 !------------
                 dzas = dzas + obst_dz(k,i,j)
             ENDDO ! * END LOOP ON K
             IF(dzas.GT.0.0_rsh)THEN
               obst_a2d(iv,i,j)     = MIN(obst_a2d(iv,i,j)/dzas,gamma)
               obst_s2d(iv,i,j)     = MIN(obst_s2d(iv,i,j)/dzas,gamma)
             ENDIF
           ENDIF ! * END TEST ON POSITION
         ENDDO ! * END LOOP ON IV
         !------------------------------
         ! Then For NoTurb, Turb and Tot
         !------------------------------
         DO iv=1,obst_nbvar
           IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
             IF(l_obst_noturb(iv)) THEN
               !-----------
               ! For NoTurb
               !-----------
               obst_a2d(obst_nbvar+1,i,j) = obst_a2d(obst_nbvar+1,i,j)+obst_a2d(iv,i,j)
               obst_s2d(obst_nbvar+1,i,j) = obst_s2d(obst_nbvar+1,i,j)+obst_s2d(iv,i,j)
             ELSE
               !---------
               ! For Turb
               !---------
               obst_a2d(obst_nbvar+2,i,j) = obst_a2d(obst_nbvar+2,i,j)+obst_a2d(iv,i,j)
               obst_s2d(obst_nbvar+2,i,j) = obst_s2d(obst_nbvar+2,i,j)+obst_s2d(iv,i,j)
             ENDIF ! * END TEST ON NOTURB
             !----------
             ! For Total
             !----------
             obst_a2d(obst_nbvar+3,i,j) = obst_a2d(obst_nbvar+3,i,j)+obst_a2d(iv,i,j)
             obst_s2d(obst_nbvar+3,i,j) = obst_s2d(obst_nbvar+3,i,j)+obst_s2d(iv,i,j)
           ENDIF ! * END TEST ON POSITION
         ENDDO ! * END LOO ON IV
         obst_a2d(obst_nbvar+1,i,j) = MIN(obst_a2d(obst_nbvar+1,i,j),gamma)
         obst_a2d(obst_nbvar+2,i,j) = MIN(obst_a2d(obst_nbvar+2,i,j),gamma)
         obst_a2d(obst_nbvar+3,i,j) = MIN(obst_a2d(obst_nbvar+3,i,j),gamma)
         obst_s2d(obst_nbvar+1,i,j) = MIN(obst_s2d(obst_nbvar+1,i,j),gamma)
         obst_s2d(obst_nbvar+2,i,j) = MIN(obst_s2d(obst_nbvar+2,i,j),gamma)
         obst_s2d(obst_nbvar+3,i,j) = MIN(obst_s2d(obst_nbvar+3,i,j),gamma)
       ENDIF ! * END TEST ON HWAT
       !---------------------------
       ! * Averaged characteristics
       !---------------------------
       spos = 0.0_rsh
       IF(hwat.GT.h0fond) THEN
         DO iv=1,obst_nbvar
           obst_height_mean(i,j) = obst_height_mean(i,j) + obst_position(iv,i,j)*obst_height(iv,i,j)
           wtmp = 0.0_rsh
           dtmp = 0.0_rsh
           sfz  = 0.0_rsh
           DO k=1,obst_kmax
             wtmp = wtmp + obst_fracz3d(iv,k,i,j) * obst_width3d(iv,k,i,j)
             dtmp = dtmp + obst_fracz3d(iv,k,i,j) * obst_dens3d(iv,k,i,j)
             sfz  = sfz  + obst_fracz3d(iv,k,i,j)
           ENDDO
           IF (sfz.GT.0.0_rsh)THEN
             wtmp = wtmp/sfz
             dtmp = dtmp/sfz
           ENDIF
           obst_dens_mean(i,j)   = obst_dens_mean(i,j)   + obst_position(iv,i,j)*dtmp
           obst_width_mean(i,j)  = obst_width_mean(i,j)  + obst_position(iv,i,j)*wtmp
           spos = spos + obst_position(iv,i,j)
         ENDDO
       ELSE !* Not enough water
         DO iv=1,obst_nbvar
           obst_dens_mean(i,j)   = obst_dens_mean(i,j)   + obst_position(iv,i,j)*obst_dens_inst(iv,i,j)
           obst_width_mean(i,j)  = obst_width_mean(i,j)  + obst_position(iv,i,j)*obst_width_inst(iv,i,j)
           obst_height_mean(i,j) = obst_height_mean(i,j) + obst_position(iv,i,j)*obst_height_inst(iv,i,j)
           spos = spos + obst_position(iv,i,j)
         ENDDO
       ENDIF ! * END TEST ON HWAT 
       IF (spos.GT.0.0_rsh)THEN
         obst_dens_mean(i,j)   = obst_dens_mean(i,j)   / spos
         obst_width_mean(i,j)  = obst_width_mean(i,j)  / spos
         obst_height_mean(i,j) = obst_height_mean(i,j) / spos
       ENDIF
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_projarea

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_obstroughness(limin,limax,ljmin,ljmax,h0,ssh)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_obstroughness  ***
   !&E
   !&E ** Purpose : Computes obstuctions bottom roughness
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2017-04-06 (F. Ganthy) Original code
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2018-04-19 (F. Ganthy) Add Abdelrhman formulation of roughness for simplified obstructions
   !&E       ! 2018-10-24 (F. Ganthy) Bug correction on averaging method for multiples obstructions
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE

   !! * Arguments
   !! * Local declaration
   INTEGER         :: iv,i,j
   REAL(KIND=rsh)  :: hwat,a,d,z0a,coef,z0tot,ctot
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   DO j=ljmin,ljmax
     DO i=limin,limax
       obst_z0obst(:,i,j) = obst_z0bed(i,j)
       hwat = h0(i,j) + ssh(i,j)
       IF(hwat.GT.h0fond) THEN
         ! ******************************************************* !
         ! *** Computation for all upward and not 3d variables *** !
         ! ******************************************************* !
         DO iv=1,obst_nbvar
           !--------------------------------
           ! *** Abdelhrman parameterization
           !--------------------------------
           IF((.NOT.l_obst_downward(iv)).AND.(.NOT.l_obst_3dobst(iv)))THEN
             IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
               IF(l_obst_abdelrough_cste(iv))THEN
                 coef = obst_c_crough_x0(iv)
               ELSE
                 coef = obst_c_crough_x1(iv) + obst_c_crough_x0(iv) * (obst_height_inst(iv,i,j)**2.0_rsh)/ &
                        ((obst_height(iv,i,j)**2.0_rsh) * obst_width_inst(iv,i,j)*obst_dens_inst(iv,i,j))
               ENDIF
               a   = obst_width_inst(iv,i,j)*obst_height_inst(iv,i,j)*(obst_height(iv,i,j)/obst_height_inst(iv,i,j))
               d   = (coef*obst_height(iv,i,j)*obst_height(iv,i,j)*obst_width_inst(iv,i,j)) / &
                     (a+coef*obst_width_inst(iv,i,j)*obst_height(iv,i,j))
               z0a = (0.5_rsh*obst_width_inst(iv,i,j)*obst_height(iv,i,j)*obst_height(iv,i,j)*a) / &
                     (a+coef*obst_width_inst(iv,i,j)*obst_height(iv,i,j))**2.0_rsh
               z0a = MAX(z0a+d,obst_z0bed(i,j))
               ! Now apply fraction
               obst_z0obst(iv,i,j) = (z0a*obst_position(iv,i,j))+((1.0_rsh-obst_position(iv,i,j))*obst_z0bed(i,j))
             ELSE
               obst_z0obst(iv,i,j)   = obst_z0bed(i,j)
             ENDIF
           ELSE
             obst_z0obst(iv,i,j)   = obst_z0bed(i,j)
           ENDIF
         ENDDO ! * END LOOP OBST_NBVAR
         ! ************************************************* !
         ! *** Total computation for NoTurb Formulations *** !
         ! ************************************************* !
         ! As a weighted-average
         !----------------------
         z0tot = 0.0_rsh
         ctot  = 0.0_rsh
         DO iv=1,obst_nbvar
           IF(l_obst_noturb(iv))THEN
             z0tot = z0tot + obst_z0obst(iv,i,j)
             ctot  = ctot  + 1.0_rsh
           ENDIF
         ENDDO
         IF(z0tot.EQ.0.0_rsh) THEN
           obst_z0obst(obst_nbvar+1,i,j) = obst_z0bed(i,j)
         ELSE
           obst_z0obst(obst_nbvar+1,i,j) = MAX(z0tot/ctot,obst_z0bed(i,j))
         ENDIF
         ! *********************************************** !
         ! *** Total computation for Turb Formulations *** !
         ! *********************************************** !
         ! As a weighted-average
         !----------------------
         z0tot = 0.0_rsh
         ctot  = 0.0_rsh
         DO iv=1,obst_nbvar
           IF(.NOT.l_obst_noturb(iv))THEN
             z0tot = z0tot + obst_z0obst(iv,i,j)
             ctot  = ctot  + 1.0_rsh
           ENDIF
         ENDDO
         IF(z0tot.EQ.0.0_rsh) THEN
           obst_z0obst(obst_nbvar+2,i,j) = obst_z0bed(i,j)
         ELSE
           obst_z0obst(obst_nbvar+2,i,j) = MAX(z0tot/ctot,obst_z0bed(i,j))
         ENDIF
         ! ************************* !
         ! *** Total computation *** !
         ! ************************* !
         ! As a weighted-average
         !----------------------
         z0tot = 0.0_rsh
         ctot  = 0.0_rsh
         DO iv=1,obst_nbvar
             z0tot = z0tot + obst_z0obst(iv,i,j)
             ctot  = ctot  + 1.0_rsh
         ENDDO
         IF(z0tot.EQ.0.0_rsh) THEN
           obst_z0obst(obst_nbvar+3,i,j) = obst_z0bed(i,j)
         ELSE
           obst_z0obst(obst_nbvar+3,i,j) = MAX(z0tot/ctot,obst_z0bed(i,j))
         ENDIF
       ENDIF ! * End Test on Hwat>H0fond
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-----------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_obstroughness

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_hydroparam(limin,limax,ljmin,ljmax)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_hydroparam  ***
   !&E
   !&E ** Purpose : Computes obstuctions parameters used for hydrodynamics
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays : obst_position, obst_height, obst_kmin, obst_kmax,
   !&E                     obst_width3d, obst_thick3d, obst_thick3d, 
   !&E                     obst_dens3d, u, v, uz, vz, h0, ssh
   !&E
   !&E ** Modified variables : obst_drag3d, z0b, obst_fu_i,
   !&E                         obst_fu_e, obst_fuz_i, obst_fuz_e, obst_fv_i, obst_fv_e,
   !&E                         obst_fvz_i, obst_fvz_e, obst_a, obst_t, obst_tau,
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2014-10-07 (F. Ganthy) Re-organisation and optimization of the routine
   !&E                                Error messages added
   !&E       ! 2014-10-16 (F. Ganthy) Some modifications
   !&E       ! 2016-03-11 (F. Ganthy) Modifications of drag coefficient computation
   !&E       ! 2016-03-11 (F. Ganthy) Add logical to choose turbulence coefficient
   !&E                                values (default or user-defined)
   !&E       ! 2016-03-11 (F. Ganthy) Add fraction of sigma layers occupied by obstructions
   !&E       ! 2016-03-21 (F. Ganthy) Better parameterization of obstructions for 2D model (l_model2d=.TRUE.)
   !&E       ! 2017-02-16 (F. Ganthy) Some modifications:
   !&E                                - Allowing multiple obstructions type in a single grid cell
   !&E                                  --> allowed multispecific computation.
   !&E                                  This imply that some tables must be allocated depending on (iv,k,i,j) or (iv,i,j)
   !&E                                  These allocations are done within obst_alloc_xyz (because tables allocated within
   !&E                                  obst_alloc_nbvar are only those read within namelist).
   !&E                                - Changes on instantaneous obstruction state variables for future coupling with Zostera growth module
   !&E                                - Differenciation of cylindric / parallelepipedic structures
   !&E                                - Taking into accounts for horizontal fractionning of obstructions (no empty grid cell)
   !&E       ! 2017-04-05 (F. Ganthy) Modification for case of 2D model
   !&E       ! 2017-04-06 (F. Ganthy) Change initialization of turbulence coefficients cmu and c2turb
   !&E       ! 2017-04-13 (F. Ganthy) Allow 3D obstructions (e.g. oyster bags)
   !&E       ! 2017-04-13 (F. Ganthy) Add simplified formulation (not using turbulence, but using roughness)
   !&E       ! 2017-10-17 (F. Ganthy) Re-organization related to simplified formulations, 3D obstructions and
   !&E                                multiple obstructions within one grid cell
   !&E       ! 2017-11-03 (F. Ganthy) Remove useless variables related to turbulence coefficients
   !&E       ! 2017-12-09 (F. Ganthy) Add some comments about units of variables
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
!    USE comvars2d,  ONLY  : h0,h0fond,z0b
!    USE comvarp2d,  ONLY  : ssh

!    USE comturb,    ONLY  : cmu,c2turb

   IMPLICIT NONE

   !! * Arguments
   !! * Local declaration
   INTEGER                           :: iv,i,j,k
   REAL(KIND=rsh)                    :: hwat,phi,ntot,clz,nclz,lz
   REAL(KIND=rsh)                    :: uzvz,fuv,tuz,tvz
   REAL(KIND=rsh),PARAMETER          :: gamma  = 1.0E-5
   !!----------------------------------------------------------------------
   !! * Executable part
   !-----------------------------------------
   ! **********************
   ! * RE-INITIALIZATION OF VARIABLES
   ! ***********************
   !-------------------
   ! Variables on (i,j)
   !-------------------
   obst_fu_i(:,:)    = 0.0_rsh
   obst_fv_i(:,:)    = 0.0_rsh
   obst_fu_e(:,:)    = 0.0_rsh
   obst_fv_e(:,:)    = 0.0_rsh
   !---------------------
   ! Variables on (k,i,j)
   !---------------------
   obst_fuz_i(:,:,:) = 0.0_rsh
   obst_fvz_i(:,:,:) = 0.0_rsh
   obst_fuz_e(:,:,:) = 0.0_rsh
   obst_fvz_e(:,:,:) = 0.0_rsh
   obst_t(:,:,:)     = 0.0_rsh
   obst_tau(:,:,:)   = 0.0_rsh
   !-----------------------
   ! Variable on (iv,k,i,j)
   !-----------------------
   obst_drag3d(:,:,:,:) = 0.0_rsh
   !*************************
   DO j=ljmin,ljmax
     DO i=limin,limax
       hwat = h0(i,j) + ssh(i,j)
       IF(hwat.GT.h0fond)THEN
         ! *************************************************** !
         ! ********** SIMPLIFIED (NoTurb) VARIABLES ********** !
         ! *************************************************** !
         IF(obst_z0obst(obst_nbvar+1,i,j) /= obst_z0bed(i,j))THEN
           z0b(i,j) = obst_z0obst(obst_nbvar+1,i,j)
         ENDIF
         ! *************************************************** !
         ! *************** TURBULENT VARIABLES *************** !
         ! *************************************************** !
         DO k=1,obst_kmax
           uzvz = SQRT(obst_uz(k,i,j)**2.0_rsh + obst_vz(k,i,j)**2.0_rsh)
           IF(uzvz.GT.0.0001_rsh)THEN
             ! *************************************************** !
             ! *********** TURBULENT VARIABLES 1ST PART ********** !
             ! *************************************************** !
             ntot = 0.0_rsh
             tuz  = 0.0_rsh
             tvz  = 0.0_rsh
             clz  = 0.0_rsh
             nclz = 0.0_rsh
             DO iv=1,obst_nbvar
               IF((obst_position(iv,i,j).GT.0.0_rsh) .AND. (.NOT.l_obst_noturb(iv)))THEN
                 IF(obst_dens3d(iv,k,i,j).GT.0.0_rsh)THEN
                   !----------------------------
                   ! * Total obstruction density
                   !----------------------------
                   ntot = ntot + obst_dens3d(iv,k,i,j)*obst_fracxy(iv,i,j)
                   !-------------------
                   ! * Drag coefficient         
                   !-------------------
                   IF(l_obst_drag_cste(iv)) THEN
                     obst_drag3d(iv,k,i,j) = obst_c_drag(iv)
                   ELSE
                     phi =  (pi/2.0_rsh) - obst_theta3d(iv,k,i,j)
                     obst_drag3d(iv,k,i,j) = gamma + ABS(phi) * (obst_c_drag(iv)-gamma)/(pi/2.0_rsh)
                   ENDIF
                     !-------------------
                     ! * Resistance force
                     !-------------------
                     ! Here, fuv has unit:
                     ! [fuv] = - * - * m * m.s-1 * m-2 * - * -
                     ! [fuv] = s-1
                     ! Here obst_fuz_e (obst_fvz_e) has unit:
                     ! [obst_fuz_e]     = [fuv] * [obst_uz]
                     ! [obst_fuz_e]     = s-1   * m.s-1
                     ! [obst_fuz_e]     = m.s-2
                     ! [obst_fuz_e*rho] = kg.m-3 * m.s-2
                     ! [obst_fuz_e*rho] = N.m-3
                     ! Here obst_fuz_i (obst_fvz_i) has unit:
                     ! [obst_fuz_i]     = [fuv]
                     ! [obst_fuz_i]     = s-1
                     !-------------------
                     fuv = 0.5_rsh * obst_drag3d(iv,k,i,j) * obst_width3d(iv,k,i,j) * uzvz * &
                           obst_dens3d(iv,k,i,j) * obst_fracxy(iv,i,j) * obst_fracz3d(iv,k,i,j)
                     obst_fuz_e(k,i,j) = obst_fuz_e(k,i,j) + (fuv * obst_uz(k,i,j))
                     obst_fvz_e(k,i,j) = obst_fvz_e(k,i,j) + (fuv * obst_vz(k,i,j))
                     obst_fuz_i(k,i,j) = obst_fuz_i(k,i,j) + fuv
                     obst_fvz_i(k,i,j) = obst_fvz_i(k,i,j) + fuv 
                     !--------------------------
                     ! * Work spent by the fluid
                     !--------------------------
                     ! Here, tuz has unit:
                     ! [tuz] = s-1 * m.s-1 * m.s-1
                     ! [tuz] = m2.s-2
                     tuz = tuz + fuv * (obst_uz(k,i,j)**2.0_rsh)
                     tvz = tvz + fuv * (obst_vz(k,i,j)**2.0_rsh)
                     !-------------------------------
                     ! * Averaging coefficient for lz
                     !-------------------------------
                     clz  = clz  + obst_c_lz(iv)
                     nclz = nclz + 1.0_rsh
                 ENDIF ! * END TEST ON OBSTRUCTION DENSITY
               ENDIF ! * END TEST ON POSITION AND TURB VARIABLE
             ENDDO ! * END LOOP ON NBVAR
             ! *************************************************** !
             ! *********** TURBULENT VARIABLES 2ND PART ********** !
             ! *************************************************** !
             IF(ntot.GT.0.0_rsh)THEN
               !--------------------------
               ! * Work spent by the fluid
               !--------------------------
               ! Here, obst_t has unit:
               ! [obst_t] = m3.s-2
               obst_t(k,i,j)    = SQRT(tuz**2.0_rsh + tvz**2.0_rsh)
               !-------------------------------------
               ! * Smallest distance between elements
               !-------------------------------------
               clz = clz/nclz
               lz  = clz * SQRT((1.0_rsh-obst_a3d(obst_nbvar+2,k,i,j)) / ntot)
               !-----------------------------------------------------
               ! * Dissipation timescale of eddies between structures
               ! ----------------------------------------------------
               IF((obst_t(k,i,j)/=0.0_rsh).AND.(c2turb/=0.0_rsh).AND.(cmu/=0.0_rsh))THEN
                 obst_tau(k,i,j) = 1.0_rsh / (1.0_rsh / (c2turb * sqrt(cmu)) * ((lz*lz) / (obst_t(k,i,j)))**(1.0_rsh/3.0_rsh))
               ENDIF
             ENDIF ! * END TEST ON OBSTRUCTION DENSITY
           ENDIF ! * END TEST ON ENOUGH VELOCITY
         ENDDO ! * END LOOP ON k
         ! *************************************************** !
         ! *********** TURBULENT VARIABLES 3RD PART ********** !
         ! *************************************************** ! 
         ! ********************************************* !
         ! ****************** 2D FORCES **************** !
         ! ********************************************* !
         !-----------------------------------
         ! * Resistance force : depth average
         !-----------------------------------
         DO k=1,obst_kmax
           obst_fu_e(i,j) = obst_fu_e(i,j) + obst_fuz_e(k,i,j) * obst_dz(k,i,j)
           obst_fv_e(i,j) = obst_fv_e(i,j) + obst_fvz_e(k,i,j) * obst_dz(k,i,j)
           obst_fu_i(i,j) = obst_fu_i(i,j) + obst_fuz_i(k,i,j) * obst_dz(k,i,j)
           obst_fv_i(i,j) = obst_fv_i(i,j) + obst_fvz_i(k,i,j) * obst_dz(k,i,j)
         ENDDO
         obst_fu_e(i,j) = obst_fu_e(i,j) / hwat
         obst_fv_e(i,j) = obst_fv_e(i,j) / hwat
         obst_fu_i(i,j) = obst_fu_i(i,j) / hwat
         obst_fv_i(i,j) = obst_fv_i(i,j) / hwat
         !---------------------------------------------
         ! * Resistance force : implicit/explicit parts
         !---------------------------------------------
         obst_fu_e(i,j) = obst_fu_e(i,j) * obst_c_exp2d
         obst_fv_e(i,j) = obst_fv_e(i,j) * obst_c_exp2d
         obst_fu_i(i,j) = obst_fu_i(i,j) * obst_c_imp2d
         obst_fv_i(i,j) = obst_fv_i(i,j) * obst_c_imp2d
       ! ********************************************* !
       ! ****************** 3D FORCES **************** !
       ! ********************************************* !
       !---------------------------------------------
       ! * Resistance force : implicit/explicit parts
       !---------------------------------------------
       DO k=1,obst_kmax
         obst_fuz_e(k,i,j) = obst_fuz_e(k,i,j) * obst_c_exp3d
         obst_fvz_e(k,i,j) = obst_fvz_e(k,i,j) * obst_c_exp3d
         obst_fuz_i(k,i,j) = obst_fuz_i(k,i,j) * obst_c_imp3d
         obst_fvz_i(k,i,j) = obst_fvz_i(k,i,j) * obst_c_imp3d
       ENDDO
       ENDIF ! * END TEST ON HWAT
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_hydroparam

   !!==========================================================================================================

   SUBROUTINE OBSTRUCTIONS_comp_bedroughness(limin,limax,ljmin,ljmax,h0,ssh)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE OBSTRUCTIONS_comp_bedroughness  ***
   !&E
   !&E ** Purpose : Computes roughness length for bottom shear stress computation
   !&E
   !&E ** Description : Here, depending on obstructions you wanted to simulate,
   !&E                  you may have to add your own parameterization
   !&E
   !&E ** Called by : obst_init and obst_update
   !&E
   !&E ** External calls :
   !&E
   !&E ** Used ij-arrays :
   !&E
   !&E ** Modified variables : 
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       ! 2012-04-20 (F. Ganthy) Original code
   !&E       ! 2018-02-05 (F. Ganthy) Change lot of calibration and organization
   !&E       ! 2018-10-24 (F. Ganthy) Bug correction on averaging method for multiples obstructions
   !&E
   !&E---------------------------------------------------------------------
   IMPLICIT NONE

   !! * Local declaration
   INTEGER :: iv,i,j
   REAL(KIND=rsh) :: hwat,z0tmp,stmp,z00,oal,oah
   REAL(KIND=rsh),PARAMETER :: epsi=1E-6
   !!----------------------------------------------------------------------
   !! * Executable part
   !-------------------------------------
   obst_z0bstress(:,:) = obst_i_z0bstress
   DO j=ljmin,ljmax
     DO i=limin,limax
       hwat = h0(i,j) + ssh(i,j)
       IF(hwat.GT.h0fond)THEN
         z0tmp = 0.0_rsh
         stmp  = 0.0_rsh
         DO iv=1,obst_nbvar
           IF((.NOT.l_obst_downward(iv)).AND.(.NOT.l_obst_3dobst(iv)).AND.(l_obst_z0bstress(iv)))THEN
             IF(obst_position(iv,i,j).GT.0.0_rsh)THEN
               IF(obst_z0bstress_option(iv).EQ.0)THEN
                 !-----------------------
                 ! Constant value is used
                 !-----------------------
                 z00 = obst_c_z0bstress(iv)
               ELSE
                 !-------------------------
                 ! Parameterization is used
                 !-------------------------
                 oah = obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_height(iv,i,j)
                 z00 = obst_c_z0bstress_x0(iv)* oah**obst_c_z0bstress_x1(iv)
                 z00 = MAX(MIN(z00,0.01_rsh),epsi)
               ENDIF ! END test on parameterization
               z0tmp = (obst_position(iv,i,j)*z00) + ((1.0_rsh-obst_position(iv,i,j))*obst_i_z0bstress)
               stmp  = stmp + 1.0_rsh
             ELSE
               z0tmp = z0tmp + obst_i_z0bstress
               stmp  = stmp + 1.0_rsh
             ENDIF           
           ELSE
             z0tmp = z0tmp + obst_i_z0bstress
             stmp  = stmp + 1.0_rsh
           ENDIF
         ENDDO ! END LOOP nbvar
         !----------
         ! Averaging
         !----------
         obst_z0bstress(i,j) = z0tmp/stmp
         obst_z0bstress(i,j) = MAX(epsi,obst_z0bstress(i,j))

       ENDIF ! * END TEST ON HWAT
     ENDDO ! * END LOOP ON i
   ENDDO ! * END LOOP ON j
   !-------------------------------------
   END SUBROUTINE OBSTRUCTIONS_comp_bedroughness

   !!==========================================================================================================

#endif
END MODULE OBSTRUCTIONS
