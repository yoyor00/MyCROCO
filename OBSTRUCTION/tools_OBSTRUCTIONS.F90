   !&E==========================================================================
   !&E                   ***  tools_OBSTRUCTIONS  ***
   !&E
   !&E ** Purpose : Contains functions tools for obstructions module
   !&E
   !&E ** Description :
   !&E
   !&E       function obsttools_abdeluz  : Computes velocity profile within or
   !&E                                     outside obstructions from depth-averaged
   !&E                                     velocity using Abdelrhman's (2003) method
   !&E
   !&E ** History :
   !&E       ! 2018-04-18 (F. Ganthy) Original code
   !&E==========================================================================
#ifdef key_OBSTRUCTIONS

   function obsttools_abdeluz(i,j,k,hwat,uv,z,z0)
   !!**********************************************************************
   !! Computes velocity profile within or outside obstructions from
   !! depth-averaged velocity using Abdelrhman's (2003) method
   !! F. Ganthy (2018-04-18)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE parameters,      ONLY: rsh
   USE comvars2d,       ONLY: h0fond
   USE comobstructions
   IMPLICIT NONE
   !! * Arguments
   INTEGER                  :: i,j,k
   REAL(KIND=rsh)           :: hwat,uv,z0,z
   REAL(KIND=rsh)           :: obsttools_abdeluz
   !! * Local declaration
   INTEGER                  :: iv
   REAL(KIND=rsh)           :: fobst,utot,divtot,utmp
   REAL(KIND=rsh)           :: ga,pa,pus,gus,uhc,d,gz0,alp,coef
   REAL(KIND=rsh),PARAMETER :: e1=EXP(1.0_rsh)
   !!----------------------------------------------------------------------
   !! * Executable part
   !!------------------
   IF(hwat.GT.h0fond)THEN
     utot   = 0.0_rsh
     divtot = 0.0_rsh
     DO iv=1,obst_nbvar
       IF((.NOT.l_obst_downward(iv)).AND.(.NOT.l_obst_3dobst(iv)).AND.(obst_height(iv,i,j).GT.0.0_rsh))THEN
         !---------------------------------------
         ! * Fraction of cell without obstruction
         !---------------------------------------
         IF(z.LE.z0)THEN
           utmp   = 0.0_rsh
         ELSE
           utmp   = uv * LOG(z/z0) / LOG(hwat/(e1*z0))
         ENDIF
         utot   = utot + (1.0_rsh-obst_position(iv,i,j))*utmp
         divtot = divtot + (1.0_rsh-obst_position(iv,i,j))
         !------------------------------------
         ! * Fraction of cell with obstruction
         !------------------------------------
         IF(obst_height(iv,i,j).LT.0.7_rsh*hwat)THEN
           ! * Fully submerged obstruction
           !------------------------------
           IF(l_obst_abdelrough_cste(iv))THEN
             coef = obst_c_crough_x0(iv)
           ELSE
             coef = obst_c_crough_x1(iv) + obst_c_crough_x0(iv) * (obst_height_inst(iv,i,j)**2.0_rsh)/ &
                    ((obst_height(iv,i,j)**2.0_rsh) * obst_width_inst(iv,i,j)*obst_dens_inst(iv,i,j))
           ENDIF
           ! Computes A
           !write(*,*) "--------"
           !write(*,*) "hw",hwat
           !write(*,*) "w",obst_width_inst(iv,i,j)
           !write(*,*) "hi",obst_height_inst(iv,i,j)
           !write(*,*) "hei",obst_height(iv,i,j)
           !write(*,*) "ac",ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j))
           !write(*,*) "cos",COS(ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j)))
           IF(obst_height(iv,i,j).GE.obst_height_inst(iv,i,j)) THEN
               ga  = obst_width_inst(iv,i,j)*MAX(obst_height_inst(iv,i,j),obst_thick_inst(iv,i,j))
           ELSE
               ga  = obst_width_inst(iv,i,j)*MAX(obst_height_inst(iv,i,j)*COS(ACOS(obst_height(iv,i,j)/obst_height_inst(iv,i,j))),obst_thick_inst(iv,i,j))
           ENDIF
           ! Computes a
           pa  = obst_width_inst(iv,i,j)*obst_dens_inst(iv,i,j)
           ! Computes u*
           pus = 0.4_rsh * (uv/LOG(hwat/(e1*z0)))
           ! Computes Z0
           gz0 = (0.5_rsh*obst_width_inst(iv,i,j)*obst_height(iv,i,j)*obst_height(iv,i,j)*ga) / &
                 (ga+coef*obst_width_inst(iv,i,j)*obst_height(iv,i,j))**2.0_rsh
           gz0 = MAX(gz0,z0)
           ! Computes d
           d   = (coef*obst_height(iv,i,j)*obst_height(iv,i,j)*obst_width_inst(iv,i,j)) / &
                 (ga+coef*obst_width_inst(iv,i,j)*obst_height(iv,i,j))
           ! Computes alpha
           alp = obst_height(iv,i,j) * (obst_c_drag(iv)*pa/(2.0_rsh*d*d))**(1.0_rsh/3.0_rsh)
           ! Computes U*
           gus = (pus*((hwat+z0)*LOG((hwat+z0)/z0) - hwat)) /                                           &
                 ((obst_height(iv,i,j)/alp)*LOG((obst_height(iv,i,j)-d+z0)/gz0)*(1.0_rsh - EXP(-alp)) + &
                 (hwat-d+z0)*LOG((hwat-d+z0)/gz0) -                                                     &
                 (obst_height(iv,i,j)-d+z0)*LOG((obst_height(iv,i,j)*d+z0)/gz0) - (hwat-obst_height(iv,i,j)))
           ! Computes u(hobst)
           uhc = gus/0.4_rsh * LOG((obst_height(iv,i,j)-d+z0)/gz0)
           ! Computes uz for corresponding obstruction variable
           IF(z.LE.obst_height(iv,i,j)) THEN
             utmp = uhc*EXP(alp*(z/obst_height(iv,i,j) - 1.0_rsh))
           ELSE
             utmp = gus/0.4_rsh * LOG((z-d+z0)/gz0)
           ENDIF
           utot   = utot   + obst_position(iv,i,j)*utmp
           divtot = divtot + obst_position(iv,i,j)
         ELSE
           ! * Emerged of constriced obstruction
           !------------------------------------
           IF(z.LE.z0)THEN
             utmp   = 0.0_rsh
           ELSE
             utmp   = uv * LOG(z/z0) / LOG(hwat/(e1*z0))
           ENDIF
           utot   = utot + obst_position(iv,i,j)*utmp
           divtot = divtot + obst_position(iv,i,j)
         ENDIF
       ELSE ! Obstruction variable incompatible with Abdelrhman
         IF(z.LE.z0)THEN
           utmp = 0.0_rsh
         ELSE
           utmp = uv * LOG(z/z0) / LOG(hwat/(e1*z0))
         ENDIF
         utot   = utot   + utmp
         divtot = divtot + 1.0_rsh
       ENDIF ! * End test on obstruction kind
     ENDDO ! * End loop on iv
     !-----------------------------
     ! * Average over all variables
     !-----------------------------
     obsttools_abdeluz = utot/divtot
   ELSE
     obsttools_abdeluz = 0.0_rsh
   ENDIF ! * End test on hawt
   !!**********************************************************************
   end function obsttools_abdeluz

   !!======================================================================

   function obsttools_z0sedim(i,j,z0sedim,hwat)
   !!**********************************************************************
   !! Computes z0sedim in presence of obstructions
   !! F. Ganthy (2019-06-12)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE parameters, ONLY  : rsh
   USE comvars2d,  ONLY  : h0,h0fond,hm
   USE comvarp2d,  ONLY  : ssh
   USE comobstructions
   IMPLICIT NONE
   !! * Arguments
   INTEGER                  :: i,j
   REAL(KIND=rsh)           :: z0sedim,hwat
   REAL(KIND=rsh)           :: obsttools_z0sedim
   !! * Local declaration
   INTEGER :: iv
   REAL(KIND=rsh) :: z0tmp,stmp,z00,oal,oah
   REAL(KIND=rsh),PARAMETER :: epsi=1E-6
   !!----------------------------------------------------------------------
   !! * Executable part
   !!------------------
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
           IF(hwat.LE.hm)THEN
             !---------------------------------------
             ! Applying correction for 2D small-depth
             !---------------------------------------
             oal = 1.0_rsh/(obst_dens_inst(iv,i,j)*obst_width_inst(iv,i,j)*obst_height_inst(iv,i,j))
             z00 = z00 * obst_c_z0bstress_x2(iv)*oal
           ENDIF ! END test on 3D or 2Dsmall depth
           z00 = MAX(MIN(z00,0.01_rsh),epsi)
         ENDIF ! END test on parameterization
         z0tmp = (obst_position(iv,i,j)*z00) + ((1.0_rsh-obst_position(iv,i,j))*z0sedim)
         stmp  = stmp + 1.0_rsh
       ELSE
         z0tmp = z0tmp + z0sedim
         stmp  = stmp + 1.0_rsh
       ENDIF
     ELSE
       z0tmp = z0tmp + z0sedim
       stmp  = stmp + 1.0_rsh
     ENDIF
   ENDDO ! END LOOP nbvar
   !----------
   ! Averaging
   !----------
   obsttools_z0sedim = z0tmp/stmp
   obsttools_z0sedim = MAX(epsi,obsttools_z0sedim)
   ! Specific value for testcase
#if defined key_casobstflume_ganthy2015_hydro
   IF((i.GE.14) .AND. (i.LE.31))THEN
     IF(obst_position(1,i,j).EQ.0.0_rsh)THEN
       ! Bare sediment
       obsttools_z0sedim = 1.0E-5
     ENDIF
   ENDIF
#elif defined key_casobstflume_ganthy2015_sedim
   IF((i.GE.5) .AND. (i.LE.11))THEN
     IF(obst_position(1,i,j).EQ.0.0_rsh)THEN
       ! Bare sediment
       obsttools_z0sedim = 1.0E-5
     ENDIF
   ENDIF
#endif
   !!**********************************************************************
   end function obsttools_z0sedim

   !!======================================================================


#endif
   !&E==========================================================================
   !&E                   *** end  OBSTRUCTIONS_tools  ***
   !&E==========================================================================
