!***************************************************************************
!***************************************************************************
!Copyright or (c) or Copr. : IFREMER
!contributor(s) : IFREMER/DYNECO/DHYSED
!
!contact Ifremer : mustang@ifremer.fr
!
!This software (MUSTANG, MUd and Sand TrAnsport modelliNG) is a Fortran F90 
!computer program whose purpose is to perform sediment transport process 
!modelling coupled to hydrodynamic models.
!Full details can be obtained on https://wwz.ifremer.fr/dyneco/MUSTANG
!
!This software is governed by the CeCILL-C license under French law and
!abiding by the rules of distribution of free software. You can use, 
!modify and/ or redistribute the software under the terms of the CeCILL-C
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************

#include "cppdefs.h"
!---------------------------------------------------------------------------
!
                     MODULE coupleur_MUSTANG
!
!---------------------------------------------------------------------------
   

#if defined MUSTANG 

   !&E==========================================================================
   !&E                   ***  MODULE  coupleur_MUSTANG  ***
   !&E
   !&E
   !&E ** Purpose : concerns coupling MUSTANG with hydro code
   !&E              
   !&E 
   !&E ** Description :
   !&E     subroutine coupl_MUSTANG_substance   ! definition number et type of variables 
   !&E     subroutine coupl_conv2MUSTANG        ! calculation of some variables needed by MUSTANG
   !&E     subroutine coupl_MUSTANG2hydro       ! transfert from MUSTANG to hydro code
   !&E
   !&E ** History :
   !&E     ! 2016-11  (B.Thouvenin) : 
   !&E
!&E===================================================================================================================
#ifdef key_MARS
#include "toolcpp.h"
#endif
#include "coupleur_define_MUSTANG.h"

   USE comMUSTANG
   USE comsubstance

   IMPLICIT NONE
   
   !! * Interface
   
   
   !! * Accessibility

   ! functions & routines of this module, called outside :
#if ! defined key_MARS 
    PUBLIC coupl_conv2MUSTANG, coupl_MUSTANG2hydro
#else
    PUBLIC coupl_conv2MUSTANG
#endif

   PRIVATE


#if ! defined key_MARS
   !! * Shared or public module variables (variables used by MUSTANG but issued from hydro model or substances module  )
   !! * for MARS_MODEL, these variables are stored in comvar.. or comsubstance which are not the same for other model

#endif
 
 CONTAINS
   !!===========================================================================

  SUBROUTINE coupl_conv2MUSTANG(ifirst,ilast,jfirst,jlast,iappel,BATHY_H0,ssh,    &
                     WATER_CONCENTRATION)

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE coupl_conv2MUSTANG  ***
   !&E
   !&E ** Purpose : transfer cv_Wat and htot et alt_cw1 computation
   !&E               only inside the domain, not at boundaries meshes
   !&E
   !&E         for MARS : ifirst=imin+2, ilast=imax-1, jfirst=jmin+2,  jlast=jmax-1
   !&E                    for interior processors : ifirst=limin, ilast=limax, jfirst=ljmin,  jlast=ljmax
   !&E                    to be inside the domain and not at boundaries :
   !&E                                              from jmin+2=jb+1 a jmax-1=jlast=jh-1 
   !&E                                          and from imin+2=id+1 a imax-1=ilast=id-1
   !&E                    ig(j)=imin+1 & id(j)=imax 
   !&E                    domain limited by the coast thanks to ig,id,jb,jh (non calculation at land)
   !&E
   !&E ** Description :  
   !&E  arguments IN : BATHY_H0,ssh,WATER_CONCENTRATION, SALINITY_MOD, TEMPERATURE_MOD
   !&E  arguments OUT: no (all variables  in comMUSTANG)
   !&E     
   !&E   variables OUT :   
   !&E       htot (total water height )
   !&E       epn_bottom (thickness of water the bottom layer)
   !&E       sal_bottom_MUSTANG,temp_bottom_MUSTANG (salinity, temperature in the water bottom layer)
   !&E       cw_bottom_MUSTANG ( concentrations in the water bottom layer)
   !&E       roswat_bot ( water density  in the water bottom layer)
   !&E     
   !&E   initial call  (iappel=0 for initialization):  
   !&E      extraction of thickness, salinity, temperature and water concentrations and densities
   !&E                 in the bottom of the water column
   !&E
   !&E   first call (iappel=1) :  
   !&E      extraction of thickness, salinity, temperature and water concentrations and densities
   !&E                 in the bottom of the water column
   !&E       calculation of total water height    
   !&E       calculation of  alt_cw1  : altitude of the computation point of Cw in the  bottom layer
   !&E     
   !&E   second call (iappel=2) :  
   !&E      extraction of thickness, salinity, temperature and water concentrations and densities
   !&E                 in the bottom of the water column
   !&E       calculation of total water height    
   !&E       si not MARS : conversion to transmit to MUSTANG the hydro variables:
   !&E                     SETTL_FLUXSUM_w2s: effective deposit flux of the particle variables during transport
   !&E    
   !&E ** Called by :  MUSTANG_init_sediment (iappel=0)
   !&E                 sed_MUSTANG_update (iappel=1) & sed_MUSTANG_deposition (iappel=2)
   !&E
   !&E ** History :
   !&E       !  2016-11  (B.Thouvenin) 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used
#ifdef key_MARS
   USE comvars2d,    ONLY : ig,id,jb,jh,fwet
   USE comvars3d,    ONLY : dsigu,sigw,sig
   USE toolgeom,     ONLY : f_dzu
   USE_MPI toolmpi,    ONLY : ex_i_rsh,ex_j_rsh
#if defined key_siggen || defined key_gencoord
   USE toolgeom,       ONLY : f_hzu
#endif
#else
#include "scalars_F90.h"
#endif
   !! * Arguments 
   INTEGER, INTENT(IN)                           :: ifirst,ilast,jfirst,jlast,iappel
   REAL(KIND=rsh),DIMENSION(ARRAY_BATHY_H0),INTENT(IN)       :: BATHY_H0                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_ELEVATION),INTENT(IN):: ssh                         
   REAL(KIND=rsh),DIMENSION(ARRAY_WATER_CONC), INTENT(IN) :: WATER_CONCENTRATION  
   !REAL(KIND=rsh),DIMENSION(ARRAY_TEMPSAL), INTENT(IN)   :: SALINITY_MOD,TEMPERATURE_MOD  
   
   !! * Local declarations
   INTEGER                            :: iv,k,i,j



   !!---------------------------------------------------------------------------
   !! * Executable part

!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j)
   DO j=jfirst,jlast
#ifdef key_MARS
      DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else
      DO i=ifirst,ilast
      ! ATTENTION : not need to calculate at boundaries meshes where MUSTANG is not applied
#endif

           ! extraction of  concentrations in the bottom of the water column
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
            sal_bottom_MUSTANG(i,j)=SALINITY_MOD(1,i,j)
            temp_bottom_MUSTANG(i,j)=TEMPERATURE_MOD(1,i,j)
            cw_bottom_MUSTANG(:,i,j)=WATER_CONCENTRATION(:,1,i,j)
#else
! CROCO vecteur au temps 1, 2 ou 3 ????
            sal_bottom_MUSTANG(i,j)=WATER_CONCENTRATION(i,j,1,1,itemp+1)
            temp_bottom_MUSTANG(i,j)=WATER_CONCENTRATION(i,j,1,1,itemp)
            cw_bottom_MUSTANG(1:nv_adv,i,j)=WATER_CONCENTRATION(i,j,1,1,itsubs1:itsubs2)
#endif
            ! thickness of the bottom water layer or altitude at the top of the bottom layer
            ! + water density in the bottom water layer
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
            epn_bottom_MUSTANG(i,j)=f_dzu(BATHY_H0(i,j),ssh(i,j),1,i,j)
        ENDIF
#else
      ! CROCO z_w connu via scalars.h
               epn_bottom_MUSTANG(i,j)=z_w(i,j,1)-z_w(i,j,0)
               roswat_bot(i,j)= rho(i,j,1)+rho0  
#endif

       ENDDO
   ENDDO
!$OMP END DO

   IF(iappel > 0 ) THEN
   ! ATTENTION : need to calculate htot at all meshes (imin+1:imax, jmin+1: jmax)
   ! at interior meshes : for ljmin-1 and ljmax+1 , BATHY_H0 and ssh are known if MPI exchange was done after change
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j)
    DO j=jfirst-1,jlast+1
      DO i=ifirst-1,ilast+1

          !  htot : total water height
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
          htot(i,j)=BATHY_H0(i,j)+ssh(i,j)
#else
        ! CROCO
          htot(i,j)=z_w(i,j,N)+h(i,j)
#endif

       ENDDO
    ENDDO
!$OMP END DO
!     MPI exchange of total heights if needed in i+1 or i-1.
!      in the case of lateral erosion and deposit slip
!      and if not MPI exchange before
!#ifdef key_MARS     
!OMPMPI barrier
!OMPMPI master
!    CALL_MPI ex_i_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,htot(liminm2:limaxp2,ljminm2:ljmaxp2))
!    CALL_MPI ex_j_rsh(-2,2,1,liminm2,limaxp2,ljminm2,ljmaxp2,htot(liminm2:limaxp2,ljminm2:ljmaxp2))
!OMPMPI end master
!OMPMPI barrier
!OMPMPI flush(htot)
!#endif

   ENDIF
   IF (iappel == 1) THEN
   ! first call before evaluation of settling velocities, erosion, consolidation, diffusion
   
!$OMP DO SCHEDULE(RUNTIME) PRIVATE(i,j)
     DO j=jfirst,jlast
#ifdef key_MARS
       DO i=MAX0(ifirst,ig(j)+1),MIN0(ilast,id(j)-1)
        IF(j.GE.jb(i)+1 .AND. j .LE. jh(i)-1) THEN
#else
       DO i=ifirst,ilast
#endif
   
          ! alt_cw1 : altitude of the computation point of Cw in the  bottom layer
          ! could be eliminated if the concentration calculation point is always in the middle of the layer
          ! but if not true, you have to calculate this altitude differently
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef key_MARS
            alt_cw1(i,j)=0.5_rsh*f_dzu(BATHY_H0(i,j),ssh(i,j),1,i,j)
#else
   ! other modele than MARS
     ! CROCO
          alt_cw1(i,j) = z_r(i,j,1)-z_w(i,j,0) 
#endif
          

#ifdef key_MARS
        ENDIF
#endif
       ENDDO
     ENDDO
!$OMP END DO

    ELSE IF (iappel==2) THEN
#if ! defined key_MARS
   ! second call before evaluation of sediment deposition :  cumulated settling flux only if no key_MARS

!$OMP DO SCHEDULE(RUNTIME)
     DO j=jfirst,jlast
      DO i=ifirst,ilast
#if defined key_noTSdiss_insed
         DO iv=1,nv_adv
#else
         DO iv=-1,nv_adv
#endif
                flx_w2s_sum(iv,i,j)=SETTL_FLUXSUM_w2s(i,j,IV_HOSTMODEL) ! flux de depot cumule effectif apres transport
         ENDDO
     ENDDO
     ENDDO
!$OMP END DO
#endif
    ENDIF

  END SUBROUTINE coupl_conv2MUSTANG      

   !!==============================================================================
#if ! defined key_MARS
  SUBROUTINE coupl_MUSTANG2hydro(ifirst,ilast,jfirst,jlast)
                                          

   !&E--------------------------------------------------------------------------
   !&E                 ***  ROUTINE coupl_MUSTANG2flx ***
   !&E
   !&E ** Purpose : transfer corflux, corfluy, flx_w2s , ws3 ,flx_s2w
   !&E               for hydro code 
   !&E
   !&E ** Description : conversion for hydro code 
   !&E  arguments OUT: no because stored in comMUSTANG
   !&E     SETTL_FLUX_w2s : deposit trends 
   !&E     corflux_SAND,corfluy_SAND : correction of horizontal flux for sands 
   !&E     EROS_FLUX_s2w : erosion flux 
   !&E     EROS_FLUX_TEMP_s2w et eros_flix_SAL : erosion flux for temperature, salinity
   !&E     
   !&E     
   !&E ** Called by :  sed_MUSTANG_update
   !&E
   !&E ** History :
   !&E       !  2016-11  (B.Thouvenin) 
   !&E
   !&E--------------------------------------------------------------------------
   !! * Modules used


   !! * Arguments 
   INTEGER, INTENT(IN)                                    :: ifirst,ilast,jfirst,jlast                  
 
   
   !! * Local declarations
   INTEGER                            :: iv,k,i,j



   !!---------------------------------------------------------------------------
   !! * Executable part

   ! exchange erosion and settling fluxes
!$OMP DO SCHEDULE(RUNTIME)

# if defined MUSTANG && defined MUSTANG_CORFLUX 

      DO j=jfirst-1,jlast+1
      DO i=ifirst,ilast
        DO iv=1,nvp
            CORFLUY_SAND(i,j+1,IV_HOSTMODEL)=corfluy(iv,i,j)
        ENDDO
      ENDDO
      ENDDO

      DO j=jfirst,jlast
      DO i=ifirst-1,ilast+1
        DO iv=1,nvp
            CORFLUX_SAND(i+1,j,IV_HOSTMODEL)=corflux(iv,i,j)
        ENDDO
      ENDDO
      ENDDO
# endif

      DO j=jfirst,jlast
      DO i=ifirst,ilast

        DO iv=1,nvp
            SETTL_FLUX_w2s(i,j,IV_HOSTMODEL)=flx_w2s(iv,i,j)
        ENDDO
        DO iv=1,nv_adv
            EROS_FLUX_s2w(i,j,IV_HOSTMODEL)=flx_s2w(iv,i,j)
        ENDDO
        ! temperature
        EROS_FLUX_s2w(i,j,ITEMP_HOSTMODEL)=flx_s2w(-1,i,j)
        ! salinity
        EROS_FLUX_s2w(i,j,ISAL_HOSTMODEL)=flx_s2w(0,i,j)

        ! no transfer of SETTL_FLUX_w2s_TEMP et SAL and for dissolved subst. because they are merged in EROS_FLUX_s2w
        ! for dissolved variables (EROS_FLUX_s2w=erosion-settling+consolidation-diffusion) 
      ENDDO
      ENDDO
!$OMP END DO


  END SUBROUTINE coupl_MUSTANG2hydro    
#endif

   !!==============================================================================
#endif

END MODULE coupleur_MUSTANG

