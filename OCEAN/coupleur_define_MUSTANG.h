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

#if defined MUSTANG

   !&E==========================================================================
   !&E                   ***  coupleur_define_MUSTANG  ***
   !&E
   !&E
   !&E ** Purpose : definitions of dimensions, variables and parameters 
   !&E               for MUSTANG         
   !&E 
   !&E ** Description : must be completed by the user
   !&E          when coupling with a hydrodynamic model
   !&E
   !&E ** History :
   !&E     ! 2018-11  (B.Thouvenin )    : creation for portability and 
   !&E                                    coupling with any hydro model 
   !&E
   !&E==========================================================================

# define sed_MUSTANG_HOST sed_MUSTANG_CROCO

!/* some specific commands used in MARS are to be commented out automatically 
!   for another hydro model (we use  define)
!   + There may be other definitions to add in order to replace compatible 
!      terms with the hydro coupled model.
!*/
# define CALL_MPI !call_MPI
# define IF_MPI !if
# define ENDIF_MPI !endif
# define PRINT_DBG !print
!/*CROCO */
# define iscreenlog stdout
# define ierrorlog stdout
# define iwarnlog stdout
# define NAME_SUBS vname(1,indxT+ntrc_salt+isubs)

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!/*         Directory where are namelists files
!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
# define REPFICNAMELIST 'FIC_NAMELIST'

!/* Spatial Grid limits definition  of loops
!/*   inside the domain - except meshes at open boundaries
!*/
# define IMIN_GRID 1
# define IMAX_GRID Lm
# define JMIN_GRID 1
# define JMAX_GRID Mm

# define IMIN_BOUCL 2
/*# define IMAX_BOUCL 202  cas TFLAT2DV avec LLm0=200 (param.h*/
# define IMAX_BOUCL 52
# define JMIN_BOUCL 2
# define JMAX_BOUCL 4

!/* dimensions table definition 
!*/
# define PROC_IN_ARRAY       GLOBAL_2D_ARRAY   
# define PROC_IN_ARRAY_m1p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m1p1  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_m2p2  GLOBAL_2D_ARRAY
# define PROC_IN_ARRAY_0p1   GLOBAL_2D_ARRAY

!/* dimensions of variables in hydro modele
!*/
# define ARRAY_EROS_FLUX_s2w PROC_IN_ARRAY,1:NT
# define ARRAY_SETTL_FLUX_w2s PROC_IN_ARRAY,1:NT
# define ARRAY_SETTL_FLUXSUM_w2s PROC_IN_ARRAY,1:NT
# define ARRAY_CORFLUX_SAND PROC_IN_ARRAY_m1p1,itsubs1:itsubs2
# define ARRAY_CORFLUY_SAND PROC_IN_ARRAY_m1p1,itsubs1:itsubs2
# define ARRAY_WATER_FLUX_INPUTS PROC_IN_ARRAY,1:NB_LAYER_WAT
# define ARRAY_BATHY_H0 GLOBAL_2D_ARRAY
# define ARRAY_WATER_ELEVATION GLOBAL_2D_ARRAY,1:4
# define ARRAY_VELOCITY_U GLOBAL_2D_ARRAY,1:4
# define ARRAY_VELOCITY_V GLOBAL_2D_ARRAY,1:4
# define ARRAY_CELL_SURF GLOBAL_2D_ARRAY
# define ARRAY_CELL_DX GLOBAL_2D_ARRAY
# define ARRAY_CELL_DY GLOBAL_2D_ARRAY
!/* dimensions of variables in hydro modele !*/
# define ARRAY_WAT_SETTL GLOBAL_2D_ARRAY,N,itsubs1:itsubs2
# define ARRAY_WATER_CONC GLOBAL_2D_ARRAY,N,3,NT
!# define ARRAY_TEMPSAL GLOBAL_2D_ARRAY,N 
# define ARRAY_DHSED GLOBAL_2D_ARRAY 
# define ARRAY_FROFON GLOBAL_2D_ARRAY 
# define ARRAY_Z0HYDRO PROC_IN_ARRAY 

#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
# define ARRAY_LATLON GLOBAL_2D_ARRAY
#endif

#define  ARRAY_morpho GLOBAL_2D_ARRAY
#define  ARRAY_h0_bedrock GLOBAL_2D_ARRAY


!/* general variable hydro , bathy, time ... defined in hydro model but using by MUSTANG
!*/
# define NUMBER_PI pi
# define NB_LAYER_WAT N
# define COORD_SIGMA sc_r
# define BATHY_H0 h
# ifdef WET_DRY
# define RESIDUAL_THICKNESS_WAT D_wetdry
# else
# define RESIDUAL_THICKNESS_WAT 0.
# endif
# define WATER_ELEVATION zeta
# define CELL_DX om_r
# define CELL_DY on_r
# define CELL_SURF surf_cell
/*# define CELL_SURF cell_surf_MUSTANG*/
# define BAROTROP_VELOCITY_U ubar
# define BAROTROP_VELOCITY_V vbar
# define TIME_STEP dt   /* in MARS :  time step declared in rlg and therefore also in MUSTANG */
# define TRANSPORT_TIME_STEP dt /* in MARS : solving equations every half time step (in rlg)*/
# define CURRENT_TIME time
# define RHOREF rho0
# define TEMPREF_LIN T0
# define SALREF_LIN S0
# define GRAVITY g
# define BOTTOM_THICK_LAYER epn_bottom
# define WAT_SETTL ws_part
# define Z0HYDRO zob
# define WATER_CONCENTRATION t  /* water concentration in hydro model (=cv_wat in MARS)*/
# define DHSED dh  
#if defined key_MUSTANG_V2 && defined key_MUSTANG_debug
# define LATITUDE latr
# define LONGITUDE lonr
#endif
#if defined key_MUSTANG_flocmod
# define WAT_CONC_ALLMUD_ijk t(i,j,k,2+imud1:nvpc) /* water concentration in hydro model (for all mud) ATTENTION order of indices */
#endif

!/* surface elevation (i,j) et courants pouvant avoir un nombre different de dimension 
!*/
# define SURF_ELEVATION_ij WATER_ELEVATION(i,j,3) 
# define COURANTU_ij BAROTROP_VELOCITY_U(i+1,j,3) 
# define COURANTV_ij BAROTROP_VELOCITY_V(i,j+1,3) 

# define COURANTV_ip1jm1 BAROTROP_VELOCITY_V(i+1,j,3) 
# define COURANTV_ip1j   BAROTROP_VELOCITY_V(i+1,j+1,3) 
# define COURANTV_im1jm1 BAROTROP_VELOCITY_V(i-1,j,3)
# define COURANTV_im1j   BAROTROP_VELOCITY_V(i-1,j+1,3)

# define COURANTU_im1jp1 BAROTROP_VELOCITY_U(i,j+1,3) 
# define COURANTU_ijp1   BAROTROP_VELOCITY_U(i+1,j+1,3) 
# define COURANTU_im1jm1 BAROTROP_VELOCITY_U(i,j-1,3) 
# define COURANTU_ijm1   BAROTROP_VELOCITY_U(i+1,j-1,3) 


!/* name of fluxes exchange between MUSTANG and hydro model 
!*/
# define EROS_FLUX_s2w flx_s2w_CROCO /* Erosion flux from sediment to water  */

/* if TEMP and SAL not inclued in other concentrations table
!#if ! defined key_MARS  
!# define EROS_FLUX_TEMP_s2w ??
!# define EROS_FLUX_SAL_s2w ??
!#endif
*/

# define SETTL_FLUX_w2s flx_w2s_CROCO  /* Tendance Flux de depot eau vers sediment (=flx_s2w in MARS) only for particulate */
# define CORFLUX_SAND corflux_CROCO /* correction transport Flux  for sands in hydro module (=corflux (iv,i,j))*/
# define CORFLUY_SAND corfluy_CROCO
# define SETTL_FLUXSUM_w2s flx_w2s_sum_CROCO /* effective deposit Flux (sum) from water to sediment in hydro model (=flx_w2s_sum in MARS)*/

/* to locate the number of variables simulated by MUSTANG in the host hydro model (used in coupleur to tranfer exchange arrays  */      
# define IV_HOSTMODEL itsubs1+iv-1
# define ITEMP_HOSTMODEL itemp
# define ISAL_HOSTMODEL itemp+1

# define WATER_FLUX_INPUT_BOTCELL phieau_CROCO(:,:,1) /* Flux d eau apporte dans la maille de fond in hydro model (=phieau in MARS) */



/* Lateral Erosion
!   in neighboring cells (could depend on grid architecture) 
!*/

# define HTOT_NEAR_E htot(i+1,j)
# define HTOT_NEAR_W htot(i-1,j)
# define HTOT_NEAR_N htot(i,j+1)
# define HTOT_NEAR_S htot(i,j-1)
# define SURF_NEAR_E CELL_SURF(i+1,j)
# define SURF_NEAR_W CELL_SURF(i-1,j)
# define SURF_NEAR_N CELL_SURF(i,j+1)
# define SURF_NEAR_S CELL_SURF(i,j-1)
# define V_NEAR_E (COURANTV_ip1jm1 + COURANTV_ip1j)
# define V_NEAR_W (COURANTV_im1jm1 + COURANTV_im1j)
# define U_NEAR_N (COURANTU_im1jp1 + COURANTU_ijp1)
# define U_NEAR_S (COURANTU_im1jm1 + COURANTU_ijm1)

!/* sliding proces of fluid mud
!   slope of the bottom in the neighboring cells (could depend on grid architecture)
!*/
!/*
# define SLOPE_W ((h0(i-1,j)-h0(i,j))/dx(i,j))
# define SLOPE_E ((h0(i+1,j)-h0(i,j))/dx(i,j))
# define SLOPE_N ((h0(i,j+1)-h0(i,j))/dy(i,j))
# define SLOPE_S ((h0(i,j-1)-h0(i,j))/dy(i,j))
!*/
# define SLOPE_W ((BATHY_H0(i-1,j)-BATHY_H0(i,j))/CELL_DX(i,j))
# define SLOPE_E ((BATHY_H0(i+1,j)-BATHY_H0(i,j))/CELL_DX(i,j))
# define SLOPE_N ((BATHY_H0(i,j+1)-BATHY_H0(i,j))/CELL_DY(i,j))
# define SLOPE_S ((BATHY_H0(i,j-1)-BATHY_H0(i,j))/CELL_DY(i,j))

#endif

