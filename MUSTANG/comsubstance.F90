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

 MODULE comsubstance
   
   !!======================================================================
   !!                   ***  MODULE comsubstance   ***
   !!
   !!  ** Purpose : declare all common variables related to substance
   !!
   !!  ** History :
   !!        ! 2018-08 (B. Thouvenin) adapted from MARS for key SUBSTANCE 
   !!======================================================================

#if defined SUBSTANCE

   implicit none

   ! default
   public

   !-------------------------------------------------------------------
   !  Definition of rsh, rlg, riosh, riolg, lchain
   !-------------------------------------------------------------------
   INTEGER,PARAMETER :: riosh = 8, riolg = 8, rlg = 8, rsh = 8
   INTEGER,PARAMETER :: lchain = 200
             
   ! ------------------------------------------------------------------
   ! Manage variable type ( sediment, particulate...)
   ! ------------------------------------------------------------------

      INTEGER           :: nv_grav  ! number of gravel variables
      INTEGER           :: nv_sand  ! number of sand variables
      INTEGER           :: nv_mud   ! number of mud variables
      INTEGER           :: nspb     ! **TODO : add each var description**
      INTEGER           :: igrav1
      INTEGER           :: igrav2
      INTEGER           :: isand1
      INTEGER           :: isand2
      INTEGER           :: imud1
      INTEGER           :: imud2
      INTEGER           :: nv_ncp
      INTEGER           :: nv_dis
      INTEGER           :: nv_sorb
      INTEGER           :: nv_fix
      INTEGER           :: nvp
      INTEGER           :: nvpc
      INTEGER           :: nv_adv
      INTEGER           :: nv_state
      INTEGER           :: nv_tot
      INTEGER           :: nv_bent

        
      INTEGER, DIMENSION(:), ALLOCATABLE           :: irkm_var_assoc ! indices des particulaires constitutives associees aux particulaires adsorbees, 

      ! ----------------------------------------------------------------------------
      ! GESTION CORRESPONDANCE ORDRE DES VARIABLES LUES ET VARIABLES DANS MODELE
      ! ----------------------------------------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE         :: irk_mod !rank (position) in the model of the variable in the array managing substances
      INTEGER, DIMENSION(:), ALLOCATABLE         :: irk_fil !rank (position) in the data file of the variable

      

   
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: cini_wat,cini_air,cobc_wat
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: typdiss
      REAL(KIND=rsh),DIMENSION(:)  ,ALLOCATABLE :: sub_flx_atm,cv_rain
      CHARACTER(LEN=lchain),DIMENSION(:) ,ALLOCATABLE :: obc_cv_name,name_var,standard_name_var
      CHARACTER(LEN=lchain),DIMENSION(:) ,ALLOCATABLE :: init_cv_name,unit_var
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE           :: ws_part
      REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE       :: ws_free_min, ws_free_max
      REAL(KIND=rsh),DIMENSION(:), PUBLIC  ,ALLOCATABLE        :: unit_modif_mudbio_N2dw
      LOGICAL, PUBLIC,DIMENSION(:),ALLOCATABLE                 :: l_subs2D

#ifdef MUSTANG
      REAL(KIND=rsh),DIMENSION(:,:)    ,ALLOCATABLE :: ws_free_para, ws_hind_para
      INTEGER,DIMENSION(:)   ,ALLOCATABLE :: ws_free_opt, ws_hind_opt
      REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE :: tocd      
      REAL(KIND=rsh), DIMENSION(:)     ,ALLOCATABLE :: cini_sed_r
      REAL(KIND=rsh),DIMENSION(:)      ,ALLOCATABLE :: diam_sed, ros
#if defined key_MUSTANG_V2 && defined key_MUSTANG_bedload
      INTEGER                          :: ibedload1, ibedload2
#endif
#ifdef key_sand2D
      LOGICAL, DIMENSION(:), ALLOCATABLE                :: l_outsandrouse
#endif

#endif

      ! ----------------------------------------------------------------------------
      ! Fixed substances in water column
      ! ----------------------------------------------------------------------------
      REAL(KIND=rsh),DIMENSION(:),ALLOCATABLE        :: cini_wat_fix
      LOGICAL,DIMENSION(:),ALLOCATABLE               :: l_out_subs_fix
      CHARACTER(LEN=lchain),DIMENSION(:),ALLOCATABLE :: init_cv_name_fix
      REAL(KIND=rsh),DIMENSION(:,:,:,:),ALLOCATABLE  :: cvfix_wat

#ifdef key_benthic
      ! ---------------------------------------------------------------------------
      ! Benthic variables
      ! ---------------------------------------------------------------------------
      CHARACTER(LEN = lchain)                            :: filespcbenthic 
      REAL(KIND = rsh), DIMENSION(:,:,:,:), ALLOCATABLE  :: cv_bent
      CHARACTER(LEN = lchain), DIMENSION(:), ALLOCATABLE :: name_var_bent,long_name_var_bent, &
                                                      standard_name_var_bent,unit_var_bent
      REAL(KIND = rsh), DIMENSION(:), ALLOCATABLE        :: cbentmax,hbentmin,hbentmax,cini_bent
      LOGICAL,DIMENSION(:), ALLOCATABLE                  :: l_bent_sedvert, l_bent_drive,l_bent_out,l_out_subs_bent
#endif /* ifdef key_benthic */

      ! ----------------------------------------------------------------------------
      ! Atmospheric fluxes for substance
      ! ----------------------------------------------------------------------------
      CHARACTER(LEN = lchain) :: file_flxatm_subs ! file name defined if atmospheric flux for substances
      LOGICAL                 :: l_subflxatm, l_subflxatm_xyt
      LOGICAL                 :: l_cvrain_readfile, l_subflxatm_readfile
      REAL(KIND = rsh)        :: sflx_sub_atm_depth

#if ! defined key_MARS
#if defined MUSTANG || defined BIOLink 
      ! ----------------------------------------------------------------------------
      !  Declaration and evaluation of surface cells if not known in hydro host model
      !    needing for MUSTANG and Bloom/oyster
      ! ----------------------------------------------------------------------------
      REAL(KIND = rsh), DIMENSION(:,:), ALLOCATABLE  :: surf_cell
#endif /* ifdef MUSTANG & BIOLink */
#endif /* ifndef MARS */

#endif /* ifdef SUBSTANCE */

   !!==============================================================================

 END MODULE comsubstance
