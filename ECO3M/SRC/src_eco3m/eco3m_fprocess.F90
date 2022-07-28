!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & V#incent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
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
!
!---------------------------------------------------------------------------
!
               MODULE mod_eco3m_fprocess
 !
 !---------------------------------------------------------------------------
 !
 ! This module is used to include all the process functions from the F_PROCESS
 ! library used in the current configuration of the biogeochemical model. 
 !
 ! It must also include the name of the light extinction function used when 
 ! the RGB key is not activated
 !
 !! \author Melika Baklouti, Vincent Faure
 !  \date 2007-07-06
 !----------------------------------------------------------------------------

    use mod_eco3m
    use mod_eco3m_user
    use mod_eco3m_id_extract
    implicit none

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------
! --------- Gross Primary production (PPB) models      -----------
!-----------------------------------------------------------------

!-- Han (2002):
#include "F_PROCESS/PPB/f_ppb_han.inc"
!-- Geider et al. (1998):
!#include "F_PROCESS/PPB/f_ppb_geid98.inc"
!
!----------------------------------------------------------------------
!---- Quota functions of growth limitation by a single nutrient ----
!----------------------------------------------------------------------

!-- Caperon & Meyer (1972):
#include "F_PROCESS/FQUOTA/f_fQPP_capmey.inc"
!-- Droop (1968):
#include "F_PROCESS/FQUOTA/f_fQPP_droop.inc"
!-- Geider et al. (1998):
#include "F_PROCESS/FQUOTA/f_fQPP_flynn.inc"
!-- Flynn (2001):
#include "F_PROCESS/FQUOTA/f_fQPP_geid.inc"

!----------------------------------------------------------------------------------------------------
!---- Quota functions of growth limitation by several nutrients (2 or more limiting nutrients) ------
!----------------------------------------------------------------------------------------------------
!
! -- Case where we the minimum is calculated on the fQ quota function:
!-- Caperon & Meyer (1972):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_capmey_min.inc"
!-- Droop (1968):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_droop_min.inc"
!-- Geider et al. (1998):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_flynn_min.inc"
!-- Flynn (2001):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_geid_min.inc"
!
! -- Case where the minimum is calculated on the Q/Qmax ratio:
!-- Caperon & Meyer (1972):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_capmey_min2.inc"
!-- Droop (1968):
!#include "F_PROCESS/FQUOTA/MULTI_LIM/f_fQPP_droop_min2.inc"
!
!---------------------------------------------------------------------------------
!---- Photoacclimatation models for the Chl dynamics (Chl is a state variable) ---
!---------------------------------------------------------------------------------
!
! Geider et al. (1998)
!#include "F_PROCESS/PChl/f_pchl_geid98.inc"
! Baklouti et al. (2006)  
#include "F_PROCESS/PChl/f_pchl_bak06.inc" 
!
!------------------------------------------------------------------------
!---- Correlations for the Chl:C ratio (Chl is not a state variable) ----
!------------------------------------------------------------------------
!
! Cloern et al. (19f_ChlC_cloern_95.inc)
#include "F_PROCESS/Chl_C/f_ChlC_cloern_95.inc"
! Smith and Tett (2000)  
!#include "F_PROCESS/Chl_C/f_ChlC_smith_tett_00.inc" 
!
!-----------------------------------------------------------------
!----     Bacterial production models (PB)                --------
!-----------------------------------------------------------------

! uptake de COD selon Monod-Michaelis-Menten: 
!-----------------------------------------
! utiliser les modeles d''uptake f_upt_monod*.inc fournis dans le 
! repertoire../../CONFIG_ECO3M/F_PROCESS/UPT/
#include "F_PROCESS/PBACT/f_gbp_monod.inc"
!-----------------------------------------------------------------
!---------   Respiration models   -------------------------------
!-----------------------------------------------------------------
#include "F_PROCESS/RESP/f_respB_vichi07.inc"
!-----------------------------------------------------------------
!-------   Release models :    excretion, exudation...  ----------
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
!---------- Uptake models     ------------------------------------
!-----------------------------------------------------------------
!
!-- Monod-Michaelis-Menten (Vmax in s-1)
#include "F_PROCESS/UPT/f_upt_monod.inc"
!-- Monod-Michaelis-Menten (Vmax in molN molC-1 s-1)
#include "F_PROCESS/UPT/f_upt_monod_C_units.inc"
!-- Monod-Michaelis-Menten (Vmax calculated by PBmaxQmax)
#include "F_PROCESS/UPT/f_upt_monod_Qmax.inc"
!
!-- Vichi et al.(2007) model for bacterial uptake/exudation of nutrients:
#include "F_PROCESS/UPT/f_uptB_nut_vichi07.inc"
!-- inhibition d''uptake nitrate par ammonium ---
#include "F_PROCESS/UPT/f_inhib_harrison.inc"
!--------------------------------------------------
!---- Quota functions limiting nutrient uptake ----
!--------------------------------------------------
!
!-- Flynn (2003):
#include "F_PROCESS/FQUOTA/f_fQupt_flynn.inc"
!-- Geider et al. (1998):
#include "F_PROCESS/FQUOTA/f_fQupt_geid.inc"
!-- Lehman et al. (1975):
#include "F_PROCESS/FQUOTA/f_fQupt_lehman.inc"
!
!---------------------------------------------------
! ----- Grazing models : case of a single prey  ----
!---------------------------------------------------
!
!-- Holling II model 
!#include "F_PROCESS/GRAZ/SINGL_PREY/f_graz_C_holling2.inc"
!
!--------------------------------------------------
! ----- Grazing models : case of several preys ---- 
!--------------------------------------------------
!
!-- Holling II model with Chesson''s (1983) preferences
!#include "F_PROCESS/GRAZ/MULT_PREYS/f_graz_C_hol2_chesson.inc"
!
!--Holling II model with Murdoch''s (1973) preferences
!#include "F_PROCESS/GRAZ/MULT_PREYS/f_graz_C_hol2_murdoch.inc"
!
!-- Evans () threshold model
!#include "F_PROCESS/GRAZ/MULT_PREYS/f_graz_C_thres_evans.inc"
!
!-- Fasham model
#include "F_PROCESS/GRAZ/MULT_PREYS/f_graz_C_Fasham90.inc"
#include "F_PROCESS/GRAZ/MULT_PREYS/f_graz_C2X.inc"
!
!--------------------------------------------------------------
! --------- Mortality models             ----------------------
!--------------------------------------------------------------
!
!-- Linear mortality law:
#include "F_PROCESS/MORT/f_mort_linear.inc"
!-- Quadratic mortality law :
#include "F_PROCESS/MORT/f_mort_quadratique.inc"
!
!-----------------------------------------------------------------
! --------- Mineralization  models       -------------------------
!-----------------------------------------------------------------
#include "F_PROCESS/REMIN/f_rem.inc"
!-----------------------------------------------------------------
! ---------  Nitrification  models       -------------------------
!-----------------------------------------------------------------
!
!-- Adapted from Vichi et al.(2007) 
#include "F_PROCESS/NITRIF/f_nitrif_vichi07diaz.inc"
!
!-----------------------------------------------------------------
!---- Mod¿le de retour du zooplancton: source de DET_C et DET_N --
!-----------------------------------------------------------------
!
!-----------------------------------------------------------------
!----------------- Temperature models ----------------------------
!-----------------------------------------------------------------
!-- Arrhenius (1889):
!#include "F_PROCESS/TEMP/f_tfunc_Arrhenius.inc"
!
!-- Q10 formulation:
#include "F_PROCESS/TEMP/f_tfunc_Q10.inc"
!
!-- Lacroix&Gregoire (2002)
#include "F_PROCESS/TEMP/f_tfunc_Lacroix02.inc"
!
!-----------------------------------------------------------------
!----- Light attenuation models in the water column           ----
!-----------------------------------------------------------------
!-- Morel (1988) model:
!#include "F_PROCESS/LIGHT/ATTENUATION/f_extinc_morel88.inc"
!-- Morel (1993) model:
#include "F_PROCESS/LIGHT/ATTENUATION/f_extinc_morel93.inc"
!-- Bricaud (1998) model :
!#include "F_PROCESS/LIGHT/ATTENUATION/f_extinc_bricaud98.inc"
!-- Levy (2005) model  
!#include "F_PROCESS/LIGHT/ATTENUATION/f_extinc_levy.inc"

END MODULE mod_eco3m_fprocess
