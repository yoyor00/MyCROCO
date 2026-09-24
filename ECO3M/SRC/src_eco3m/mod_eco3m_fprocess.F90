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
    ! This file is therefore one of the few files that the user can modify.
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
    !----------------------------------------------------------------------------
    !-----------------------------UPTAKE DE NUTRIMENTS---------------------------
    !----------------------------------------------------------------------------
    !FONCTION DE MONOD:
#include "F_PROCESS/UPT/f_upt_monod_cell.inc"
    !----------------------------------------------------------------------------
    !-----------------------INHIBITION UPTAKE PAR UN NUTRIMENT-------------------
#include "F_PROCESS/UPT/f_inhib_Frost_Franzen.inc"
#include "F_PROCESS/UPT/f_cp_inhib_Frost_Franzen.inc"
    !----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !------------------------------MORTALITE NATURELLE----------------------------
    !-----------------------------------------------------------------------------
    !FONCTION DE MORTALITE:
#include "F_PROCESS/MORT/f_mort_quadratique_cell.inc"
#include "F_PROCESS/MORT/f_mort_linear_cell.inc"
#include "F_PROCESS/MORT/f_mort_linear_implicit.inc"
    !-----------------------------------------------------------------------------
    !----------------------LIMITATION PAR UN QUOTA INTRACELLULAIRE----------------
    !-----------------------------------------------------------------------------
    !FONCTION DE GEIDER:
#include "F_PROCESS/FQUOTA/f_fQPPX_geid_2020.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPX_geid_2021.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPX_geid_2019.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPX_geid_split_2018_hind3D.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPX_geid_split_2018_ssbug.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPX_geid_split_2019.inc"
#include "F_PROCESS/FQUOTA/f_fQPPX_geid_split_2020.inc"
!#include "F_PROCESS/FQUOTA/f_fQRESP_X_geid.inc"
#include "F_PROCESS/FQUOTA/f_fQRESPX_geid_2020.inc"
!#include "F_PROCESS/FQUOTA/f_fQRESPX_geid_2021.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPC_geid_2018_hind3D.inc"
#include "F_PROCESS/FQUOTA/f_fQPPC_geid_2020.inc"
!#include "F_PROCESS/FQUOTA/f_fQRESPC_geid_2018_hind3D.inc"
!#include "F_PROCESS/FQUOTA/f_fQPPC_geid_2019_new.inc"
!#include "F_PROCESS/FQUOTA/f_fQRESP_C_geid.inc"
!#include "F_PROCESS/FQUOTA/f_fQRESPC_geid_2019_new.inc"
#include "F_PROCESS/FQUOTA/f_fQRESPC_geid_2020.inc"
    !-----------------------------------------------------------------------------
    !--------------------------------PHOTOSYNTHESE--------------------------------
    !-----------------------------------------------------------------------------
    !FONCTION DE HAN:
#include "F_PROCESS/PPB/f_ppb_han.inc"
    !-----------------------------------------------------------------------------
    !---------------------------------CROISSANCE----------------------------------
    !-----------------------------------------------------------------------------
    !FONCTION DE CROISSANCE:
!#include "F_PROCESS/GROWTH/f_growth_cell_droop.inc"
#include "F_PROCESS/GROWTH/f_growth_cell_droop_2020.inc"
#include "F_PROCESS/GROWTH/f_growth_cell_droop_conv.inc"
!#include "F_PROCESS/GROWTH/f_grow_cell_droop_conv_NP.inc"
    !-----------------------------------------------------------------------------
    !-------------------- LIMITATION PAR DES QUOTA INTRACELLULAIRES---------------
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !----------------------NITRIFICATION------------------------------------------
#include "F_PROCESS/NITRIF/f_nitrif_linear.inc"
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !---------------------------------RESPIRATION---------------------------------
    !-----------------------------------------------------------------------------
    !FONCTION DE THINGSTAD(COUT DE LA CROISSANCE):
#include "F_PROCESS/RESP/f_resp.inc"
#include "F_PROCESS/RESP/f_resp_maint.inc"
    !-----------------------------------------------------------------------------
    !----------------------PRODUCTION DE CHLOROPHYLLE-----------------------------
    !-----------------------------------------------------------------------------
!#include "F_PROCESS/PChl/f_lchl_bak06_2018.inc"
!#include "F_PROCESS/PChl/f_lchl_bak06_2019.inc"
#include "F_PROCESS/PChl/f_lchl_bak06_2020.inc"
#include "F_PROCESS/PChl/f_lchl_bak06_2021.inc"
!#include "F_PROCESS/PChl/f_Chl_bak2021.inc"
!#include "F_PROCESS/PChl/f_lchl_bak06_2019b.inc"
!#include "F_PROCESS/FQUOTA/f_fQChl.inc"
!#include "F_PROCESS/FQUOTA/f_fQChl_2019.inc"
#include "F_PROCESS/FQUOTA/f_fQChl_2020.inc"
    !-----------------------------------------------------------------------------
    !------------------------FONCTION DE DIAZOTROPHIE-----------------------------
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !-------------------------REMINERALISATION DE LA MO---------------------------
    !-----------------------------------------------------------------------------
#include "F_PROCESS/REMIN/f_hydrolysis.inc"
#include "F_PROCESS/REMIN/f_hydrolysis_XC.inc"
!#include "F_PROCESS/REMIN/f_rem_DOM.inc"
#include "F_PROCESS/REMIN/f_rem_DOM_baretta95.inc"
    !-----------------------------------------------------------------------------
    !--------------------------PREDATION DU ZOOPLANCTON---------------------------
    !-----------------------------------------------------------------------------
#include "F_PROCESS/GRAZ/M_PROIES/f_graz_KOOIJ_KTW_cell.inc"
#include "F_PROCESS/GRAZ/M_PROIES/f_graz_hol2_mpreys_X.inc"
    !-----------------------------------------------------------------------------
#include "F_PROCESS/LIGHT/ATTENUATION/f_extinc_morel88.inc"
    !-----------------------------------------------------------------------------
END MODULE mod_eco3m_fprocess
