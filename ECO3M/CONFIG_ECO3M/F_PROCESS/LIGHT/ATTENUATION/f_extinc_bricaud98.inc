!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
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
!--------------------------------------------------------------------
!
       function f_extinc_bricaud98(Kw,Chl_tot,dz)
!    
!-------------------------------------------------------------------------
! This function calculates the attenuation coefficient of the incident light by diffusion
! according to Bricaud''s (1998) model. This model takes into account chlorophyll as
! well as CDOM absorption and also includes dispersion effects (scattering)
!
! Correlation for Kbio compiles 1166 spectra relative to different oceanic case 1 waters 
!
! Important: this correlation is rigorously valid only for the wavelength lambda =  440 nm
! ----------
!
! Variables and parameters and their units:
!-----------------------------------------
!  Kw  : diffuse attenuation coefficient in pure marine water, m-1
!  Kbio: diffuse attenuation coefficient due to biological particles, m-1
!
! Ref:  Bricaud, A. et al. (1998) J. Geophys. Res., 103(C13):31,033-31,044
!
! Melika Baklouti 08/07/2008
!-------------------------------------------------------------------------
!
! USE COUPLEUR_PHY_BIO
! USE VAR_GLOBAL
! USE MOD_FCHAIN

 Implicit None

!-- function:
 Real(8):: f_extinc_bricaud98

!-- Arguments
 real(8) :: Kw  
 real(8),optional :: Chl_tot,dz

!-- Local variables:
 real(8) :: Kd, Kbio

!-- Diffusion attenuation coefficient due to biological particles
!   for the wavelength 440 nm. 
!   In the following relation Chl_tot is in g m^-3 :
!
    Kbio = 0.0378 * (Chl_tot*1.d3)**0.627
   
!-- Diffusion attenuation coefficient :
   Kd = Kw + Kbio
 
!-- Attenuation term when Chl is in g m-3
   f_extinc_bricaud98 = exp(-Kd * dz) 
 
end function f_extinc_bricaud98
!--------------------------------------------------------------------
