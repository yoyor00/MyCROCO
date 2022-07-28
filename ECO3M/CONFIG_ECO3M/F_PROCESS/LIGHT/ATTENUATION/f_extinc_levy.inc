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
!--------------------------------------------------------------------
!
       function f_extinc_levy(Chl_tot,dz,Kw_r,Kw_b,Kbio_r,Kbio_b,lr,lb)
!    
!-------------------------------------------------------------------------
! This function calculates light attenuation in the water column 
! according to Levy''s(2005) model that considers 2 wavelengths
! 
! Variables and parameters and their units:
!-----------------------------------------
!  Kw_r : absorption of pure water at  676 nm (rouge) en m-1 (Pope et Fry 1997 ), typical val. 0.45
!  Kw_b : absorption of pure water at 440 nm (bleu) en m-1 (Pope et Fry 1997 ), typical val. 0.00635
!  Kbio_r : absorption coefficient in red  (m2*mg-1) , typical val. 0.0181 (Bricaud et al. 2004)
!  Kbio_b : absorption coefficient in blue (m2*mg-1), typical val. 0.0654!0.0383 (Bricaud et al. 2004)
!  lr  : exponant for phytoplankton absorption in red (-), typical val. 0.776 (Bricaud et al. 2004)
!  lb  : exposant for phytoplankton absorption in blue (-), typical val. 0.728 or 0.651 (Bricaud et al. 2004)
!  Chl_tot : Chl. conc. in g/m3
!  dz      : distance between two nodes in the vertical direction, m
!
! Ref:  Lévy, M. et al. (2005) JGR, 110, C07S21:1-28
!
! Melika Baklouti 19/07/2007
!-------------------------------------------------------------------------
!
! USE COUPLEUR_PHY_BIO
! USE VAR_GLOBAL
! USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8):: f_extinc_levy

!-- arguments
 real(8) :: Chl_tot,dz
 real(8) :: Kw_r,Kw_b,Kbio_r,Kbio_b,lr,lb

 
!-- relation pour Chl exprime en mg m-3
   f_extinc_levy = 0.5*(exp(-dz*(Kw_r+ Kbio_r*(Chl_tot*1.d3)**lr))+ exp(-dz*(Kw_b+Kbio_b*(Chl_tot*1.d3)**lb)))    
! 
end function f_extinc_levy
!--------------------------------------------------------------------
