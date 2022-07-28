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
       function f_extinc_morel88(Kw,Chl_tot,dz)
!    
!-------------------------------------------------------------------------
! Cette fonction calcule le coefficient d''attenuation de la lumière descendante par diffusion
! selon le modèle de Morel (1988) qui tient compte de l''absorption par la chlorophylle a et la
! pheotinine a, et qui inclue également les effets de dispersion (scattering)
!
!
! La corrélation pour Kbio est intégrée sur le spectre des longueurs d''onde utilisables pour
! la photosynthèse
!
! Variables et parametres et leurs unites:
!-----------------------------------------
!  Kw  : coefficient d''attenuation par diffusion de l''eau de mer pure, m-1
!  Kbio: coefficient d''attenuation par diffusion du à autre chose que l''eau (et notamment
!         aux particules biologiques), m-1
!
! Attention: validité de la loi: pour (Chl_tot*1.d3) comprise entre 0.03 et 30 mgChla / m3
!----------  d''après Morel page 10,757
!
! Ref:  Morel, A. et al. (1988) J. Geophys. Res., 93(C9):10,749:10,768
!
! Melika Baklouti 19/07/2007
!-------------------------------------------------------------------------
!
! USE COUPLEUR_PHY_BIO
! USE VAR_GLOBAL
! USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8):: f_extinc_morel88

!-- arguments
 real(8) :: Kw  
 real(8),optional :: Chl_tot,dz

!-- variables locales:
 real(8) :: Kd, Kbio

!-- coefficient d''attenuation par diffusion du fait des particules biologiques:
   Kbio =  0.0518 * (Chl_tot*1.d3)**(1-0.572)  
  
!-- coefficient d''attenuation par diffusion :
   Kd = Kw + Kbio
 
!-- relation pour Chl exprime en mg m-3
   f_extinc_morel88 = exp(-Kd * dz) 
 
end function f_extinc_morel88
!--------------------------------------------------------------------
