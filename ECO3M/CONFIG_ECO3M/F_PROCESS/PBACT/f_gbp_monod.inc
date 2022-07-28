!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
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
    function f_gbp_monod(ili,jcol,inut,k_demi,mu_bact_NR)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule l''absorption/assimilation (uptake) d''un nutriment 
! d''indice inut par le modele de Monod-Michaelis-Menten 
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! k_demi: cste de demi saturation, (mol X)/m3
! mu_bact_NR : taux de croissance max en NR conditions., en  s^{-1}
! inut : parametre donnant l'indice de la concentration
!         en nutriment limitant  
!         
!
! f_gbp_monod : vitesse d''absorption, en (mol X).m-3 s^{-1}
!
! ATTENTION: 
!    * ili doit etre l'indice d'une des concentrations (C ou N ou P...)
!      de l'organisme consommant du nutriment
!    * il faut fournir l'indice inut de la concentration en nutriment limitant
!      car il peut differer de jcol 
!
! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
!      Michaelis, L. & Menten, M. L. (1913), 49:333-369
!
! F. Diaz 19/12/2007
!--------------------------------------------------------------------
 Implicit None

 Real(8) :: f_gbp_monod(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: inut,k_demi,mu_bact_NR
 Integer :: ili,jcol,ik

!-- variables locales
 Integer :: iinut,iorg_bact,iX
 Character(L_VAR)  :: elmt_nut

!-- Organisme de la biomasse :
 iorg_bact = var(ili)%idorg 

!-- indice de la conc. en nutriments:
 iinut=nint(inut)

!-- element X apporte par le nutriment 
 elmt_nut = trim(adjustl(VAR(iinut)%elmt))  

!--concentration en X du plancton 
 iX= f_idorg2id(iorg_bact,elmt_nut)

  f_gbp_monod(:,:,:) =  mu_bact_NR * var(iinut)%conc(:,:,:)/(var(iinut)%conc(:,:,:) + k_demi) * var(iX)%conc(:,:,:)
! do ik=1,nscp_bac
!   if (PBB_NR(ik)%idorg == iorg_bact) then
!    PBB_NR(ik)%val(:,:,:) = mu_bact_NR * var(iinut)%conc(:,:,:)/(var(iinut)%conc(:,:,:) + k_demi) * var(iX)%conc(:,:,:)
!
!    exit
!	 endif
!enddo
 
end function f_gbp_monod
