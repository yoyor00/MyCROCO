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
    function f_upt_monod_Qmax(ili,jcol,inut,k_demi,Qmax,mu_max)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule l''absorption/assimilation (uptake) d''un nutriment 
! d''indice inut par le modele de Monod-Michaelis-Menten en posant que Vmax = mu_max . Qmax
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! inut   : paramètre donnant l'indice de la concentration
!          en nutriment limitant  (lorsqu'il differe de jcol)
! k_demi : cste de demi saturation, (mol X)/m3
! Qmax   : quota intracellulaire maximum, en  mol X / (mol C)
! mu_max : argument facultatif qui représente:
!              - la vitesse SPECIFIQUE de production bactérienne max pour les bactéries (s^{-1})
!              - la vitesse SPECIFIQUE mu_PPB_NR de production primaire en absence de limitation par les nuts
!                (mais obtenue à l'irradiance courante E et non pas à Emax) pour le phytoplancton 
!                Cette vitesse est calculée dans la fonction donnant la production primaire à chaque pas de temps
!                de calcul de la biologie 
!         
!
! f_upt_monod_Qmax : vitesse d''absorption, en (mol X).m-3 s^{-1}
!
! ATTENTION: 
!    * ili doit être l'indice d'une des concentrations (C ou N ou P...)
!      de l'organisme consommant du nutriment
!    *  il faut fournir la concentration  du nutriment consommé via l'indice inut 
!
! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
!      Michaelis, L. & Menten, M. L. (1913), 49:333-369
!      Geider, R. et al. (1998) Limnol. & Oceanogr. 43(4):679-694
!
! Melika Baklouti 6/07/2007
!--------------------------------------------------------------------
 Implicit None

 Real(8) :: f_upt_monod_Qmax(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: k_demi,Qmax,inut
 Real(8),optional ::mu_max
 Integer :: ili,jcol

!-- variables locales
 Integer :: iinut,iorg_plancton,iC,i
 Real(8) :: Vmax(nx_min:nx_max,ny_min:ny_max,nz_max)
 Character(L_VAR)     :: chain

!-- organisme de la biomasse :
 iorg_plancton = var(ili)%idorg 

!-- indice de la variable contenant la conc. en nut:
 iinut = nint(inut)
  
!-- indice de la conc en carbone de la biomasse:
   chain = 'C'
   iC = f_idorg2id(iorg_plancton,chain)

 if (present (mu_max)) then
!-- Vitesse max d'absorption en s-1 x (mol X)/(mol C) =  molX/molC/s
   Vmax = mu_max * Qmax
 else
!-- Vitesse max d'absorption en s^-1 x (mol X)/(mol C) =  molX/m3/s
! NB: le vecteur mu_PPB_NR a ete calcule avant d'appeler cette fonction normalement
!
   do i=1,nscp_phy
     if (mu_PPB_NR(i)%idorg == iorg_plancton) then
        Vmax(:,:,:)= mu_PPB_NR(i)%val(:,:,:) * Qmax
        exit
     endif
   enddo
 endif

!-- vitesse d'absorption en (mol X) m-3 s-1 :
  f_upt_monod_Qmax =  Vmax * var(iinut)%conc/(var(iinut)%conc + k_demi) * var(iC)%conc

end function f_upt_monod_Qmax
!--------------------------------------------------------------------
