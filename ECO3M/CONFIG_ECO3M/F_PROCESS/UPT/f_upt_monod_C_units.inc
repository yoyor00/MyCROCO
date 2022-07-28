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
    function f_upt_monod_C_units(ili,jcol,inut,k_demi,Vmax)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule l''absorption d''un nutriment d''indice
! inut par le modele de Monod-Michaelis-Menten 
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! inut : parametre donnant l'indice de la concentration
!         intracellulaire en element X (lorsqu'il differe de ili), où
!         X represente le nutriment (ex: X=N ou P,...)
! k_demi: cste de demi saturation, (mol X)/m3
! Vmax : vitesse CARBONE-SPECIFIQUE d'absorption max., en  mol X (mol C)^{-1} s^{-1}
!
! f_upt_monod_C_units : vitesse d''absorption, en (mol X).m-3 s^{-1}
!
! ATTENTION: 
!    * si ili ne représente pas la conc. intracellulaire en élément X 
!      (X étant l'élément limitant absorbé), il faut fournir cette concentration
!      via l'indice inut 
!
! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
!      Michaelis, L. & Menten, M. L. (1913), 49:333-369
!
! Melika Baklouti 6/07/2007
!--------------------------------------------------------------------
 Implicit None

 Real(8) :: f_upt_monod_C_units(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: inut,k_demi,Vmax
 Integer :: ili,jcol

!-- variables locales
 Integer :: iinut,iorg_plancton,iC,id_plancton_X
 Character(L_VAR)    :: chain

!-- organisme de la biomasse :
 iorg_plancton = var(ili)%idorg 
 chain = 'C'
 
!-- indice de la concentration en carbone de la biomasse
 iC = f_idorg2id(iorg_plancton,chain)

!-- indice de la concentration en element nutritif: 
   iinut=nint(inut)

   f_upt_monod_C_units =  Vmax * var(iinut)%conc/(var(iinut)%conc + k_demi) * var(iC)%conc


end function f_upt_monod_C_units
!--------------------------------------------------------------------
