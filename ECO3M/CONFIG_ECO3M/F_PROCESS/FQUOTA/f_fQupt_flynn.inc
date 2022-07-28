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
       function f_fQupt_flynn(ili,jcol,indX,Qmin,Qmax,K_Q)
!    
!-------------------------------------------------------------------------
! Cette fonction calcule le facteur adimensionnel (fonction
! de quotas) qui limite l'ASSIMILATION du nutriment X
! dont la conc. interne  est donnee par  la variable d''indice indX
!
! La fonction de quota est donnée par le modèle de Flynn et al. (2003) 
!
! ATTENTION:  On laisse cette fonction délibérément devenir négative pour 
!             traduire les situations où l'organisme évacue plus de X qu'il n'en  
!             n'absorbe. Par ailleurs, cette fonction ne peut excéder la valeur 1
!            (cf Baklouti et al., 2006 pour plus de détails sur les modèles d'exudation)
!
!
! Variables et parametres et leurs unites:
!-----------------------------------------
! indX: indice de la conc. interne en nutriment X dont la disponibilité est limitante
! Qmin: valeur minimale du quota intracellulaire de X/C (mol X (mol C)^-1)
! Qmax: valeur maximale du quota intracellulaire de X/C (mol X (mol C)^-1)
! K_Q : paramètre de courbure de la fonction, sans dimension
!
!
! Ref:  Flynn, K. J. et al.(2003) Prog. Oceanogr., 56:249-279
!       Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
!
! Melika Baklouti 04/09/2007
!-------------------------------------------------------------------------
 Implicit None

!-- fonction:
 Real(8):: f_fQupt_flynn(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: Qmin,Qmax,indX,K_Q
 Integer :: ili,jcol,iindX

!--variables locales:
 Integer :: iiC,i,j,k
 Real(8) :: Q,hQ
 Integer:: identite
 Character(L_VAR)::chain

 iindX = nint(indX)
 identite = var(iindX)%idorg
 chain = 'C'
 iiC = f_idorg2id(identite,chain)
 
!-- Calcul du quota cellulaire:
do i=nx_min,nx_max
  do j=ny_min,ny_max
    do k=1,nz_max
    Q = var(iindX)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-15)! quota cellulaire X:C
    hQ =  (1 - Q/Qmax) / (1 - Q/Qmax + K_Q)

    f_fQupt_flynn(i,j,k) =  min(1.d0,hQ)
   enddo
  enddo
enddo

end function f_fQupt_flynn
