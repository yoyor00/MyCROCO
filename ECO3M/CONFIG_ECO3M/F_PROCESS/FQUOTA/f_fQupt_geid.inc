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
!---------------------------------------------------------------------------
!
       function f_fQupt_geid(ili,jcol,indX,Qmin,Qmax,n,valmin)
!    
!---------------------------------------------------------------------------
! Cette fonction calcule le facteur adimensionnel (fonction
! de quotas) qui limite l'ASSIMILATION du nutriment X
! dont la conc. interne est donnee par la variable d''indice indX
!
! La fonction de quota est donnée par le modèle de Geider et al. (1998) 
!
! ATTENTION:  On laisse cette fonction délibérément devenir négative pour 
!             traduire les situations où l'organisme évacue plus de X qu'il n'en  
!             n'absorbe. Par ailleurs, cette fonction ne peut excéder la valeur 1
!            (cf Baklouti et al., 2006 pour plus de détails sur les modèles d'exudation)
!
! NB: cette fonction est plutot adaptée au nitrate mais pas
!     vraiment aux autres nutriments (cf Geider et al., 1998)
!
! Variables et parametres et leurs unites:
!-----------------------------------------
! indX: indice de la conc. interne en nutriment X dont la disponibilité est limitante
! Qmin: valeur minimale du quota intracellulaire de X/C (mol X (mol C)^-1)
! Qmax: valeur maximale du quota intracellulaire de X/C (mol X (mol C)^-1)
! n   : facteur de forme de courbe (=0.06 dans Geider et al., 1998)
!
! Ref:  Geider, R. J. et al.(1998) Limnol. Oceanogr., 43:679-694
!       Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
!
! Melika Baklouti 04/09/2007
!-------------------------------------------------------------------------
 Implicit None

!-- fonction:
 Real(8):: f_fQupt_geid(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: Qmin, Qmax,indX,n
 Real(8),optional :: valmin
 Integer :: ili,jcol

!--variables locales:
 Integer :: iiC,i,j,k
 Real(8) :: Q,hQ
 Integer :: iorg_X,iindX
 Character(L_VAR)::chain
 
 iindX = nint(indX)
 
 !-- indice de l'organisme :
 iorg_X = var(iindx)%idorg

 chain = 'C'
 !-- conc en C de l'organisme:
 iiC = f_idorg2id(iorg_X,chain)

!-- Calcul du quota cellulaire:
do i=nx_min,nx_max
  do j=ny_min,ny_max
    do k=1,nz_max
    Q = var(iindX)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-15)! quota cellulaire X:C
    hQ = ((Qmax-Q) / (Qmax - Qmin))
    if (hQ >= 0.d0) then 
         f_fQupt_geid(i,j,k) =  min(1.d0,hQ**n)
    else
         f_fQupt_geid(i,j,k) =  - (abs(hQ))**n
    endif
    if (present(valmin)) f_fQupt_geid(i,j,k) = max(valmin,f_fQupt_geid(i,j,k))
   enddo
  enddo
enddo

end function f_fQupt_geid
!--------------------------------------------------------------------
