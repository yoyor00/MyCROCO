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
function f_hydrolysis(ili,jcol,indX,TT,indBACC,QminC,QmaxC)
    !
    !---------------------------------------------------------------------------
    ! Cette fonction calcule la vitesse d'hydrolyse de la MOP detritique par les bacteries
    ! a partir d'un temps de turnover fixe
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indX  : indice de la conc. en MOP (POC,PON,...)
    ! indBACC: indice de la biomasse en C des bacteries
    ! TT = temps de turnover en jours (cette valeur peut devenir extremement elevee 
    ! pour la matiere refractaire)
    !
    !
    ! R. MAURIAC  07/09/2009
    ! Derniere revision:  M. Baklouti 29/01/2018
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_hydrolysis(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: TT,indX,indBACC,QminC,QmaxC
    Integer :: ili,jcol

    !--variables locales:
    Integer :: i,j,k,iiMOP,iiCell,iiMOPC,iiC,identite
    real(8) :: nn
    real(8) :: TTsec,fQ,Q
    character(L_VAR) :: chain,chain2,chain3

    !-- Indice de la concentration en MOP
    iiMOP = nint(indX)
    iiC = nint(indBACC)

    ! Indice de  la concentration en cellule des bacteries
    chain = 'cell'
    identite = var(iiC)%idorg
    iiCell = f_idorg2id(identite,chain)

    ! Indice de  la concentration en C de la matiere particulaire detritique
    chain = 'C'
    identite = var(iiMOP)%idorg
    iiMOPC = f_idorg2id(identite,chain)

    nn = 0.2!0.06!0.6 ! 1.0 !0.6!0.06

    ! -- Calcul du temps de turnover en second
    TTsec = TT * 24 * 3600

    !-- Calcul de la fonction de remineralisation:
 do k=1,nz_max
   do j=ny_min,ny_max
     do i=nx_min,nx_max
        Q = var(iiC)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
   ! Carbon quota function for bacteria
        fQ = ((QmaxC - Q)/(QmaxC-QminC))
        if (fQ >= 0.0e0) fQ = fQ**nn
                
        fQ = max(0.0e0,fQ)
        fQ = min(1.0e0,fQ)
        f_hydrolysis(i,j,k) = 0.0d0
!-- seuil min en bacteries                
        if (var(iiCell)%conc(i,j,k) .gt. 1000   &
                .and. (var(iiMOPC)%conc(i,j,k) > (var(iiMOPC)%conc(i,j,k)/TTsec)*fQ *dt_bio)) then
            f_hydrolysis(i,j,k) = (var(iiMOP)%conc(i,j,k)/TTsec)*fQ   
        endif
     enddo
   enddo
 enddo

end function f_hydrolysis
!
