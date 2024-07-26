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
function f_hydrolysis_XC(ili,jcol,indX,TT,indBAC,Qmin,Qmax)
    !
    !---------------------------------------------------------------------------
    ! Cette fonction calcule la vitesse d'hydrolyse de la MOP detritique par les bacteries
    ! a partir d'un temps de turnover fixe
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indX: indice de la conc. en MOP (POC,PON,...)
    ! TT = temps de turnover en jours (cette valeur peut devenir extremement elevee
    ! pour la matiere refractaire)
    !
    !
    ! Refs: Polimene et al. 2006
    !
    ! R. MAURIAC  07/09/2009
    ! Derniere revision:  M. Baklouti 11/11/2016
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_hydrolysis_XC(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: TT,indX,indBAC,Qmin,Qmax
    Integer :: ili,jcol,iiX,identite

    !--variables locales:
    Integer :: i,j,k,iiMOP,iiCell,iiC
    real(8) :: nn
    real(8) :: TTsec,fQ,QC,XC,XCmin,XCmax
    character(L_VAR) :: chain,chain2,chain3

    !-- Indice de la concentration en MOP
    iiMOP = nint(indX)
    iiX = nint(indBAC)

    ! Indice de  la concentration en cellule des bacteries
    chain = 'cell'
    identite = var(iiX)%idorg
    iiCell = f_idorg2id(identite,chain)
    chain = 'C'
    iiC = f_idorg2id(identite,chain)

    nn = 0.6 ! pour sim173 au lieu de 0.3 jusqu'a present, car comparaison avec f_hydrol qui utilise nn=0.6
!-- MB: bug corrige le 05/07/2018
    XCmin =1.d0/Qmax
    XCmax = 1.d0/Qmin
!-- MB: bug corrige le 05/07/2018

    ! -- Calcul du temps de turnover en second
    TTsec = TT * 24 * 3600

    !-- Calcul de la fonction dd'hdrolyse:
 do k=1,nz_max
   do j=ny_min,ny_max
     do i=nx_min,nx_max
          if (iiX ==iiC) then 
             QC = var(iiC)%conc(i,j,k)/(var(iicell)%conc(i,j,k) + 1.d-80)
           ! Calcul de la fonction de quota
             fQ = (Qmax -QC)/(Qmax-Qmin)
          else
             XC = var(iiX)%conc(i,j,k)/(var(iiC)%conc(i,j,k) + 1.d-80)
           ! Quota function of the C:X ratio in bacteria 
             fQ = ((XC - XCmin)/(XCmax-XCmin))
           endif
           if (fQ >= 0.0e0) fQ = fQ**nn
           fQ = max(0.0e0,fQ)
           fQ = min(1.0e0,fQ)
!-- seuil min en bacteries                
           if (var(iiCell)%conc(i,j,k) > 1000 &
                .and. (var(iiMOP)%conc(i,j,k) > (var(iiMOP)%conc(i,j,k)/TTsec)*fQ *dt_bio)) then
              f_hydrolysis_XC(i,j,k) = (var(iiMOP)%conc(i,j,k)/TTsec)*fQ   
            else
              f_hydrolysis_XC(i,j,k) = 0.0e0
            endif
         enddo
    enddo
 enddo

end function f_hydrolysis_XC
!