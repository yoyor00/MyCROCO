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
function f_hydrol_CNP_2025_clean(ili,jcol,indMOP,TT,indBACC,QminC,QmaxC,QminN,QmaxN,QminP,QmaxP)
    !
    !---------------------------------------------------------------------------
    ! This function calcules the rate of particulate detrital organic matter (DET) by bacteria
    ! based on a fixed turnover time
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indMOP  : index of DET concentration
    ! indBACC: index of bacteria C biomass 
    ! TT = turnover time of DET in days (cette valeur peut devenir extremement elevee 
    ! pour la matiere refractaire)
    !
    !
    ! R. Mauriac  07/09/2009
    ! Last revision:  M. Baklouti 09/07/2025
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_hydrol_CNP_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: TT,indMOP,indBACC,QminC,QmaxC
    Real(8),optional :: QminN,QmaxN,QminP,QmaxP
    Integer :: ili,jcol

    !--local variables:
    integer :: i,j,k,iiMOP
    integer :: identite,iiC,iiCell,iiN,iiP
    real(8) :: nn,nnC,nnX
    real(8) :: TTsec,QC,QN,QP,fQC,fQN,fQP,fQ
    character(L_VAR) :: chain,chain2,chain3

    !-- Index of  concentration of detrital particulate organic matter
    iiMOP = nint(indMOP)
    !-- Index of  cabon concentration of bacteria
    iiC = nint(indBACC)

    ! Indices of bacterial concentrations and abundance
    chain = 'cell'
    identite = var(iiC)%idorg
    iiCell = f_idorg2id(identite,chain)
    chain = 'N'
    iiN = f_idorg2id(identite,chain)
    chain = 'P'
    iiP = f_idorg2id(identite,chain)


    nnC = 0.2
    nnX = 2

    ! -- Calcul du temps de turnover en second
    TTsec = TT * 24 * 3600

    !-- Calcul de la fonction de remineralisation:
 do k=1,nz_max
   do j=ny_min,ny_max
     do i=nx_min,nx_max
        QC = var(iiC)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
   ! Carbon quota function for bacteria
        fQC = ((QmaxC - QC)/(QmaxC-QminC))
   ! Nitrogen quota function for bacteria
        QN = var(iiN)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
        fQN = ((QmaxN - QN)/(QmaxN-QminN))
   ! Phosphorous quota function for bacteria
        QP = var(iiP)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
        fQP = ((QmaxP - QP)/(QmaxP-QminP))
                
        fQC = max(0.0e0,fQC)
        fQC = min(1.0e0,fQC)
        fQN = max(0.0e0,fQN)
        fQN = min(1.0e0,fQN)
        fQP = max(0.0e0,fQP)
        fQP = min(1.0e0,fQP)

        fQ = max(fQC**nnC,fQN**nnX,fQP**nnX)

        f_hydrol_CNP_2025_clean(i,j,k) = 0.0d0

!-- test on  the presence of bacteria and on DET availability                 
        if (var(iiCell)%conc(i,j,k) .gt. 100  .and. &
                    var(iiMOP)%conc(i,j,k) > 1.1 * var(iiMOP)%conc(i,j,k)/TTsec*fQ *dt_bio ) then

            f_hydrol_CNP_2025_clean(i,j,k) = var(iiMOP)%conc(i,j,k)/TTsec*fQ   

        endif
     enddo
   enddo
 enddo

end function f_hydrol_CNP_2025_clean
!
