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
function f_fQPPCX_rem_2026(ili,jcol,indX,QCmin,QCmax,QXmin,QXmax)
    !
    !--------------------------------------------------------------------
    ! This function calculates the quota function for remineralization as a funtion
    ! of the carbon cell quota 
    !
    ! This function is based on the formulation provided in Geider et al. (1968) 
    !
    ! Variables,  parameters and their units :
    !-----------------------------------------
    ! indC   : index of carbon concentration 
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Melika Baklouti 01/07/2011; Last version July 2025
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_fQPPCX_rem_2026(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indX,QCmin,QCmax,QXmin,QXmax
    Integer :: ili,jcol

    !-- local variables
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiX
    Real(8) :: QC,fQC,QX,fQX,hQC
    Character(L_VAR)::chain2
    Real:: nn

    iiX = nint(indX)
    identite = var(iiX)%idorg
    chain2 = 'cell'
    iiCell = f_idorg2id(identite,chain2)
    chain2 = 'C'
    iiC = f_idorg2id(identite,chain2)


    nn = 0.06


!-- Calcul of the C quota function and the part of mineralized flux  :

    do k = 1, nz_max
        do j = ny_min, ny_max
            do i = nx_min, nx_max

            !-- carbon cell quota :
                QC = var(iiC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                QX = var(iiX)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                fQC = (QCmax - QC) / (QCmax - QCmin)
                fQX = (QXmax - QX) / (QXmax - QXmin)
                fQC = min(1.d0,fQC)
                fQC = max(0.d0,fQC)
                fQX = min(1.d0,fQX)
                fQX = max(0.d0,fQX)

		hQC = 1 - fQC**nn

                f_fQPPCX_rem_2026(i,j,k) = min(hQC,fQX**4)
            enddo
        enddo
    enddo



end function f_fQPPCX_rem_2026
!--------------------------------------------------------------------
