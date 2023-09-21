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
function f_fQPPC_geid_2020(ili,jcol,indC,Qmin,Qmax,CP_min,CP_max,CN_min,CN_max)
    !
    !--------------------------------------------------------------------
    ! This function calculates the quota function for C exudation as a funtion
    ! of the carbon cell quota and the intracellular C:N and C:P ratios
    !
    ! This function is based on the formulation provided in Geider et al. (1968) 
    !
    ! Variables,  parameters and their units :
    !-----------------------------------------
    ! indC   : index of carbon concentration 
    ! CP_min : C:P minimum ratio  
    ! CP_max : C:P maximum ratio
    ! CN_min : C:N minimum ratio
    ! CN_max : C:N maximum ratio
    !
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Melika Baklouti 01/07/2011; Last version May 2020
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_fQPPC_geid_2020(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,Qmin,Qmax
    Real(8), optional :: CP_min,CP_max,CN_min,CN_max
    Integer :: ili,jcol

    !-- local variables
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiN,iiP
    Real(8) :: QNC,hQNC,fQNC,QPC,hQPC,fQPC,QC,hQC,fQC,hQ,NC_min,NC_max,PC_min,PC_max,nn
    Character(L_VAR)::chain2,chain3,chain4
    Real:: nn_PHYL, nn_PHYS,nn_hetero

    iiC = nint(indC)
    identite = var(iiC)%idorg
    chain2 = 'cell'
    chain3 ='N'
    chain4 ='P'
    iiCell = f_idorg2id(identite,chain2)


    nn_PHYL = 0.3
    nn_PHYS = 0.3
    nn_hetero = 0.06  ! heterotrophs and autotrophs have a different strategy for C

    if (var(iiC)%scomp == "phyl" .or. var(iiC)%scomp == "PHYL" .or.  var(iiC)%scomp == "Phyl") then
        nn = nn_PHYL
    elseif (var(iiC)%scomp == "phys" .or. var(iiC)%scomp == "PHYS" .or.  var(iiC)%scomp == "Phys") then
        nn = nn_PHYS
    else
        nn = nn_hetero
    endif

    !-- initialization:
    hQNC = 0.d0
    hQPC = 0.d0

    !-- N:C min and N:C max ratios:
    if (present (CN_min)) then
        iiN  = f_idorg2id(identite,chain3)
        !-- conversion de C/X en X/C
        NC_min = 1.d0/CN_max
        NC_max = 1.d0/CN_min
    endif

    !-- P:C min and P:C max ratios:
    if (present (CP_min)) then
        iiP  = f_idorg2id(identite,chain4)
        !-- conversion de C/X en X/C
        PC_min = 1.d0/CP_max
        PC_max = 1.d0/CP_min
    endif

!-- Calcul of the quota function for C exudation  :

if (present (CN_min) .and. present(CP_min)) then
    do k = 1, nz_max
        do j = ny_min, ny_max
            do i = nx_min, nx_max

            !-- carbon cell quota :
                QC = var(iiC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                fQC = (Qmax - QC) / (Qmax - Qmin)
                fQC = min(1.d0,fQC)
                fQC = max(0.d0,fQC)
                fQC = fQC**nn

            !-- exudation function (0 <= hQC <= 1)
                hQC = 1.d0 - fQC

            !-- test on C cell quota :

                if (QC <= Qmin) then
                    f_fQPPC_geid_2020(i,j,k) = 0.0d0
                elseif (QC >= Qmax) then
                    f_fQPPC_geid_2020(i,j,k) = hQC
                else !-- C cell quota is in the range [Qmin;Qmax]

                    QNC = var(iiN)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+ 1.d-80) ! quota N:C

                  !-- Quota function for net C uptake based on N:C ratio
                    fQNC = (QNC - NC_min) / (NC_max - NC_min)

                    fQNC = min(1.d0,fQNC)
                    fQNC = max(0.d0,fQNC)
                    fQNC = fQNC**nn
                    !-- Quota function for exudation
                    hQNC = 1.d0 - fQNC

                    QPC = var(iiP)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-80)   ! quota P:C

                  !-- Quota function of net uptake of C based on P:C ratio
                    fQPC = (QPC - PC_min) / (PC_max - PC_min)

                    fQPC = min(1.d0,fQPC)
                    fQPC = max(0.d0,fQPC)
                    fQPC = fQPC**nn

                  !-- Quota function for exudation
                    hQPC = 1.d0 - fQPC
                    hQ = max(hQC,hQNC,hQPC)
                    f_fQPPC_geid_2020(i,j,k) = hQ
                endif
            enddo
        enddo
    enddo

else !-- case where the model is only expressed in carbon :
    do k = 1, nz_max
        do j = ny_min, ny_max
            do i = nx_min, nx_max

            !-- carbon cell quota :
                QC = var(iiC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                fQC = (Qmax - QC) / (Qmax - Qmin)
                fQC = min(1.d0,fQC)
                fQC = max(0.d0,fQC)
                fQC = fQC**nn

            !-- exudation function (0 <= hQC <= 1)
                hQC = 1.d0 - fQC
                f_fQPPC_geid_2020(i,j,k) = hQC

            enddo
        enddo
    enddo
endif

end function f_fQPPC_geid_2020
!--------------------------------------------------------------------
