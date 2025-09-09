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
function f_fQChl_2020(ili,jcol,indC,Qmin,Qmax,N_min_ratio,N_max_ratio)
    !
    !--------------------------------------------------------------------
    ! Cette fonction calcule le facteur d'exudation de carbone en fonction
    ! des contraintes sur les rapports C:N et C:P
    !
    ! La fonction de quota est donnée par le modèle de Geider et al. (1968) 
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indC        : indice de la concentration en carbone
    ! N_min_ratio : rapport C:N minimum
    ! N_max_ratio : rappport C:N maximum
    !
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Melika Baklouti 01/07/2011
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_fQChl_2020(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,Qmin,Qmax,N_min_ratio,N_max_ratio
    Integer :: ili,jcol

    !--variables locales:
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiN,iiP
    Real(8) :: Q1,fQ1,Q3,fQ3,NC_min,NC_max
    Real    :: nn,nn_PHYL,nn_PHYS
    Character(L_VAR)::chain2,chain3,chain4

    iiC = nint(indC)
    identite = var(iiC)%idorg
    chain2 = 'cell'
    chain3 ='N'
    iiCell = f_idorg2id(identite,chain2)


    nn_PHYL = 4!6!0.06
    nn_PHYS = 0.06
    if  (var(ili)%scomp == "phyl" .or. var(ili)%scomp == "PHYL" .or.  var(ili)%scomp == "Phyl") then
       nn = nn_PHYL
    elseif (var(ili)%scomp == "phys" .or. var(ili)%scomp == "PHYS" .or.  var(ili)%scomp == "Phys") then
       nn = nn_PHYS
    endif
    !-- initialisation:

    !-- Calcul de la fonction de limitation :

        iiN  = f_idorg2id(identite,chain3)
        !-- conversion de C/X en X/C
        NC_min = 1.d0/N_max_ratio
        NC_max = 1.d0/N_min_ratio

    !-- quota cellulaire en C :
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                !-- test sur le contenu cellulaire en C:

                            Q1 = var(iiN)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+ 1.d-80) ! quota N:C
                            ! Calcul des fQ d'uptake de C
                                fQ1 = ( Q1 - NC_min) / (NC_max - NC_min)
                                fQ1 = max(0.d0,fQ1)
                                fQ1 = min(1.d0,fQ1)
                                f_fQChl_2020(i,j,k) = fQ1**nn
            enddo
        enddo
    enddo

    !WRITE(*,*) 'fin FUNCTION f_fQChl_2020 ',ili,jcol
end function f_fQChl_2020
!--------------------------------------------------------------------
