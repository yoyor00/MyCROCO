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
function f_fQRESPC_geid_2020(ili,jcol,indC,Qmin,Qmax,CP_min,CP_max,CN_min,CN_max)
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
    ! CP_min : rapport C:P minimum
    ! CP_max : rapport C:P maximum
    ! CN_min : rapport C:N minimum
    ! CN_max : rappport C:N maximum
    !
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Melika Baklouti 01/07/2011
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_fQRESPC_geid_2020(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,Qmin,Qmax
    Real(8), optional :: CP_min,CP_max,CN_min,CN_max
    Integer :: ili,jcol

    !--variables locales:
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiN,iiP
    Real(8) :: QNC,hQNC,fQNC,QPC,hQPC,fQPC,QC,hQC,fQC,hQ,NC_min,NC_max,PC_min,PC_max,nn
    Character(L_VAR)::chain2,chain3,chain4


    !-- Calcul de la fonction de limitation :

    f_fQRESPC_geid_2020 = 1.d0 - f_fQPPC_geid_2020(ili,jcol,indC,Qmin,Qmax,CP_min,CP_max,CN_min,CN_max)

    !WRITE(*,*) 'fin FUNCTION f_fQRESPC_geid_2020 '
end function f_fQRESPC_geid_2020
!--------------------------------------------------------------------
