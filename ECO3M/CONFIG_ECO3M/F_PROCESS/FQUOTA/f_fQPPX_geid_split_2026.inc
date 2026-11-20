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
function f_fQPPX_geid_split_2026(ili,jcol,indX,Qmin,Qmax,CXmin,CXmax,QminC,QmaxC)
    !
    !--------------------------------------------------------------------
    ! Cette fonction permet de contraidre les rapport C/N et C/P et N/P entre le valeurs 
    ! min_ratio et max_ratio
    !
    ! La fonction de quota est donnée par le modèle de Geider et al. (1968) 
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indX: indice de la conc. interne en nutriment X dont le rapport est a contraindre
    !
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Romain MAuriac 4/09/2007
    ! Melika Baklouti Avril 2017
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_fQPPX_geid_split_2026(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: CNmin,CNmax,indX,Qmin,Qmax
    Real(8) :: CXmin,CXmax
    Real(8), optional :: QminC,QmaxC
    Integer :: ili,jcol


 f_fQPPX_geid_split_2026 = 1.0d0 - f_hQPPX_geid_split_2026(ili,jcol,indX,Qmin,Qmax,CXmin,CXmax,QminC,QmaxC)
                

end function f_fQPPX_geid_split_2026
!--------------------------------------------------------------------
