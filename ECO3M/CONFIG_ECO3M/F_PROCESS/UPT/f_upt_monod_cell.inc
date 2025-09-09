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
function f_upt_monod_cell(ili,jcol,ibio,inut,k_demi,Vmax)
    !
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule l''absorption/assimilation (uptake) d''un nutriment 
    ! d''indice inut par le modele de Monod-Michaelis-Menten 
    !
    ! Variables et parametres et leurs unites:
    !----------------------------------------- 
    ! ibio  : indice d'une des concentrations (C ou N ou P...)
    !         de l'organisme consommant du nutriment
    ! inut  : parametre donnant l'indice de la concentration
    !         en nutriment limitant (specifie car il peut différer de jcol) 
    ! k_demi: cste de demi saturation, (mol X)/l
    ! Vmax  : vitesse SPECIFIQUE d'absorption max., en  molX.cell-1.s^{-1}
    !
    ! f_upt_monod : vitesse d''absorption, en (mol X).l-1 s^{-1}
    !
    ! ATTENTION: 
    !    * Vmax doit etre imperativement en mol X.cell-1.s-1
    !
    ! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
    !      Michaelis, L. & Menten, M. L. (1913), 49:333-369
    !
    ! Melika Baklouti 01/07/2011
    !--------------------------------------------------------------------
    !
    !-- Variables globales
    ! Use VAR_GLOBAL

    Implicit None

    Real(8) :: f_upt_monod_cell(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: ibio,inut,k_demi,Vmax
    Integer :: ili,jcol

    !-- variables locales
    Integer :: iibio,iinut,iorg_plancton,iX,iCell
    Character(L_VAR)  :: cell
    Integer::i,j,k
    !logical :: isnan

    !-- indice de la conc. en nutriments:
    iinut=nint(inut)

    !-- Organisme de la biomasse planctonique:
    iibio = nint(ibio)
    iorg_plancton = var(iibio)%idorg 


    !--concentration en cell du plancton:
    cell = 'cell'
    iCell= f_idorg2id(iorg_plancton,cell)

    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                if (var(iinut)%conc(i,j,k) > Vmax*var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + k_demi) &
                                                              * var(iCell)%conc(i,j,k) * dt_bio) then
                  f_upt_monod_cell(i,j,k) = &
                     (Vmax * var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + k_demi)) * var(iCell)%conc(i,j,k)
                else 
                  f_upt_monod_cell(i,j,k) = 0.d0
                endif
           enddo
      enddo
  enddo
      
end function f_upt_monod_cell
!--------------------------------------------------------------------
