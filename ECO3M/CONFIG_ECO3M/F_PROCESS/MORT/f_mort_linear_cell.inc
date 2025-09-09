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
function f_mort_linear_cell(ili,jcol,kmort)
    !      
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule la vitesse de mortalite du zoo, phyto,
    ! ou celle des bactéries par lyse virale selon une loi LINEAIRE
    !
    !
    ! Variables et parametres et leurs unites:
    !----------------------------------------- 
    !
    ! kmort : cste de vitesse de la mortalite (s^-1)
    !
    ! ATTENTION: c'est la variable d'indice ili qui disparait par mortalité
    !
    !
    ! Melika Baklouti 09/07/2007
    ! Melika Baklouti 29/05/2017
    !--------------------------------------------------------------------


    Implicit None

    !-- fonction:
    Real(8)             :: f_mort_linear_cell(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)             :: kmort !-- cste spécifique de mortalité, s^-1
    Integer             :: ili,jcol,i,j,k,iiCell,iorg
    Character(L_CHAIN)  :: chain
    
!--variables locales:
    chain = 'cell'
    iorg = var(ili)%idorg
    iiCell = f_idorg2id(iorg,chain)


    !-- Vitesse de mortalite naturelle (mol.X .m-3 . s-1):
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                if (var(iiCell)%conc(i,j,k) > var(iiCell)%conc(i,j,k) * kmort * 20 * 60) then
                    f_mort_linear_cell(i,j,k)= var(ili)%conc(i,j,k) * kmort
                else
                    f_mort_linear_cell(i,j,k) = 0.0
                endif
            enddo
        enddo
    enddo

    ! write(*,*) 'fin de la fonction f_mort_linear_cell'

End function f_mort_linear_cell
!--------------------------------------------------------------------
