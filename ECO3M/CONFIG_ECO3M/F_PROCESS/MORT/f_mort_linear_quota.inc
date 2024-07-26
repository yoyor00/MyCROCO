!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
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
function f_mort_linear_quota(ili,jcol,kmort,QminC,QmaxC)
    !      
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule la vitesse de mortalite du zoo, phyto,
    ! ou celle des bact�ries par lyse virale selon une loi LINEAIRE
    ! tant que le quota cellulaire en carbone est superieur a QminC
    ! En de�a, la constante de vitesse augment de maniere quadratique
    !
    ! Variables et parametres et leurs unites:
    !----------------------------------------- 
    !
    ! kmort : cste de vitesse de la mortalite (s^-1)
    !
    ! ATTENTION: c'est la variable d'indice ili qui disparait par mortalit�
    !
    !
    ! Melika Baklouti 29/05/2017
    !--------------------------------------------------------------------


    Implicit None

    !-- fonction:
    Real(8)             :: f_mort_linear_quota(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)             :: kmort !-- cste sp�cifique de mortalit�, s^-1
    Real(8)             :: QminC,QmaxC !-- quotas cellulaires en carbone, mol C m^-3
    Integer             :: ili,jcol

    !--variables locales:
    Character(L_CHAIN)  :: chain
    Integer:: i,j,k,iiCell,iorg,iiC    
    Real(8):: kmort_quota
    Real(8):: QCrel(nx_min:nx_max,ny_min:ny_max,nz_max) 
 
    iorg = var(ili)%idorg
    chain = 'cell'
    iiCell = f_idorg2id(iorg,chain)
    chain = 'C'
    iiC = f_idorg2id(iorg,chain)


    !-- Vitesse de mortalite naturelle (mol.X .m-3 . s-1):
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                 QCrel(i,j,k) = (var(iiC)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+ 1.d-80) - QminC)/(QmaxC-QminC)
                 if (QCrel(i,j,k) >  0) then
                    kmort_quota = kmort
                 else
                    kmort_quota = kmort + 9 * kmort * (2*QCrel(i,j,k))**2 ! vitesse de mortalite x 10 quand QCrel = - 0.5
                 endif

                if (var(iiCell)%conc(i,j,k) > var(iiCell)%conc(i,j,k) * kmort_quota * 20 * 60) then
                    f_mort_linear_quota(i,j,k)= var(ili)%conc(i,j,k) * kmort_quota
                else
                    f_mort_linear_quota(i,j,k) = 0.0d0
                endif
            enddo
        enddo
    enddo

    ! write(*,*) 'fin de la fonction f_mort_linear_quota'

End function f_mort_linear_quota
!--------------------------------------------------------------------