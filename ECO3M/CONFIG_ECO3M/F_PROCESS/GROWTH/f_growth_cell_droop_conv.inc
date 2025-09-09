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
function f_growth_cell_droop_conv(ili,jcol,mumax,indcell,indX1,Qmin1,Qmax1,iconvC,iconv)
    !
    !---------------------------------------------------------------------------
    ! Fonction pour convertir la croissance des copepodes en carbone et autre unites
    !
    ! Melika Baklouti 10/04/2012
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_growth_cell_droop_conv(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: Qmin1,Qmax1,indX1,mumax,indcell
    Real(8),optional :: iconvC,iconv
    Integer :: ili,jcol

    !--variables locales:
    Integer :: iorg1,iorgCell
    Integer :: iindX1,iiCell,iCell_prey,iorg_prey
    Integer :: iiconvC,iiconv
    Integer :: i,j,k
    Real(8) :: Q1,fQ1
    Real(8) :: fQ
    Character(L_VAR)    :: chain


    !-- Indice de la concentration en cellules du predateur
    iiCell = nint(indcell)
    iorgCell =var(iiCell)%idorg

    if (present(iconvC) .and. present (iconv)) then
        iiconv = nint(iconv)
        iiconvC = nint(iconvC)
        iorg_prey =var(iiconv)%idorg
        chain='cell'
        iCell_prey = f_idorg2id(iorg_prey,chain)
    endif
    

    !-- Indice de la concentration du premiere element (obligatoire)
    iindX1 = nint(indX1)
    iorg1=var(iindX1)%idorg

    !-- test sur les organismes :
    if (iorgCell /= iorg1 ) then
        write(*,*) 'pb avec f_growth_cell_droop_conv'
        STOP
    endif

    !-- Calcul de la fonction de limitation :
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                fQ1=0.d0

                Q1 = var(iindX1)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calcul des fQ de croissance
                if (Q1 >= 0.d0)   fQ1 = 1.d0 - Qmin1/(Q1+1.d-80) 


                !MB 11/04/2011 Ajout de limites 
                fQ1 = min(fQ1,1.d0)

                fQ1 = max(0.d0,fQ1)


                !-- calcul du quota minimum:
                fQ =  fQ1
                !-- calcul du taux de croissance brute (cell.l-1.s-1)
                f_growth_cell_droop_conv(i,j,k) =  mumax*fQ*var(iiCell)%conc(i,j,k)

                !-- conversion en carbone :
                f_growth_cell_droop_conv(i,j,k) =  f_growth_cell_droop_conv(i,j,k)* 0.85* Qmin1 
                if (present(iconvC) .and. present (iconv)) then

                    !-- conversion supplementaire (si unites de iconv et iconvC differentes):
                    if (iiconvC /= iiconv) then
                        f_growth_cell_droop_conv(i,j,k) =  &
                            f_growth_cell_droop_conv(i,j,k) * VAR(iiconv)%conc(i,j,k) / (VAR(iiconvC)%conc(i,j,k) + 1.d-80)
                    endif

                endif

                !-- mise en place de seuils (on s'assure que la concentration est suffisante dans la maille pour un pas de temps)  :

                if (VAR(iiconv)%conc(i,j,k) <  f_growth_cell_droop_conv(i,j,k) * dt_bio ) then 
                        f_growth_cell_droop_conv(i,j,k) = 0.d0
                endif

            enddo
        enddo
    enddo
end function f_growth_cell_droop_conv
!--------------------------------------------------------------------
