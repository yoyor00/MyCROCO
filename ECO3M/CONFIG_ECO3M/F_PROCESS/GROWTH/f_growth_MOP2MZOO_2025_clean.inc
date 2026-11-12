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
function f_growth_MOP2MZOO_2025_clean(ili,jcol,indcell,indC,QC_MIN,QC_MAX,indDETC)
    !
    !---------------------------------------------------------------------------
    ! Fonction pour convertir une croissance specifique (en s-1) en une croissance en biomasse (mol X m-3 s-1) 
    !
    ! Melika Baklouti 16/10/2024
    !-------------------------------------------------------------------------
    USE mod_eco3m_user
    Implicit None

    !-- fonction:
    Real(8):: f_growth_MOP2MZOO_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,indcell,indDETC,QC_MIN,QC_MAX
    Integer :: ili,jcol

    !--variables locales:
    Real(8) :: QC,fQ
    Integer :: iorgX,iorgCell
    Integer :: iindC,iiCell,iindDETC
    Integer :: i,j,k,l,ll,nn
    Character(L_VAR)    :: chain


    nn = 0.6

    !-- Indice de la concentration en individus 
    iiCell = nint(indcell)
    iorgCell =var(iiCell)%idorg


    !-- Indice de la concentration de l'élément C
    iindC = nint(indC)
    iorgX=var(iindC)%idorg

    !-- Indice de la concentration en DET
    iindDETC = nint(indDETC)

    !-- Vérification  :
    if (iorgCell /= iorgX ) then
        write(*,*) 'pb avec f_growth_MOP2MZOO_2025_clean'
        STOP
    endif

    !-- recherche de l'indice de mu_growth qui correspond à iiCell
    do l = 1,nscp_cell
        if (mu_growth(l)%idorg == iiCell) then
           ll = l
        endif
    enddo
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                QC = var(iindC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                fQ = (QC_MAX - QC)/(QC_MAX - QC_MIN)
                fQ = min(1.d0,fQ)
                fQ = max(0.d0,fQ)

            !-- consumption of ten percent of the structure weight per day : 
                f_growth_MOP2MZOO_2025_clean(i,j,k) =  0.8*mu_growth(ll)%val(i,j,k)* var(iindC)%conc(i,j,k) * fQ**nn
                !f_growth_MOP2MZOO_2025_clean(i,j,k) =  mu_growth(ll)%val(i,j,k)* var(iindC)%conc(i,j,k) * fQ**nn
!                f_growth_MOP2MZOO_2025_clean(i,j,k) =   var(iindC)%conc(i,j,k) *0.001  /3600.0/24.0    * fQ**nn
                !-- mise en place de seuils (on s'assure que la concentration est suffisante dans la maille pour un pas de temps)  :
                   if (VAR(iindDETC)%conc(i,j,k) <  1.1 * f_growth_MOP2MZOO_2025_clean(i,j,k) * dt_bio ) then 
                        f_growth_MOP2MZOO_2025_clean(i,j,k) = 0.d0
                   endif
            enddo
        enddo
    enddo
end function f_growth_MOP2MZOO_2025_clean
!--------------------------------------------------------------------
