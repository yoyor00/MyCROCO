!***************************************************************************
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
!--------------------------------------------------------------------
!
function f_ppb_han(ili,jcol,a_chl,phi_max_C,tau,&
        sig_ps2,k_d,k_r,iphyC)
    !
    ! Cette fonction calcule la vitesse de production primaire BRUTE de la variable
    ! d'etat dont la concentration est donnee par la variable d'indice ili
    !
    ! variables parametres et leurs unites: 
    !--------------------------------------
    ! iphyC     : concentration de phyto en carbone (mol C/l). Si absente, ili utilise par defaut
    ! a_chl     : coeff. d''absorption spec. de la Chl, m^2 / (g Chl)
    ! phi_max_C : rendement quantique max., (mol C)/ J
    ! tau       : temps de turnover des electrons, s
    ! sig_ps2   : section efficace des PSII, m^2 / Joules
    ! kd_Han    : constante d'endommagement des PSII, adim ds modele de Han
    ! kr        : constante de r�paration  des PSII, s^-1 
    ! f_ppb_han : en molC/l/s
    !
    ! cette fonction permet aussi de remplir le tableau muphy_NR
    ! muphy_NR : vitesse specifique de production primaire en absence de lim. 
    !            par les nutriments, en s^-1
    !
    ! Ref:  Han, B. P. (2002) J. Theor. Biol. 214:519-527
    !
    ! Melika Baklouti 6/07/2007
    !-------------------------------------------------------------------------
    !-- Variables globales

    Implicit None

    !-- arguments
    Real(8) :: a_chl,phi_max_C,tau,sig_ps2,k_d,k_r
    Real(8),optional:: iphyC
    Real(8) :: f_ppb_han(nx_min:nx_max,ny_min:ny_max,nz_max)
    Integer :: ili,jcol, ii

    !-- variables locales:
    Integer:: iorg,iichl,iiC,iiphyC
    Real(8)  :: po(nx_min:nx_max,ny_min:ny_max,nz_max),Chl(nx_min:nx_max,ny_min:ny_max,nz_max)
    Character(L_VAR)::elmt

    !-- pourcentage de PSII ouverts:
    po = 1.e0/(1.d0 + sig_ps2 &
        * E_PARZ * tau + k_d/k_r * (sig_ps2 * E_PARZ)**2 *tau)

    !-- organisme:
    if (present(iphyC)) then
      iiphyC = nint(iphyC)
      iorg = var(iiphyC)%idorg
    else
      iorg = var(ili)%idorg
      iiphyC = ili
    endif

    !-- cas d'un rapport CHL:C constant:
    if (.not. CHL_C_BOOL) then
        Chl = CHL_C0 * var(iiphyC)%conc
    elseif ((CHL_C_BOOL) .and. (.not. CHL_C_FUNC)) then

        !-- Indice correspondant a la conc. chlorophylienne du phyto:
        elmt ='Chl'
        iichl = f_idorg2id (iorg,elmt)
        Chl = var(iichl)%conc
    else !-- cas ou CHL:C est calcule par une fction analytique:
        elmt ='C'
        iiC = f_idorg2id (iorg,elmt)
        Chl = CHL_C(1)%val * var(iiC)%conc
    endif

    !-- Vitesse de production primaire (en molC/m3/s):
    f_ppb_han =  a_chl * phi_max_C * po * E_PARZ * Chl


    !-- Sauvegarde de la vitesse specifique de croissance en absence
    !-- de limitation par les nutriments (Nutrient Replete), en s^-1, il faut:
    ! Important: verifier que l'allocation dyn. de mu_PPB_NR
    ! dans la subroutine ALLOC_VAR_Eco3M du fichier coupleur.F90
    ! est operationnelle
    !
    do ii = 1, nscp_phy
        if (mu_PPB_NR(ii)%idorg == iorg) then
            mu_PPB_NR(ii)%val = &
                a_chl*phi_max_C*po*E_PARZ*Chl/(var(iiphyC)%conc+1.0d-80)
            exit
        endif
    enddo

end function f_ppb_han