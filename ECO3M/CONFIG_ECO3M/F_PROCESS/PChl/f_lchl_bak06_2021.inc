!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr
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
!encouraged to load and testhe software''s suitability as regards their
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
!
function f_lchl_bak06_2021(ili,jcol,thetaN_max,alpha,phimax,Q_MIN_PHY_C,Q_MAX_PHY_C)
    !
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule la variation de chlorophylle en fonction
    ! de la production primaire par un modele presente dans
    ! (Baklouti et al, 2006) et utilisant le coefficient
    ! rho_Chl tel que defini dans (Geider, 1998).
    !
    ! /!\ modifications par rapport a la fonction f_lchl_bak06 /!\
    ! cette fonction calcule uniquement la limitation de la production de chl
    ! par le rapport Chl/N et l'irradiance
    ! elle doit etre multipliee par les fonctions d'uptake des differentes sources
    ! d'azote de l'organisme (uptake et diazotrophie eventuelle) dans le fichier
    ! config.ini
    !
    !
    ! variables et parametres et leurs unites:
    !----------------------------------------
    ! thetaN_max  = Chl:N maximum, gChl/molN
    ! alpha_star = pente a l''origine de PP^* vs E, molC.m^2/(gChl.J)
    ! rho_chl    = efficacite de photosynthese, gChl/molN
    ! Vm         = vitesse max d'absorption du nutriment azote
    !               son unite est d�termin�e par l'indice iXorC
    !               si iXorC correspond � la conc. de la biomasse azot�e <--> Vm est en s^{-1}
    !               si iXorC correspond � la conc. de la biomasse cellulaire <--> Vm est en molN. (cell)^{-1} s^{-1}
    !               si iXorC  vaut 0 Vm doit etre calcule selon la formule PBmax.Qmax et la valeur
    !               de Qmax est donn�e dans le param�tre Vm
    ! inut       = indice de la concentration du nutriment azot� absorb
    ! K          = constante de demi saturation de l'absorption de nutriment, molX m^{-3}
    ! ili        = contient l'indice de la concentration en chlorophylle de l'autotrophe
    !
    ! f_lchl_bak06_2021 = Vitesse de production de chlorophylle par mol de N absorbee, gChl/molN
    !
    !-- ordres de grandeur: thetaN_max= 0.3 * 14 gChl/molN (Geider, 1998)
    !
    ! Refs:   Baklouti et al. (2006) Progr. in Oceanogr., 71:1-33
    !         Geider et al. (1998) Limnol. Oceanog., 43(4):679-694
    !
    ! Revisions:
    ! Melika Baklouti 27/06/2008
    ! Clement ALDEBERT 09/07/2012
    ! Remi Pages, Melika Baklouti 20/04/2018
    !------------------------------------------------------------------
    !-- Variables globales

    Implicit None

    !-- arguments
    Real(8) :: thetaN_max,alpha_star,alpha,phimax,Q_MIN_PHY_C,Q_MAX_PHY_C
    Real(8) :: f_lchl_bak06_2021(nx_min:nx_max,ny_min:ny_max,nz_max)
    Integer :: ili,jcol

    !-- variables locales:
    Character(L_CHAIN) :: chaine
    integer :: iorg
    integer :: i,j,k, ii,iiphy
    real    :: nn
    integer :: iiC,iiChl,iiN, iicell
    real(8) :: QC,fratio, CHL_C_MAX, thetaN_max_local,fCHL_C
    real(8) :: rhochl(nx_min:nx_max,ny_min:ny_max,nz_max), theta(nx_min:nx_max,ny_min:ny_max,nz_max)
    real(8) :: thetaN(nx_min:nx_max,ny_min:ny_max,nz_max)
    real :: nn_PHYS,nn_PHYL

    !-- Initialisations:
    alpha_star = alpha * phimax
    nn_PHYS = 0.06! 0.6!0.06
    nn_PHYL = 0.06

    !-- Definition d'un rapport Chl:C max different pour PHYS et PHYL
    
    if  (var(ili)%scomp == "phyl" .or. var(ili)%scomp == "PHYL" .or.  var(ili)%scomp == "Phyl") then
!       CHL_C_MAX = 1.2      ! (C:Chl) must be > 10  gC/g Chl to synthetise new Chl
        nn=nn_PHYL
    elseif (var(ili)%scomp == "phys" .or. var(ili)%scomp == "PHYS" .or.  var(ili)%scomp == "Phys") then
!       CHL_C_MAX = 0.6      ! (C:Chl) must be > 20  gC/g Chl to synthetise new Chl
        nn=nn_PHYS
    else
        STOP "problem in function f_lchl_bak06_2021"
    endif

    !-- Recherche du numero d'organisme de l'autotrophe:
    iorg = var(ili)%idorg

    !-- de l'indice de sa conc en C :
    chaine = 'C'
    iiC = f_idorg2id (iorg,chaine)

    !-- de l'indice de sa conc en N :
    chaine = 'N'
    iiN = f_idorg2id (iorg,chaine)

    !-- et de celle de sa conc en Chl :
    chaine = 'Chl'
    iiChl = f_idorg2id (iorg,chaine)
    !-- et de celle de sa conc en Chl :
    chaine = 'cell'
    iicell = f_idorg2id (iorg,chaine)

   
    !--rapport CHL:C (g Chl/mol C):
    if ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .FALSE.) ) then
        theta = VAR(iiChl)%conc / (VAR(iiC)%conc + 1d-80)
        thetaN =  VAR(iiChl)%conc / (VAR(iiN)%conc + 1d-80)
    elseif ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .TRUE.) ) then
        theta = CHL_C(1)%val
        thetaN =  VAR(iiChl)%conc / (VAR(iiN)%conc+1d-80)
    elseif ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .TRUE.)) then
        theta = CHL_C(1)%val
        thetaN = CHL_C(1)%val * VAR(iiC)%conc / (VAR(iiN)%conc+1d-80)
    else
        theta = CHL_C0
        thetaN = CHL_C0 * VAR(iiC)%conc / (VAR(iiN)%conc+1d-80)
    endif

!--  Find the index of phytoplankton among all phytoplancton PFTs
    do ii=1,nscp_phy
        if (mu_PPB_NR(ii)%idorg == iorg) then
            iiphy = ii
            exit
        endif
    enddo


!-- Chlorophyll production rate (in gChl/molN):

 if ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .FALSE.)) then
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                     E_PARZ(i,j,k) = max(E_PARZ(i,j,k),0.0d0)

!                     CHL_C_MAX =  0.6 * exp(-0.5 * E_PARZ(i,j,k)) + 0.6 ! SIM 492 
!                     CHL_C_MAX =  0.55 * exp(-0.25 * E_PARZ(i,j,k)) + 0.6 ! SIM 465
!                     CHL_C_MAX =  0.6 * exp(-0.2 * E_PARZ(i,j,k)) + 0.6 ! SIM 417
!                     CHL_C_MAX =  0.7 * exp(-0.5 * E_PARZ(i,j,k)) + 0.6
                     CHL_C_MAX =  0.6 * exp(-0.3 * E_PARZ(i,j,k)) + 0.6 ! SIM 381

!                     thetaN_max_local = CHL_C_MAX * 3.53 ! 3.53 = (C:N)min (celui des proteines)             
!                     thetaN_max_local = min(thetaN_max_local,thetaN_max)
                     thetaN_max_local = thetaN_max
                      
    !-- Calcul du coefficient ro_chl de l'expression de Geider (1998)
                     rhochl(i,j,k) = mu_PPB_NR(iiphy)%val(i,j,k) / (alpha_star * E_PARZ(i,j,k) * theta(i,j,k) +1d-80) !-- egal a po finalement
                     rhochl(i,j,k)  = min(1.0d0, rhochl(i,j,k) )
                     rhochl(i,j,k)  = max(0.0d0, rhochl(i,j,k) )

                     if (thetaN(i,j,k) <= thetaN_max_local .and. theta(i,j,k) < CHL_C_MAX .and. E_PARZ(i,j,k) >= 0.5) then 

                          !-- Regulation by C quota in phytoplankton
                          QC = VAR(iiC)%conc(i,j,k) / (VAR(iicell)%conc(i,j,k) + 1.d-80)
                          fratio= (Q_MAX_PHY_C-QC)/(Q_MAX_PHY_C-Q_MIN_PHY_C)
  			  fratio = min(1.0d0,fratio)
                          fratio = max(0.0d0,fratio)
      			  
                          !-- Regulation by Chl:C ratio in phytoplankton
!                          fCHL_C = (CHL_C_MAX - theta(i,j,k))/(CHL_C_MAX - 0.12)
!                          fCHL_C = (1.d0 - theta(i,j,k)/CHL_C_MAX)/(1.d0-theta(i,j,k)/CHL_C_MAX + 0.05)
!                          fCHL_C = (1.d0 - thetaN(i,j,k)/thetaN_max_local)/(1.d0-thetaN(i,j,k)/thetaN_max_local + 0.05)
                         fCHL_C = ((thetaN_max_local - thetaN(i,j,k) )/(thetaN_max_local - 0.42))
                          fCHL_C = min(1.0d0,fCHL_C)
                          fCHL_C = max(0.0d0,fCHL_C)
                          fCHL_C = fCHL_C**0.9!0.06

                          f_lchl_bak06_2021(i,j,k) =  thetaN_max_local * rhochl(i,j,k) * fCHL_C * fratio**nn
                    else
                          f_lchl_bak06_2021(i,j,k) =  0.d0
                    endif
!                    IF (ISNAN(f_lchl_bak06_2021(i,j,k))) THEN
!                       WRITE(*,*) 'FUNCTION f_lchl_bak06_2021 contains NaN',i,j,k
!                       STOP
!                   ENDIF
            enddo
        enddo
    enddo
 endif
 if ((CHL_C_BOOL .EQV. .FALSE.) .AND. (CHL_C_FUNC .EQV. .FALSE.)) then
      VAR(iiChl)%conc = CHL_C0 * VAR(iiC)%conc
  endif


end function f_lchl_bak06_2021
