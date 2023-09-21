!***************************************************************************
!***************************************************************************
!Copyright or ï¿œ or Copr. CNRS/IRD/Universite de la Mediterranee
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
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
function f_graz_KOOIJ_KTW_cell(ili,jcol,ipred,iprey,iprey_X,nbprey,Beta,&
        iCell1,iC1,Imax1,F1,pref1,iCell2,iC2,Imax2,F2,pref2,icell3,iC3,Imax3,F3,pref3,iCell4,iC4,Imax4,F4,pref4)
    !
    !--------------------------------------------------------------------
    ! Cette fonction calcule la vitesse de predation d'un predateur
    ! en presence de PLUSIEURS PROIES potentielles
    ! Reponse fonctionnelle: 
    !     * de type classe 3 multiple functional responses (cf. Gentleman et al. 2003 DSR 50: 2847-2875)
    !     * tient compte des interferences entre proies
    !     * preferences ponderees,i.e. variables et dependantes de la concentration en proies (permet la coexistences des proies)
    !     * maximal feeding 
    !     * active prey switching 
    !
    ! Attention:
    ! ----------
    ! 
    ! * cette fonction s''applique pour un FLUX  impliquant
    !   un predateur dont la concentration est donnee par VAR(ipred)%conc et
    !   une proie dont la concentration est donnee par VAR(iprey)%conc
    !
    ! * Le calcul se base sur les concentrations en C du predateur et des proies,
    !  et on obtient une vitesse BRUTE de grazing Graz en [mol C] m-3 s^{-1}
    !  qui ne tient pas compte des pertes dues a  l'egestion, et au sloppy feeding
    !
    ! Variables et parametres et leurs unites:
    !----------------------------------------- 
    ! pref_i             : preference fixe du predateur pour la proie i, (sans unite [0,1])
    ! nbprey             : nb de proies susceptibles d''etre broutees par le predateur (entre 1 et 4)
    ! cminpre            : concentration en dessous de laquelle il n'y a pas predation (en abondances)
    ! Im_i               : vitesse d'ingestion maximum  (en s-1)
    ! iCell_i            : indice des proies potentielles exprimees en abondances 
    ! iC_i               : indice des proies potentielles exprimees en C 
    ! Beta		     : constante qui determine la forme de la "feeding probability Q" =2 dans le papier de référence
    !                     beta = 1: reponse de type II; beta = 2: reponse de type III
    ! f_graz_KOOIJ_KTW_cell    : vitesse de consommation de la proie iprey (en ind m-3)
    !                          par le predateur ipred (en ind m-3)
    !
    !
    ! Ref: Vallina et al. Progress in Oceanography 120 (2014),93-109
    !
    ! Frederic Diaz, Melika Baklouti  02/05/2014
    !--------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8)             :: f_graz_KOOIJ_KTW_cell(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)             :: ipred,iprey,iprey_X,nbprey,iCell1,iC1,pref1,Beta,F1,Imax1
    Integer             :: ili,jcol
    Real(8), optional   :: iCell2,iCell3,iCell4, iC2,iC3,iC4
    Real(8), optional   :: pref2,pref3,pref4 
    Real(8), optional   :: F2,F3,F4,Imax2,Imax3,Imax4 
    !--variables locales:
    Integer              :: int_nbprey,iiprey,iXprey,iipred,pos_iprey,iscp_pred,iscp_prey
    Integer              :: iCell_pred,iCell_prey,iorg_prey,iorg_pred
    Integer              :: i,j,k,ii,mm
    Integer              :: indCell(nint(nbprey)) ! indice de la conc en abondance des proies
    Integer              :: indC(nint(nbprey)) ! indice de la conc en carbone des proies
    Real(8)              :: pref(nint(nbprey)),vfilt(nint(nbprey)),Imax(nint(nbprey))
    Real(8)              :: ress(nx_min:nx_max,ny_min:ny_max,nz_max),ress_v(nx_min:nx_max,ny_min:ny_max,nz_max)
    Real(8)              :: pref_v(nint(nbprey),nx_min:nx_max,ny_min:ny_max,nz_max)
    Real(8)                 :: sum_selec
    ! Character(L_VAR)     :: elmt
    Character(L_CHAIN)   :: chaine

    !-- conversion en entier:
    int_nbprey = nint(nbprey)
    iipred = nint(ipred)
    iiprey = nint(iprey)
    iXprey = nint(iprey_X)


    !-- initialisations:
    indCell = 0
    indC = 0
    pref = 0.d0
    f_graz_KOOIJ_KTW_cell=0.d0

    !
    !-- premiere proie (obligatoire)
    indCell(1) = nint(iCell1)
    indC(1) = nint(iC1)
    pref(1)= pref1
    vfilt(1) = F1
    Imax(1) = Imax1

    Select case (int_nbprey)
    case (2)
        indCell(2) = nint(iCell2)
        indC(2) = nint(iC2)
        pref(2)= pref2
        vfilt(2) = F2
        Imax(2) = Imax2

    case (3)
        indCell(2) = nint(iCell2)
        indC(2) = nint(iC2)
        pref(2)= pref2 
        vfilt(2) = F2
        Imax(2) = Imax2

        indCell(3) = nint(iCell3)
        indC(3) = nint(iC3)
        pref(3)= pref3
        vfilt(3) = F3
        Imax(3) = Imax3

    case (4)
        indCell(2) = nint(iCell2)
        indC(2) = nint(iC2)
        pref(2)= pref2
        vfilt(2) = F2
        Imax(2) = Imax2

        indCell(3) = nint(iCell3)
        indC(3) = nint(iC3)
        pref(3)= pref3 
        vfilt(3) = F3
        Imax(3) = Imax3

        indCell(4) = nint(iCell4)
        indC(4) = nint(iC4)
        pref(4)= pref4
        vfilt(4) = F4
        Imax(4) = Imax4

    end select

    !-- Initialisation du vecteur des preferences variables (selectivites):
    pref_v=0.d0
    !-- Initialisation de la ressource totale (calculee avec pref fixes et variables): 
    ress = 0.d0
    ress_v = 0.d0
    
    !-- Calcul en chaque point de l'espace, des selectivites modulant les preferences fixes par la conc de proies (KTW):
    do k = 1, nz_max 
        do j = ny_min, ny_max
            do i = nx_min, nx_max
                sum_selec = 0.d0
                do ii = 1, int_nbprey
                    !-- Indice correspondant a la conc. en C de la proie (la modulation par la conc de proies se fait par les biomasses en C car il y 
                    !   a trop d'ecarts entre les abondances l'effet de la modulation serait quasi nul):
                    sum_selec = sum_selec + pref(ii) * var(indC(ii))%conc(i,j,k)
                enddo

                do ii=1,int_nbprey
                    pref_v(ii,i,j,k) = pref(ii) * var(indC(ii))%conc(i,j,k) / (sum_selec+1.d-80)
                    ! Ressource totale :
                    ress(i,j,k)    = ress(i,j,k) + vfilt(ii)/Imax(ii) * pref(ii) * var(indCell(ii))%conc(i,j,k)
                    ress_v(i,j,k)  = ress_v(i,j,k) + vfilt(ii)/Imax(ii) * pref_v(ii,i,j,k) * var(indCell(ii))%conc(i,j,k)
                enddo
            enddo
        enddo
    enddo


    !-- numeros organisme de la proie et du predateur:
    iorg_prey = VAR(iiprey)%idorg
    iorg_pred = VAR(iipred)%idorg

    !-- Indice correspondant a la conc. en cellule de la proie:
    chaine = 'cell'
    iCell_prey = f_idorg2id(iorg_prey,chaine)
    iCell_pred = f_idorg2id(iorg_pred,chaine)


    !-- Creation et sauvegarde du tableau des vitesses specifiques de grazing du predateur i
    !   sur la proie j (Gamma_j x G) en s^-1 :

    !-- initialisations:
    iscp_pred = 0
    iscp_prey = 0

    !-- recherche de indice de organisme du predateur (zoo) et de une proie parmi le zoo:
    do mm = 1,nscp_zoo
        if (iscp_zoo(mm)==VAR(iipred)%idorg) then
            iscp_pred = mm
        elseif (iscp_zoo(mm)==VAR(iiprey)%idorg) then
            iscp_prey = mm
        endif
    enddo

    !-- si non trouvee au sein du zoo, recherche de indice de organisme de la proie dans 
    !   le vecteur iscp_phy puis iscp_bact:

    if( iscp_prey == 0) then
        if (allocated (iscp_phy)) then
            do mm = 1,nscp_phy
                if (iscp_phy(mm)==VAR(iiprey)%idorg) then
                    iscp_prey = mm + nscp_zoo
                    exit
                endif
            enddo
        endif
    endif
    if (iscp_prey ==0) then
        if (allocated (iscp_bac)) then
            do mm = 1,nscp_bac
                if (iscp_bac(mm)==VAR(iiprey)%idorg) then
                    iscp_prey = mm +nscp_zoo + nscp_phy
                    exit
                endif
            enddo
        endif
    endif


    !--recherche de la position pos_iprey (entre 1 et nbprey) relatif a 
    !  la conc en cellules de la proie iprey dans le vecteur indCell:
    pos_iprey=0

    do i=1,int_nbprey
        if (indCell(i)==iCell_prey) then
            pos_iprey = i
            exit
        endif
    enddo
    if (pos_iprey ==0) then
        write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_KOOIJ_KTW_cell',ili,jcol
        read(*,*)
    endif

    !-- vitesse specifique de grazing de la proie iscp_prey par le predateur iscp_pred: 
    do k=1,nz_max
      do j=ny_min,ny_max
        do i=nx_min,nx_max
           mu_graz(iscp_pred,iscp_prey)%val(i,j,k) =  &
             vfilt(pos_iprey) * pref_v(pos_iprey,i,j,k) * var(iCell_prey)%conc(i,j,k) /(ress_v(i,j,k)+1d-80) * &
             (ress(i,j,k))**beta/(1.d0 + (ress(i,j,k))**beta) 

 !--vitesse de grazing brute (Gamma_i x G x Zc ) en  proies m^-3 s-1
        if (var(iCell_prey)%conc(i,j,k)  > mu_graz(iscp_pred,iscp_prey)%val(i,j,k) * var(iCell_pred)%conc(i,j,k) * dt_bio) then
            if (iXprey == iCell_prey) then
              f_graz_KOOIJ_KTW_cell(i,j,k) = mu_graz(iscp_pred,iscp_prey)%val(i,j,k) * var(iCell_pred)%conc(i,j,k)                  ! mu x Z_Cell
            else
      !-- vitesse de grazing brute dans l'unite de iXprey / m^3 /s
              f_graz_KOOIJ_KTW_cell(i,j,k) = mu_graz(iscp_pred,iscp_prey)%val(i,j,k) * var(iCell_pred)%conc(i,j,k)  &                ! mu x Z_Cell
                   * (var(iXprey)%conc(i,j,k)/(var(iCell_prey)%conc(i,j,k)+1.d-80))                  ! mu x ZX
            endif
        else
              f_graz_KOOIJ_KTW_cell(i,j,k) = 0.d0
              mu_graz(iscp_pred,iscp_prey)%val(i,j,k) = 0.d0
        endif
       enddo
     enddo
    enddo

    ! write(*,*) 'fin de la fonction f_graz_KOOIJ_KTW_cell'

End function f_graz_KOOIJ_KTW_cell
!--------------------------------------------------------------------
