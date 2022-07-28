!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Universite de la Mediterraneee
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
   function f_graz_C_hol2_murdoch(ili,jcol,ipred,iprey,Gmax,nbprey,&
     ind1,Imax1,vF1,ind2,Imax2,vF2,ind3,Imax3,vF3,ind4,Imax4,vF4,ind5,Imax5,vF5,&
     ind6,Imax6,vF6,ind7,Imax7,vF7,ind8,Imax8,vF8)
!
!--------------------------------------------------------------------
! Cette fonction calcule la vitesse de predation d'un predateur
! dont la concentration est donnee par VAR(ipred)%conc sur une
! proie dont la concentration est donnee par VAR(iprey)%conc
! en presence de PLUSIEURS PROIES potentielles selon le modele de
! Holling II adapte au cas de  plusieurs proies.
!
! La reponse fonctionnelle est: 
!     * basee sur la reponse fonctionnelle de type Holling II (etablie pour une proie)
!     * tient compte des interferences entre proies
!     * preferences fixes calculees selon le modele de Murdoch (1973)
!
! Remarques IMPORTANTES :
! ---------------------
!
! * Le calcul se base sur les concentrations en carbone (C) du predateur et des proies,
!  et on obtient une vitesse BRUTE de grazing Graz en [mol C] m-3 s^{-1}
! * Si les biomasses ne sont pas exprimees en C dans le modele, tout autre element
!   ou nb de cellules peut etre utilise  a condition de respecter l'homogeneite.
! * Si le modele est multi-element, les flux de predation pour les elements autres que l'element
!   "principal" doivent etre calcules avec la fonction f_graz_C_hol2_murdoch_X 
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! Imaxi               : Vitesse d'ingestion max du predateur de la proie i, en [pred]^{-1} [proie]^{-1} s^{-1}
! nbprey              : nb de proies susceptibles d''etre broutees par le predateur
! vFi                 : vitesse de filtration ou plus generalement vitesse de capture des proies en [pred]^{-1} s^{-1}
! Gmax                : vitesse specifique maximum de grazing  s^{-1}
! indi                : indice des proies potentielles exprimees dans l'unite de biomasse adequate (generalement C)
! f_graz_C_hol2_murdoch  : vitesse de consommation de la proie iprey (exprimee generalement en mol C m-3)
!                       par le predateur ipred (exprime en mol C m-3), en mol C m^-3 s-1
!
!
! Ref: Holling, C. (1959)  Canadian Entomologist, 91:385-398
!      Murdoch, W. W. (1973)  J. Appl. Ecol., 10:335-342
!      Tian, R. C. (2006) Ecol. Model., 193:363-386
!
! Melika Baklouti 29/01/2008
!--------------------------------------------------------------------
 USE DEF_TYPE
 USE MOD_FCHAIN
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL

 Implicit None

!-- fonction:
 Real(8)             :: f_graz_C_hol2_murdoch(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8)             :: ipred,iprey,Gmax,nbprey,ind1,Imax1,vF1
 Integer             :: ili,jcol
 Real(8), optional   :: ind2,ind3,ind4,ind5,ind6,ind7,ind8
 Real(8), optional   :: Imax2,Imax3,Imax4,Imax5,Imax6,Imax7,Imax8
 Real(8), optional   :: vF2,vF3,vF4,vF5,vF6,vF7,vF8
 
!--variables locales:
 Integer              :: int_nbprey,iiprey,iipred,vfilt_iprey,iscp_pred,iscp_prey
 Integer              :: iC_pred,iC_prey,iorg_prey,iorg_pred
 Integer              :: i,j,k,ii,jj,kk,ir,mm,istat
 Integer, Allocatable :: indC(:) ! indice de la conc en carbone des proies
 Real(8) ,Allocatable :: v_Filt(:,:,:,:),Imax(:)
 Real(8), Allocatable :: denom_G(:,:,:)
 Character(L_VAR)     :: elmt
 Character(L_CHAIN) :: chaine

!-- conversion en entier:
 int_nbprey = nint(nbprey)
 iipred = nint(ipred)
 iiprey = nint(iprey)

 
!-- Allocations dynamiques:
Allocate (v_Filt(int_nbprey,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de v_Filt dans la fonction v_Filt'

Allocate (indC(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de indC dans la fonction f_graz_hol2'

Allocate (Imax(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de Imax dans la fonction f_graz_hol2'


Allocate (denom_G(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de denom_G dans la fonction f_graz_hol2'


!-- initialisations:
indC = 0
Imax = 0.d0
!
!-- premiere proie (obligatoire)
 indC(1) = nint(ind1)
 Imax(1)= Imax1
 v_Filt(1,:,:,:)= vf1
 
 Select case (int_nbprey)

 case (2)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2
 
 case (3)
 indC(2) = nint(ind2)
 Imax(2)= Imax2 
 v_Filt(2,:,:,:)= vf2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 v_Filt(3,:,:,:)= vf3
 
 case (4)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2
 
 indC(3) = nint(ind3)
 Imax(3)= Imax3 
 v_Filt(3,:,:,:)= vf3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4
 v_Filt(4,:,:,:)= vf4
 
 case (5)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 v_Filt(3,:,:,:)= vf3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 v_Filt(4,:,:,:)= vf4
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 v_Filt(5,:,:,:)= vf5

 case (6)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 v_Filt(3,:,:,:)= vf3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 v_Filt(4,:,:,:)= vf4
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 v_Filt(5,:,:,:)= vf5

indC(6) = nint(ind6)
 Imax(6)= Imax6 
 v_Filt(6,:,:,:)= vf6

 case (7)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 v_Filt(3,:,:,:)= vf3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 v_Filt(4,:,:,:)= vf4
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 v_Filt(5,:,:,:)= vf5
 
 indC(6) = nint(ind6)
 Imax(6)= Imax6 
 v_Filt(6,:,:,:)= vf6

 indC(7) = nint(ind7)
 Imax(7)= Imax7 
 v_Filt(7,:,:,:)= vf7

 case (8)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 v_Filt(2,:,:,:)= vf2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 v_Filt(3,:,:,:)= vf3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 v_Filt(4,:,:,:)= vf4
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 v_Filt(5,:,:,:)= vf5
 
indC(6) = nint(ind6)
Imax(6)= Imax6 
 v_Filt(6,:,:,:)= vf6

indC(7) = nint(ind7)
Imax(7)= Imax7 
 v_Filt(7,:,:,:)= vf7

indC(8) = nint(ind8)
Imax(8)= Imax8 
v_Filt(8,:,:,:)= vf8
 end select


!-- numeros d'organisme de la proie et du predateur:
iorg_prey = VAR(iiprey)%idorg
iorg_pred = VAR(iipred)%idorg

!-- Indice correspondant a la conc. en carbone de la proie:
 chaine = 'C'
 iC_prey = f_idorg2id(iorg_prey,chaine)
 iC_pred = f_idorg2id(iorg_pred,chaine)


!-- Iteration sur l'ensemble des proies concernees pour le calcul du denominateur:
 denom_G = 1.d0 !-- car denominateur de G = 1+ ...

!-- Ressource totale exprimee en C:
!-- modele Holling II:

 do ii=1,int_nbprey
   denom_G(:,:,:) = denom_G(:,:,:) + v_Filt(ii,:,:,:)/Imax(ii) * var(indC(ii))%conc(:,:,:)
 enddo

!-- Creation et sauvegarde du tableau des vitesses specifiques de grazing du predateur i
!   sur la proie j (G_i) en s^-1 :

!-- initialisations:
 iscp_pred = 0
 iscp_prey = 0

!-- recherche de l''indice de l organisme du predateur (zoo) et d''une proie parmi le zoo
!   dans le tableau iscp_zoo :
 do mm = 1,nscp_zoo
   if (iscp_zoo(mm)==VAR(iipred)%idorg) then
     iscp_pred = mm
   elseif (iscp_zoo(mm)==VAR(iiprey)%idorg) then
     iscp_prey = mm
  endif
 enddo

!-- si non trouve au sein du zoo, recherche de l''indice de l organisme de la proie dans 
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


!--recherche du numero vfilt_iprey (entre 1 et nbprey) relatif a 
!  la proie dans le vecteur v_filt:
 vfilt_iprey=0

 do i=1,int_nbprey
   if (indC(i)==iC_prey) then
       vfilt_iprey = i
       exit
   endif
 enddo
 if (vfilt_iprey ==0) then
   write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_C_hol2_murdoch'
   read(*,*)
 endif

!-- vitesse specifique de grazing de la proie iscp_prey par le predateur iscp_pred: 

 mu_graz(iscp_pred,iscp_prey)%val(:,:,:) =  v_filt(vfilt_iprey,:,:,:) * var(iC_prey)%conc(:,:,:)/ denom_G(:,:,:)     ! G (en s^-1)


!-- vitesse de grazing brute (G x Zc ) en mol C m^-3 s-1

f_graz_C_hol2_murdoch(:,:,:) =   mu_graz(iscp_pred,iscp_prey)%val(:,:,:)   * var(iC_pred)%conc(:,:,:)       ! mu x Z_C


 Deallocate (v_Filt,indC,Imax,denom_G)

 !write(*,*) 'fin de la fonction f_graz_C_hol2_murdoch',ili,jcol,vfilt_iprey,f_graz_C_hol2_murdoch(1:1:1)
 
 End function f_graz_C_hol2_murdoch
!--------------------------------------------------------------------
