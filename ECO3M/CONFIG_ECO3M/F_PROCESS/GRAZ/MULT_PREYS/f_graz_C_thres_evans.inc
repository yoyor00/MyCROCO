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
   function f_graz_C_thres_evans(ili,jcol,ipred,iprey,Fmax,Imax,seuil,nbprey,&
     ind1,vF1,ind2,vF2,ind3,vF3,ind4,vF4,ind5,vF5,ind6,vF6,ind7,vF7,ind8,vF8)
!
!--------------------------------------------------------------------
! Cette fonction calcule la vitesse de predation d'un predateur
! en presence de PLUSIEURS PROIES potentielles, a l'aide d'une fonction a seuil.
! Elle est adaptée de la formulation d'Evans (1988) pour se ramener à une expression
! de type Holling II, ce qui explique le fait que le seuil soit adimensionnel et
! non pas homogène à une concentration de proies. Cf Gentleman et al. (2003) pour
! plus de détails.
!
! Reponse fonctionnelle: 
!     * de classe I à seuil
!     * tient compte des interferences entre proies
!     * preferences fixes 
!
! Attention:
! ----------
! 
! * cette fonction s''applique pour un FLUX  impliquant
!   un prédateur dont la concentration est donnee par VAR(ipred)%conc et 
!   une proie dont la concentration est donnee par VAR(iprey)%conc
!
! * Le calcul se base sur les concentrations en C du prédateur et des proies,
!  et on obtient une vitesse BRUTE de grazing Graz en [mol C] m-3 s^{-1}
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! Imax                : Vitesse d'ingestion max du predateur (identique pour toutes les proies), en [pred]^{-1} [proie]^{-1} s^{-1}
! nbprey              : nb de proies susceptibles d''être broutées par le prédateur
! Fmax                : vitesse maximum de filtration des proies en [pred]^{-1} s^{-1}
! indi                : indice des proies potentielles exprimées en C
! prefi		      : preference (fixe) pour la proie i
! seuil		      : adimensionnel
! f_graz_C_thres_evans  : vitesse de consommation de la proie iprey (exprimée en mol C m-3)
!                       par le prédateur ipred (exprimé en mol C m-3)
!
!
! Ref: Evans, G.T. (1988)  Limnol. & Oceanogr., 33(5):1027-1036
!      Lancelot, C. et al.  (2000)  Deep Sea Research I, 47:1621-1662
!      Gentleman, W. et al. (2003)  Deep Sea Research I, 50:2847-2875
!
! Melika Baklouti 24/04/2008
!--------------------------------------------------------------------
 USE DEF_TYPE
 USE MOD_FCHAIN
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL

 Implicit None

!-- fonction:
 Real(8)             :: f_graz_C_thres_evans(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8)             :: ipred,iprey,Fmax,nbprey,ind1,Imax,vF1
 Integer             :: ili,jcol
 Real(8), optional   :: ind2,ind3,ind4,ind5,ind6,ind7,ind8
 Real(8), optional   :: vF2,vF3,vF4,vF5,vF6,vF7,vF8
 
!--variables locales:
 Integer              :: int_nbprey,iiprey,iipred,vFilt_iprey,iscp_pred,iscp_prey
 Integer              :: iC_pred,iC_prey,iorg_prey,iorg_pred
 Integer              :: i,j,k,ii,jj,kk,ir,mm,istat
 Integer, Allocatable :: indC(:) ! indice de la conc en carbone des proies
 Real(8) ,Allocatable :: v_Filt(:,:,:,:)
 Real(8), Allocatable :: coeff_G(:,:,:)
 Real(8)              :: seuil
 Character(L_VAR)     :: elmt
 Character(L_CHAIN) :: chaine

!-- conversion en entier:
 int_nbprey = nint(nbprey)
 iipred = nint(ipred)
 iiprey = nint(iprey)

 
!-- Allocations dynamiques:

Allocate (v_Filt(int_nbprey,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de v_filt dans la fonction f_graz_hol2'

Allocate (indC(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de indC dans la fonction f_graz_hol2'

Allocate (coeff_G(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de coeff_G dans la fonction f_graz_hol2'


!-- initialisations:
indC = 0
f_graz_C_thres_evans=0.d0

!
!-- premiere proie (obligatoire)
 indC(1) = nint(ind1)
 v_Filt(1,:,:,:)= vF1
 
 Select case (int_nbprey)
 case (2)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2
 
 case (3)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2 

 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3
 
 case (4)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2
 
 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3 
 
 indC(4) = nint(ind4)
 v_Filt(4,:,:,:)= vF4
 
 case (5)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2

 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3
 
 indC(4) = nint(ind4)
 v_Filt(4,:,:,:)= vF4 
 
 indC(5) = nint(ind5)
 v_Filt(5,:,:,:)= vF5 

 case (6)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2

 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3
 
 indC(4) = nint(ind4)
 v_Filt(4,:,:,:)= vF4 
 
 indC(5) = nint(ind5)
 v_Filt(5,:,:,:)= vF5 

indC(6) = nint(ind6)
 v_Filt(6,:,:,:)= vF6 

 case (7)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2

 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3
 
 indC(4) = nint(ind4)
 v_Filt(4,:,:,:)= vF4 
 
 indC(5) = nint(ind5)
 v_Filt(5,:,:,:)= vF5 
 
 indC(6) = nint(ind6)
 v_Filt(6,:,:,:)= vF6 

 indC(7) = nint(ind7)
 v_Filt(7,:,:,:)= vF7 

 case (8)
 indC(2) = nint(ind2)
 v_Filt(2,:,:,:)= vF2

 indC(3) = nint(ind3)
 v_Filt(3,:,:,:)= vF3
 
 indC(4) = nint(ind4)
 v_Filt(4,:,:,:)= vF4 
 
 indC(5) = nint(ind5)
 v_Filt(5,:,:,:)= vF5 
 
indC(6) = nint(ind6)
v_Filt(6,:,:,:)= vF6 

indC(7) = nint(ind7)
v_Filt(7,:,:,:)= vF7 

indC(8) = nint(ind8)
v_Filt(8,:,:,:)= vF8 
 end select

!------------- Calcul des vitesses de filtration ----------------------


!-- Alternatives possibles pour la vitesse de filtration:

!  v_Filt(:,:,:,:) = 0.d0
! do ii= 1, int_nbprey
!  do i=nx_min,nx_max
!    do j=ny_min,ny_max
!     do k=1,nz_max
!1/ Cas d'une vitesse variant avec la concentration :
!      v_Filt(ii,i,j,k) =  max(1.5d-4, Fmax * (1 - 1.04* (sum_prey(i,j,k))**2/((sum_prey(i,j,k))**2 + 0.0063*5)))
!      if (v_Filt(ii,i,j,k) > Fmax) v_Filt(ii,i,j,k) = Fmax
!2/ Cas d'une vitesse constante :
!      v_Filt(ii,i,j,k) =   Fmax 
!   enddo
!  enddo
! enddo
!enddo

!-- numeros d'organisme de la proie et du predateur:
 iorg_prey = VAR(iiprey)%idorg
 iorg_pred = VAR(iipred)%idorg

!-- Indice correspondant a la conc. en carbone de la proie:
 chaine = 'C'
 iC_prey = f_idorg2id(iorg_prey,chaine)
 iC_pred = f_idorg2id(iorg_pred,chaine)


!-- Iteration sur l'ensemble des proies concernees pour le calcul de 
!   la ressource totale pour ce zoo:
 coeff_G = 0.d0

!-- Ressource totale exprimée en C:
!-- modele Holling II:

 do ii=1,int_nbprey
    coeff_G(:,:,:) = coeff_G(:,:,:) + v_Filt(ii,:,:,:)/Imax * var(indC(ii))%conc(:,:,:) 
 enddo

!-- Création et sauvegarde du tableau des vitesses specifiques de grazing du prédateur i
!   sur la proie j (Gamma_j x G) en s^-1 :

!-- initialisations:
 iscp_pred = 0
 iscp_prey = 0

!-- recherche de l''indice de l organisme du predateur (zoo) et d''une proie parmi le zoo:
 do mm = 1,nscp_zoo
   if (iscp_zoo(mm)==VAR(iipred)%idorg) then
     iscp_pred = mm
   elseif (iscp_zoo(mm)==VAR(iiprey)%idorg) then
     iscp_prey = mm
  endif
 enddo

!-- si non trouvé au sein du zoo, recherche de l''indice de l organisme de la proie dans 
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


!--recherche du numero vFilt_iprey (entre 1 et nbprey) relatif a 
!  la conc en carbone de la proie dans le vecteur vFilt:
 vFilt_iprey=0

 do i=1,int_nbprey
   if (indC(i)==iC_prey) then
       vFilt_iprey = i
       exit
   endif
 enddo
 if (vFilt_iprey ==0) then
   write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_C_thres_evans'
   read(*,*)
 endif

!-- vitesse specifique de grazing de la proie iscp_prey par le predateur iscp_pred: 

  do i=nx_min,nx_max
    do j=ny_min,ny_max
     do k=1,nz_max
       if (coeff_G(i,j,k) > seuil) then
          mu_graz(iscp_pred,iscp_prey)%val(i,j,k) =  (coeff_G(i,j,k) - seuil) / (1. + coeff_G(i,j,k) - seuil)  &
                            * v_Filt(vFilt_iprey,i,j,k)  * var(iC_prey)%conc(i,j,k)/(coeff_G(i,j,k)+1.d-20) 
       else
         mu_graz(iscp_pred,iscp_prey)%val(i,j,k) = 0.d0
       endif
     enddo
   enddo
  enddo
 
!-- vitesse de grazing brute (Gamma_i x G x Zc ) en mol C m^-3 s-1

f_graz_C_thres_evans(:,:,:) =   mu_graz(iscp_pred,iscp_prey)%val(:,:,:)   * var(iC_pred)%conc(:,:,:)                  ! mu x Z_C

 Deallocate (v_Filt,indC,coeff_G)

! write(*,*) 'fin de la fonction f_graz_C_thres_evans'
 
 End function f_graz_C_thres_evans
!--------------------------------------------------------------------
