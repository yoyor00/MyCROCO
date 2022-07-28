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
   function f_graz_C_hol2_chesson(ili,jcol,ipred,iprey,Fmax,nbprey,&
     ind1,Imax1,ind2,Imax2,ind3,Imax3,ind4,Imax4,ind5,Imax5,ind6,Imax6,ind7,Imax7,ind8,Imax8)
!
!--------------------------------------------------------------------
! Cette fonction calcule la vitesse de predation d'un predateur
! en presence de PLUSIEURS PROIES potentielles
! Reponse fonctionnelle: 
!     * de type Holling II
!     * tient compte des interferences entre proies
!     * preferences fixes calculees selon le modele de Chesson (1983)
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
!  qui tient compte des pertes dues a  l'egestion, et au sloppy feeding
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! Imaxi               : Vitesse d'ingestion max du predateur de la proie i, (en [pred]^{-1} [proie]^{-1} s^{-1})
! nbprey              : nb de proies susceptibles d''etre broutees par le predateur (entre 1 et 8)
! Fi                  : vitesse de filtration des proies (en [pred]^{-1} s^{-1})
! Fmax                : vitesse maximum de filtration des proies (en [pred]^{-1} s^{-1})
! indi                : indice des proies potentielles exprimees en C
! f_graz_C_hol2_chesson  : vitesse de consommation de la proie iprey (en mol C m-3)
!                       par le predateur ipred (en mol C m-3)
!
!
! Ref: Holling, C. (1959)  Canadian Entomologist, 91:385-398
!      Chesson, J. (1983)  Ecology, 64:1297-1304
!
! Melika Baklouti 4/06/2008
!--------------------------------------------------------------------
 USE DEF_TYPE
 USE MOD_FCHAIN
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL

 Implicit None

!-- fonction:
 Real(8)             :: f_graz_C_hol2_chesson(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8)             :: ipred,iprey,Fmax,nbprey,ind1,Imax1
 Integer             :: ili,jcol
 Real(8), optional   :: ind2,ind3,ind4,ind5,ind6,ind7,ind8
 Real(8), optional   :: Imax2,Imax3,Imax4,Imax5,Imax6,Imax7,Imax8
 
!--variables locales:
 Integer              :: int_nbprey,iiprey,iipred,vpref_iprey,iscp_pred,iscp_prey
 Integer              :: iC_pred,iC_prey,iorg_prey,iorg_pred
 Integer              :: i,j,k,ii,jj,kk,ir,mm,istat
 Integer, Allocatable :: indC(:) ! indice de la conc en carbone des proies
 Real(8) ,Allocatable :: v_Filt(:,:,:,:),v_pref(:,:,:,:),sum_filt(:,:,:),Imax(:)
 Real(8), Allocatable :: ress(:,:,:),num_G(:,:,:),denom_G(:,:,:),sum_prey(:,:,:)
 Character(L_VAR)     :: elmt
 Character(L_CHAIN) :: chaine

!-- conversion en entier:
 int_nbprey = nint(nbprey)
 iipred = nint(ipred)
 iiprey = nint(iprey)

 
!-- Allocations dynamiques:
Allocate (v_pref(int_nbprey,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de v_pref dans la fonction f_graz_hol2'

Allocate (v_Filt(int_nbprey,nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de v_pref dans la fonction v_Filt'

Allocate (sum_filt(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de sum_filt dans la fonction f_graz_hol2'

Allocate (sum_prey(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de sum_prey dans la fonction f_graz_hol2'

Allocate (indC(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de indC dans la fonction f_graz_hol2'

Allocate (Imax(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de Imax dans la fonction f_graz_hol2'

Allocate (ress(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de ress dans la fonction f_graz_hol2'

Allocate (num_G(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de num_G dans la fonction f_graz_hol2'

Allocate (denom_G(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de denom_G dans la fonction f_graz_hol2'


!-- initialisations:
indC = 0
Imax = 0.d0
f_graz_C_hol2_chesson=0.d0

!
!-- premiere proie (obligatoire)
 indC(1) = nint(ind1)
 Imax(1)= Imax1
 
 Select case (int_nbprey)
 case (2)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 
 case (3)
 indC(2) = nint(ind2)
 Imax(2)= Imax2 

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 
 case (4)
 indC(2) = nint(ind2)
 Imax(2)= Imax2
 
 indC(3) = nint(ind3)
 Imax(3)= Imax3 
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4
 
 case (5)
 indC(2) = nint(ind2)
 Imax(2)= Imax2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 

 case (6)
 indC(2) = nint(ind2)
 Imax(2)= Imax2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 

indC(6) = nint(ind6)
 Imax(6)= Imax6 

 case (7)
 indC(2) = nint(ind2)
 Imax(2)= Imax2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 
 indC(6) = nint(ind6)
 Imax(6)= Imax6 

 indC(7) = nint(ind7)
 Imax(7)= Imax7 

 case (8)
 indC(2) = nint(ind2)
 Imax(2)= Imax2

 indC(3) = nint(ind3)
 Imax(3)= Imax3
 
 indC(4) = nint(ind4)
 Imax(4)= Imax4 
 
 indC(5) = nint(ind5)
 Imax(5)= Imax5 
 
indC(6) = nint(ind6)
Imax(6)= Imax6 

indC(7) = nint(ind7)
Imax(7)= Imax7 

indC(8) = nint(ind8)
Imax(8)= Imax8 
 end select

!------------- Calcul des vitesses de filtration ----------------------

 v_Filt(:,:,:,:) = 0.d0
 sum_prey = 0.d0

! Calcul de la conc. totale de proies:
 do ii=1,int_nbprey
  sum_prey(:,:,:) = sum_prey(:,:,:) + var(indC(ii))%conc(:,:,:) 
 enddo

 do ii=1,int_nbprey
  do i=nx_min,nx_max
    do j=ny_min,ny_max
     do k=1,nz_max
!1/ Cas d'une vitesse variant avec la concentration (introduire la correlation desiree):
!      v_Filt(ii,i,j,k) =  max(1.5d-4, Fmax * (1 - 1.04* (sum_prey(i,j,k))**2/((sum_prey(i,j,k))**2 + 0.0063*5)))
!      if (v_Filt(ii,i,j,k) > Fmax) v_Filt(ii,i,j,k) = Fmax
!2/ Cas d'une vitesse constante :
      v_Filt(ii,i,j,k) =   Fmax 
   enddo
  enddo
 enddo
 enddo
!------------  Fin du calcul des vitesses de filtration ---------------

!-- Somme des vitesses de filtration des organismes presents (clearance rates):
 sum_filt = 0.d0
 do ii=1,int_nbprey
  do i=nx_min,nx_max
    do j=ny_min,ny_max
     do k=1,nz_max
            sum_filt(i,j,k) = sum_filt(i,j,k) + v_Filt(ii,i,j,k)
   enddo
  enddo
 enddo
enddo

!-- Initialisation du vecteur des preferences :
v_pref=0.d0

!-- Calcul des preferences:
 do ii=1,int_nbprey
  do i=nx_min,nx_max
    do j=ny_min,ny_max
     do k=1,nz_max
          if (sum_filt(i,j,k)>1.d-20) then
              v_pref(ii,i,j,k) = v_Filt(ii,i,j,k)/sum_filt(i,j,k)
         endif
     enddo
   enddo
  enddo
enddo

!-- numeros d'organisme de la proie et du predateur:
iorg_prey = VAR(iiprey)%idorg
iorg_pred = VAR(iipred)%idorg

!-- Indice correspondant a la conc. en carbone de la proie:
 chaine = 'C'
 iC_prey = f_idorg2id(iorg_prey,chaine)
 iC_pred = f_idorg2id(iorg_pred,chaine)


!-- Iteration sur l'ensemble des proies concernees pour le calcul de 
!   la ressource totale pour ce zoo:
 ress = 0.d0
 denom_G = 1.d0 !-- car denominateur de G = 1+ ...
 num_G= 0.d0

!-- Ressource totale exprimee en C:
!-- modele Holling II:

 do ii=1,int_nbprey
    ress(:,:,:)    = ress(:,:,:) + v_pref(ii,:,:,:) * var(indC(ii))%conc(:,:,:) !--ressource
    denom_G(:,:,:) = denom_G(:,:,:) + v_Filt(ii,:,:,:)/Imax(ii) * (var(indC(ii))%conc(:,:,:) )
    num_G(:,:,:)   = num_G(:,:,:) + v_Filt(ii,:,:,:) * (var(indC(ii))%conc(:,:,:))
 enddo

!-- Creation et sauvegarde du tableau des vitesses specifiques de grazing du préedateur i
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

!-- si non trouvee au sein du zoo, recherche de l''indice de l organisme de la proie dans 
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


!--recherche du numero vpref_iprey (entre 1 et nbprey) relatif a 
!  la conc en carbone de la proie dans le vecteur vpref:
 vpref_iprey=0

 do i=1,int_nbprey
   if (indC(i)==iC_prey) then
       vpref_iprey = i
       exit
   endif
 enddo
 if (vpref_iprey ==0) then
   write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_C_hol2_chesson'
   read(*,*)
 endif

!-- vitesse specifique de grazing de la proie iscp_prey par le predateur iscp_pred: 

 mu_graz(iscp_pred,iscp_prey)%val(:,:,:) =  (v_pref(vpref_iprey,:,:,:) * var(iC_prey)%conc(:,:,:))/(ress(:,:,:)+1.d-20) &  ! Gamma_i 
                                                           * num_G(:,:,:)/ denom_G(:,:,:)         ! G


!-- vitesse de grazing brute (Gamma_i x G x Zc ) en mol C m^-3 s-1

f_graz_C_hol2_chesson(:,:,:) =   mu_graz(iscp_pred,iscp_prey)%val(:,:,:)   * var(iC_pred)%conc(:,:,:)                  ! mu x Z_C


 Deallocate (v_pref,v_Filt,indC,sum_filt,Imax,ress,num_G,denom_G,sum_prey)

! write(*,*) 'fin de la fonction f_graz_C_hol2_chesson'
 
 End function f_graz_C_hol2_chesson
!--------------------------------------------------------------------
