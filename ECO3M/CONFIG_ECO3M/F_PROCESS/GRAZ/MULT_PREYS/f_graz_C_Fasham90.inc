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
!
   function f_graz_C_Fasham90(ili,jcol,ipred,iprey,Gmax,kzoo,cminprey,nbprey,&
     ind1,Pref1,ind2,Pref2,ind3,Pref3,ind4,Pref4,ind5,Pref5,ind6,Pref6,ind7,Pref7,ind8,Pref8)
!
!--------------------------------------------------------------------
! Cette fonction calcule la vitesse de predation d'un predateur
! en presence de PLUSIEURS PROIES potentielles
! Reponse fonctionnelle: 
!     * de type classe 3 multiple functional responses (cf. Gentleman et al. 2003 DSR 50: 2847-2875)
!     * tient compte des interferences entre proies
!     * preferences ponderees,i.e. variables et dependantes de la concentration en proies (permet la coexistences des proies)
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
! Prefi               : preference fixe du predateur pour la proie i, (sans unite [0,1])
! nbprey              : nb de proies susceptibles d''etre broutees par le predateur (entre 1 et 8)
! cminprey            : concentration en dessous de laquelle il n''y a pas pr�dation (en mmolC m-3)
! Gmax                : vitesse maximum de broutage (en s-1)
! indi                : indice des proies potentielles exprimees en C
! kzoo                : constante de demi-saturation du grazing (en mol C m-3)
! f_graz_C_Fasham90  : vitesse de consommation de la proie iprey (en mol C m-3)
!                       par le predateur ipred (en mol C m-3)
!
!
! Ref: Fasham et al. (1990) Journal of Marine Systems, 48: 591-539
!
! Frederic Diaz 30/09/2008
!--------------------------------------------------------------------
 Implicit None

!-- fonction:
 Real(8)             :: f_graz_C_Fasham90(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8)             :: ipred,iprey,Gmax,nbprey,ind1,Pref1,kzoo,cminprey
 Integer             :: ili,jcol
 Real(8), optional   :: ind2,ind3,ind4,ind5,ind6,ind7,ind8
 Real(8), optional   :: Pref2,Pref3,Pref4,Pref5,Pref6,Pref7,Pref8
 
!--variables locales:
 Integer              :: int_nbprey,iiprey,iipred,vpref_iprey,iscp_pred,iscp_prey
 Integer              :: iC_pred,iC_prey,iorg_prey,iorg_pred
 Integer              :: i,j,k,ii,jj,kk,ir,mm,istat
 Integer, Allocatable :: indC(:) ! indice de la conc en carbone des proies
 Real(8), Allocatable :: ress(:,:,:),denom_G(:,:,:),Pref(:)
 Character(L_VAR)     :: elmt
 Character(L_CHAIN)   :: chaine

!-- conversion en entier:
 int_nbprey = nint(nbprey)
 iipred = nint(ipred)
 iiprey = nint(iprey)

 
!-- Allocations dynamiques:

Allocate (indC(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de indC dans la fonction f_graz_Fasham'

Allocate (Pref(int_nbprey),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de Pref dans la fonction f_graz_Fasham'

Allocate (ress(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de num_G dans la fonction f_graz_Fasham'

Allocate (denom_G(nx_min:nx_max,ny_min:ny_max,nz_max),STAT=istat)
if (istat /=0) write(*,*) 'probleme d''allocation de denom_G dans la fonction f_graz_Fasham'


!-- initialisations:
indC = 0
Pref = 0.d0
f_graz_C_Fasham90=0.d0

!
!-- premiere proie (obligatoire)
 indC(1) = nint(ind1)
 Pref(1)= Pref1
 
 Select case (int_nbprey)
 case (2)
 indC(2) = nint(ind2)
 Pref(2)= Pref2
 
 case (3)
 indC(2) = nint(ind2)
 Pref(2)= Pref2 

 indC(3) = nint(ind3)
 Pref(3)= Pref3
 
 case (4)
 indC(2) = nint(ind2)
 Pref(2)= Pref2
 
 indC(3) = nint(ind3)
 Pref(3)= Pref3 
 
 indC(4) = nint(ind4)
 Pref(4)= Pref4
 
 case (5)
 indC(2) = nint(ind2)
 Pref(2)= Pref2

 indC(3) = nint(ind3)
 Pref(3)= Pref3
 
 indC(4) = nint(ind4)
 Pref(4)= Pref4 
 
 indC(5) = nint(ind5)
 Pref(5)= Pref5 

 case (6)
 indC(2) = nint(ind2)
 Pref(2)= Pref2

 indC(3) = nint(ind3)
 Pref(3)= Pref3
 
 indC(4) = nint(ind4)
 Pref(4)= Pref4 
 
 indC(5) = nint(ind5)
 Pref(5)= Pref5 

 indC(6) = nint(ind6)
 Pref(6)= Pref6 

 case (7)
 indC(2) = nint(ind2)
 Pref(2)= Pref2

 indC(3) = nint(ind3)
 Pref(3)= Pref3
 
 indC(4) = nint(ind4)
 Pref(4)= Pref4 
 
 indC(5) = nint(ind5)
 Pref(5)= Pref5 
 
 indC(6) = nint(ind6)
 Pref(6)= Pref6 

 indC(7) = nint(ind7)
 Pref(7)= Pref7 

 case (8)
 indC(2) = nint(ind2)
 Pref(2)= Pref2

 indC(3) = nint(ind3)
 Pref(3)= Pref3
 
 indC(4) = nint(ind4)
 Pref(4)= Pref4 
 
 indC(5) = nint(ind5)
 Pref(5)= Pref5 
 
 indC(6) = nint(ind6)
 Pref(6)= Pref6 

 indC(7) = nint(ind7)
 Pref(7)= Pref7 

 indC(8) = nint(ind8)
 Pref(8)= Pref8 
 end select

!-- numeros organisme de la proie et du predateur:
iorg_prey = VAR(iiprey)%idorg
iorg_pred = VAR(iipred)%idorg

!-- Indice correspondant a la conc. en carbone de la proie:
 chaine = 'C'
 iC_prey = f_idorg2id(iorg_prey,chaine)
 iC_pred = f_idorg2id(iorg_pred,chaine)


!-- Iteration sur l''ensemble des proies concernees pour le calcul de 
!   la ressource totale pour ce zoo:
 ress = 0.d0
 denom_G = 0.d0 

!-- Ressource totale exprimee en C:
!-- modele Holling II:

 do ii=1,int_nbprey
    ress(:,:,:)    = ress(:,:,:) + Pref(ii) * var(indC(ii))%conc(:,:,:) !--ressource
    denom_G(:,:,:) = denom_G(:,:,:) + Pref(ii) * (var(indC(ii))%conc(:,:,:) *var(indC(ii))%conc(:,:,:) )
    
 enddo
!write(*,*) 'coucou'
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
   write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_C_Fasham90'
   read(*,*)
 endif

!-- vitesse specifique de grazing de la proie iscp_prey par le predateur iscp_pred: 
 do i=nx_min,nx_max
 do j=ny_min,ny_max
  do k=1,nz_max
     if (var(iC_prey)%conc(i,j,k).ge.cminprey) then
 mu_graz(iscp_pred,iscp_prey)%val(i,j,k) =  Gmax*(Pref(vpref_iprey) * (var(iC_prey)%conc(i,j,k))**2.)/(kzoo*(ress(i,j,k))+ &   
                                                           denom_G(i,j,k))       
!write(*,*) 'Pref(vpref_iprey)', Pref(vpref_iprey),mu_graz(iscp_pred,iscp_prey)%val(:,:,:), IndC(vpref_iprey) 
f_graz_C_Fasham90(i,j,k) =   mu_graz(iscp_pred,iscp_prey)%val(i,j,k)   * var(iC_pred)%conc(i,j,k)                  ! mu x Z_C
else 
mu_graz(iscp_pred,iscp_prey)%val(i,j,k) = 0.

!-- vitesse de grazing brute (Gamma_i x G x Zc ) en mol C m^-3 s-1

f_graz_C_Fasham90(i,j,k) =   0.               ! mu x Z_C

endif
enddo
enddo
enddo

 Deallocate (indC,Pref,ress,denom_G)

! write(*,*) 'fin de la fonction f_graz_C_Fasham90'
 
 End function f_graz_C_Fasham90
!--------------------------------------------------------------------
