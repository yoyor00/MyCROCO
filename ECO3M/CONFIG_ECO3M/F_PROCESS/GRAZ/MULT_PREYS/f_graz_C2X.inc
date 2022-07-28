!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Universite de la Mediterranee
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
            function f_graz_C2X(ili,jcol,ipred,iprey)
!
!--------------------------------------------------------------------
! Cette fonction est a utiliser paralllelement a une des fonctions f_graz_C_* dans laquelle
! le flux de predation exprime en carbone par unite de volume et de temps a ete calcule 
! mu_graz 
! 
! Cette fonction permet d''exprimer la consommation d'une proie exprimee en un element X par un predateur
! dont la concentration est aussi exprimee en X, en presence d''une ou plusieurs proies potentielles
!
! La fonction f_graz_C2X se borne donc a convertir le flux de predation en carbone 
! en une unite de biomasse differente
!
! Attention:
! ----------
!
! * Le calcul se base sur les concentrations en C du predateur et des proies,
!  et on obtient une vitesse BRUTE de grazing  en [mol X] m-3 s^{-1}
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! ipred    : indice de la variable d''etat representant la concentration X du predateur
! iprey    : indice de la variable d''etat representant la concentration X de la proie
! mu_graz  : ce tableau est une variable globale. Il contient les vitesses specifiques (en s^-1) de grazing des predateurs.
!            Ce tableau est rempli lors de l''appel de la fonction  f_graz_C_* (e.g. f_graz_C_hol2_murdoch)
!            l''element (i,j) de ce tableau represente la vitesse de grazing du predateur VAR(i) sur la proie VAR(j)
!
! f_graz_C2X  : vitesse de consommation (en mol X m^-3 s-1) de la proie iprey par le predateur ipred
!
!
! Melika Baklouti 29/01/2008
!--------------------------------------------------------------------
 Implicit None

!-- fonction:
 Real(8)             :: f_graz_C2X(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Integer             :: ili,jcol
 Real(8)             :: ipred,iprey
!--variables locales:
 Integer             :: iX_pred,iX_prey,iC_pred,iC_prey,iscp_prey,iscp_pred,mm
 Integer             :: iorg_prey,iorg_pred,iorg_ili,iorg_jcol
 Character(L_VAR)    :: elmt
 Character(L_VAR)    :: chain

!-- indices des conc en element X du predateur et de la proie
iX_pred=nint(ipred)
iX_prey=nint(iprey)

!--recherche des numeros d''organisme de la proie et du predateur :
 iorg_prey= var(iX_prey)%idorg
 iorg_pred= var(iX_pred)%idorg


!-- concentration en C de la proie et du predateur:
 chain='C'
 iC_prey = f_idorg2id(iorg_prey,chain)
 iC_pred = f_idorg2id(iorg_pred,chain)


!-- Rappel: mu_graz = vitesse specifique de grazing du predateur (Gamma x G) en s^-1
!          calcule pour chaque couple pred/proie dans la fonction f_graz_C_*

!-- Initialisations:
 iscp_pred = 0
 iscp_prey = 0

!-- Recherche de l''indice de l organisme du predateur (zoo):
 do mm = 1,nscp_zoo
   if (iscp_zoo(mm)==iorg_pred) then
     iscp_pred = mm
   elseif (iscp_zoo(mm)==iorg_prey) then
     iscp_prey = mm
  endif
 enddo
!-- Si la proie n''a pas ete trouve dans le zoo, on recherche l''indice de l organisme 
!   de la proie dans le vecteur iscp_phy puis au sein de iscp_bact:

if( iscp_prey == 0) then
 if (allocated (iscp_phy)) then
  do mm = 1,nscp_phy
   if (iscp_phy(mm)==iorg_prey) then
     iscp_prey = mm + nscp_zoo
     exit
  endif
 enddo
 endif
endif
 if (iscp_prey ==0) then
  if (allocated (iscp_bac)) then
   do mm = 1,nscp_bac
     if (iscp_bac(mm)==iorg_prey) then
        iscp_prey = mm +nscp_zoo + nscp_phy
        exit
     endif
   enddo
  endif
 endif

!-- vitesse de grazing en mol X m^-3 s-1

        f_graz_C2X(:,:,:) =   mu_graz(iscp_pred,iscp_prey)%val(:,:,:)   * var(iC_pred)%conc(:,:,:)    &
               * var(iX_prey)%conc(:,:,:) /  (var(iC_prey)%conc(:,:,:)  + 1.d-20)           !mu x Z_C x (X/C)_prey

 !write(*,*) 'Fin de la fonction f_graz_C2X'
 
 End function f_graz_C2X
!--------------------------------------------------------------------
