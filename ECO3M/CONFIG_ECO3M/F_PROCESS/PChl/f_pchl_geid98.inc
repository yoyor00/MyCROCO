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
!
   function f_pchl_geid98(ili,jcol,theta_max,alpha_star,inut1,Vm1,K1,iXorC1, &
                                          inut2,Vm2,K2,iXorC2,inut3,Vm3,K3,iXorC3) 
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la variation de chlorophylle en fonction
! de la production primaire par un modele de (Geider et al., 1998). 
!
! variables et parametres et leurs unites:
!----------------------------------------
! theta_max  = Chl:N maximum, gChl/molN
! alpha_star = pente a l''origine de PP^* vs E, molC.m^2/(gChl.J)
! rho_chl    = efficacite de photosynthese, gChl/molN
! Vm         = vitesse max d'absorption du nutriment azote
!               son unite est déterminée par l'indice iXorC
!               si iXorC correspond à la conc. de la biomasse azotée <--> Vm est en s^{-1}
!               si iXorC correspond à la conc. de la biomasse carbonée <--> Vm est en molN. (mol C)^{-1} s^{-1}
!               si iXorC  vaut 0 Vm doit etre calcule selon la formule PBmax.Qmax et la valeur
!               de Qmax est donnée dans le paramètre Vm
! inut       = indice de la concentration du nutriment azoté absorbé
! K          = constante de demi saturation de l'absorption de nutriment, molX m^{-3}
! ili        = contient l'indice de la concentration en chlorophylle de l'autotrophe
!
! f_pchl_geid98 = Vitesse de production de chlorophylle, gChl/m3/s
!
!-- ordres de grandeur: theta_max= 0.3 * 14 gChl/molN (Geider, 98)
!
! Refs: Geider et al. (1998) Limnol. Oceanog., 43(4):679-694 
!
! Melika Baklouti 27/06/2008
!------------------------------------------------------------------
!-- Variables globales
 Use DEF_TYPE
 Use VAR_GLOBAL
 Use COUPLEUR_PHY_BIO
 Use MOD_FCHAIN

 Implicit None

!-- arguments
 Real(8) :: theta_max,alpha_star
 Real(8) :: inut1,Vm1,K1,iXorC1
 Real(8), optional :: inut2,Vm2,K2,iXorC2,inut3,Vm3,K3,iXorC3
 Real(8) :: f_pchl_geid98(nx_min:nx_max,ny_min:ny_max,nz_max)
 Integer :: ili,jcol

!-- variables locales:
 Character(L_CHAIN) :: chaine,chaine2
 integer:: iorg,iiXorC1,iiXorC2,iiXorC3,i
 integer:: iinut1,iinut2,iinut3
 integer:: iiC,iiChl
 real(8),Allocatable :: rhochl(:,:,:), theta(:,:,:)
 real(8),Allocatable :: vit_upt_N(:,:,:),Vmax(:,:,:)
 real(8),Allocatable :: VN1(:,:,:),VN2(:,:,:),VN3(:,:,:)
 real(8) :: Qmax

 Allocate (rhochl(nx_min:nx_max,ny_min:ny_max,nz_max),theta(nx_min:nx_max,ny_min:ny_max,nz_max))
 Allocate (vit_upt_N(nx_min:nx_max,ny_min:ny_max,nz_max),Vmax(nx_min:nx_max,ny_min:ny_max,nz_max))
 Allocate (VN1(nx_min:nx_max,ny_min:ny_max,nz_max))
 if (present(inut2)) Allocate (VN2(nx_min:nx_max,ny_min:ny_max,nz_max))
 if (present(inut3)) Allocate (VN3(nx_min:nx_max,ny_min:ny_max,nz_max))


!-- Initialisations:
  vit_upt_N = 0.d0
  Vmax = 0.d0

!-- Conversions en entiers:
iinut1 = nint(inut1)
iiXorC1 = nint(iXorC1)

if (present (inut2)) then
  iinut2 = nint(inut2)
  iiXorC2 = nint(iXorC2)
  if (present (inut3)) then 
    iinut3 = nint(inut3)
    iiXorC3 = nint(iXorC3)
  endif
endif

!-- Recherche du numero d'organisme de l'autotrophe: 
 iorg = var(ili)%idorg

!-- de l'indice de sa conc en C :
 chaine = 'N'
 iiC = f_idorg2id (iorg,chaine)

!-- et de celle de sa conc en Chl :
 chaine = 'Chl'
 iiChl = f_idorg2id (iorg,chaine)

!--rapport CHL:C (g Chl/mol C):
 theta = VAR(iiChl)%conc / (VAR(iiC)%conc + 1.d-15) 

!-- Calcul du coefficient ro_chl de l'expression de Geider (1998)
 do i=1,nscp_phy
   if (mu_PPB_NR(i)%idorg == iorg) then
!-- rhoChl (g Chl / mol N):
      rhochl = theta_max * mu_PPB_NR(i)%val(:,:,:) / (alpha_star * E_PARZ(:,:,:) * theta(:,:,:) + 1.d-15) !--1d-15 pour denom non nul
      Vmax(:,:,:)= mu_PPB_NR(i)%val(:,:,:) 
      exit
   endif
 enddo

if (iiXorC1 == 0) then
! calcul de Vm (en  molX/m3/s) selon PBmaxQmax et la valeur de Qmax est transmise par le parametre Vm1
   Qmax = Vm1
   Vmax = Vmax * Qmax * VAR(iiC)%conc ! conversion de la PPB en (mol X) m^{-3} s^{-1}
else
   Vmax = Vm1 * VAR(iiXorC1)%conc  ! conversion de la PPB en (mol X) m^{-3} s^{-1}
endif
   
VN1 = Vmax * VAR(iinut1)%conc / (K1 + VAR(iinut1)%conc)
vit_upt_N =  VN1 

!-- Test sur la presence d'autres nutriments azotes:
if (present (inut2)) then
  if (iiXorC2 == 0) then
   Qmax = Vm2
   Vmax = Vmax * Qmax * VAR(iiC)%conc ! conversion de la PPB en mol X m-3 s-1
  else
   Vmax  = Vm2  * VAR(iiXorC2)%conc  ! conversion de Vmax en mol X m-3 s-1
  endif

  VN2 = Vmax * VAR(iinut2)%conc / (K2 + VAR(iinut2)%conc)
  vit_upt_N = vit_upt_N + VN2

  if (present (inut3)) then 
    if (iiXorC3 == 0) then
      Qmax = Vm3
      Vmax = Vmax * Qmax * VAR(iiC)%conc ! conversion de la PPB en mol X m-3 s-1
    else
      Vmax  = Vm3  * VAR(iiXorC3)%conc  ! conversion de Vmax en mol X m-3 s-1
    endif

    VN3 = Vmax * VAR(iinut3)%conc / (K3 + VAR(iinut3)%conc)
    vit_upt_N = vit_upt_N + VN3
  endif
endif

!-- Vitesse de production de chlorophylle (en gChl/m3/s):
 f_pchl_geid98 =  rhochl *  vit_upt_N

Deallocate (rhochl,theta,vit_upt_N,Vmax,VN1)
if (present(inut2)) DeAllocate (VN2)
if (present(inut3)) DeAllocate (VN3)

!-- Affichage eventuel à l'écran (pour tests):
! write(*,*) 'tps,f_pchl_geid98 =',tps/3600,f_pchl_geid98

 end function f_pchl_geid98
