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
   function f_graz_C_holling2(ili,jcol,ipred,iprey,F,Imax)
!
!--------------------------------------------------------------------
! This function calculates the grazing rate according to Holling formulation.
!
! Attention:
! ----------
! 
! * This function applies for un FLUX  implying one predator and one prey
!   the concentration of which are respectively given by VAR(ipred)%conc
!   and VAR(iprey)%conc
!
! * The calculation is based on Carbon concentrations of predator and preys,
!   and provides  a grazing rate in  [mol C] m-3 s^{-1}
!
! Variables and  parameters and their units:
!------------------------------------------- 
! Imax               : Maximum ingestion rate of the prey by predator, (in [pred]^{-1} [prey]^{-1} s^{-1})
! F                  : Filtration rate (in [pred]^{-1} s^{-1})
! f_graz_C_holling2  : Grazing rate of the prey iprey by predator ipred (in mol C m^{-3} s^{-1})
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

!-- function:
 Real(8)             :: f_graz_C_holling2(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- Arguments
 Real(8)             :: ipred,iprey,F,Imax
 Integer             :: ili,jcol
 
!-- Local variables :
 Integer              :: iiprey,iipred,iscp_pred,iscp_prey
 Integer              :: iC_pred,iC_prey,iorg_prey,iorg_pred
 Integer              :: i,j,k,ii,jj,kk,ir,mm,istat
 Character(L_VAR)     :: elmt
 Character(L_CHAIN)   :: chaine

!-- conversion in integer:
 iipred = nint(ipred)
 iiprey = nint(iprey)


!-- initialization:
Imax = 0.d0
f_graz_C_holling2=0.d0

!
!-- Organism of the predator and the prey :
iorg_prey = VAR(iiprey)%idorg
iorg_pred = VAR(iipred)%idorg

!-- Search of the carbon concentrations of predator and prey :
 chaine = 'C'
 iC_prey = f_idorg2id(iorg_prey,chaine)
 iC_pred = f_idorg2id(iorg_pred,chaine)

!-- Search of the position of the predator and prey in table grazphi :
!-- initialisations:
 iscp_pred = 0
 iscp_prey = 0

!-- Search of predator and prey inside the zooplankton compartment:
 do mm = 1,nscp_zoo
   if (iscp_zoo(mm)==VAR(iipred)%idorg) then
     iscp_pred = mm
   elseif (iscp_zoo(mm)==VAR(iiprey)%idorg) then
     iscp_prey = mm
  endif
 enddo

!-- If the prey is not in the zooplankton compartment, search inside the phytoplankton
!    compartment :

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
!-- If the prey is not in the zooplankton compartment, search inside the bacterial
!    compartment :
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


!-- Specific grazing rate of prey  iscp_prey by predator iscp_pred: 

 mu_graz(iscp_pred,iscp_prey)%val(:,:,:) =  Imax /(1 + Imax/(F * var(iC_prey)%conc(:,:,:))) 


!-- Grazing rate  (G x Zc ) in mol C m^-3 s-1

f_graz_C_holling2(:,:,:) =   mu_graz(iscp_pred,iscp_prey)%val(:,:,:)   * var(iC_pred)%conc(:,:,:)                  ! mu x Z_C


! write(*,*) 'End of function f_graz_C_holling2'
 
 End function f_graz_C_holling2
!--------------------------------------------------------------------
