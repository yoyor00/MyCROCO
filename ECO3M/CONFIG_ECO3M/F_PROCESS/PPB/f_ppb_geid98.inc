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
 function f_ppb_geid98(ili,jcol,Pcref,Tref,AE,alpha,Qmin,Qmax)
!
!--------------------------------------------------------------------
!
! This function calculates de photosynthesis rate of phytoplankton VAR(ili) according
! to Geider et al. (1998) formulation. Growth is either limited by nitrogen
! availability (when Qmin and Qmax are provided) or proceeds in nutrient replete
! conditions
!
! 
! Variables and parameters and their units:
!----------------------------------------- 
! ili   : index of the phytoplankton carbon in vector VAR
! Pcref : C-specific photosynthesis rate at temperature Tref 
! Tref  : Reference temperature (in K or C according to TEMP_BIO)
! AE    : Slope of the linear region of the Arrhenius plot (-)
! Qmin  : minimum N/C ratio in phytoplankton  (mol N . mol C^-1)
! Qmax  : maximum N/C ratio in phytoplankton  (mol N . mol C^-1)
! alpha : nitial slope of the  P vs E curve (mol C. g Chl^-1 . µmol photons^-1 . m^2)
!
! Refs: Geider et al. (1998) Limnol. Oceanogr. 43(4):679-694
!
! Melika Baklouti 7/07/2008
!--------------------------------------------------------------------
!
! Global variables 
Use  DEF_TYPE
Use COUPLEUR_PHY_BIO
USE  VAR_GLOBAL

 Implicit None
!-- function
Real(8) :: f_ppb_geid98(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- Arguments
integer    :: ili,jcol
real(8)    :: Pcref,Tref,AE,alpha
real(8), optional :: Qmin,Qmax

!-- Local variables :
integer    :: i,j,k,iiC
real(8)    :: Pcmax    ! Maximum C-specific photosynthesis rate  (s^-1)
real(8)    :: Chl      ! Phytoplankton chlorophyll (g.m^-3);
real(8)    :: Q        ! nitrogen quoat in phytoplankton (molN molC^-1)
real(8)    :: theta    ! Chla/C  ratio (g.mol^-1) 
integer    :: iorg,iichl,iiN
Character(L_VAR)::elmt

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

!-- index of phytoplankton organism:
 iorg = VAR(ili)%idorg

!-- case of a constant CHL:C ratio:
if (.not. CHL_C_BOOL) then
        Chl = CHL_C0 * var(ili)%conc(i,j,k)
        theta = CHL_C0 !-- Chl:C ratio
elseif (CHL_C_BOOL  .and.  .not. CHL_C_FUNC ) then
!-- Search for phytoplankton chlorophyll index
   elmt ='Chl'
   iichl = f_idorg2id (iorg,elmt)
   Chl = var(iichl)%conc(i,j,k)
   theta = Chl/(var(ili)%conc(i,j,k)+1.d-15)
else
   elmt ='C'
   iiC = f_idorg2id (iorg,elmt)
   Chl = CHL_C(1)%val * var(iiC)%conc
endif


if (present (Qmin)) then
!-- Search for phytoplankton nitrogen index :
   elmt ='N'
   iiN = f_idorg2id (iorg,elmt)

!-- nitrogen ratio in phytoplankton:
Q = var(iiN)%conc(i,j,k)/(var(ili)%conc(i,j,k)+1.d-15)
endif
 
if (present (Qmin)) then
!-- Maximum C_specific photosynthesis rate :
 Pcmax  = Pcref * exp(AE * (1./TEMP_BIO(i,j,k) - 1./Tref))*(Q-Qmin)/(Qmax-Qmin)

!--  C_specific photosynthesis rate (Pcphot in Geider et al.,98)
 mu_PPB(iorg)%val(i,j,k) = Pcmax *(1.-exp(-alpha*theta*E_PARZ(i,j,k)/Pcmax))

!-- photosynthesis rate (molC m-3 s-1)
 f_ppb_geid98(i,j,k) = mu_PPB(iorg)%val(i,j,k)*VAR(ili)%conc(i,j,k)

else !-- nutrient replete rate

 Pcmax  = Pcref * exp(AE * (1./TEMP_BIO(i,j,k) - 1./Tref))

!--  C_specific photosynthesis rate (Pcphot in Geider et al.,98)
 mu_PPB_NR(iorg)%val(i,j,k) = Pcmax *(1.-exp(-alpha*theta*E_PARZ(i,j,k)/Pcmax))

!-- photosynthesis rate (molC m-3 s-1)
 f_ppb_geid98(i,j,k) = mu_PPB_NR(iorg)%val(i,j,k)*VAR(ili)%conc(i,j,k)
endif

    enddo
  enddo
enddo

End function f_ppb_geid98
