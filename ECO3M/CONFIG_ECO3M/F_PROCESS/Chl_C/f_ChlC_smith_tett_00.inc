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
!
      function f_ChlC_smith_tett_00(Qmin,Qmax,ChlNmin,ChlNmax)
!      
!--------------------------------------------------------------------
!
! Cette fonction calcule le rapport Chl:C selon la formulation
! de  Smith & tett (2000) 
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
!
! Qmin : ratio N/C minimum pour le phytoplancton considéré (mol N . mol C^-1)
! Qmax : ratio N/C maximum pour le phytoplancton considéré (mol N . mol C^-1)
! ChlNmin : ratio Chl/N minimum pour le phytoplancton considéré (g Chl . mol N^-1)
! ChlNmax : ratio Chl/N  maximum pour le phytoplancton considéré (g Chl . mol N^-1)
! 
! Refs: Smith and Tett (2000) Journal of Marine Systems 26:1-36
!
! Christel Pinazo et Vincent Faure 19/02/2008
!--------------------------------------------------------------------
Use VAR_GLOBAL
Use DEF_TYPE
Use COUPLEUR_PHY_BIO
Use MOD_FCHAIN
implicit none

!-- fonction:
 Real(8)  :: f_ChlC_smith_tett_00(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- Arguments
real(8)  :: Qmin,Qmax,ChlNmin,ChlNmax

!-- Variables locales

integer    :: i,j,k,ii
real(8)    :: theta    ! Rapport Chla/C (g.mol^-1) 
real(8)    :: Q        ! Rapport N/C (mol N .mol C^-1) 
real(8)    :: N        ! Concentration en azote phytoplanctonique (mol.m^-3) 
real(8)    :: C        ! Concentration en carbone phytoplanctonique (mol.m^-3) 
real(8), allocatable :: sum(:,:,:)
integer    :: iorg,iin,iic
Character(L_VAR)::elmt

Allocate (sum(nx_min:nx_max,ny_min:ny_max,nz_max))

do ii = 1, nscp_phy
   iorg = iscp_phy(ii)
!-- Indice correspondant a la conc. en carbone du phyto:
   elmt ='C'
   iic = f_idorg2id (iorg,elmt)

!-- Indice correspondant a la conc. en azote du phyto:
   elmt ='N'
   iin = f_idorg2id (iorg,elmt)

do i=nx_min,nx_max
  do j=ny_min, ny_max
    do k=1,nz_max

   C = var(iic)%conc(i,j,k)
   N = var(iin)%conc(i,j,k)

   Q = N / (C + 1.d-20)
   f_ChlC_smith_tett_00(i,j,k) = Q*(ChlNmin+(ChlNmax-ChlNmin)*(Q-Qmin)/(Qmax-Qmin))
   sum = sum + f_ChlC_smith_tett_00(i,j,k)
    end do
  end do
 end do
!-- On prend une valeur moyenne sur l'ensembles des autotrophes:
  CHL_C(1)%val(:,:,:) = sum(:,:,:) / nscp_phy
enddo

End function f_ChlC_smith_tett_00

