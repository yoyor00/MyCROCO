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
        function f_uptB_nut_vichi07(ili,jcol,inut,upt_max,Qbact,kinut)
!
!--------------------------------------------------------------------
!
! This function calculates the absorption/exudation of ammonium and phosphate by bacteria
! according to Vichi et al. (2007) model. Note that the bacterial concentration is 
! necessarily in variable VAR(ili)%conc.
! 
! variables and units:
!---------------------
! k       : half-saturation constant, in molX.m-3
! Qbact   : Optimal X/C quota of bacteria, in molX/molC
! upt_max : Max. specific rate of uptake/exudation, in s-1,
! inut    : index of the concentration of the limiting nutrient

! Refs: Vichi M., N. Pinardi, S. Masina (2007) Journal of Marine Systems, 64: 89-109
!
! Frederic Diaz 15/07/2008
!----------------------------------------------------------------------
 Implicit None

!-- function:
real(8)::f_uptB_nut_vichi07(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
Real(8) :: upt_max,Qbact,kinut,inut
Integer :: ili,jcol

!-- local variables :
Integer :: i,j,k,iinut,iC,iX,iorg_plancton
character(L_CHAIN)    :: chain
Character(L_VAR)      :: elmt_nut
real(8)               :: Q,func

!-- index of nutrient concentration:
iinut=nint(inut)

!-- Plankton organisme :
 iorg_plancton = var(ili)%idorg 

!-- element X given by the nutrient : 
 elmt_nut = trim(adjustl(VAR(iinut)%elmt))  

!-- plankton concentration in X 
 iX= f_idorg2id(iorg_plancton,elmt_nut) 

!-- plankton concentration in C
 chain='C'
 iC= f_idorg2id(iorg_plancton,chain) 



do i=nx_min,nx_max
  do j=nx_min,ny_max
    do k=1,nz_max
       Q=var(iX)%conc(i,j,k)/(var(iC)%conc(i,j,k)+1.d-15)
 
       if (Q > Qbact) then
            func=-1.d0
       elseif (Q < Qbact) then
            func= -1.*var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + kinut)
       else 
            func= 0.d0
       endif 

       f_uptB_nut_vichi07(i,j,k) = func * upt_max * (Q-Qbact) * var(iC)%conc(i,j,k)
  
    enddo 
  enddo
enddo   
    


!write(*,*)'f_uptB_nut_vichi07=',f_uptB_nut_vichi07(:,:,:)

end function f_uptB_nut_vichi07
