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
      function f_ChlC_cloern_95(inut,KX)
!      
!--------------------------------------------------------------------
!
! This function calculates the Chl:C ratio according to Cloern et al.(1995) 
! formulation
!
! Variables and parametres and their units:
!------------------------------------------ 
!
! inut:               index of the most limiting nutrient X
! KX  :               half-saturation constant for X uptake
! f_ChlC_cloern_95:   Chl:C ratio in gChl (molC)^{-1}
!
! Refs: Cloern et al. (1995) Limnol. Oceanogr., 40(7):1313-1321
!
! Christel Pinazo et Vincent Faure 15/07/2008
!--------------------------------------------------------------------
implicit none

!-- function:
 Real(8)  :: f_ChlC_cloern_95(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- Arguments
real(8)  :: inut,KX

!-- Local variables:
integer    :: iinut,i,j,k

iinut= nint(inut)

!-- Chl:C in (g Chl) (gC)^{-1}:

do i= nx_min,nx_max
  do j= ny_min,ny_max
    do k = 1, nz_max
     f_ChlC_cloern_95(i,j,k) = 0.003 + 0.0154 * exp(0.050*TEMP_BIO(i,j,k))*exp(-0.059*E_PARZ(i,j,k)*0.4) &
                  * var(iinut)%conc(i,j,k)/(KX+var(iinut)%conc(i,j,k))
    enddo
  enddo
enddo

!-- conversion in gChl (molC)^{-1}:
    f_ChlC_cloern_95 = f_ChlC_cloern_95 * 12.0

!    do i=1,nscp_phy
!      CHL_C(i)%val(:,:,:) = f_ChlC_cloern_95(:,:,:)
      CHL_C(1)%val(:,:,:) = f_ChlC_cloern_95(:,:,:)
!    enddo
!write(*,*) 'fin de f_ChlC_cloern_95'

End function f_ChlC_cloern_95

