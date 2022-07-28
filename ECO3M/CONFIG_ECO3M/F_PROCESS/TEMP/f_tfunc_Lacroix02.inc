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
!---------------------------------------------------------------------------
!
       function f_tfunc_Lacroix02(ili,jcol,beta,t_opt,t_let)
!    
!---------------------------------------------------------------------------
! Cette fonction calcule l''effet de la température sur les processus physiologiques
! selon une courbe en cloche. Fonction variant entre zero et l''unite. 
!
!
! Variables et parametres et leurs unites:
!-----------------------------------------
!beta: facteur de forme agissant sur la largeur de la courbe de temperature autour
!      de l''optimum
!t_opt: temperature optimale de developpement du phytoplancton, (en C)
!t_let: temperature minimale de croissance, (en C)
!
!Ref:  Lacroix et Gregoire (2002) JMS, 37:229-258
!
!Frederic Diaz 13/09/2007
!-------------------------------------------------------------------------

Implicit None  

!-- fonction
Real(8) :: f_tfunc_Lacroix02(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
Real(8) :: beta,t_opt,t_let
Integer :: ili,jcol

!-- Local variables:
Real(8),Allocatable:: tfunc(:,:,:)

Allocate(tfunc(nx_min:nx_max,ny_min:ny_max,nz_max))

tfunc = (2.*(1.-beta)*(TEMP_BIO(:,:,:)-t_let)/(t_opt-t_let))/(((TEMP_BIO(:,:,:)-t_let)/(t_opt-t_let))**2.&
   +2.*(-beta)*(TEMP_BIO(:,:,:)-t_let)/(t_opt-t_let)+1.)

f_tfunc_Lacroix02(:,:,:)=MAX(tfunc,0.d0)

Deallocate (tfunc)

end function f_tfunc_Lacroix02
!--------------------------------------------------------------------
