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
!------------------------------------------------------------------------------
!
!
       function f_fQPP_droop_min2(ili,jcol,indX1,Qmin1,Qmax1,indX2,Qmin2,Qmax2,&
                                                  indX3,Qmin3,Qmax3,indx4,Qmin4,Qmax4)
!
!------------------------------------------------------------------------------
! Cette fonction calcule le facteur adimensionnel (fonction
! de quotas) qui limite la CROISSANCE d'un organisme du fait de 2 (ou +) elements
! chimiques X1 et X2 dont les conc. internes sont données par les variables
! d''indice indX1 et indX2
!
! La limitation par plusieurs nutriments repose ici sur la fonction de quota
! de Droop (1968) et sur la loi du minimum:
!               f_fQPP_droop_min_2 = 1- Qmin/Q_i where Q_i/Qmax_i = min_j(Q_j/Qmax_j)
!
! ATTENTION:  On laisse cette fonction délibérément devenir négative pour 
!             traduire les situations où l'organisme évacue plus de carbone qu'il n'en  
!             n'absorbe. Par ailleurs, cette fonction ne peut excéder la valeur 1
!            (cf Baklouti et al., 2006 pour plus de détails sur les modèles d'exudation)
!
! Variables et parametres et leurs unites:
!-----------------------------------------
! indX: indice de la conc. interne en nutriment X dont la disponibilité est limitante
! Qmin: valeur minimale du quota intracellulaire de X/C (mol X (mol C)^-1)
!
! Refs: Droop, M. R.        (1968) JMBA of the United kingdom  48:689-733
!       Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
!
! Frederic Diaz, Melika Baklouti 11/07/2008
!-------------------------------------------------------------------------
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL
 USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8):: f_fQPP_droop_min2(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: Qmin1, Qmax1,indX1,indX2,Qmin2,Qmax2
 Real(8) , optional :: Qmin3,Qmax3, indX3,indX4,Qmin4,Qmax4
 Integer :: ili,jcol

!--variables locales:
 Integer :: iiC,iorg1,iorg2,iorg3,iorg4
 Integer :: iindX1,iindX2,iindx3,iindx4
 Integer :: nbnut, istat
 Integer :: i,j,k
 Real(8), parameter :: epsilon = 1.d-10
 Real(8)  :: Q1,fQ1,Q2,fQ2
 Real(8)  :: Q3,fQ3,Q4,fQ4,fQ
 Character(L_VAR)::chain

 
iindX1 = nint(indX1)
iindX2 = nint(indX2)
iorg1 = var(iindX1)%idorg
iorg2 = var(iindX2)%idorg

!-- initialization of the nutrient numbers
   nbnut = 2

if (present(indx3))then
   iindx3 = nint(indX3)
   iorg3 = var(iindX3)%idorg
   nbnut = 3
   if (present(indx4))then
     iindx4 = nint(indX4)
     iorg4 = var(iindX4)%idorg
     nbnut = 4
   endif
endif

!-- test on organisms :
 Select case (nbnut)
  case(2)
   if (iorg1 /= iorg2 ) then
    write(*,*) 'pb avec f_fQPP_droop_min.inc pour le flux var(',ili,') vers var(',jcol,')'
    STOP
  endif
  case(3)
    if (iorg1 /= iorg2 .or. iorg1 /= iorg3 ) then
     write(*,*) 'pb avec f_fQPP_droop_min.inc pour le flux var(',ili,') vers var(',jcol,')'
     STOP
   endif
  case(4)
   if (iorg1 /= iorg2 .or. iorg1 /= iorg3 .or. iorg1 /= iorg4) then
     write(*,*) 'pb avec f_fQPP_droop_min.inc pour le flux var(',ili,') vers var(',jcol,')'
     STOP
   endif
 end select

!-- Carbon concentration:
 chain = 'C'
 iiC = f_idorg2id(iorg1,chain)

!-- Calcul de la fonction de limitation :
do i=nx_min,nx_max
  do j=ny_min,ny_max
    do k=1,nz_max
           Q1 = var(iindX1)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
           Q2 = var(iindX2)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
	   fQ =  min(Q1/Qmax1,Q2/Qmax2) 
	   
       if (nbnut > 2) then
           Q3 = var(iindX3)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
	   fQ =  min(Q1/Qmax1,Q2/Qmax2,Q3/Qmax3) 
	   
       if (nbnut > 3) then
           Q4 = var(iindX4)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
	   fQ =  min(Q1/Qmax1,Q2/Qmax2,Q3/Qmax3,Q4/Qmax4)
         endif
       endif

     Select case (nbnut)
     
       case(2)
          if (abs(fQ -  Q1/Qmax1) < epsilon) then 
	    fQ1 = 1.- Qmin1/(Q1+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ1) 
	    
	  else
	    fQ2 = 1.- Qmin2/(Q2+1.d-20)
	    f_fQPP_droop_min2 =  min(1.d0, fQ2)
	  endif
	  
       case(3)
          if (abs(fQ - Q1/Qmax1) < epsilon) then 
	    fQ1 = 1.- Qmin1/(Q1+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ1) 
	    
	  elseif (abs(fQ- Q2/Qmax2) < epsilon) then
	    fQ2 = 1.- Qmin2/(Q2+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ2)
	    
	  else
	    fQ3 = 1.- Qmin3/(Q3+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ3)
	  endif 
	  
       case(4)
          if (abs(fQ- Q1/Qmax1) < epsilon) then 
	    fQ1 = 1.- Qmin1/(Q1+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ1) 
	    
	  elseif (abs(fQ - Q2/Qmax2) < epsilon) then
	    fQ2 = 1.- Qmin2/(Q2+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ2)
	    
	  elseif (abs(fQ- Q3/Qmax3) < epsilon) then
	    fQ3 = 1.- Qmin3/(Q3+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ3)
	  
	  else
	    fQ4 = 1.- Qmin4/(Q4+1.d-20)
	    f_fQPP_droop_min2(i,j,k) =  min(1.d0, fQ4)
	  endif
     end select	   
    enddo
  enddo
enddo


! write(*,*) 'End of function  f_fQPP_droop_min2'

end function f_fQPP_droop_min2
!--------------------------------------------------------------------
