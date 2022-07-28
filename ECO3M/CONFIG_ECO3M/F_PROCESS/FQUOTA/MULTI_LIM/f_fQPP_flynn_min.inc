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
!--------------------------------------------------------------------------
!
   function f_fQPP_flynn_min(ili,jcol,indX1,Qmin1,Qmax1,K_Q1,indX2,Qmin2,Qmax2,K_Q2,&
                                          indX3,Qmin3,Qmax3,K_Q3,indX4,Qmin4,Qmax4,K_Q4)
!
!--------------------------------------------------------------------------
! Cette fonction calcule le facteur adimensionnel (fonction
! de quotas) qui limite la CROISSANCE d'un organisme du fait de 2 (ou +) elements
! chimiques X1 et X2 dont les conc. internes sont données par les variables
! d''indice indX1 et indX2
!
! La limitation par plusieurs nutriments repose ici sur la fonction de quota
! de Flynn et al.(2001) et sur la loi du minimum:
!
!   f_fQPP_flynn_min = min_i ((Q_i - Qmin)/(Q_i - Qmin + K_Q) * (Qmax - Qmin + K_Qi) / (Qmax - Qmin) )
!
! ATTENTION:  On laisse cette fonction delibérément devenir négative pour 
!             traduire les situations où l'organisme évacue plus de carbone qu'il n'en  
!             n'absorbe. Par ailleurs, cette fonction ne peut excéder la valeur 1
!            (cf Baklouti et al., 2006 pour plus de détails sur les modèles d'exudation)
!
! Variables et parametres et leurs unites:
!-----------------------------------------
! indX: indice de la conc. interne en nutriment X dont la disponibilité est limitante
! Qmin: valeur minimale du quota intracellulaire de X/C (mol X (mol C)^-1)
! Qmax: valeur maximale du quota intracellulaire de X/C (mol X (mol C)^-1) 
! K_Q : parametre de courbure de la fonction, dimension d'un quota 
!       intracellulaire (mol X (mol C)^-1)
!
! Refs: Flynn et al. (2001) Limnol. & Oceanogr. 46(7):1794-1802
!       Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
!
! Melika Baklouti 27/06/2008
!-------------------------------------------------------------------------
 USE COUPLEUR_PHY_BIO
 USE VAR_GLOBAL
 USE MOD_FCHAIN

 Implicit None

!-- fonction:
 Real(8) :: f_fQPP_flynn_min(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
 Real(8) :: Qmin1,Qmax1,indX1,K_Q1,Qmin2,Qmax2,indX2,K_Q2
 Real(8),optional :: Qmin3,Qmax3,indX3,K_Q3,Qmin4,Qmax4,indX4,K_Q4
 Integer :: ili,jcol

!--variables locales:
 Integer :: iiC,iorg1,iorg2,iorg3,iorg4
 Integer :: nbnut
 Integer :: i,j,k
 Integer :: iindX1,iindX2,iindx3,iindx4
 Real(8) :: Q1,fQ1,Q2,fQ2
 Real(8) :: Q3,fQ3,Q4,fQ4
 Character(L_VAR)::chain

 iindX1 = nint(indX1)
 iindX2 = nint(indX2)
 iorg1=var(iindX1)%idorg
 iorg2=var(iindX2)%idorg

 nbnut = 2
 if (present(indx3)) then
   iindx3 = nint(indx3)
   iorg3=var(iindX3)%idorg
   nbnut = 3
   if (present(indx4)) then
     iindx4 = nint(indx4)
     iorg4=var(iindX4)%idorg
     nbnut = 4
   endif
 endif
   

!-- Test sur les organismes:
 Select case (nbnut)
  case(2)
   if (iorg1 /= iorg2 ) then
    write(*,*) 'pb avec f_fQPP_flynn_min.inc pour le flux var(',ili,') vers var(',jcol,')'
    STOP
  endif
  case(3)
    if (iorg1 /= iorg2 .or. iorg1 /= iorg3 ) then
     write(*,*) 'pb avec f_fQPP_flynn_min.inc pour le flux var(',ili,') vers var(',jcol,')'
     STOP
   endif
  case(4)
   if (iorg1 /= iorg2 .or. iorg1 /= iorg3 .or. iorg1 /= iorg4) then
     write(*,*) 'pb avec f_fQPP_flynn_min.inc pour le flux var(',ili,') vers var(',jcol,')'
     STOP
   endif
 end select
 
 chain = 'C'
 iiC = f_idorg2id(iorg1,chain)
 

!-- Calcul de la fonction de limitation :
do i=nx_min,nx_max
  do j=ny_min,ny_max
    do k=1,nz_max
        Q1 = var(iindX1)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
        Q2 = var(iindX2)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
        if (nbnut > 2) then
           Q3 = var(iindX3)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
           if (nbnut > 3) then
             Q4 = var(iindX4)%conc(i,j,k) / (var(iiC)%conc(i,j,k) + 1.d-15)
           endif
        endif
           
        fQ1 = (Q1 - Qmin1) / (Q1 - Qmin1 + K_Q1) * (Qmax1 - Qmin1 + K_Q1) / (Qmax1 - Qmin1)
        fQ2 = (Q2 - Qmin2) / (Q2 - Qmin2 + K_Q2) * (Qmax2 - Qmin2 + K_Q2) / (Qmax2 - Qmin2)

        if (nbnut > 2) then
          fQ3 = (Q3 - Qmin3) / (Q3 - Qmin3 + K_Q3) * (Qmax3 - Qmin3 + K_Q3) / (Qmax3 - Qmin3)
          if (nbnut > 3) then
            fQ4 = (Q4 - Qmin4) / (Q4 - Qmin4 + K_Q4) * (Qmax4 - Qmin4 + K_Q4) / (Qmax4 - Qmin4)
          endif
        endif

!-- calcul du quota minimum:
           Select case (nbnut)
           case(2)
             f_fQPP_flynn_min(i,j,k) =  min(1.d0, fQ1, fQ2) 
           case(3)
             f_fQPP_flynn_min(i,j,k) =  min(1.d0, fQ1, fQ2 , fQ3)
           case(4)
             f_fQPP_flynn_min(i,j,k) =  min(1.d0, fQ1, fQ2, fQ3, fQ4)
           end select
    enddo
  enddo
enddo


end function f_fQPP_flynn_min
!--------------------------------------------------------------------
