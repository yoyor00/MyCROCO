!--------------------------------------------------------------------
! Derniere modification: 01/07/2008
!
 function f_respB_vichi07(ili,jcol,rbact,bge,QTEN)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la respiration bacterienne a partir du
! taux de croissance nette bacterienne
!
! Vichi M., N. Pinardi, S. Masina 2007. Journal of Marine Systems, 64: 89-109
! variables et unites:
!---------------------
! rbact: taux de respiration basal bacteries (s-1)
! bge: efficacite de croissance bacteries (sans unite)
!
 Implicit None

!-- arguments
Real(8) :: rbact,bge,QTEN
Real(8) :: f_respB_vichi07(nx_min:nx_max,ny_min:ny_max,nz_max)
Integer :: ili,jcol,ik

!-- variables locales:
 Integer:: iorg_bact
 
!-- organisme:
 iorg_bact = var(ili)%idorg
 
do ik=1,nscp_bac
 
	if (PBB_NR(ik)%idorg == iorg_bact) then
	
!-- vitesse de respiration bacterienne (en mmolC.m-3.s-1):

		  f_respB_vichi07 =  rbact *QTEN**((TEMP_BIO(:,:,:)-10.)/10.)*var(ili)%conc + (1.-bge)*PBB_NR(ik)%val(:,:,:)
	
	exit

	endif
 
enddo

end function f_respB_vichi07
