!--------------------------------------------------------------------
!
         function f_nitrif_vichi07diaz(ili,jcol,tx_nit,QTEN)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la nitrification (regeneration du nitrate)
! adaptee de Vichi et al. 2007; ajout d une parametrisation sur la lumiere 
! conformement aux observations de Ward (2000).
!
! variables et unites:
!---------------------
! tx_nit : taux max de nitrification , en s-1
! QTEN : coefficient empirique pour la nitrification, en °C-1
! ili: indice de la concentration en nitrate
! jcol: indice de la concentration en ammonium
! 
! Refs: Vichi et al. (2007) JMS 64: 89-109
!       Ward (2000) In: Kirchman DL., Ed.,
! Microbial Ecology of the Oceans, Wiley-Liss, New York, pp 427-453
!
! Frederic Diaz 10/07/2008
!---------------------------------------------------------------------

 Implicit None

!-- fonction
real(8):: f_nitrif_vichi07diaz(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
Real(8) :: tx_nit,QTEN
Integer :: ili,jcol

! 3/12/2018 MB: corrrection bug de la division possible par 0
f_nitrif_vichi07diaz(:,:,:)=  tx_nit * var(jcol)%conc(:,:,:)* QTEN**((TEMP_BIO(:,:,:)-10.)/10.)*&
           (1.-E_PARZ(:,:,:)/(E_PAR(nx_max,ny_max)+1d-30))
	   
end function f_nitrif_vichi07diaz
