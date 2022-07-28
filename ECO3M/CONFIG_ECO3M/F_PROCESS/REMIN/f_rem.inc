!--------------------------------------------------------------------
!
            function f_rem(ili,jcol,tx_rem)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule la reminéralisation des détritus (ds et dl en don), de chl
!
!
!Cette fonction calcule aussi la reminéralisation du frustule de Si des
!diatomées en acide silicique
!
! variables et unites:
!---------------------
! tx_rem : taux de reminéralisation maxi en s-1
!
 Implicit None
 
!-- fonction
real(8) ::f_rem(nx_min:nx_max,ny_min:ny_max,nz_max),tfunc(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
Real(8) :: tx_rem
Integer :: ili,jcol


tfunc(:,:,:)=2.95**((TEMP_BIO(:,:,:)-11.0)/10.0)

f_rem(:,:,:)=  tx_rem * tfunc * var(ili)%conc(:,:,:)


end function f_rem
