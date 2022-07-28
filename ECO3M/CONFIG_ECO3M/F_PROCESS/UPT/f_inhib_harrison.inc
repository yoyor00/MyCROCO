!--------------------------------------------------------------------
!Dernière modification: 11/10/2007
!
     function f_inhib_harrison(ili,jcol,ki,im)
!
!--------------------------------------------------------------------
!
! Cette fonction calcule le term d'inhibition par ammonium au cours de uptake de nitrate 
!Harrison et al. 1996. Limnology Oceanogr....
!
! variables et unites:
!---------------------
! ki : constante d'inhibition de uptake nitrate par ammonium,
! im : paramètre d'inhibition 

 Implicit None
 
!-- fonction
real(8) :: f_inhib_harrison(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- arguments
Real(8) :: ki,im
Integer :: iorg,ili,jcol,id

!-- variables locales:
 character(L_CHAIN)::chain1,chain2
   
 chain1 =  'NH4'
 chain2 =  'N'
  
id=f_scomp2id(chain1,chain2)

f_inhib_harrison(:,:,:) =  1.-im * var(id)%conc(:,:,:)/(var(id)%conc(:,:,:) + ki)
       

end function f_inhib_harrison
