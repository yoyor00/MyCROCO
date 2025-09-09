!--------------------------------------------------------------------
!
function f_nitrif_linear(ili,jcol,iNH4,tx_nit)
    !
    !--------------------------------------------------------------------
    !
    ! Cette fonction calcule la nitrification (regeneration du nitrate)
    ! conformement aux observations de Ward (2000).
    !
    ! variables et unites:
    !---------------------
    ! tx_nit : taux max de nitrification , en s-1
    ! iNH4: indice de la concentration en ammonium
    !
    ! Refs: !       Ward (2000) In: Kirchman DL., Ed.,
    ! Microbial Ecology of the Oceans, Wiley-Liss, New York, pp 427-453
    !
    ! M. Baklouti 07/05/2012
    !---------------------------------------------------------------------
    ! Variables globales

    Implicit None

    !-- fonction
    real(8):: f_nitrif_linear(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: tx_nit,iNH4
    Integer :: ili,jcol,iiNH4
    !-- local variables
    Integer :: i,j,k

    iiNH4=nint(iNH4)

    f_nitrif_linear = 0.d0

     do k = 1,nz_max
        do j = ny_min,ny_max
            do i = nx_min,nx_max
               if (VAR(iiNH4)%conc(i,j,k) > tx_nit * VAR(iiNH4)%conc(i,j,k) * dt_bio) then
	          f_nitrif_linear(i,j,k) = tx_nit * VAR(iiNH4)%conc(i,j,k)
               endif
            enddo
       enddo
     enddo

end function f_nitrif_linear

