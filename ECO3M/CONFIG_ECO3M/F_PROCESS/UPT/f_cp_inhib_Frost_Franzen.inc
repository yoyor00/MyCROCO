!--------------------------------------------------------------------
!Derniere modification: 12/07/2011
!
function f_cp_inhib_Frost_Franzen(ili,jcol,inut,kn)
    !
    !--------------------------------------------------------------------
    !
    ! Cette fonction compense le terme d'inhibition par ammonium au cours de uptake de nitrate 
    ! Frost and Franzen 1992. MEPS 83: 291-303
    !
    ! variables et unites:
    !---------------------
    ! kn : constante d'inhibition de l'uptake de nitrate ou DON par ammonium, ou de DOP par PO4, molX m-3
    !
    ! Romain Mauriac
    !--------------------------------------------------------------------------------------------

    ! Variables globales

    Implicit None

    !-- fonction
    real(8) :: f_cp_inhib_Frost_Franzen(nx_min:nx_max,ny_min:ny_max,1:nz_max)

    !-- arguments
    Real(8) :: kn,inut
    Integer :: ili,jcol

    !-- variables locales:
    integer :: id,i,j,k

    id=nint(inut)
     
    f_cp_inhib_Frost_Franzen = 0.d0
    do k = 1,nz_max
        do j = ny_min,ny_max
            do i = nx_min,nx_max
            if (var(id)%conc(i,j,k) > 1d-10) then
              f_cp_inhib_Frost_Franzen(i,j,k) = 1.d0/(1.d0 + kn/(var(id)%conc(i,j,k) + 1.d-80))
!              f_cp_inhib_Frost_Franzen(i,j,k) = 0.d0
            endif
            enddo
        enddo
    enddo

    !do i = nx_min,nx_max
    !  do j = ny_min,ny_max
    !    do k = 1,nz_max
    !IF (ISNAN(f_cp_inhib_Frost_Franzen(i,j,k))) THEN
    !WRITE(*,*) 'FUNCTION f_cp_inhib_Frost_Franzen contains NaN',i,j,k
    !STOP
    !ENDIF
    !enddo
    !enddo
    !enddo

end function f_cp_inhib_Frost_Franzen
