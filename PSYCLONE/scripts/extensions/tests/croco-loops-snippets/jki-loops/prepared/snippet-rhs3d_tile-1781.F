subroutine snippet()
  integer*4 :: itrc
  integer*4 :: nt
  integer*4 :: i
  integer*4 :: j
  integer*4 :: k
  integer*4 :: nadv
  integer*4 :: iend
  integer*4 :: lm
  integer*4 :: n
  integer*4 :: jstr
  integer*4 :: jend
  integer*4 :: istr
  integer*4 :: mm
  integer*4 :: iic
  integer*4 :: ntstart
  integer*4 :: cff
  integer*4 :: dt
  integer*4 :: nnew
  integer*4 :: nstp
  integer*4 :: gamma
  integer*4 :: cff1
  integer*4 :: cff2
  integer*4 :: indx
  integer*4 :: jstrv
  integer*4 :: cdt
  integer*4 :: istru
  integer*4 :: ru
  integer*4 :: epsil
  integer*4 :: itemp
  integer*4 :: akt
  integer*4 :: nrhs
  integer*4 :: cfr
  integer*4 :: g
  integer*4 :: grho
  integer*4 :: onefifth
  integer*4 :: halfgrho
  integer*4 :: onetwelfth
  integer*4 :: jendr
  integer*4 :: istrr
  integer*4 :: iendr
  integer*4 :: knew
  integer*4 :: jstrr
  integer*4 :: aa
  integer*4 :: akv
  integer*4 :: cc
  integer*4, dimension(10,10) :: work
  integer*4, dimension(10,10) :: fe
  integer*4, dimension(10,10) :: pm
  integer*4, dimension(10,10) :: pn
  integer*4, dimension(10,10) :: hz_bak
  integer*4, dimension(10,10) :: fc
  integer*4, dimension(10,10) :: cf
  integer*4, dimension(10,10) :: dz
  integer*4, dimension(10,10) :: dr
  integer*4, dimension(10,10) :: on_u
  integer*4, dimension(10,10) :: rufrc
  integer*4, dimension(10,10) :: rvfrc
  integer*4, dimension(10,10) :: bc
  integer*4, dimension(10,10,10) :: fx
  integer*4, dimension(10,10,10) :: huon
  integer*4, dimension(10,10,10) :: hvom
  integer*4, dimension(10,10,10) :: hz
  integer*4, dimension(10,10,10) :: hz_half
  integer*4, dimension(10,10,10) :: rv
  integer*4, dimension(10,10,10) :: we
  integer*4, dimension(10,10,10) :: u
  integer*4, dimension(10,10,10) :: stflx
  integer*4, dimension(10,10,10) :: btflx
  integer*4, dimension(10,10,10) :: z_r
  integer*4, dimension(10,10,10) :: dc
  integer*4, dimension(10,10,10) :: vfx_3d
  integer*4, dimension(10,10,10) :: vfe_3d
  integer*4, dimension(10,10,10) :: ufx_3d
  integer*4, dimension(10,10,10) :: ufe_3d
  integer*4, dimension(10,10,10) :: rho
  integer*4, dimension(10,10,10) :: p
  integer*4, dimension(10,10,10) :: z_w
  integer*4, dimension(10,10,10) :: om_v
  integer*4, dimension(10,10,10) :: du_avg1
  integer*4, dimension(10,10,10) :: du_avg2
  integer*4, dimension(10,10,10) :: ubar
  integer*4, dimension(10,10,10) :: dv_avg1
  integer*4, dimension(10,10,10) :: dv_avg2
  integer*4, dimension(10,10,10) :: vbar
  integer*4, dimension(10,10,10) :: romega
  integer*4, dimension(10,10,10,10) :: v
  integer*4, dimension(10,10,10,10,10) :: t

  do j = jstr, jend, 1
    do k = 1, n, 1
      do i = istru, iend, 1
        dc(i,k) = ufx_3d(i,j,k)
      enddo
    enddo
    do i = istru, iend, 1
      fc(i,0) = 0.
      cf(i,0) = 0.
    enddo
    do k = 1, n - 1, +1
      do i = istru, iend, 1
        cff = 1. / (2. * dc(i,k + 1) + dc(i,k) * (2. - fc(i,k - 1)))
        fc(i,k) = cff * dc(i,k + 1)
        cf(i,k) = cff * (6. * (u(i,j,k + 1,nrhs) - u(i,j,k,nrhs)) - dc(i,k) * cf(i,k - 1))
      enddo
    enddo
    do i = istru, iend, 1
      cf(i,n) = 0.
    enddo
    do k = n - 1, 1, -1
      do i = istru, iend, 1
        cf(i,k) = cf(i,k) - fc(i,k) * cf(i,k + 1)
      enddo
    enddo
    do k = 1, n - 1, 1
      do i = istru, iend, 1
        fc(i,k) = ufe_3d(i,j,k) * (u(i,j,k,nrhs) + dc(i,k) * (0.33333333333333 * cf(i,k) + 0.16666666666667 * cf(i,k - 1)))
      enddo
    enddo
    do i = istru, iend, 1
      fc(i,n) = 0.
      fc(i,0) = 0.
    enddo
    do k = 1, n, 1
      do i = istru, iend, 1
        ru(i,j,k) = ru(i,j,k) - fc(i,k) + fc(i,k - 1)
      enddo
    enddo
  enddo

end subroutine snippet
