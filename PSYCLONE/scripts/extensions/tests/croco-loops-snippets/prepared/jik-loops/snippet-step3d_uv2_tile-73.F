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
    do i = istru, iend, 1
      fc(i,0) = 0.
      fc(i,n) = 0.
    enddo
    do k = 1, n - 1, 1
      do i = istru, iend, 1
        fc(i,k) = -dt * (akv(i,j,k) + akv(i - 1,j,k)) / (z_r(i,j,k + 1) + z_r(i - 1,j,k + 1) - z_r(i,j,k) - z_r(i - 1,j,k))
      enddo
    enddo
    do k = 1, n, 1
      do i = istru, iend, 1
        dc(i,k) = u(i,j,k,nnew)
        bc(i,k) = 0.5 * (hz(i,j,k) + hz(i - 1,j,k)) - fc(i,k) - fc(i,k - 1)
      enddo
    enddo
    do i = istru, iend, 1
      cff = 1. / bc(i,1)
      cf(i,1) = cff * fc(i,1)
      dc(i,1) = cff * dc(i,1)
    enddo
    do k = 2, n - 1, 1
      do i = istru, iend, 1
        aa = fc(i,k - 1)
        cc = fc(i,k)
        cff = 1. / (bc(i,k) - aa * cf(i,k - 1))
        cf(i,k) = cff * cc
        dc(i,k) = cff * (dc(i,k) - aa * dc(i,k - 1))
      enddo
    enddo
    do i = istru, iend, 1
      aa = fc(i,n - 1)
      dc(i,n) = (dc(i,n) - aa * dc(i,n - 1)) / (bc(i,n) - aa * cf(i,n - 1))
      cf(i,0) = 0.5 * (hz(i,j,n) + hz(i - 1,j,n))
      dc(i,0) = cf(i,0) * dc(i,n)
    enddo
    do k = n - 1, 1, -1
      do i = istru, iend, 1
        dc(i,k) = dc(i,k) - cf(i,k) * dc(i,k + 1)
        cff = 0.5 * (hz(i,j,k) + hz(i - 1,j,k))
        cf(i,0) = cf(i,0) + cff
        dc(i,0) = dc(i,0) + cff * dc(i,k)
      enddo
    enddo
    do i = istru, iend, 1
      dc(i,0) = (dc(i,0) * on_u(i,j) - du_avg1(i,j,nnew)) / (cf(i,0) * on_u(i,j))
    enddo
    do k = 1, n, 1
      do i = istru, iend, 1
        u(i,j,k,nnew) = dc(i,k) - dc(i,0)
      enddo
    enddo
    do k = 1, n, 1
      do i = istru, iend, 1
      enddo
    enddo
  enddo

end subroutine snippet
