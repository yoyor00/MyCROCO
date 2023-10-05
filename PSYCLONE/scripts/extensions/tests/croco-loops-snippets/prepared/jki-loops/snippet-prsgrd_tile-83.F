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
  integer*4 :: nnew
  integer*4 :: nstp
  integer*4 :: indx
  integer*4 :: jstrv
  integer*4 :: istru
  integer*4 :: itemp
  integer*4 :: nrhs
  integer*4 :: jendr
  integer*4 :: istrr
  integer*4 :: iendr
  integer*4 :: knew
  integer*4 :: jstrr
  real*4 :: cff
  real*4 :: dt
  real*4 :: gamma
  real*4 :: cff1
  real*4 :: cff2
  real*4 :: cdt
  real*4 :: epsil
  real*4 :: cfr
  real*4 :: g
  real*4 :: grho
  real*4 :: onefifth
  real*4 :: halfgrho
  real*4 :: onetwelfth
  real*4 :: aa
  real*4 :: cc
  real*4, dimension(0:10) :: fc1d
  real*4, dimension(0:10) :: cf1d
  real*4, dimension(0:10) :: dc1d
  real*4, dimension(0:10) :: bc1d
  real*4, dimension(0:10) :: dz1d
  real*4, dimension(0:10) :: dr1d
  real*4, dimension(0:10,0:10) :: fx
  real*4, dimension(0:10,0:10) :: work
  real*4, dimension(0:10,0:10) :: fe
  real*4, dimension(0:10,0:10) :: pm
  real*4, dimension(0:10,0:10) :: pn
  real*4, dimension(0:10,0:10) :: fc
  real*4, dimension(0:10,0:10) :: cf
  real*4, dimension(0:10,0:10) :: dz
  real*4, dimension(0:10,0:10) :: dr
  real*4, dimension(0:10,0:10) :: on_u
  real*4, dimension(0:10,0:10) :: rufrc
  real*4, dimension(0:10,0:10) :: rvfrc
  real*4, dimension(0:10,0:10) :: bc
  real*4, dimension(0:10,0:10) :: dc
  real*4, dimension(0:10,0:10) :: om_v
  real*4, dimension(0:10,0:10) :: dv_avg2
  real*4, dimension(0:10,0:10) :: du_avg2
  real*4, dimension(0:10,0:10,0:10) :: huon
  real*4, dimension(0:10,0:10,0:10) :: hvom
  real*4, dimension(0:10,0:10,0:10) :: hz
  real*4, dimension(0:10,0:10,0:10) :: hz_half
  real*4, dimension(0:10,0:10,0:10) :: rv
  real*4, dimension(0:10,0:10,0:10) :: we
  real*4, dimension(0:10,0:10,0:10) :: stflx
  real*4, dimension(0:10,0:10,0:10) :: btflx
  real*4, dimension(0:10,0:10,0:10) :: z_r
  real*4, dimension(0:10,0:10,0:10) :: vfx_3d
  real*4, dimension(0:10,0:10,0:10) :: vfe_3d
  real*4, dimension(0:10,0:10,0:10) :: ufx_3d
  real*4, dimension(0:10,0:10,0:10) :: ufe_3d
  real*4, dimension(0:10,0:10,0:10) :: rho
  real*4, dimension(0:10,0:10,0:10) :: p
  real*4, dimension(0:10,0:10,0:10) :: z_w
  real*4, dimension(0:10,0:10,0:10) :: du_avg1
  real*4, dimension(0:10,0:10,0:10) :: ubar
  real*4, dimension(0:10,0:10,0:10) :: dv_avg1
  real*4, dimension(0:10,0:10,0:10) :: vbar
  real*4, dimension(0:10,0:10,0:10) :: romega
  real*4, dimension(0:10,0:10,0:10) :: ru
  real*4, dimension(0:10,0:10,0:10) :: work_3d
  real*4, dimension(0:10,0:10,0:10) :: hz_bak
  real*4, dimension(0:10,0:10,0:10) :: fx_3d
  real*4, dimension(0:10,0:10,0:10) :: fe_3d
  real*4, dimension(0:10,0:10,0:10) :: akv
  real*4, dimension(0:10,0:10,0:10,0:10) :: v
  real*4, dimension(0:10,0:10,0:10,0:10) :: u
  real*4, dimension(0:10,0:10,0:10,0:10) :: akt
  real*4, dimension(0:10,0:10,0:10,0:10,0:10) :: t

  do j = jstrv - 1, jend, 1
    do k = 1, n - 1, 1
      do i = istru - 1, iend, 1
        dz(i,k) = z_r(i,j,k + 1) - z_r(i,j,k)
        dr(i,k) = rho(i,j,k + 1) - rho(i,j,k)
      enddo
    enddo
    do i = istru - 1, iend, 1
      dr(i,n) = dr(i,n - 1)
      dr(i,0) = dr(i,1)
      dz(i,n) = dz(i,n - 1)
      dz(i,0) = dz(i,1)
    enddo
    do k = n, 1, -1
      do i = istru - 1, iend, 1
        cff = 2. * dz(i,k) * dz(i,k - 1)
        dz(i,k) = cff / (dz(i,k) + dz(i,k - 1))
        cfr = 2. * dr(i,k) * dr(i,k - 1)
        if (cfr > epsil) then
          dr(i,k) = cfr / (dr(i,k) + dr(i,k - 1))
        else
          dr(i,k) = 0.
        end if
      enddo
    enddo
    do i = istru - 1, iend, 1
      p(i,j,n) = g * z_w(i,j,n) + grho * (rho(i,j,n) + 0.5 * (rho(i,j,n) - rho(i,j,n - 1)) * (z_w(i,j,n) - z_r(i,j,n)) / (z_r(i,j,n) - z_r(i,j,n - 1))) * (z_w(i,j,n) - z_r(i,j,n))
    enddo
    do k = n - 1, 1, -1
      do i = istru - 1, iend, 1
        p(i,j,k) = halfgrho * ((rho(i,j,k + 1) + rho(i,j,k)) * (z_r(i,j,k + 1) - z_r(i,j,k)) - onefifth * ((dr(i,k + 1) - dr(i,k)) * (z_r(i,j,k + 1) - z_r(i,j,k) - onetwelfth * (dz(i,k + 1) + dz(i,k))) - (dz(i,k + 1) - dz(i,k)) * (rho(i,j,k + 1) - rho(i,j,k) - onetwelfth * (dr(i,k + 1) + dr(i,k)))))
      enddo
    enddo
  enddo

end subroutine snippet
