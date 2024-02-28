subroutine snippet()
  integer, parameter :: nsa = 28
  integer, parameter :: n2d = 100
  integer, parameter :: n3d = 100
  integer, parameter :: npp = 10
  integer, parameter :: trd = 1
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
  REAL :: A2d(N2d, NSA, 0 : NPP - 1)
  REAL :: A3d(N3d, 8, 0 : NPP - 1)
  COMMON /private_scratch/ A2d, A3d

  do j = jstrv - 1, jend, 1
    do i = istru - 1, iend, 1
      do k = 1, n - 1, 1
        dr1d(k) = rho(i,j,k + 1) - rho(i,j,k)
        dz1d(k) = z_r(i,j,k + 1) - z_r(i,j,k)
      enddo
      dz1d(0) = dz1d(1)
      dz1d(n) = dz1d(n - 1)
      dr1d(0) = dr1d(1)
      dr1d(n) = dr1d(n - 1)
      do k = n, 1, -1
        if (cfr > epsil) then
          dr1d(k) = cfr / (dr1d(k) + dr1d(k - 1))
        else
          dr1d(k) = 0.
        end if
        cfr = 2. * dr1d(k) * dr1d(k - 1)
        dz1d(k) = cff / (dz1d(k) + dz1d(k - 1))
        cff = 2. * dz1d(k) * dz1d(k - 1)
      enddo
      p(i,j,n) = g * z_w(i,j,n) + grho * (rho(i,j,n) + 0.5 * (rho(i,j,n) - rho(i,j,n - 1)) * (z_w(i,j,n) - z_r(i,j,n)) / (z_r(i,j,n) - z_r(i,j,n - 1))) * (z_w(i,j,n) - z_r(i,j,n))
      do k = n - 1, 1, -1
        p(i,j,k) = halfgrho * ((rho(i,j,k + 1) + rho(i,j,k)) * (z_r(i,j,k + 1) - z_r(i,j,k)) - onefifth * ((dr1d(k + 1) - dr1d(k)) * (z_r(i,j,k + 1) - z_r(i,j,k) - onetwelfth * (dz1d(k + 1) + dz1d(k))) - (dz1d(k + 1) - dz1d(k)) * (rho(i,j,k + 1) - rho(i,j,k) - onetwelfth * (dr1d(k + 1) + dr1d(k)))))
      enddo
    enddo
  enddo

end subroutine snippet
