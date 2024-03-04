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
  integer*4 :: jmin
  integer*4 :: jmax
  integer*4 :: imin
  integer*4 :: imax
  integer*4 :: kstp
  integer*4 :: kbak
  integer*4 :: kold
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
  COMMON /private_scratch/ A2d
  COMMON /private_scratch/ A3d

  n = 10
  jstr = n - 1
  jstrr = jstr - 1
  jstrv = jstr + 1
  istru = jstr + 1
  istrr = jstr + 1
  jend = n - 1
  jendr = jend + 1
  istr = n - 1
  iend = n - 1
  iendr = iend + 1
  jmin = jstr - 2
  jmax = jend + 1
  imin = istr - 2
  imax = iend + 1
  knew = 1
  nstp = 1
  nnew = 1
  lm = n - 1
  mm = n - 1
  nadv = nstp
  itrc = 1
  itemp = 1
  kstp = 1
  kbak = 2
  kold = 3
  nrhs = 1
  nt = 10
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      do k = 1, n, 1
        dc1d(k) = vfx_3d(i,j,k)
      enddo
      fc1d(0) = 0.
      cf1d(0) = 0.
      do k = 1, n - 1, +1
        cff = 1. / (2. * dc1d(k + 1) + dc1d(k) * (2. - fc1d(k - 1)))
        fc1d(k) = cff * dc1d(k + 1)
        cf1d(k) = cff * (6. * (v(i,j,k + 1,nrhs) - v(i,j,k,nrhs)) - dc1d(k) * cf1d(k - 1))
      enddo
      cf1d(n) = 0.
      do k = n - 1, 1, -1
        cf1d(k) = cf1d(k) - fc1d(k) * cf1d(k + 1)
      enddo
      do k = 1, n - 1, 1
        fc1d(k) = vfe_3d(i,j,k) * (v(i,j,k,nrhs) + dc1d(k) * (0.33333333333333 * cf1d(k) + 0.16666666666667 * cf1d(k - 1)))
      enddo
      fc1d(n) = 0.
      fc1d(0) = 0.
      do k = 1, n, 1
        rv(i,j,k) = rv(i,j,k) - fc1d(k) + fc1d(k - 1)
      enddo
    enddo
  enddo

end subroutine snippet
