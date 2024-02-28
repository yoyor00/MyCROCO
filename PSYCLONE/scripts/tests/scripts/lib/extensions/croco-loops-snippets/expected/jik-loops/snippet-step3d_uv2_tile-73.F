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

  do j = jstr, jend, 1
    do i = istru, iend, 1
      fc1d(0) = 0.
      fc1d(n) = 0.
      do k = 1, n - 1, 1
        fc1d(k) = -dt * (akv(i,j,k) + akv(i - 1,j,k)) / (z_r(i,j,k + 1) + z_r(i - 1,j,k + 1) - z_r(i,j,k) - z_r(i - 1,j,k))
      enddo
      do k = 1, n, 1
        dc1d(k) = u(i,j,k,nnew)
        bc1d(k) = 0.5 * (hz(i,j,k) + hz(i - 1,j,k)) - fc1d(k) - fc1d(k - 1)
      enddo
      cff = 1. / bc1d(1)
      cf1d(1) = cff * fc1d(1)
      dc1d(1) = cff * dc1d(1)
      do k = 2, n - 1, 1
        aa = fc1d(k - 1)
        cc = fc1d(k)
        cff = 1. / (bc1d(k) - aa * cf1d(k - 1))
        cf1d(k) = cff * cc
        dc1d(k) = cff * (dc1d(k) - aa * dc1d(k - 1))
      enddo
      aa = fc1d(n - 1)
      dc1d(n) = (dc1d(n) - aa * dc1d(n - 1)) / (bc1d(n) - aa * cf1d(n - 1))
      cf1d(0) = 0.5 * (hz(i,j,n) + hz(i - 1,j,n))
      dc1d(0) = cf1d(0) * dc1d(n)
      do k = n - 1, 1, -1
        dc1d(k) = dc1d(k) - cf1d(k) * dc1d(k + 1)
        cff = 0.5 * (hz(i,j,k) + hz(i - 1,j,k))
        cf1d(0) = cf1d(0) + cff
        dc1d(0) = dc1d(0) + cff * dc1d(k)
      enddo
      dc1d(0) = (dc1d(0) * on_u(i,j) - du_avg1(i,j,nnew)) / (cf1d(0) * on_u(i,j))
      do k = 1, n, 1
        u(i,j,k,nnew) = dc1d(k) - dc1d(0)
      enddo
      do k = 1, n, 1
      enddo
    enddo
  enddo

end subroutine snippet
