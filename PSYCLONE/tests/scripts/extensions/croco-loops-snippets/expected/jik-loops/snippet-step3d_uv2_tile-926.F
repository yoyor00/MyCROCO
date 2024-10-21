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
  real*4, dimension(0:11) :: fc1d
  real*4, dimension(0:11) :: cf1d
  real*4, dimension(0:11) :: dc1d
  real*4, dimension(0:11) :: bc1d
  real*4, dimension(0:11) :: dz1d
  real*4, dimension(0:11) :: dr1d
  real*4, dimension(0:11,0:11) :: fx
  real*4, dimension(0:11,0:11) :: work
  real*4, dimension(0:11,0:11) :: fe
  real*4, dimension(0:11,0:11) :: pm
  real*4, dimension(0:11,0:11) :: pn
  real*4, dimension(0:11,0:11) :: fc
  real*4, dimension(0:11,0:11) :: cf
  real*4, dimension(0:11,0:11) :: dz
  real*4, dimension(0:11,0:11) :: dr
  real*4, dimension(0:11,0:11) :: on_u
  real*4, dimension(0:11,0:11) :: rufrc
  real*4, dimension(0:11,0:11) :: rvfrc
  real*4, dimension(0:11,0:11) :: bc
  real*4, dimension(0:11,0:11) :: dc
  real*4, dimension(0:11,0:11) :: om_v
  real*4, dimension(0:11,0:11) :: dv_avg2
  real*4, dimension(0:11,0:11) :: du_avg2
  real*4, dimension(0:11,0:11,0:11) :: huon
  real*4, dimension(0:11,0:11,0:11) :: hvom
  real*4, dimension(0:11,0:11,0:11) :: hz
  real*4, dimension(0:11,0:11,0:11) :: hz_half
  real*4, dimension(0:11,0:11,0:11) :: rv
  real*4, dimension(0:11,0:11,0:11) :: we
  real*4, dimension(0:11,0:11,0:11) :: stflx
  real*4, dimension(0:11,0:11,0:11) :: btflx
  real*4, dimension(0:11,0:11,0:11) :: z_r
  real*4, dimension(0:11,0:11,0:11) :: vfx_3d
  real*4, dimension(0:11,0:11,0:11) :: vfe_3d
  real*4, dimension(0:11,0:11,0:11) :: ufx_3d
  real*4, dimension(0:11,0:11,0:11) :: ufe_3d
  real*4, dimension(0:11,0:11,0:11) :: rho
  real*4, dimension(0:11,0:11,0:11) :: p
  real*4, dimension(0:11,0:11,0:11) :: z_w
  real*4, dimension(0:11,0:11,0:11) :: du_avg1
  real*4, dimension(0:11,0:11,0:11) :: ubar
  real*4, dimension(0:11,0:11,0:11) :: dv_avg1
  real*4, dimension(0:11,0:11,0:11) :: vbar
  real*4, dimension(0:11,0:11,0:11) :: romega
  real*4, dimension(0:11,0:11,0:11) :: ru
  real*4, dimension(0:11,0:11,0:11) :: work_3d
  real*4, dimension(0:11,0:11,0:11) :: hz_bak
  real*4, dimension(0:11,0:11,0:11) :: fx_3d
  real*4, dimension(0:11,0:11,0:11) :: fe_3d
  real*4, dimension(0:11,0:11,0:11) :: akv
  real*4, dimension(0:11,0:11,0:11,0:11) :: v
  real*4, dimension(0:11,0:11,0:11,0:11) :: u
  real*4, dimension(0:11,0:11,0:11,0:11) :: akt
  real*4, dimension(0:11,0:11,0:11,0:11,0:11) :: t
  REAL :: A2d(N2d, NSA, 0 : NPP - 1)
  REAL :: A3d(N3d, 8, 0 : NPP - 1)
  COMMON /private_scratch/ A2d
  COMMON /private_scratch/ A3d

  n = 10 - 1
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
  nt = 9
  indx = 3 - nstp
  do j = jstr, jendr, 1
    do i = istrr, iendr, 1
      dc1d(0) = 0.
      cf1d(0) = 0.
      fc1d(0) = 0.
      do k = 1, n, +1
        dc1d(k) = 0.5 * (hz(i,j,k) + hz(i,j - 1,k)) * om_v(i,j)
        dc1d(0) = dc1d(0) + dc1d(k)
        cf1d(0) = cf1d(0) + dc1d(k) * v(i,j,k,nnew)
      enddo
      dc1d(0) = 1. / dc1d(0)
      cf1d(0) = dc1d(0) * (cf1d(0) - dv_avg1(i,j,nnew))
      vbar(i,j,knew) = dc1d(0) * dv_avg1(i,j,nnew)
      do k = n, 1, -1
        v(i,j,k,nnew) = v(i,j,k,nnew) - cf1d(0)
        fc1d(k) = 0.75 * hvom(i,j,k) + 0.125 * dc1d(k) * (v(i,j,k,nstp) + v(i,j,k,nnew))
        fc1d(0) = fc1d(0) + fc1d(k)
      enddo
      fc1d(0) = dc1d(0) * (fc1d(0) - dv_avg2(i,j))
      do k = 1, n, +1
        hvom(i,j,k) = fc1d(k) - dc1d(k) * fc1d(0)
      enddo
    enddo
  enddo

end subroutine snippet
