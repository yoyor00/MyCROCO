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

  do itrc = 1, nt, 1
    do k = 1, n, 1
      do j = jstr, jend, 1
        do i = MAX(istr - 1, 1), MIN(iend + 2, lm + 1), 1
          fx(i,j) = t(i,j,k,nadv,itrc) - t(i - 1,j,k,nadv,itrc)
        enddo
      enddo
      if (istr == 1) then
        do j = jstr, jend, 1
          fx(0,j) = fx(1,j)
        enddo
      end if
      if (iend == lm) then
        do j = jstr, jend, 1
          fx(lm + 2,j) = fx(lm + 1,j)
        enddo
      end if
      do j = jstr, jend, 1
        do i = istr - 1, iend + 1, 1
          work(i,j) = 0.5 * (fx(i + 1,j) + fx(i,j))
        enddo
      enddo
      do j = jstr, jend, 1
        do i = istr, iend + 1, 1
          fx(i,j) = 0.5 * (t(i,j,k,nadv,itrc) + t(i - 1,j,k,nadv,itrc) - 0.333333333333 * (work(i,j) - work(i - 1,j))) * huon(i,j,k)
        enddo
      enddo
      do j = MAX(jstr - 1, 1), MIN(jend + 2, mm + 1), 1
        do i = istr, iend, 1
          fe(i,j) = t(i,j,k,nadv,itrc) - t(i,j - 1,k,nadv,itrc)
        enddo
      enddo
      if (jstr == 1) then
        do i = istr, iend, 1
          fe(i,0) = fe(i,1)
        enddo
      end if
      if (jend == mm) then
        do i = istr, iend, 1
          fe(i,mm + 2) = fe(i,mm + 1)
        enddo
      end if
      do j = jstr - 1, jend + 1, 1
        do i = istr, iend, 1
          work(i,j) = 0.5 * (fe(i,j + 1) + fe(i,j))
        enddo
      enddo
      do j = jstr, jend + 1, 1
        do i = istr, iend, 1
          fe(i,j) = 0.5 * (t(i,j,k,nadv,itrc) + t(i,j - 1,k,nadv,itrc) - 0.333333333333 * (work(i,j) - work(i,j - 1))) * hvom(i,j,k)
        enddo
      enddo
      if (iic == ntstart) then
        cff = 0.5 * dt
        do j = jstr, jend, 1
          do i = istr, iend, 1
            t(i,j,k,nnew,itrc) = hz(i,j,k) * t(i,j,k,nstp,itrc) - cff * pm(i,j) * pn(i,j) * (fx(i + 1,j) - fx(i,j) + fe(i,j + 1) - fe(i,j))
          enddo
        enddo
      else
        cff = (1. - gamma) * dt
        cff1 = 0.5 + gamma
        cff2 = 0.5 - gamma
        do j = jstr, jend, 1
          do i = istr, iend, 1
            t(i,j,k,nnew,itrc) = cff1 * hz(i,j,k) * t(i,j,k,nstp,itrc) + cff2 * hz_bak(i,j,k) * t(i,j,k,indx,itrc) - cff * pm(i,j) * pn(i,j) * (fx(i + 1,j) - fx(i,j) + fe(i,j + 1) - fe(i,j))
          enddo
        enddo
      end if
    enddo
  enddo

end subroutine snippet
