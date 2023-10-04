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
