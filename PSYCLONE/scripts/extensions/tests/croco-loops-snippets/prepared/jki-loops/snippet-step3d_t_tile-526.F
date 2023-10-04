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
    do itrc = 1, nt, 1
      do k = 1, n - 1, 1
        do i = istr, iend, 1
          fc(i,k) = t(i,j,k + 1,nadv,itrc) - t(i,j,k,nadv,itrc)
        enddo
      enddo
      do i = istr, iend, 1
        fc(i,0) = fc(i,1)
        fc(i,n) = fc(i,n - 1)
      enddo
      do k = 1, n, 1
        do i = istr, iend, 1
          cff = 2. * fc(i,k) * fc(i,k - 1)
          if (cff > epsil) then
            cf(i,k) = cff / (fc(i,k) + fc(i,k - 1))
          else
            cf(i,k) = 0.
          end if
        enddo
      enddo
      do k = 1, n - 1, 1
        do i = istr, iend, 1
          fc(i,k) = 0.5 * (t(i,j,k,nadv,itrc) + t(i,j,k + 1,nadv,itrc) - 0.333333333333 * (cf(i,k + 1) - cf(i,k))) * we(i,j,k)
        enddo
      enddo
      do i = istr, iend, 1
        fc(i,0) = 0.
        fc(i,n) = 0.
        cf(i,0) = dt * pm(i,j) * pn(i,j)
      enddo
      do k = 1, n, 1
        do i = istr, iend, 1
          t(i,j,k,nnew,itrc) = t(i,j,k,nnew,itrc) - cf(i,0) * (fc(i,k) - fc(i,k - 1))
        enddo
      enddo
      do i = istr, iend, 1
        fc(i,n) = dt * stflx(i,j,itrc)
        fc(i,0) = -dt * btflx(i,j,itrc)
      enddo
      if (itrc == itemp) then
        do k = 1, n - 1, 1
          do i = istr, iend, 1
            fc(i,k) = 0.
          enddo
        enddo
      end if
      if (itrc == itemp) then
        do k = 1, n, 1
          do i = istr, iend, 1
            t(i,j,k,nnew,itrc) = t(i,j,k,nnew,itrc) + fc(i,k) - fc(i,k - 1)
          enddo
        enddo
      end if
      indx = MIN(itrc, itemp)
      do i = istr, iend, 1
        fc(i,1) = dt * akt(i,j,1,indx) / (z_r(i,j,2) - z_r(i,j,1))
        cff = 1. / (hz(i,j,1) + fc(i,1))
        cf(i,1) = cff * fc(i,1)
        dc(i,1) = cff * t(i,j,1,nnew,itrc)
      enddo
      do k = 2, n - 1, +1
        do i = istr, iend, 1
          fc(i,k) = dt * akt(i,j,k,indx) / (z_r(i,j,k + 1) - z_r(i,j,k))
          cff = 1. / (hz(i,j,k) + fc(i,k) + fc(i,k - 1) * (1. - cf(i,k - 1)))
          cf(i,k) = cff * fc(i,k)
          dc(i,k) = cff * (t(i,j,k,nnew,itrc) + fc(i,k - 1) * dc(i,k - 1))
        enddo
      enddo
      do i = istr, iend, 1
        t(i,j,n,nnew,itrc) = (t(i,j,n,nnew,itrc) + fc(i,n - 1) * dc(i,n - 1)) / (hz(i,j,n) + fc(i,n - 1) * (1. - cf(i,n - 1)))
      enddo
      do k = n - 1, 1, -1
        do i = istr, iend, 1
          t(i,j,k,nnew,itrc) = dc(i,k) + cf(i,k) * t(i,j,k + 1,nnew,itrc)
        enddo
      enddo
    enddo
  enddo

end subroutine snippet
