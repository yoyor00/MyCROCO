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
