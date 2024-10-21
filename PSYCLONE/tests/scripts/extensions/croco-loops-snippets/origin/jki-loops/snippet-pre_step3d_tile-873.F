do j = jstr, jend, 1
  do k = 1, n, 1
    do i = istr, iend, 1
      dc(i,k) = 1. / hz_half(i,j,k)
    enddo
  enddo
  do i = istr, iend, 1
    dc(i,0) = cdt * pn(i,j) * pm(i,j)
  enddo
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
        t(i,j,k,nnew,itrc) = dc(i,k) * (t(i,j,k,nnew,itrc) - dc(i,0) * (fc(i,k) - fc(i,k - 1)))
      enddo
    enddo
  enddo
enddo
