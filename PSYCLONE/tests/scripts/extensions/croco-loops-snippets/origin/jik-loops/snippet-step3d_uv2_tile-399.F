do j = jstrv, jend, 1
  do i = istr, iend, 1
    fc(i,0) = 0.
    fc(i,n) = 0.
  enddo
  do k = 1, n - 1, 1
    do i = istr, iend, 1
      fc(i,k) = -dt * (akv(i,j,k) + akv(i,j - 1,k)) / (z_r(i,j,k + 1) + z_r(i,j - 1,k + 1) - z_r(i,j,k) - z_r(i,j - 1,k))
    enddo
  enddo
  do k = 1, n, 1
    do i = istr, iend, 1
      dc(i,k) = v(i,j,k,nnew)
      bc(i,k) = 0.5 * (hz(i,j,k) + hz(i,j - 1,k)) - fc(i,k) - fc(i,k - 1)
    enddo
  enddo
  do i = istr, iend, 1
    cff = 1. / bc(i,1)
    cf(i,1) = cff * fc(i,1)
    dc(i,1) = cff * dc(i,1)
  enddo
  do k = 2, n - 1, 1
    do i = istr, iend, 1
      aa = fc(i,k - 1)
      cc = fc(i,k)
      cff = 1. / (bc(i,k) - aa * cf(i,k - 1))
      cf(i,k) = cff * cc
      dc(i,k) = cff * (dc(i,k) - aa * dc(i,k - 1))
    enddo
  enddo
  do i = istr, iend, 1
    aa = fc(i,n - 1)
    dc(i,n) = (dc(i,n) - aa * dc(i,n - 1)) / (bc(i,n) - aa * cf(i,n - 1))
    cf(i,0) = 0.5 * (hz(i,j,n) + hz(i,j - 1,n))
    dc(i,0) = cf(i,0) * dc(i,n)
  enddo
  do k = n - 1, 1, -1
    do i = istr, iend, 1
      dc(i,k) = dc(i,k) - cf(i,k) * dc(i,k + 1)
      cff = 0.5 * (hz(i,j,k) + hz(i,j - 1,k))
      cf(i,0) = cf(i,0) + cff
      dc(i,0) = dc(i,0) + cff * dc(i,k)
    enddo
  enddo
  do i = istr, iend, 1
    dc(i,0) = (dc(i,0) * om_v(i,j) - dv_avg1(i,j,nnew)) / (cf(i,0) * om_v(i,j))
  enddo
  do k = 1, n, 1
    do i = istr, iend, 1
      v(i,j,k,nnew) = dc(i,k) - dc(i,0)
    enddo
  enddo
  do k = 1, n, 1
    do i = istr, iend, 1
    enddo
  enddo
enddo
