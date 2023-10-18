do j = jstrv, jend, 1
  do k = 1, n, 1
    do i = istr, iend, 1
      dc(i,k) = vfx_3d(i,j,k)
    enddo
  enddo
  do i = istr, iend, 1
    fc(i,0) = 0.
    cf(i,0) = 0.
  enddo
  do k = 1, n - 1, +1
    do i = istr, iend, 1
      cff = 1. / (2. * dc(i,k + 1) + dc(i,k) * (2. - fc(i,k - 1)))
      fc(i,k) = cff * dc(i,k + 1)
      cf(i,k) = cff * (6. * (v(i,j,k + 1,nrhs) - v(i,j,k,nrhs)) - dc(i,k) * cf(i,k - 1))
    enddo
  enddo
  do i = istr, iend, 1
    cf(i,n) = 0.
  enddo
  do k = n - 1, 1, -1
    do i = istr, iend, 1
      cf(i,k) = cf(i,k) - fc(i,k) * cf(i,k + 1)
    enddo
  enddo
  do k = 1, n - 1, 1
    do i = istr, iend, 1
      fc(i,k) = vfe_3d(i,j,k) * (v(i,j,k,nrhs) + dc(i,k) * (0.33333333333333 * cf(i,k) + 0.16666666666667 * cf(i,k - 1)))
    enddo
  enddo
  do i = istr, iend, 1
    fc(i,n) = 0.
    fc(i,0) = 0.
  enddo
  do k = 1, n, 1
    do i = istr, iend, 1
      rv(i,j,k) = rv(i,j,k) - fc(i,k) + fc(i,k - 1)
    enddo
  enddo
enddo
