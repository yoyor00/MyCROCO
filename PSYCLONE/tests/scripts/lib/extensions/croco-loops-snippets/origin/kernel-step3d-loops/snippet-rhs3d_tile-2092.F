do k = 1, n, 1
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      vfe_3d(i,j,k) = 0.5625 * (we(i,j,k) + we(i,j - 1,k)) - 0.0625 * (we(i,j + 1,k) + we(i,j - 2,k))
    enddo
  enddo
enddo
