do k = 1, n, 1
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ufe_3d(i,j,k) = 0.5625 * (we(i,j,k) + we(i - 1,j,k)) - 0.0625 * (we(i + 1,j,k) + we(i - 2,j,k))
    enddo
  enddo
enddo
