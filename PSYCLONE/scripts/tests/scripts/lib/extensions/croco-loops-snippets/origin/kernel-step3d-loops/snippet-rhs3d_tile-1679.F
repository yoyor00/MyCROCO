do k = 1, n, 1
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ufx_3d(i,j,k) = 0.5625 * (hz(i,j,k) + hz(i - 1,j,k)) - 0.0625 * (hz(i + 1,j,k) + hz(i - 2,j,k))
    enddo
  enddo
enddo
