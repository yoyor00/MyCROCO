do k = 1, n, 1
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      vfx_3d(i,j,k) = 0.5625 * (hz(i,j,k) + hz(i,j - 1,k)) - 0.0625 * (hz(i,j + 1,k) + hz(i,j - 2,k))
    enddo
  enddo
enddo
