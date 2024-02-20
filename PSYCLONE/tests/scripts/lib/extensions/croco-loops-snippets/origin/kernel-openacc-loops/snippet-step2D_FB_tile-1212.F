do j = jstrv, jend, 1
  do i = istr, iend, 1
    rvbar(i,j) = rvbar(i,j) - 0.25 * (vfe(i,j) + vfe(i,j - 1))
  enddo
enddo
