do j = jstr, jend, 1
  do i = istr, iend, 1
    rubar(i,j) = rubar(i,j) - ufx(i,j) + ufx(i - 1,j) - ufe(i,j + 1) + ufe(i,j)
    rvbar(i,j) = rvbar(i,j) - vfx(i + 1,j) + vfx(i,j) - vfe(i,j) + vfe(i,j - 1)
  enddo
enddo
