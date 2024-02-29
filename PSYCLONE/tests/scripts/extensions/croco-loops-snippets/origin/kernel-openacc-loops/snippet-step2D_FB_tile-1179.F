do j = jstr, jend, 1
  do i = istru, iend, 1
    rubar(i,j) = rubar(i,j) + 0.25 * (ufx(i,j) + ufx(i - 1,j))
  enddo
enddo
