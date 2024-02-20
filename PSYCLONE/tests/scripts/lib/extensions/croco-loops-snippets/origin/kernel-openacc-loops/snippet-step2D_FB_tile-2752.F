do j = jstr, jend, 1
  do i = istr, iend, 1
    vmax = MAX(vmax, ABS(ubar(i,j,knew)))
  enddo
enddo
