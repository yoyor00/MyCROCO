do j = jstrv, jend, 1
  do i = istr, iend, 1
    work(i,j) = pm_v(i,j) * pn_v(i,j)
  enddo
enddo
