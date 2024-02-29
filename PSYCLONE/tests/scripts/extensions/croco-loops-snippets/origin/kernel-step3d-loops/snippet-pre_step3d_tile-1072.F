do j = jstr, jend, 1
  do i = istru, iend, 1
    work(i,j) = pm_u(i,j) * pn_u(i,j)
  enddo
enddo
