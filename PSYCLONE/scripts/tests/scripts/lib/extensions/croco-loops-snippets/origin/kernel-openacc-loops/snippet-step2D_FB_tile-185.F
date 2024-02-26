do j = jmin, jmax, 1
  do i = imin, imax, 1
    drhs(i,j) = cff1 * zeta(i,j,kstp) + cff2 * zeta(i,j,kbak) + cff3 * zeta(i,j,kold) + h(i,j)
  enddo
enddo
