do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    duon(i,j) = zeta(i,j,kstp) + h(i,j)
  enddo
enddo
