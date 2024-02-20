do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    dnew(i,j) = zeta_new(i,j) + h(i,j)
    zeta(i,j,knew) = zeta_new(i,j)
  enddo
enddo
