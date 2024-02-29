do j = jstr - 1, jend, 1
  do i = istr, iend, 1
    vfe(i,j) = 0.25 * (dvom(i,j) + dvom(i,j + 1)) * (vrhs(i,j) + vrhs(i,j + 1))
    ufe(i,j + 1) = 0.25 * (dvom(i,j + 1) + dvom(i - 1,j + 1)) * (urhs(i,j + 1) + urhs(i,j))
  enddo
enddo
