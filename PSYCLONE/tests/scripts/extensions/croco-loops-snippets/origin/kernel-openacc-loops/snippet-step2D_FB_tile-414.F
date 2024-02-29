do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    zeta_new(i,j) = zeta(i,j,kstp) + dtfast * pm(i,j) * pn(i,j) * (duon(i,j) - duon(i + 1,j) + dvom(i,j) - dvom(i,j + 1))
  enddo
enddo
