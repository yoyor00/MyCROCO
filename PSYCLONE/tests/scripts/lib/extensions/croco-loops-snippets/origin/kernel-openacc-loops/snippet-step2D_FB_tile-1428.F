do j = jstr, jend, 1
  do i = istru, iend, 1
    rubar(i,j) = rubar(i,j) - rdrg * ubar(i,j,kstp) * om_u(i,j) * on_u(i,j)
  enddo
enddo
