do j = jstr, jend, 1
  do i = istru, iend, 1
    cff = 0.25 * (vbar(i,j,kstp) + vbar(i,j + 1,kstp) + vbar(i - 1,j,kstp) + vbar(i - 1,j + 1,kstp))
    rubar(i,j) = rubar(i,j) - ubar(i,j,kstp) * (rdrg + rdrg2 * SQRT(ubar(i,j,kstp) * ubar(i,j,kstp) + cff * cff)) * om_u(i,j) * on_u(i,j)
  enddo
enddo
