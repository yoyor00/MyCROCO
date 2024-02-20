do j = jstrv, jend, 1
  do i = istr, iend, 1
    cff = 0.25 * (ubar(i,j,kstp) + ubar(i + 1,j,kstp) + ubar(i,j - 1,kstp) + ubar(i + 1,j - 1,kstp))
    rvbar(i,j) = rvbar(i,j) - vbar(i,j,kstp) * (rdrg + rdrg2 * SQRT(cff * cff + vbar(i,j,kstp) * vbar(i,j,kstp))) * om_v(i,j) * on_v(i,j)
  enddo
enddo
