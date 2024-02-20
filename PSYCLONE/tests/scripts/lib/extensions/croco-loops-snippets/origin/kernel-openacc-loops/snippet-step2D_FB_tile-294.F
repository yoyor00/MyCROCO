do j = jmin + 1, jmax, 1
  do i = istr - 1, iend + 1, 1
    vrhs(i,j) = cff1 * vbar(i,j,kstp) + cff2 * vbar(i,j,kbak) + cff3 * vbar(i,j,kold)
    dvom(i,j) = 0.5 * (drhs(i,j) + drhs(i,j - 1)) * om_v(i,j) * vrhs(i,j)
  enddo
enddo
