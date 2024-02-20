do j = jstrv, jend, 1
  do i = istr, iend, 1
    rvbar(i,j) = rvbar(i,j) - rdrg * vbar(i,j,kstp) * om_v(i,j) * on_v(i,j)
  enddo
enddo
