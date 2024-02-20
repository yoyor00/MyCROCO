do j = jstr, jend, 1
  do i = istr, iend, 1
    cf(i,j) = om_r(i,j) * on_r(i,j) / dt
  enddo
enddo
