do j = jstrv, jend, 1
  dv_avg1(iend + 1,j,nnew) = dv_avg1(iend + 1,j,nnew) + cff1 * (dnew(iend + 1,j) + dnew(iend + 1,j - 1)) * vbar(iend + 1,j,knew) * om_v(iend + 1,j)
enddo
