do j = jstrv, jend, 1
  dv_avg1(istr - 1,j,nnew) = dv_avg1(istr - 1,j,nnew) + cff1 * (dnew(istr - 1,j) + dnew(istr - 1,j - 1)) * vbar(istr - 1,j,knew) * om_v(istr - 1,j)
enddo
