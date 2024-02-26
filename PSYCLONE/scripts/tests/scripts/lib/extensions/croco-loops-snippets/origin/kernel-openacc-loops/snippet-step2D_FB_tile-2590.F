do i = istrr, iendr, 1
  dv_avg1(i,jstrv - 1,nnew) = dv_avg1(i,jstrv - 1,nnew) + cff1 * (dnew(i,jstrv - 1) + dnew(i,jstrv - 2)) * vbar(i,jstrv - 1,knew) * om_v(i,jstrv - 1)
enddo
