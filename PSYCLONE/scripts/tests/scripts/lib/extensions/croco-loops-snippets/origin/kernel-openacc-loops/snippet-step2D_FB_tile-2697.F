do i = istrr, iendr, 1
  dv_avg1(i,jend + 1,nnew) = dv_avg1(i,jend + 1,nnew) + cff1 * (dnew(i,jend + 1) + dnew(i,jend)) * vbar(i,jend + 1,knew) * om_v(i,jend + 1)
enddo
