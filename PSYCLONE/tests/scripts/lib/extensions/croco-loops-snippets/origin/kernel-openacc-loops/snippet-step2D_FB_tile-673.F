do j = jstrr, jendr, 1
  do i = istrr, iendr, 1
    zt_avg1(i,j) = zt_avg1(i,j) + cff1 * zeta(i,j,knew)
    du_avg2(i,j) = du_avg2(i,j) + cff2 * duon(i,j)
    dv_avg2(i,j) = dv_avg2(i,j) + cff2 * dvom(i,j)
  enddo
enddo
