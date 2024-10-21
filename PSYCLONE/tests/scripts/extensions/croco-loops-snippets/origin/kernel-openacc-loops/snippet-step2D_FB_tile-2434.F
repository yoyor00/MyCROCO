do j = jstrr, jendr, 1
  du_avg1(iend + 1,j,nnew) = du_avg1(iend + 1,j,nnew) + cff1 * (dnew(iend + 1,j) + dnew(iend,j)) * ubar(iend + 1,j,knew) * on_u(iend + 1,j)
enddo
