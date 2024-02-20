do j = jstrr, jendr, 1
  du_avg1(istru - 1,j,nnew) = du_avg1(istru - 1,j,nnew) + cff1 * (dnew(istru - 1,j) + dnew(istru - 2,j)) * ubar(istru - 1,j,knew) * on_u(istru - 1,j)
enddo
