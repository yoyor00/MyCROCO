do i = istru, iend, 1
  du_avg1(i,jstr - 1,nnew) = du_avg1(i,jstr - 1,nnew) + cff1 * (dnew(i,jstr - 1) + dnew(i - 1,jstr - 1)) * ubar(i,jstr - 1,knew) * on_u(i,jstr - 1)
enddo
