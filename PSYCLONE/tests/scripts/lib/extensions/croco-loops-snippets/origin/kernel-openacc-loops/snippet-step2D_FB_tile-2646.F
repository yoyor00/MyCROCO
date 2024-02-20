do i = istru, iend, 1
  du_avg1(i,jend + 1,nnew) = du_avg1(i,jend + 1,nnew) + cff1 * (dnew(i,jend + 1) + dnew(i - 1,jend + 1)) * ubar(i,jend + 1,knew) * on_u(i,jend + 1)
enddo
