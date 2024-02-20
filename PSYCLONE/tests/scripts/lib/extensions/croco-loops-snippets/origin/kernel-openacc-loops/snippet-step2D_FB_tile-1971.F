do j = jstr, jend, 1
  do i = istru, iend, 1
    dunew = (duon(i,j) + duon(i - 1,j)) * ubar(i,j,kstp) + cff * (pm(i,j) + pm(i - 1,j)) * (pn(i,j) + pn(i - 1,j)) * (rubar(i,j) + rufrc(i,j))
    ubar(i,j,knew) = dunew / (dnew(i,j) + dnew(i - 1,j))
    du_avg1(i,j,nnew) = du_avg1(i,j,nnew) + cff1 * on_u(i,j) * dunew
  enddo
enddo
