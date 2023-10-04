do k = 1, n, 1
  do j = jstrv - 1, jend, 1
    do i = istru - 1, iend, 1
      hz_half(i,j,k) = cff1 * hz(i,j,k) + cff2 * hz_bak(i,j,k) - cff * pm(i,j) * pn(i,j) * (huon(i + 1,j,k) - huon(i,j,k) + hvom(i,j + 1,k) - hvom(i,j,k) + we(i,j,k) - we(i,j,k - 1))
    enddo
  enddo
enddo
