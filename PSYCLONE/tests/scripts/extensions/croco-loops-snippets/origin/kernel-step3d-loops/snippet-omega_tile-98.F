do k = 1, n, 1
  do j = jstr, jend, 1
    do i = istr, iend, 1
      romega(i,j,k) = -huon(i + 1,j,k) + huon(i,j,k) - hvom(i,j + 1,k) + hvom(i,j,k) - (hz(i,j,k) - hz_bak(i,j,k)) * cf(i,j)
    enddo
  enddo
enddo
