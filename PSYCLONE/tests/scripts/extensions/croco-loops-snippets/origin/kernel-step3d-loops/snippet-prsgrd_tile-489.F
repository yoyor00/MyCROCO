do k = n, 1, -1
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ru(i,j,k) = 0.5 * (hz(i,j,k) + hz(i - 1,j,k)) * on_u(i,j) * (p(i - 1,j,k) - p(i,j,k))
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      rv(i,j,k) = 0.5 * (hz(i,j,k) + hz(i,j - 1,k)) * om_v(i,j) * (p(i,j - 1,k) - p(i,j,k))
    enddo
  enddo
enddo
