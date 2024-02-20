do k = levsfrc, n, 1
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ru(i,j,k) = ru(i,j,k) + ufx(i,j) * (hz(i,j,k) + hz(i - 1,j,k))
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      rv(i,j,k) = rv(i,j,k) + vfe(i,j) * (hz(i,j,k) + hz(i,j - 1,k))
    enddo
  enddo
enddo
