do k = n, levsfrc, -1
  do j = jstrv - 1, jend, 1
    do i = istru - 1, iend, 1
      ufe(i,j) = ufe(i,j) + hz(i,j,k)
    enddo
  enddo
enddo
