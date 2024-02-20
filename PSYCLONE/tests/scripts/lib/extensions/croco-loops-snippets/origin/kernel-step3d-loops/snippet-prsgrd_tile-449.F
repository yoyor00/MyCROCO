do k = n - 1, 1, -1
  do j = jstrv - 1, jend, 1
    do i = istru - 1, iend, 1
      p(i,j,k) = p(i,j,k) + p(i,j,k + 1)
    enddo
  enddo
enddo
