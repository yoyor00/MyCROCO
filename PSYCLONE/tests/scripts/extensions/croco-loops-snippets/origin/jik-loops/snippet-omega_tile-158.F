do j = jstr, jend, 1
  do i = istr, iend, 1
    we(i,j,0) = 0.
    fc(i,0) = 0.
    dc(i,0) = 0.
  enddo
  do k = 1, n, +1
    do i = istr, iend, 1
      fc(i,k) = fc(i,k - 1) + romega(i,j,k)
      dc(i,k) = dc(i,k - 1) + hz(i,j,k)
    enddo
  enddo
  do i = istr, iend, 1
    dc(i,n) = fc(i,n) / dc(i,n)
  enddo
  do k = 1, n - 1, +1
    do i = istr, iend, 1
      we(i,j,k) = fc(i,k) - dc(i,n) * dc(i,k)
    enddo
  enddo
  do i = istr, iend, 1
    we(i,j,n) = 0.
  enddo
enddo
