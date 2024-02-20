do j = jstr, jend, 1
  do i = istru, iend, 1
    rufrc(i,j) = ru(i,j,1)
  enddo
  do k = 2, n, 1
    do i = istru, iend, 1
      rufrc(i,j) = rufrc(i,j) + ru(i,j,k)
    enddo
  enddo
  if (j >= jstrv) then
    do i = istr, iend, 1
      rvfrc(i,j) = rv(i,j,1)
    enddo
    do k = 2, n, 1
      do i = istr, iend, 1
        rvfrc(i,j) = rvfrc(i,j) + rv(i,j,k)
      enddo
    enddo
  end if
enddo
