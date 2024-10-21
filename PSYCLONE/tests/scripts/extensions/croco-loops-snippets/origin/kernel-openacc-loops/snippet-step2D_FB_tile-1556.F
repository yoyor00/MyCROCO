do j = jstr, jend, 1
  do i = istru, iend, 1
    cff = rufrc(i,j) - rubar(i,j)
    rufrc(i,j) = cff1 * cff + cff2 * rufrc_bak(i,j,3 - nstp) + cff3 * rufrc_bak(i,j,nstp)
    rufrc_bak(i,j,nstp) = cff
  enddo
enddo
