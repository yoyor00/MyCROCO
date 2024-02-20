do j = jstrv, jend, 1
  do i = istr, iend, 1
    cff = rvfrc(i,j) - rvbar(i,j)
    rvfrc(i,j) = cff1 * cff + cff2 * rvfrc_bak(i,j,3 - nstp) + cff3 * rvfrc_bak(i,j,nstp)
    rvfrc_bak(i,j,nstp) = cff
  enddo
enddo
