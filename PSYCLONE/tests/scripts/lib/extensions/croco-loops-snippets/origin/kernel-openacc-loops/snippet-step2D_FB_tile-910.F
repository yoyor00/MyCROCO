do j = jstr, jend, 1
  do i = istr - 1, iend, 1
    ufx(i,j) = 0.25 * (duon(i,j) + duon(i + 1,j)) * (urhs(i,j) + urhs(i + 1,j))
    vfx(i + 1,j) = 0.25 * (duon(i + 1,j) + duon(i + 1,j - 1)) * (vrhs(i + 1,j) + vrhs(i,j))
  enddo
enddo
