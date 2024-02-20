do j = jstrv, jend, 1
  do i = istr, iend, 1
    vfe(i,j) = svstr(i,j) * 4. / ((ufe(i,j - 1) + ufe(i,j)) * (pm(i,j - 1) + pm(i,j)) * (pn(i,j - 1) + pn(i,j)))
  enddo
enddo
