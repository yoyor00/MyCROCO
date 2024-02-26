do j = jstr, jend, 1
  do i = istru, iend, 1
    ufx(i,j) = bustr(i,j) * 4. / ((ufe(i - 1,j) + ufe(i,j)) * (pm(i - 1,j) + pm(i,j)) * (pn(i - 1,j) + pn(i,j)))
  enddo
enddo
