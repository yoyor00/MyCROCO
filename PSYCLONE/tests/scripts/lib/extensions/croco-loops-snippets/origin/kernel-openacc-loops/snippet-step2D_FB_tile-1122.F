do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    cff = drhs(i,j) * fomn(i,j)
    ufx(i,j) = cff * (vrhs(i,j) + vrhs(i,j + 1))
    vfe(i,j) = cff * (urhs(i,j) + urhs(i + 1,j))
  enddo
enddo
