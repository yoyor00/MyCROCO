do k = 1, n, 1
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      cff = 2. / (hz_half(i,j,k) + hz_half(i,j - 1,k))
      v(i,j,k,nnew) = (cff1 * v(i,j,k,nstp) * 0.5 * (hz(i,j,k) + hz(i,j - 1,k)) + cff2 * v(i,j,k,indx) * 0.5 * (hz_bak(i,j,k) + hz_bak(i,j - 1,k)) + cdt * work(i,j) * rv(i,j,k)) * cff
      v(i,j,k,indx) = v(i,j,k,nstp) * 0.5 * (hz(i,j,k) + hz(i,j - 1,k))
    enddo
  enddo
enddo
