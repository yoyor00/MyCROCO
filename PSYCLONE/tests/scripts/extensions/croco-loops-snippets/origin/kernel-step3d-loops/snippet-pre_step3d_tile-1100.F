do k = 1, n, 1
  do j = jstr, jend, 1
    do i = istru, iend, 1
      cff = 2. / (hz_half(i,j,k) + hz_half(i - 1,j,k))
      u(i,j,k,nnew) = (u(i,j,k,nstp) * 0.5 * (hz(i,j,k) + hz(i - 1,j,k)) + cdt * work(i,j) * ru(i,j,k)) * cff
      u(i,j,k,indx) = u(i,j,k,nstp) * 0.5 * (hz(i,j,k) + hz(i - 1,j,k))
    enddo
  enddo
enddo
