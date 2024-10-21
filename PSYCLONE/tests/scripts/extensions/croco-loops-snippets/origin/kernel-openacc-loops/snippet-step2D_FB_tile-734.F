do j = jstr, jend, 1
  do i = istr, iend, 1
    rubar(i,j) = cff * on_u(i,j) * ((h(i - 1,j) + h(i,j)) * (ufe(i - 1,j) - ufe(i,j)) + vfe(i - 1,j) - vfe(i,j) + (h(i - 1,j) - h(i,j)) * (vfx(i - 1,j) + vfx(i,j) + 0.333333333333 * (rhoa(i - 1,j) - rhoa(i,j)) * (ufx(i - 1,j) - ufx(i,j))))
    rvbar(i,j) = cff * om_v(i,j) * ((h(i,j - 1) + h(i,j)) * (ufe(i,j - 1) - ufe(i,j)) + vfe(i,j - 1) - vfe(i,j) + (h(i,j - 1) - h(i,j)) * (vfx(i,j - 1) + vfx(i,j) + 0.333333333333 * (rhoa(i,j - 1) - rhoa(i,j)) * (ufx(i,j - 1) - ufx(i,j))))
  enddo
enddo
