do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    ufx(i,j) = cff0 * zeta_new(i,j) + cff1 * zeta(i,j,kstp) + cff2 * zeta(i,j,kbak) + cff3 * zeta(i,j,kold)
    ufe(i,j) = (1. + rhos(i,j)) * ufx(i,j)
    vfe(i,j) = ufe(i,j) * ufx(i,j)
    vfx(i,j) = ufx(i,j) * (rhos(i,j) - rhoa(i,j))
  enddo
enddo
