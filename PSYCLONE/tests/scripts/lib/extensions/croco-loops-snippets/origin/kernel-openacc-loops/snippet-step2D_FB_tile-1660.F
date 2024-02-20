do j = jstrv - 1, jend, 1
  do i = istru - 1, iend, 1
    ufx(i,j) = zeta_new(i,j) - zeta(i,j,kstp)
    ufe(i,j) = (1. + rhos(i,j)) * ufx(i,j)
    vfe(i,j) = ufe(i,j) * (zeta_new(i,j) + zeta(i,j,kstp))
    vfx(i,j) = ufx(i,j) * (rhos(i,j) - rhoa(i,j))
  enddo
enddo
