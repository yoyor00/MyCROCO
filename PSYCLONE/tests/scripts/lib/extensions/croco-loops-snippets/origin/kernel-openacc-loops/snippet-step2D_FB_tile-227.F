do j = jstr - 1, jend + 1, 1
  do i = imin + 1, imax, 1
    urhs(i,j) = cff1 * ubar(i,j,kstp) + cff2 * ubar(i,j,kbak) + cff3 * ubar(i,j,kold)
    duon(i,j) = 0.5 * (drhs(i,j) + drhs(i - 1,j)) * on_u(i,j) * urhs(i,j)
  enddo
enddo
