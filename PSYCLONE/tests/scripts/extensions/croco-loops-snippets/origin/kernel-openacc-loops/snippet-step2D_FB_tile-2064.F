do j = jstrv, jend, 1
  do i = istr, iend, 1
    dvnew = (duon(i,j) + duon(i,j - 1)) * vbar(i,j,kstp) + cff * (pm(i,j) + pm(i,j - 1)) * (pn(i,j) + pn(i,j - 1)) * (rvbar(i,j) + rvfrc(i,j))
    vbar(i,j,knew) = dvnew / (dnew(i,j) + dnew(i,j - 1))
    dv_avg1(i,j,nnew) = dv_avg1(i,j,nnew) + cff1 * om_v(i,j) * dvnew
  enddo
enddo
