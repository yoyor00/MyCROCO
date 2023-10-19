do j = jstr, jendr, 1
  do i = istrr, iendr, 1
    dc(i,0) = 0.
    cf(i,0) = 0.
    fc(i,0) = 0.
  enddo
  do k = 1, n, +1
    do i = istrr, iendr, 1
      dc(i,k) = 0.5 * (hz(i,j,k) + hz(i,j - 1,k)) * om_v(i,j)
      dc(i,0) = dc(i,0) + dc(i,k)
      cf(i,0) = cf(i,0) + dc(i,k) * v(i,j,k,nnew)
    enddo
  enddo
  do i = istrr, iendr, 1
    dc(i,0) = 1. / dc(i,0)
    cf(i,0) = dc(i,0) * (cf(i,0) - dv_avg1(i,j,nnew))
    vbar(i,j,knew) = dc(i,0) * dv_avg1(i,j,nnew)
  enddo
  do k = n, 1, -1
    do i = istrr, iendr, 1
      v(i,j,k,nnew) = v(i,j,k,nnew) - cf(i,0)
      fc(i,k) = 0.75 * hvom(i,j,k) + 0.125 * dc(i,k) * (v(i,j,k,nstp) + v(i,j,k,nnew))
      fc(i,0) = fc(i,0) + fc(i,k)
    enddo
  enddo
  do i = istrr, iendr, 1
    fc(i,0) = dc(i,0) * (fc(i,0) - dv_avg2(i,j))
  enddo
  do k = 1, n, +1
    do i = istrr, iendr, 1
      hvom(i,j,k) = fc(i,k) - dc(i,k) * fc(i,0)
    enddo
  enddo
enddo
