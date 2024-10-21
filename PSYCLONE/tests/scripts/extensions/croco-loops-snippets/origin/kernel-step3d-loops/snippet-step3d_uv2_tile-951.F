do j = jstrr, jendr, 1
  do i = istr, iendr, 1
    dc(i,0) = 0.
    cf(i,0) = 0.
    fc(i,0) = 0.
  enddo
  do k = 1, n, +1
    do i = istr, iendr, 1
      dc(i,k) = 0.5 * (hz(i,j,k) + hz(i - 1,j,k)) * on_u(i,j)
      dc(i,0) = dc(i,0) + dc(i,k)
      cf(i,0) = cf(i,0) + dc(i,k) * u(i,j,k,nnew)
    enddo
  enddo
  do i = istr, iendr, 1
    dc(i,0) = 1. / dc(i,0)
    cf(i,0) = dc(i,0) * (cf(i,0) - du_avg1(i,j,nnew))
    ubar(i,j,knew) = dc(i,0) * du_avg1(i,j,nnew)
  enddo
  do k = n, 1, -1
    do i = istr, iendr, 1
      u(i,j,k,nnew) = u(i,j,k,nnew) - cf(i,0)
      fc(i,k) = 0.75 * huon(i,j,k) + 0.125 * dc(i,k) * (u(i,j,k,nstp) + u(i,j,k,nnew))
      fc(i,0) = fc(i,0) + fc(i,k)
    enddo
  enddo
  do i = istr, iendr, 1
    fc(i,0) = dc(i,0) * (fc(i,0) - du_avg2(i,j))
  enddo
  do k = 1, n, +1
    do i = istr, iendr, 1
      huon(i,j,k) = fc(i,k) - dc(i,k) * fc(i,0)
    enddo
  enddo
enddo
