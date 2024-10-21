do j = jstrv - 1, jend, 1
  do k = 1, n - 1, 1
    do i = istru - 1, iend, 1
      dz(i,k) = z_r(i,j,k + 1) - z_r(i,j,k)
      dr(i,k) = rho(i,j,k + 1) - rho(i,j,k)
    enddo
  enddo
  do i = istru - 1, iend, 1
    dr(i,n) = dr(i,n - 1)
    dr(i,0) = dr(i,1)
    dz(i,n) = dz(i,n - 1)
    dz(i,0) = dz(i,1)
  enddo
  do k = n, 1, -1
    do i = istru - 1, iend, 1
      cff = 2. * dz(i,k) * dz(i,k - 1)
      dz(i,k) = cff / (dz(i,k) + dz(i,k - 1))
      cfr = 2. * dr(i,k) * dr(i,k - 1)
      if (cfr > epsil) then
        dr(i,k) = cfr / (dr(i,k) + dr(i,k - 1))
      else
        dr(i,k) = 0.
      end if
    enddo
  enddo
  do i = istru - 1, iend, 1
    p(i,j,n) = g * z_w(i,j,n) + grho * (rho(i,j,n) + 0.5 * (rho(i,j,n) - rho(i,j,n - 1)) * (z_w(i,j,n) - z_r(i,j,n)) / (z_r(i,j,n) - z_r(i,j,n - 1))) * (z_w(i,j,n) - z_r(i,j,n))
  enddo
  do k = n - 1, 1, -1
    do i = istru - 1, iend, 1
      p(i,j,k) = halfgrho * ((rho(i,j,k + 1) + rho(i,j,k)) * (z_r(i,j,k + 1) - z_r(i,j,k)) - onefifth * ((dr(i,k + 1) - dr(i,k)) * (z_r(i,j,k + 1) - z_r(i,j,k) - onetwelfth * (dz(i,k + 1) + dz(i,k))) - (dz(i,k + 1) - dz(i,k)) * (rho(i,j,k + 1) - rho(i,j,k) - onetwelfth * (dr(i,k + 1) + dr(i,k)))))
    enddo
  enddo
enddo
