do k = 1, n, 1
  do j = jstrv - 1, jend, 1
    do i = istru - 1, iend, 1
      cff = 0.5 * hz(i,j,k) * fomn(i,j)
      ufx_3d(i,j,k) = cff * (v(i,j,k,nrhs) + v(i,j + 1,k,nrhs))
      vfe_3d(i,j,k) = cff * (u(i,j,k,nrhs) + u(i + 1,j,k,nrhs))
    enddo
  enddo
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ru(i,j,k) = ru(i,j,k) + 0.5 * (ufx_3d(i,j,k) + ufx_3d(i - 1,j,k))
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      rv(i,j,k) = rv(i,j,k) - 0.5 * (vfe_3d(i,j,k) + vfe_3d(i,j - 1,k))
    enddo
  enddo
  do j = jstr, jend, 1
    do i = MAX(istru - 1, 2), MIN(iend + 1, lm), 1
      wrk1_3d(i,j,k) = u(i - 1,j,k,nrhs) - 2. * u(i,j,k,nrhs) + u(i + 1,j,k,nrhs)
      wrk2_3d(i,j,k) = huon(i - 1,j,k) - 2. * huon(i,j,k) + huon(i + 1,j,k)
    enddo
  enddo
  if (istr == 1) then
    do j = jstr, jend, 1
      wrk1_3d(1,j,k) = wrk1_3d(2,j,k)
      wrk2_3d(1,j,k) = wrk2_3d(2,j,k)
    enddo
  end if
  if (iend == lm) then
    do j = jstr, jend, 1
      wrk1_3d(lm + 1,j,k) = wrk1_3d(lm,j,k)
      wrk2_3d(lm + 1,j,k) = wrk2_3d(lm,j,k)
    enddo
  end if
  do j = jstr, jend, 1
    do i = istru - 1, iend, 1
      cffx = u(i,j,k,nrhs) + u(i + 1,j,k,nrhs)
      if (cffx > 0.) then
        curvx = wrk1_3d(i,j,k)
      else
        curvx = wrk1_3d(i + 1,j,k)
      end if
      ufx_3d(i,j,k) = 0.25 * (cffx + gamma * curvx) * (huon(i,j,k) + huon(i + 1,j,k) - 0.125 * (wrk2_3d(i,j,k) + wrk2_3d(i + 1,j,k)))
    enddo
  enddo
  do j = MAX(jstrv - 1, 2), MIN(jend + 1, mm), 1
    do i = istr, iend, 1
      wrk1_3d(i,j,k) = v(i,j - 1,k,nrhs) - 2. * v(i,j,k,nrhs) + v(i,j + 1,k,nrhs)
      wrk2_3d(i,j,k) = hvom(i,j - 1,k) - 2. * hvom(i,j,k) + hvom(i,j + 1,k)
    enddo
  enddo
  if (jstr == 1) then
    do i = istr, iend, 1
      wrk1_3d(i,1,k) = wrk1_3d(i,2,k)
      wrk2_3d(i,1,k) = wrk2_3d(i,2,k)
    enddo
  end if
  if (jend == mm) then
    do i = istr, iend, 1
      wrk1_3d(i,mm + 1,k) = wrk1_3d(i,mm,k)
      wrk2_3d(i,mm + 1,k) = wrk2_3d(i,mm,k)
    enddo
  end if
  do j = jstrv - 1, jend, 1
    do i = istr, iend, 1
      cffe = v(i,j,k,nrhs) + v(i,j + 1,k,nrhs)
      if (cffe > 0.) then
        curve = wrk1_3d(i,j,k)
      else
        curve = wrk1_3d(i,j + 1,k)
      end if
      vfe_3d(i,j,k) = 0.25 * (cffe + gamma * curve) * (hvom(i,j,k) + hvom(i,j + 1,k) - 0.125 * (wrk2_3d(i,j,k) + wrk2_3d(i,j + 1,k)))
    enddo
  enddo
  do j = MAX(jstr - 1, 1), MIN(jend + 1, mm), 1
    do i = istru, iend, 1
      wrk1_3d(i,j,k) = u(i,j - 1,k,nrhs) - u(i,j,k,nrhs) + (u(i,j + 1,k,nrhs) - u(i,j,k,nrhs))
    enddo
  enddo
  if (jstr == 1) then
    do i = istru, iend, 1
      wrk1_3d(i,0,k) = wrk1_3d(i,1,k)
    enddo
  end if
  if (jend == mm) then
    do i = istru, iend, 1
      wrk1_3d(i,mm + 1,k) = wrk1_3d(i,mm,k)
    enddo
  end if
  do j = jstr, jend + 1, 1
    do i = istru - 1, iend, 1
      wrk2_3d(i,j,k) = hvom(i - 1,j,k) - 2. * hvom(i,j,k) + hvom(i + 1,j,k)
    enddo
  enddo
  do j = jstr, jend + 1, 1
    do i = istru, iend, 1
      cffx = u(i,j,k,nrhs) + u(i,j - 1,k,nrhs)
      cffe = hvom(i,j,k) + hvom(i - 1,j,k)
      if (cffe > 0.) then
        curvx = wrk1_3d(i,j - 1,k)
      else
        curvx = wrk1_3d(i,j,k)
      end if
      ufe_3d(i,j,k) = 0.25 * (cffx + gamma * curvx) * (cffe - 0.125 * (wrk2_3d(i,j,k) + wrk2_3d(i - 1,j,k)))
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = MAX(istr - 1, 1), MIN(iend + 1, lm), 1
      wrk1_3d(i,j,k) = v(i - 1,j,k,nrhs) - v(i,j,k,nrhs) + (v(i + 1,j,k,nrhs) - v(i,j,k,nrhs))
    enddo
  enddo
  if (istr == 1) then
    do j = jstrv, jend, 1
      wrk1_3d(0,j,k) = wrk1_3d(1,j,k)
    enddo
  end if
  if (iend == lm) then
    do j = jstrv, jend, 1
      wrk1_3d(lm + 1,j,k) = wrk1_3d(lm,j,k)
    enddo
  end if
  do j = jstrv - 1, jend, 1
    do i = istr, iend + 1, 1
      wrk2_3d(i,j,k) = huon(i,j - 1,k) - 2. * huon(i,j,k) + huon(i,j + 1,k)
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = istr, iend + 1, 1
      cffe = v(i,j,k,nrhs) + v(i - 1,j,k,nrhs)
      cffx = huon(i,j,k) + huon(i,j - 1,k)
      if (cffx > 0.) then
        curve = wrk1_3d(i - 1,j,k)
      else
        curve = wrk1_3d(i,j,k)
      end if
      vfx_3d(i,j,k) = 0.25 * (cffe + gamma * curve) * (cffx - 0.125 * (wrk2_3d(i,j,k) + wrk2_3d(i,j - 1,k)))
    enddo
  enddo
  do j = jstr, jend, 1
    do i = istru, iend, 1
      ru(i,j,k) = ru(i,j,k) - ufx_3d(i,j,k) + ufx_3d(i - 1,j,k) - ufe_3d(i,j + 1,k) + ufe_3d(i,j,k)
    enddo
  enddo
  do j = jstrv, jend, 1
    do i = istr, iend, 1
      rv(i,j,k) = rv(i,j,k) - vfx_3d(i + 1,j,k) + vfx_3d(i,j,k) - vfe_3d(i,j,k) + vfe_3d(i,j - 1,k)
    enddo
  enddo
enddo
