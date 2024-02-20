do itrc = 1, nt, 1
  do k = 1, n, 1
    do j = jstr, jend, 1
      do i = MAX(istr - 1, 1), MIN(iend + 2, lm + 1), 1
        fx(i,j) = t(i,j,k,nadv,itrc) - t(i - 1,j,k,nadv,itrc)
      enddo
    enddo
    if (istr == 1) then
      do j = jstr, jend, 1
        fx(0,j) = fx(1,j)
      enddo
    end if
    if (iend == lm) then
      do j = jstr, jend, 1
        fx(lm + 2,j) = fx(lm + 1,j)
      enddo
    end if
    do j = jstr, jend, 1
      do i = istr - 1, iend + 1, 1
        work(i,j) = 0.5 * (fx(i + 1,j) + fx(i,j))
      enddo
    enddo
    do j = jstr, jend, 1
      do i = istr, iend + 1, 1
        fx(i,j) = 0.5 * (t(i,j,k,nadv,itrc) + t(i - 1,j,k,nadv,itrc) - 0.333333333333 * (work(i,j) - work(i - 1,j))) * huon(i,j,k)
      enddo
    enddo
    do j = MAX(jstr - 1, 1), MIN(jend + 2, mm + 1), 1
      do i = istr, iend, 1
        fe(i,j) = t(i,j,k,nadv,itrc) - t(i,j - 1,k,nadv,itrc)
      enddo
    enddo
    if (jstr == 1) then
      do i = istr, iend, 1
        fe(i,0) = fe(i,1)
      enddo
    end if
    if (jend == mm) then
      do i = istr, iend, 1
        fe(i,mm + 2) = fe(i,mm + 1)
      enddo
    end if
    do j = jstr - 1, jend + 1, 1
      do i = istr, iend, 1
        work(i,j) = 0.5 * (fe(i,j + 1) + fe(i,j))
      enddo
    enddo
    do j = jstr, jend + 1, 1
      do i = istr, iend, 1
        fe(i,j) = 0.5 * (t(i,j,k,nadv,itrc) + t(i,j - 1,k,nadv,itrc) - 0.333333333333 * (work(i,j) - work(i,j - 1))) * hvom(i,j,k)
      enddo
    enddo
    if (iic == ntstart) then
      cff = 0.5 * dt
      do j = jstr, jend, 1
        do i = istr, iend, 1
          t(i,j,k,nnew,itrc) = hz(i,j,k) * t(i,j,k,nstp,itrc) - cff * pm(i,j) * pn(i,j) * (fx(i + 1,j) - fx(i,j) + fe(i,j + 1) - fe(i,j))
        enddo
      enddo
    else
      cff = (1. - gamma) * dt
      cff1 = 0.5 + gamma
      cff2 = 0.5 - gamma
      do j = jstr, jend, 1
        do i = istr, iend, 1
          t(i,j,k,nnew,itrc) = cff1 * hz(i,j,k) * t(i,j,k,nstp,itrc) + cff2 * hz_bak(i,j,k) * t(i,j,k,indx,itrc) - cff * pm(i,j) * pn(i,j) * (fx(i + 1,j) - fx(i,j) + fe(i,j + 1) - fe(i,j))
        enddo
      enddo
    end if
  enddo
enddo
