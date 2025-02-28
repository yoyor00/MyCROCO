       if (mynode==0) write(6,*) "CAUTION: Mask is modified."
        do i=-2,Lm+3+padd_X
        do j=-2,Mm+3+padd_E
        if (  ( xr(i,j).ge.-0.015.and.xr(i,j).le.0.815
     &    .and.yr(i,j).ge.-1.405.and.yr(i,j).le.-1.355)
     &    .or. ( xr(i,j).ge.-0.085.and.xr(i,j).le.0.205
     &    .and.yr(i,j).ge.2.945.and.yr(i,j).le.3.205)
     &    .or. ( xr(i,j).ge.0.105.and.xr(i,j).le.0.195
     &    .and.yr(i,j).ge.3.195.and.yr(i,j).le.6.285)
     &     ) then
         rmask(i,j)=0
        endif
        enddo
        enddo
