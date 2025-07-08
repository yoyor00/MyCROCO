       water_lev=0.000762   ! Water level difference between
                             ! Atlantic and Mediterranean (m)
        Sref=1003.3-1000.-R0 ! EXP52 ! 1004.9-1000.-R0 ! sigma-t - R0 
        drho=19.             ! rho_med - rho_atl
                             ! e.g. 1024.5 - 1005.2 = 19.3 for EXP88
        ! ===========================================================================
        ! -------------------------- Linear seawater EOS: ---------------------------
        !  sigma-t = rho_{atl_or_med}-1000 = R0-Tcoef(T-T0)+Scoef(S_{atl_or_med}-S0)
        ! ---------------------------------------------------------------------------
        ! --------------------- With the previous parameters: -----------------------
        ! - [Assume there is no tempature and S0 = 35 PSU and Scoef = 0.78 PSU^{-1} -
        ! --------- (cf. Vallis - Atmospheric and Oceanic Fluid Dynamics)] ----------
        ! ---------------------------------------------------------------------------
        ! ----------- S_atl = (sigma-t       - R0)/Scoef + S0 = 3.2051  PSU ---------
        ! ----------- S_med = (sigma-t +drho - R0)/Scoef + S0 = 27.9487 PSU ---------
        ! ---------------------------------------------------------------------------
        ! --- Output netcdf will give sigma-t instead of rho. Thus, the actual rho --
        ! -------- is given by sigma-t + 1000., i.e. rho+1000. in the netcdf --------
        ! ===========================================================================
#    ifdef MPI
        HR_grid=.true.
        ! --- Define [x7, x1, x2, xgate]
        ! ------ and [y5, y2, y7, y1   ]:
        if (HR_grid) then
           ! --- HR:
           !x_vect = [0.9344, 4.41485, 9.425, 10.00]
           x_vect = [0.94, 4.41, 9.425, 10.00]
           !y_vect = [3.69507, 5.675, 6.58523, 13.005]
           y_vect = [3.7, 5.675, 6.565, 13.005]
        else
           ! --- BR:
           !x_vect = [0.850, 4.30, 9.50, 10.00]
           !y_vect = [3.775, 5.80, 7.00, 13.15]
           !slope=(y_vect(3)-y_vect(1))/
     !&            (x_vect(1)-x_vect(2))
        endif
        ! --- CS(0,0) ==> grid_cntrd = true
        grid_cntrd=.true.
        if (grid_cntrd) then
          x_vect=x_vect-9.21
          y_vect(1)=y_vect(1)-5.10
          y_vect(2:4)=y_vect(2:4)-4.80
        endif
        if (HR_grid) then
           slope=(y_vect(3)-y_vect(1))/
     &           (x_vect(1)-x_vect(2))
        endif
        ! --- Parameters for "linspace":
        num_val=1000
        start_val=x_vect(1) !0.60-9.21
        end_val=x_vect(2)   !4.45-9.21
        step = (end_val-start_val)/(num_val-1)
        ! --- adj parameter: adjust the
        ! --- straight line to the wall
        if (HR_grid) then
           adj=0. !-0.0083 !-0.009 !-0.023
        else
           if (grid_cntrd) then
              adj=-0.30 
           else
              adj=0.23
           endif
        endif
        ! =========================
        ! --- Update ==> Real gate:
        ! =========================
        brutal_LEx=.false. ! True:  straight door with depth --> brutal 
                           ! False: inclined door with depth --> smooth
        ! -----------------------------------------------------
        ! --- Note: real gate does not match with the edges ---
        ! ------- of the coast. Thus, I preliminary use -------
        ! ----------- Python to extrapolate the real ----------
        ! ------------ position before updating it ------------
        ! -----------------------------------------------------
        ! --- Real position:
        !x_gate = [ 0.92,  1.08] ! x-coordinates: [x_south;x_north]
        !y_gate = [-0.68, -0.19] ! y-coordinates: [y_south;y_north]
        ! --- Extrapolated position:
        x_gate = [ 0.88  ,  1.12  ]
        y_gate = [-0.8025, -0.0675]
        if (brutal_LEx) then
           ! --- Parameters for the door:
           slope_g=(y_gate(1)-y_gate(2))/
     &             (x_gate(1)-x_gate(2))
           num_g=1000
           start_g=x_gate(1)
           end_g=x_gate(2)
           step_g=(end_g-start_g)/(num_g-1)
        else
           z_gate = [0., 0.] ! z-coordinate of the previous points
                             ! Only useful for an inclined door.
           ! ====================================================
           ! ---- Definition of two direction vectors in the ----
           ! ------- plane of the door + a normal vector --------
           ! ====================================================
           ! --- The first vector is given using the position ---
           ! ----------- of the gate. Assuming that, ------------
           ! ---------- A = (x_gate(1), y_gate(1)) and ----------
           ! ---------- B = (x_gate(2), y_gate(2)) then, --------
           ! ----------- vec{u} = (x_gate(2)-x_gate(1), ---------
           ! --------------------- y_gate(2)-y_gate(1), 0) ------
           ! ----- Now, define a third point located at the -----
           ! --- bottom of the door and noted C = (Cx,Cy,-1). ---
           ! -- In order to obtain an inclined door, Cx and Cy --
           ! ------ coordinates must be different from the ------ 
           ! ---- coordinates of the door. Moreover, vec{AC} ----
           ! ----- is a direction vector if the dot product ----- 
           ! --- vec{AB} . vec{AC} .ne. 1. The simpliest way ----
           ! --- to determine vec{AC} coordinates, is to use ----
           ! - the definition of an orthogonal vector, so that --
           ! --------- vec{AB} . vec{AC} = 0. Doing this, ------- 
           ! --------- we find a second direction vector: -------
           ! ----------- vec{v} = (y_gate(1)-y_gate(2), ---------
           ! --------------------- x_gate(2)-x_gate(1),-1). -----
           ! ---- Finally, the normal vector to the plane of ----
           ! -- the door is given by the cross product between --
           ! --- vec{u} and vec{v}: vec{n} = vec{u} x vec{v} = --
           ! ---- (y_gate(1)-y_gate(2), x_gate(2)-x_gate(1), ---- 
           ! - (x_gate(2)-x_gate(1))^2-(y_gate(2)-y_gate(1))^2) -
           ! -- Thus, the equation of the plane of the door is -- 
           ! -- simply nx*x+ny*y+nz*z+d=0. As A belongs to the -- 
           ! - plane, d=-nx*Ax-ny*Ay, since Az = 0. It follows -- 
           ! -------- that the equation of the plane is: --------  
           !  [y_gate(1)-y_gate(2)]*x + [x_gate(2)-x_gate(1)]*y 
           !                         + [(x_gate(2)-x_gate(1))**2. 
           !                       - (y_gate(2)-y_gate(1))**2.]*z 
           !                    - [(y_gate(1)-y_gate(2))*x_gate(1)
           !                  + (x_gate(2)-x_gate(1))*y_gate(1)]=0
           ! ----------------------------------------------------
           ! ----- Important: note that to modify the door ------ 
           ! -- inclination, simply modify the coordinates of ---
           ! -- the C-point. For example, to increase (reduce) --
           ! --- the inclination, multiply (divide) Cx and Cy ---
           ! --------- by two: C_new = (Cx*2,Cy*2,-1). ----------
           ! ---- This requires to modify the vector vec{AC} ----
           ! --- vec{AC} = (Cx-Ax,Cy-Ay,-1). As the previous ----
           ! ---- Cx-Ax = y_gate(1)-y_gate(2), then the new ----- 
           ! ---- Cx-Ax = 2*y_gate(1)-2*y_gate(2)+Ax, for an ----
           ! ---------- increase of the inclination. ------------
           ! ====================================================
           ! ----- To play with the inclination of the door -----
           ! ====================================================
           incl = -10.0 !1. ! Coeff to in(de)crease the door inclination
           Cx = y_gate(1)-y_gate(2)+x_gate(1) ! x-coordinate of C
           Cy = x_gate(2)-x_gate(1)+y_gate(1) ! y-coordinate of C
           ! ====================================================
           ! - i_lck [i-coordinate of the 1st direction vector, -
           ! -------- i_coordinate of the 2nd direction vector, -
           ! -------- i_coordinate of the normal vector] --------
           ! ====================================================
           x_lck = [x_gate(2)-x_gate(1),y_gate(1)-y_gate(2),
     &                                 y_gate(1)-y_gate(2)]
           y_lck = [y_gate(2)-y_gate(1),x_gate(2)-x_gate(1),
     &                                 x_gate(2)-x_gate(1)]
           z_lck = [0., -1., (x_gate(2)-x_gate(1))*(incl*Cy-y_gate(1))-(
     &                        y_gate(2)-y_gate(1))*(incl*Cx-x_gate(1))]
        endif
#    endif

#    ifdef CENTRIFUGE
        do j=JR_RANGE
          do i=IR_RANGE
              zeta(i,j,:) = (f(i,j)/2.)**2*(ray(i,j)**2)/(2.0*g)
          enddo
        enddo
#    endif

#    ifndef LEVEL_DIFF
        flag_BG = 1.
#    endif

#    ifdef TEMPERATURE
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              if(xr(i,j).le.0.5*xl) then
                t(i,j,k,1,itemp)=22.0
              else
                t(i,j,k,1,itemp)=18.0
              endif
              t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            enddo
          enddo
        enddo
#    endif /* TEMPERATURE */
#    ifdef SALINITY
      ! --- MPI ==> Think local
#     ifdef MPI
        do k=1,N
          do j=JR_RANGE
            do i=IR_RANGE
              t(i,j,k,1,isalt)=Sref/Scoef + S0
              t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
              if ( yr(i,j) >= y_vect(4)   .or. ( yr(i,j) >= y_vect(2)  .and.
     &             xr(i,j) >= x_vect(3) ) .or.   yr(i,j) <= y_vect(1) ) then
     !&                                   .or.   xr(i,j) >= x_vect(4) )then 
     ! --- [Uncomment the above line if you don't want the tilted gate] ---
                t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                   flag_BG(i,j  ) = 1
                   zeta(i,j,1) = zeta(i,j,1)+water_lev
                   zeta(i,j,2) = zeta(i,j,1)
                endif
              endif
              ! --- Handle the tilted wall:
              allocate( x_dr(num_val), y_dr(num_val) )
              do l=1,num_val
                x_dr(l)=start_val+(l-1)*step
              enddo
              y_dr=slope*(x_dr-x_vect(1))+y_vect(3)
              do l=1,num_val
                if ( xr(i,j) <= x_dr(l) .and. 
     &               yr(i,j) <= y_dr(l)+adj ) then
                  t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                  t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                  if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                     flag_BG(i,j  ) = 1
                     zeta(i,j,1) = zeta(i,j,1)+water_lev
                     zeta(i,j,2) = zeta(i,j,1)
                  endif
                endif
              enddo
              deallocate( x_dr, y_dr )
              ! ----------------------------------------------
              ! --- Handle the bays of mediterranean water ---
              ! ----------------------------------------------
              ! --- First bay:
              if ( xr(i,j) > -0.305 .and. yr(i,j) < -1.175 ) then
                 t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                 t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                 if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                    flag_BG(i,j  ) = 1
                    zeta(i,j,1) = zeta(i,j,1)+water_lev
                    zeta(i,j,2) = zeta(i,j,1)
                 endif
              endif
              ! --- Second bay:
              if ( xr(i,j) > 0.165 .and. yr(i,j) < 3.075  .and.
     &             yr(i,j) > 2.985                       ) then
                 t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                 t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                 if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                    flag_BG(i,j  ) = 1
                    zeta(i,j,1) = zeta(i,j,1)+water_lev
                    zeta(i,j,2) = zeta(i,j,1)
                 endif
              endif
              ! --- Handle the tilted gate:
              ! -------------------------------------------
              ! --- [Comment the following lines if you ---
              ! ------- don't want the tilted gate] -------
              ! -------------------------------------------
              ! --- Straight door: brutal LEx:
              if (brutal_LEx) then
                 allocate( x_drg(num_g), y_drg(num_g) )
                 do l=1,num_g
                    x_drg(l)=start_g+(l-1)*step_g
                 enddo
                 if ( xr(i,j) >= end_g ) then
                   t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                   t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                   if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                      flag_BG(i,j  ) = 1
                      zeta(i,j,1) = zeta(i,j,1)+water_lev
                      zeta(i,j,2) = zeta(i,j,1)
                   endif
                 endif
                 y_drg=slope_g*(x_drg-x_gate(1))+y_gate(1)
                 do l=1,num_g
                   if ( xr(i,j) >= x_drg(l) .and.
     &                  yr(i,j) <= y_drg(l) ) then
                      t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                      t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                      if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                         flag_BG(i,j  ) = 1
                         zeta(i,j,1) = zeta(i,j,1)+water_lev
                         zeta(i,j,2) = zeta(i,j,1)
                      endif
                   endif
                 enddo
                 deallocate( x_drg, y_drg )
              ! ===========================================
              ! --- Inclined door with depth: smooth LEx --
              ! ===========================================
              ! --- Define the distance from a point of ---
              ! --- the door to the plane of the door: ----
              ! ---- If A = (x_gate(1),y_gate(1)) then ----
              ! ------------ the distance is: -------------
              ! --- d_lck = |vec{n} . vec{MA}| / ||vec{n}||
              ! -------------------------------------------
              ! ===========================================
              else
                 if ( xr(i,j) >= x_gate(2) ) then
                    t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                    t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                    if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                       flag_BG(i,j  ) = 1
                       zeta(i,j,1) = zeta(i,j,1)+water_lev
                       zeta(i,j,2) = zeta(i,j,1)
                    endif
                 endif
                 d_lck=( (x_gate(1)-xr(i,j))   *x_lck(3)
     &                  +(y_gate(1)-yr(i,j))   *y_lck(3)
     &                  +(z_gate(1)-z_r(i,j,k))*z_lck(3) )
     &              /(x_lck(3)**2.+y_lck(3)**2.+z_lck(3)**2.)
                 ! --- Box to "safely" change the salinity 
                 ! --- only in the strait:
                 if ( xr(i,j) .ge. x_gate(1)-1.8 .and. 
     &                xr(i,j) .le. x_gate(2)+0.8 .and.
     &                yr(i,j) .ge. y_gate(1)-0.5 .and.
     &                yr(i,j) .le. y_gate(2)+0.5 ) then
                    ! --- Based on the distance and on the mask, 
                    ! --- change the salinity of mediterranean water:
                    if (d_lck .le. 0. .and. d_lck .ge. -0.05 .and. 
     &                                           rmask(i,j) .ne. 0) then
                       t(i,j,k,1,isalt)=(Sref+drho*(1.+20.*d_lck))/Scoef
     &                                  + S0
                       t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                    elseif (d_lck .gt. 0. .and. rmask(i,j) .ne. 0) then
                       t(i,j,k,1,isalt)=(Sref+drho)/Scoef + S0
                       t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
                       if ((k .eq. N) .and. (flag_BG(i,j).eq.0)) then
                          flag_BG(i,j  ) = 1
                          zeta(i,j,1) = zeta(i,j,1)+water_lev
                          zeta(i,j,2) = zeta(i,j,1)
                       endif
                    endif
                 endif
              endif
            enddo
          enddo
        enddo
      ! --- No MPI ==> Vectorise
#     else
        t(:,:,:,1,isalt)=(Sref+drho)/Scoef + S0 
        t(:,:,:,2,isalt)=t(:,:,:,1,isalt)
        ! ============================
        ! --- Implement lock exchange:
        ! ============================
        allocate( i_inter(1), j_inter(1) )
        ! --- Define [x7, x1, x2, xgate]
        ! ------ and [y5, y2, y7, y1   ]:
        x_vect = [0.850, 4.30, 9.50, 10.00]
        y_vect = [3.775, 5.80, 7.00, 13.05]
        ! --- Construct i_vect = [i1, i2, i7, ig] 
        ! ----------and j_vect = [j1, j2, j5, j7]:
        do i=1,4
          i_inter=minloc(xr(:,1), mask=xr(:,1)>=x_vect(i))
          i_vect(i)=i_inter(1)
          j_inter=minloc(yr(1,:), mask=yr(1,:)>=y_vect(i))
          j_vect(i)=j_inter(1)
        enddo
        ! --- End of construction ==> release memory:
        deallocate( i_inter, j_inter )
        ! --- Create a straight line along the wall:
        slope=(y_vect(1)-y_vect(3))/(x_vect(2)-x_vect(1))
        lth=size(xr(1,i_vect(1):i_vect(2)))
        allocate( x_dr(lth), y_dr(lth) )
        x_dr(1)=x_vect(1)
        do i=1,lth-2
          x_dr(i+1)=x_dr(i)+(x_vect(2)-x_vect(1))/lth
        enddo
        x_dr(lth)=x_vect(2)
        y_dr=slope*(x_dr-x_vect(1))+y_vect(3)
        ! --- Water distribution:
        t(i_vect(2):i_vect(3),j_vect(1):j_vect(4),:,1,isalt) = 
     &  Sref/Scoef + S0
        t(i_vect(2):i_vect(3),j_vect(1):j_vect(4),:,2,isalt) = 
     &  t(i_vect(2):i_vect(3),j_vect(1):j_vect(4),:,1,isalt)
        ! ---
        t(i_vect(3):i_vect(4),j_vect(1):j_vect(2),:,1,isalt) = 
     &  Sref/Scoef + S0
        t(i_vect(3):i_vect(4),j_vect(1):j_vect(2),:,2,isalt) =
     &  t(i_vect(3):i_vect(4),j_vect(1):j_vect(2),:,1,isalt)
        ! ---
        t(1:i_vect(1),:,:,1,isalt) = Sref/Scoef + S0
        t(1:i_vect(1),:,:,2,isalt) =
     &  t(1:i_vect(1),:,:,1,isalt)
        ! ---
        allocate( id(1) )
        do i=1,size(x_dr)
          id=minloc(yr(1,:), mask=yr(1,:)>=y_dr(i))
          t(i_vect(1)+i,id(1)-1:size(yr(1,:)),:,1,isalt) = 
     &    Sref/Scoef + S0
          t(i_vect(1)+i,id(1)-1:size(yr(1,:)),:,2,isalt) =
     &    t(i_vect(1)+i,id(1)-1:size(yr(1,:)),:,1,isalt)
        enddo
        ! --- End of distribution ==> release memory:
        deallocate( x_dr, y_dr, id )
#     endif
#    endif
        
#    ifdef BATHY_SLOPE       
      myslope2=0.
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef MASKING
          if (umask(i,j).gt.0.) then
# endif
             myslope2(i,j)=sqrt(( (zeta(i,j,1)-zeta(i-1,j,1))
     &    *pm(i,j))**2)
# ifdef MASKING
          endif
# endif
        enddo
      enddo

      do j=JstrV,Jend
        do i=Istr,Iend
# ifdef MASKING
          if (vmask(i,j).gt.0.) then
# endif
             myslope2(i,j)=sqrt(
     &       myslope2(i,j)**2+( (zeta(i,j,1)-zeta(i,j-1,1))
     &    *pn(i,j))**2)
# ifdef MASKING
          endif
# endif
        enddo
      enddo
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        myslope2(-2,-2))
#      endif
!$acc end kernels
