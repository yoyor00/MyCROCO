!
!
! **********************************************************************
! **********************************************************************
       subroutine BIOSEDIM
!
! Appele par biologie.f a chaque passage: calcul des flux de sedimenta-
! tion
!
! ======================================================================
!
        Use comrunmod
        Use comdynmod
        Use mod_varphy_coupl
! Added by Alain:
	Use mod_eco3m
#ifdef MIGR
	REAL(8) :: ind_maxP,ind_maxL,rad,temp1,temp2,ind_PL,diffZ
	REAL(8) :: above_down,above_up,below_down,below_up
	REAL(8) :: grad_above,grad_below

	ind_PL = 1
        ! Find where PHYL is the largest...
	ind_maxP = 1
	do jtr=1,jptract
          if (VAR(jtr)%scomp == 'PHYL') then
		!write(*,*) "elmt is ", VAR(ivar)%elmt
             if (VAR(jtr)%elmt == 'cell') then
		!write(*,*) "elmt is ", VAR(ivar)%elmt
	        temp1 = 0.0
                do jk = 1,nzt
                   temp2 = VAR(jtr)%conc(1,1,jk)
		   !write(*,*) "conc at depth ", DEPT(jk), " m is ",VAR(ivar)%conc(1,1,jk)
                   if (temp2.GT.temp1) then
                     temp1 = temp2
                     ind_maxP = jk
                   endif
                enddo
              endif
           endif
        enddo

        ! Find where E_PARZ is at the isolume needed...
        ind_maxL = -1
        do jk = 1,nzt
           !write(*,*) "E_PARZ at ", DEPT(jk), " m is ",E_PARZ(1,1,jk)
           if (E_PARZ(1,1,jk).LT.0.2919) then
              if (ind_maxL == -1) then
                 ind_maxL = jk
              endif
           endif
        enddo

	! If isolume is deepest, set up swimming +/- 50 m around
        ! If PHYL is deepest, set up swimming +/- 20 m around
	if (ind_maxL.GT.ind_maxP) then
           ind_PL = ind_maxL
           rad = 50
        else
           ind_PL = ind_maxP
           rad = 20
        endif
	
	! Calculate the profile of Zooplankton velocities... this is necessary
	! because all the speeds need to be known for each calculation of flux...
	do jk=1,nzt
	    diffZ = (DEPT(ind_PL)-DEPT(jk))
	    if (diffZ.GT.rad) then
		VITZOO(jk) = 0.03
	    elseif (diffZ.LT.(-rad)) then
		VITZOO(jk) = -0.03
	    else
		VITZOO(jk) = 0.03*(diffZ/rad)
	    endif
	enddo
!	write(*,*) "ind_PL is", ind_PL 
!	write(*,*) "ind_maxL is ", ind_maxL
!	write(*,*) "ind_maxP is ", ind_maxP
!	write(*,*) "rad is ", rad
#endif

!        write(*,*) "Sed speed ", vitsed
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
#ifdef MIGR
          if (VAR(jtr)%scomp == 'COP') then
	     ! Re-set the SEDFLUX...
	     do jk=1,nzt
!		if (VAR(jtr)%elmt == 'cell') then
!		   write(*,*) "COPCell at ", DEPW(jk), " m is ",TENEUR(jk,jtr)
!		endif
		SEDFLUX(jk,jtr) = 0
	     enddo

	     ! TRY EXPLICIT GRAD EVALUATION : TRIAL
	     ! Calculate surface flux (in/out btwn 1st,2nd cell
!	     grad_below = (TENEUR(1,jtr)-TENEUR(2,jtr))/DE3T(1)
!	     grad_above = 0
!	     if (VITZOO(1).GE.0) then
!		SEDFLUX(1,jtr) = -grad_below*VITZOO(1)	     
!	     else
!		SEDFLUX(1,jtr) = grad_below*VITZOO(1)
!	     endif
!
!	     do jk=2,nzt-1
!		grad_above = (TENEUR(jk-1,jtr)-TENEUR(jk,jtr))/DE3T(jk-1)
!		grad_below = (TENEUR(jk,jtr)-TENEUR(jk+1,jtr))/DE3T(jk)
!		SEDFLUX(jk,jtr) = grad_above*VITZOO(jk-1) - grad_below*VITZOO(jk)
!	     enddo
!
!	     ! get bottom cell
!	     grad_above = (TENEUR(nzt-1,jtr)-TENEUR(nzt,jtr))/DE3T(nzt-1)
!	     grad_below = 0.0
!	     SEDFLUX(nzt,jtr) = grad_above*VITZOO(nzt-1)
	
	     above_up = 0.0  ! Upwards flux from current grid to above cell
	     above_down = 0.0	! Downward flux from above cell to current cell
	     below_up = 0.0	! Upward flux from deeper cell to current cell
	     below_down = 0.0	! Downward flux from current cell to deeper cell

	     ! Calculate flux from the surface...
	     ! If 'upwelling' at surface, ignore it, else calc flux out...
	     if (VITZOO(1).GE.0) then
		below_down = 1.0
	     else
		below_up = 1.0
	     endif
	     SEDFLUX(1,jtr) = -below_up*VITZOO(1)*TENEUR(2,jtr)/DE3T(1) &
	&		      - below_down*VITZOO(1)*TENEUR(1,jtr)/DE3T(1)
!	     SEDFLUX(1,jtr) = -below_up*VITZOO(1)*(TENEUR(2,jtr)+TENEUR(1,jtr))/2 &
!	&                     - below_down*VITZOO(1)*(TENEUR(2,jtr)+TENEUR(1,jtr))/2
!	     SEDFLUX(1,jtr) = -below_up*VITZOO(1)*TENEUR(2,jtr)/DE3T(1) - &
!	&			below_down*VITZOO(1)*TENEUR(1,jtr)/DE3T(1)

	     do jk=2,nzt-1
		! Re-set the coefficients for switching...
		above_up = 0.0
		above_down = 0.0
		below_up = 0.0
		below_down = 0.0
		if (VITZOO(jk-1).GE.0) then
		    above_down = 1.0
		else
		    above_up = 1.0
		endif
		if (VITZOO(jk).GE.0) then
		    below_down = 1.0
		else
		    below_up = 1.0
		endif
		! FLUXES ARE: a) input from above, b) output to above,
		! c) output to below, d) intput from below
		SEDFLUX(jk,jtr) = above_down*VITZOO(jk-1)*TENEUR(jk-1,jtr)/DE3T(jk)&
	&		+ above_up*VITZOO(jk-1)*TENEUR(jk,jtr)/DE3T(jk) &
	&		- below_down*VITZOO(jk)*TENEUR(jk,jtr)/DE3T(jk) &
	&		- below_up*VITZOO(jk)*TENEUR(jk+1,jtr)/DE3T(jk)	

!		SEDFLUX(jk,jtr) = above_down*VITZOO(jk-1)*(TENEUR(jk-1,jtr)+TENEUR(jk,jtr))/2 + &
!	&		above_up*VITZOO(jk-1)*(TENEUR(jk-1,jtr)+TENEUR(jk,jtr))/2 - &
!	&               below_down*VITZOO(jk)*(TENEUR(jk,jtr)+TENEUR(jk+1,jtr))/2 - &
!	&		below_up*VITZOO(jk)*(TENEUR(jk+1,jtr)+TENEUR(jk,jtr))/2

!		SEDFLUX(jk,jtr) = above_down*VITZOO(jk-1)*TENEUR(jk-1,jtr)/DE3T(jk-1) + &
!	&			above_up*VITZOO(jk-1)*TENEUR(jk,jtr)/DE3T(jk-1) - &
!	&			below_down*VITZOO(jk)*TENEUR(jk,jtr)/DE3T(jk) - &
!	&			below_up*VITZOO(jk)*TENEUR(jk+1,jtr)/DE3T(jk)

!	        SEDFLUX(1,jtr) = -VITZOO(1)*TENEUR(1,jtr)/DE3T(1)
	     enddo

	     below_up = 0.0
	     below_down = 0.0
             if (VITZOO(nzt-1).GE.0) then
		 above_down = 1.0
	      else
		 above_up = 1.0
             endif
	     SEDFLUX(nzt,jtr) = above_up*VITZOO(nzt-1)*TENEUR(nzt,jrt)/DE3T(jk-1)&
	&       + above_down*VITZOO(nzt-1)*TENEUR(nzt-1,jtr)/DE3T(jk-1)
!	     SEDFLUX(nzt,jtr) = above_up*VITZOO(nzt-1)*(TENEUR(nzt,jtr)+TENEUR(nzt-1,jtr))/2 + &
!	&		above_down*VITZOO(nzt-1)*(TENEUR(nzt,jtr)+TENEUR(nzt-1,jtr))/2

!	     SEDFLUX(nzt,jtr) = above_up*VITZOO(nzt-1)*TENEUR(nzt,jtr)/DE3T(jk-1) + &
!	&			above_down*VITZOO(nzt-1)*TENEUR(nzt-1,jtr)/DE3T(jk-1)

	     ! Bottom cell

	  else   ! NOT A COPEPOD
             
             ! Surface flux...
	     SEDFLUX(1,jtr)=-VITSED(jtr)*TENEUR(1,jtr)/DE3T(1)
             do jk=2,nzt-2
		SEDFLUX(jk,jtr)=VITSED(jtr)*(TENEUR(jk-1,jtr)-TENEUR(jk,jtr))&
	&              /DE3T(jk)
	     enddo
	     SEDFLUX(nzt-1,jtr)=VITSED(jtr)*(TENEUR(nzt-2,jtr))&
	&           /DE3T(nzt-1)
	     SEDFLUX(nzt,jtr)=0.0

	     
          endif  ! END COPEPOD IF             
#else
	   ! THIS MEANS THERE IS NO MIGRATION!
	   SEDFLUX(1,jtr)=-VITSED(jtr)*TENEUR(1,jtr)/DE3T(1)
	   do jk=2,nzt-2
              SEDFLUX(jk,jtr)=VITSED(jtr)*(TENEUR(jk-1,jtr)-TENEUR(jk,jtr))&
	&            /DE3T(jk)
           enddo
	   SEDFLUX(nzt-1,jtr)=VITSED(jtr)*(TENEUR(nzt-2,jtr))&
	&         /DE3T(nzt-1)
	   SEDFLUX(nzt,jtr)=0.0

#endif
	enddo
       
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt
           SEDBIO(jk,jtr)=SEDBIO(jk,jtr)+SEDFLUX(jk,jtr)
         enddo
       enddo
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt-1
           TENEUR(jk,jtr)=TENEUR(jk,jtr)+dts*SEDFLUX(jk,jtr)
         enddo
       enddo
!
       return
       end
