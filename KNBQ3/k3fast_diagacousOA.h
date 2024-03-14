! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_diagacousOA.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#    if defined ONLINE_ANALYSIS && defined NBQ
C$OMP BARRIER
!C$OMP MASTER
!      call check_tab3d(rho_nbq,'step3d_fastrho_nbq(a)','r')
!      call check_tab3d(z_w,'z_w (a)','r')
!      call check_tab2d(xr,'xr (a)','r')
!      call check_tab2d(yr,'yr (a)','r')
!C$OMP END MASTER
#    endif   
    
      ic=1
      ivc=1
      ! Name OA variable
      !print *,"DIAG SACOUS Var name =",tgvnam_oa(tv_oa(ic))
     
!$acc kernels if(compute_on_device) default(present)
        if (JstrV.le.Jend) then
          jvar1=JstrV-2
          jvar2=Jend+1
        else
          jvar1=JstrV-1
          jvar2=Jend
        endif
!        write(47,*) 42
        do j=jvar1,jvar2
         do i=IstrU-1,Iend
          do k=-N_sl+1,N
             mvoa1(i,j,k)=max(mvoa1(i,j,k), 
     &  abs(var3d_oa( i,j,k, tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ))))
             mvoa2(i,j,k)=max(mvoa2(i,j,k), 
     &  abs(var3d_oa( i,j,k, tvar_oa( tc_oa(ic+1), tvc_oa(ic+1), ivc ))))
          enddo
         enddo
        enddo
!$acc end kernels     
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_diagacousOA.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
