! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_AGRIF1.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
      if (.NOT.Agrif_Root()) THEN
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Zt_avg3(i,j,iif)=zeta(i,j,knew)
          enddo
        enddo 
        do j=JstrR,JendR
          do i=Istr,IendR
            DU_avg3(i,j,iif) = 0.5*(h(i  ,j)+zeta(i  ,j,knew)+
     &                              h(i-1,j)+zeta(i-1,j,knew)) 
     &                               *on_u(i,j)*ubar(i,j,knew)
          enddo
        enddo 
        do j=Jstr,JendR
          do i=IstrR,IendR
            DV_avg3(i,j,iif) = 0.5*(h(i,j  )+zeta(i,j  ,knew)+
     &                              h(i,j-1)+zeta(i,j-1,knew)) 
     &                               *om_v(i,j)*vbar(i,j,knew)
          enddo
        enddo
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
        call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                          Zt_avg3(START_2D_ARRAY,iif))
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg3(START_2D_ARRAY,iif))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg3(START_2D_ARRAY,iif))
#  endif
#  ifdef RVTK_DEBUG_ADVANCED
C$OMP BARRIER
C$OMP MASTER
        call check_tab2d(Zt_avg3(:,:,iif),'Zt_avg3 step3d_fast','r'
      &    ,ondevice=.TRUE.)
        call check_tab2d(DU_avg3(:,:,iif),'DU_avg3 step3d_fast','u'
      &    ,ondevice=.TRUE.)
        call check_tab2d(DV_avg3(:,:,iif),'DV_avg3 step3d_fast','v'
      &    ,ondevice=.TRUE.)
C$OMP END MASTER
#  endif   
      endif
#  ifdef AGRIF_CONSERV_VOL
      if (iif.eq.nfast) then
        if (agrif_root()) then
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
            enddo
          enddo
        else
          do j=JstrR,JendR
            do i=IstrR,IendR
              DU_avg1(i,j,5) = dt * DU_avg2(i,j)
              DV_avg1(i,j,5) = dt * DV_avg2(i,j)
              DU_avg1(i,j,4) = DU_avg1(i,j,4) + DU_avg1(i,j,5)
              DV_avg1(i,j,4) = DV_avg1(i,j,4) + DV_avg1(i,j,5)
            enddo
          enddo       
        endif
      endif
#  endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_AGRIF1.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
