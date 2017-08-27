
!
!-----------------------------------------------------------------------
!
! Get filtered rhs terms (M2FILTER_NONE)
! and multiply by dx*dy to get units of rho*Hz*dx*dy*ru
!
!-----------------------------------------------------------------------
!     
        if (LAST_2D_STEP) then
!
!-----------------------------
! Finalizes rho-avg1 variables
!-----------------------------   
!     
# ifdef NBQ_MASS
          do j=Jstr,Jend
            do i=Istr,Iend
              rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
            enddo
          enddo 
          do k=1,N
            do j=Jstr,Jend
              do i=Istr,Iend
#  ifdef NBQ_IJK
                rho_nbq_avg1(i,j,k)=1.d0+
     &         (rho_nbq(i,j,k)*Hzr_half_nbq_inv(i,j,k)+rho(i,j,k)/rho0 )   !XXX
#  else
                l_nbq=ijk2lq_nh(i,j,k)
                rho_nbq_avg1(i,j,k)=1.d0+
     &               (rhp_nbq_a(l_nbq)+rho(i,j,k)/rho0  )
#  endif  /* NBQ_IJK */
     !           rho_nbq_avg1(i,j,k)=1.d0
     !           rho_nbq_avg2(i,j,k)=rho_nbq_avg2(i,j,k)/real(nfast)
     !           rho_nbq_avg2(i,j,k)=1.d0+
     !&               (rho_nbq_avg2(i,j,k)+rho(i,j,k)/rho0  )
              enddo
            enddo 
          enddo
# endif /* NBQ_MASS */
!          rhobar_nbq_avg1=1.
!          rho_nbq_avg1=1.
!
!-----------------------------
! Finalizes Ru_avg2 variables
!-----------------------------   

          do k=1,N
            do j=Jstr,Jend
              do i=IstrU,Iend
     !    !     ru_nbq_avg1(i,j,k)=ru_nbq_ext(i,j,k)
                ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
# ifdef NBQ_ZETAW
     &             -ruext_nbq_2d_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
# else
     &             -ruext_nbq_2d_old(i,j)*(Hz_half(i-1,j,k)+Hz_half(i,j,k))
# endif
# ifdef NBQ_IJK
                ru_nbq_avg2(i,j,k)=
     &             ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt
     &             -ru_int_nbq(i,j,k)-(ruext_nbq_2d_sum(i,j)/nfast)*
# ifdef NBQ_ZETAW
     &             (Hz(i,j,k)+Hz(i-1,j,k)))*on_u(i,j)*om_u(i,j)
# else
     &             (Hz_half(i,j,k)+Hz_half(i-1,j,k)))*on_u(i,j)*om_u(i,j)
#endif
# else
                 l_nbq=ijk2lmom_nh(i,j,k,1)
                 ru_nbq_avg2(i,j,k)=
     &             ((qdm_nbq_a(l_nbq)-ru_nbq_avg2(i,j,k))/dt
     &             -ru_int_nbq(i,j,k)-(ruext_nbq_2d_sum(i,j)/nfast)*
     &             (Hz(i,j,k)+Hz(i-1,j,k)))*on_u(i,j)*om_u(i,j)
# endif


              enddo
            enddo 
          enddo          

          do k=1,N
            do j=JstrV,Jend
              do i=Istr,Iend             
          !     rv_nbq_avg1(i,j,k)=rv_nbq_ext(i,j,k) 
                rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
# ifdef NBQ_ZETAW
     &            -rvext_nbq_2d_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
# else
     &            -rvext_nbq_2d_old(i,j)*(Hz_half(i,j-1,k)+Hz_half(i,j,k))
# endif

# if defined NBQ_IJK
                rv_nbq_avg2(i,j,k)=
     &             ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt
     &             -rv_int_nbq(i,j,k)-(rvext_nbq_2d_sum(i,j)/nfast)*
# ifdef NBQ_ZETAW
     &             (Hz(i,j,k)+Hz(i,j-1,k)))*on_v(i,j)*om_v(i,j)   
# else         
     &             (Hz_half(i,j,k)+Hz_half(i,j-1,k)))*on_v(i,j)*om_v(i,j)         
# endif        
# else
                l_nbq=ijk2lmom_nh(i,j,k,2)
                rv_nbq_avg2(i,j,k)=
     &             ((qdm_nbq_a(l_nbq)-rv_nbq_avg2(i,j,k))/dt
     &             -rv_int_nbq(i,j,k)-(rvext_nbq_2d_sum(i,j)/nfast)*
     &             (Hz(i,j,k)+Hz(i,j-1,k)))*on_v(i,j)*om_v(i,j)
# endif
              enddo
            enddo 
          enddo

          do k=1,N   ! Francis: switch from k=0 to k=1 here
            do j=Jstr,Jend
              do i=Istr,Iend
         !      rw_nbq_avg1(i,j,k)=rw_nbq_ext(i,j,k)
# if defined NBQ_IJK
                rw_nbq_avg2(i,j,k)=
     &             ((qdmw_nbq(i,j,k)-rw_nbq_avg2(i,j,k))/dt
     &             -rw_int_nbq(i,j,k))*on_r(i,j)*om_r(i,j)
# else
                l_nbq=ijk2lmom_nh(i,j,k,3)
                 rw_nbq_avg2(i,j,k)=
     &              ((qdm_nbq_a(l_nbq)-rw_nbq_avg2(i,j,k))/dt
     &              -rw_int_nbq(i,j,k))*on_r(i,j)*om_r(i,j)
# endif
              enddo
            enddo 
          enddo   

        endif

