      do j=JU_RANGE          
        do i=IU_RANGE
          DC3(i,j,0)=0.
          CF3(i,j,0)=0.
          FC3(i,j,0)=0.
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            DC3(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC3(i,j,0)=DC3(i,j,0)+DC3(i,j,k)
            CF3(i,j,0)=CF3(i,j,0)+DC3(i,j,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=IU_RANGE
          DC3(i,j,0)=1./DC3(i,j,0)
          CF3(i,j,0)=DC3(i,j,0)*(CF3(i,j,0)-DU_avg1(i,j,nnew))
#  ifdef MRL_WCI
     &                                        +ust2d(i,j)
#  endif
          ubar(i,j,knew)=DC3(i,j,0)*DU_avg1(i,j,nnew)
#  ifdef MRL_WCI
     &                                        -ust2d(i,j)
#   ifdef MASKING
          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
     &                +ust2d(i,j)*(umask(i,j)-1.0)
#   endif
#  endif
#  ifdef WET_DRY
          ubar(i,j,knew)=ubar(i,j,knew)*umask_wet(i,j)
#  endif
        enddo
        do k=N,1,-1
          do i=IU_RANGE
            u(i,j,k,nnew)=(u(i,j,k,nnew)-CF3(i,j,0))
#  ifdef MASKING
     &                                   *umask(i,j)
#   ifdef MRL_WCI
     &                  +ust(i,j,k)*(umask(i,j)-1.0)
#   endif
#  endif
#  ifdef WET_DRY
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
#   ifdef MRL_WCI
            ust(i,j,k)=ust(i,j,k)*umask_wet(i,j)
#   endif
#  endif
            workr(i,j,k) = DC3(i,j,k)*(u(i,j,k,nstp)+u(i,j,k,nnew)
#  ifdef MRL_WCI
     &                                          +2.0*ust(i,j,k)
#  endif
     &                                                         ) 
          enddo
        enddo
      enddo
!
! Exchange MPI domain or periodic boundaries
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   if defined LERAY_FILTER_9PTS
      call exchange_u3d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   elif defined LERAY_FILTER_7PTS
      call exchange_u3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   else
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   endif
#  endif    
!
! Flux spatial filtering
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
# ifdef LERAY_FILTER_9PTS
            if(u_fwidth_array(i,j).eq.9) then
              work2d(i,j) = sum( workr(i-4:i+4,j-4:j+4,k)
     &                      * filter_weights(-4:4,-4:4) )
     &                   * inv_weight_sum9
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
            if(u_fwidth_array(i,j).eq.7) then
              work2d(i,j) = sum( workr(i-3:i+3,j-3:j+3,k)
     &                      * filter_weights(-3:3,-3:3) )
     &                   * inv_weight_sum7
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
            if(u_fwidth_array(i,j).eq.5) then
              work2d(i,j) = sum( workr(i-2:i+2,j-2:j+2,k)
     &                      * filter_weights(-2:2,-2:2) )
     &                   * inv_weight_sum5
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
            if(u_fwidth_array(i,j).eq.3) then
              work2d(i,j) = sum( workr(i-1:i+1,j-1:j+1,k)
     &                      * filter_weights(-1:1,-1:1) )
     &                   * inv_weight_sum3
            endif
# endif
            if(u_fwidth_array(i,j).lt.3) then
              work2d(i,j) = workr(i,j,k)
            endif
          enddo
        enddo
        workr(IstrU:Iend,Jstr:Jend,k) = work2d(IstrU:Iend,Jstr:Jend)
      enddo
!      
! Exchange MPI domain or periodic boundaries
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#  endif 
!
      do j=JU_RANGE
        do k=N,1,-1
          do i=IU_RANGE
#  define EPSIL 0.125
#  define DELTA 0.75
            FC3(i,j,k)=DELTA*Huon(i,j,k) + EPSIL*workr(i,j,k)
            FC3(i,j,0)=FC3(i,j,0)+FC3(i,j,k)
          enddo
        enddo
        do i=IU_RANGE
          FC3(i,j,0)=DC3(i,j,0)*(FC3(i,j,0)-DU_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            Huon(i,j,k)=FC3(i,j,k)-DC3(i,j,k)*FC3(i,j,0)
          enddo
        enddo
      enddo

      do j=JU_RANGE
