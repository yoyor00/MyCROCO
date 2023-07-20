      enddo

      do j=JV_RANGE    
        do i=IV_RANGE
          DC3(i,j,0)=0.
          CF3(i,j,0)=0.
          FC3(i,j,0)=0.
        enddo
        do k=1,N,+1
          do i=IV_RANGE
            DC3(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
            DC3(i,j,0)=DC3(i,j,0)+DC3(i,j,k)
            CF3(i,j,0)=CF3(i,j,0)+DC3(i,j,k)*v(i,j,k,nnew)
          enddo
        enddo
        do i=IV_RANGE
          DC3(i,j,0)=1./DC3(i,j,0)
          CF3(i,j,0)=DC3(i,j,0)*(CF3(i,j,0)-DV_avg1(i,j,nnew))
#  ifdef MRL_WCI
     &                                        +vst2d(i,j)
#  endif
          vbar(i,j,knew)=DC3(i,j,0)*DV_avg1(i,j,nnew)
#  ifdef MRL_WCI
     &                                        -vst2d(i,j)
#   ifdef MASKING
          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
     &                +vst2d(i,j)*(vmask(i,j)-1.0)
#   endif
#  endif
#  ifdef WET_DRY
          vbar(i,j,knew)=vbar(i,j,knew)*vmask_wet(i,j)
#  endif
        enddo
        do k=N,1,-1
          do i=IV_RANGE
            v(i,j,k,nnew)=(v(i,j,k,nnew)-CF3(i,j,0))
#  ifdef MASKING
     &                                   *vmask(i,j)
#   ifdef MRL_WCI
     &                  +vst(i,j,k)*(vmask(i,j)-1.0)
#   endif
#  endif
#  ifdef WET_DRY
            v(i,j,k,nnew)=v(i,j,k,nnew)*vmask_wet(i,j)
#   ifdef MRL_WCI
            vst(i,j,k)=vst(i,j,k)*vmask_wet(i,j)
#   endif
#  endif
            workr(i,j,k) = DC3(i,j,k)*(v(i,j,k,nstp)+v(i,j,k,nnew)
#  ifdef MRL_WCI
     &                                          +2.0*vst(i,j,k)
#  endif
     &                                                         ) 
          enddo
        enddo
      enddo
!
! Exchange MPI domain or periodic boundaries
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
#   if defined LERAY_FILTER_9PTS
      call exchange_v3d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   elif defined LERAY_FILTER_7PTS
      call exchange_v3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   else
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#   endif
#  endif    
!
! Flux spatial filtering
      do k=1,N
        do j=JstrV,Jend
          do i=Istr,Iend
# ifdef LERAY_FILTER_9PTS
            if(v_fwidth_array(i,j).eq.9) then
              work2d(i,j) = sum( workr(i-4:i+4,j-4:j+4,k)
     &                      * filter_weights(-4:4,-4:4) )
     &                   * inv_weight_sum9
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
            if(v_fwidth_array(i,j).eq.7) then
              work2d(i,j) = sum( workr(i-3:i+3,j-3:j+3,k)
     &                      * filter_weights(-3:3,-3:3) )
     &                   * inv_weight_sum7
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
            if(v_fwidth_array(i,j).eq.5) then
              work2d(i,j) = sum( workr(i-2:i+2,j-2:j+2,k)
     &                      * filter_weights(-2:2,-2:2) )
     &                   * inv_weight_sum5
            endif
# endif
# if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
            if(v_fwidth_array(i,j).eq.3) then
              work2d(i,j) = sum( workr(i-1:i+1,j-1:j+1,k)
     &                      * filter_weights(-1:1,-1:1) )
     &                   * inv_weight_sum3
            endif
# endif
            if(v_fwidth_array(i,j).lt.3) then
              work2d(i,j) = workr(i,j,k)
            endif
          enddo
        enddo
        workr(Istr:Iend,JstrV:Jend,k) = work2d(Istr:Iend,JstrV:Jend)
      enddo
!      
! Exchange MPI domain or periodic boundaries
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                             workr(START_2D_ARRAY,1))
#  endif 
!
      do j=JV_RANGE
        do k=N,1,-1
          do i=IV_RANGE
            FC3(i,j,k)=DELTA*Hvom(i,j,k) + EPSIL*workr(i,j,k)
#  undef DELTA
#  undef EPSIL
            FC3(i,j,0)=FC3(i,j,0)+FC3(i,j,k)
          enddo
        enddo
        do i=IV_RANGE
          FC3(i,j,0)=DC3(i,j,0)*(FC3(i,j,0)-DV_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=IV_RANGE
            Hvom(i,j,k)=FC3(i,j,k)-DC3(i,j,k)*FC3(i,j,0)
          enddo
        enddo
      enddo

      do j=JV_RANGE
