#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#    if defined LERAY_FILTER_9PTS
      call exchange_u2d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                    work2d(START_2D_ARRAY))
      call exchange_v2d_4pts_tile (Istr,Iend,Jstr,Jend,
     &                    work2d2(START_2D_ARRAY))
#    elif defined LERAY_FILTER_7PTS
      call exchange_u2d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                    work2d(START_2D_ARRAY))
      call exchange_v2d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                    work2d2(START_2D_ARRAY))
#    else
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   work2d(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   work2d2(START_2D_ARRAY))
#    endif      
#   endif
      do j=JU_EXT_RANGE
        do i=imin+1,imax
#   ifdef LERAY_FILTER_9PTS
          if(u_fwidth_array(i,j).eq.9) then
            urhs_smooth(i,j) = sum( work2d(i-4:i+4,j-4:j+4)
     &                         * filter_weights(-4:4,-4:4) )
     &                         * inv_weight_sum9
          endif
#   endif          
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
          if(u_fwidth_array(i,j).eq.7) then 
            urhs_smooth(i,j) = sum( work2d(i-3:i+3,j-3:j+3)
     &                         * filter_weights(-3:3,-3:3) )
     &                         * inv_weight_sum7
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
          if(u_fwidth_array(i,j).eq.5) then 
            urhs_smooth(i,j) = sum( work2d(i-2:i+2,j-2:j+2)
     &                         * filter_weights(-2:2,-2:2) )
     &                         * inv_weight_sum5
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
          if(u_fwidth_array(i,j).eq.3) then 
            urhs_smooth(i,j) = sum( work2d(i-1:i+1,j-1:j+1)
     &                         * filter_weights(-1:1,-1:1) )
     &                         * inv_weight_sum3
          endif
#   endif
          if(u_fwidth_array(i,j).lt.3) then
              urhs_smooth(i,j) = work2d(i,j)
          endif
        enddo
      enddo
!
      do j=jmin+1,jmax
        do i=IV_EXT_RANGE
#   ifdef LERAY_FILTER_9PTS
          if(v_fwidth_array(i,j).eq.9) then
            vrhs_smooth(i,j) = sum( work2d2(i-4:i+4,j-4:j+4)
     &                         * filter_weights(-4:4,-4:4) )
     &                         * inv_weight_sum9
          endif
#   endif          
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS
          if(v_fwidth_array(i,j).eq.7) then 
            vrhs_smooth(i,j) = sum( work2d2(i-3:i+3,j-3:j+3)
     &                         * filter_weights(-3:3,-3:3) )
     &                         * inv_weight_sum7
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS
          if(v_fwidth_array(i,j).eq.5) then 
            vrhs_smooth(i,j) = sum( work2d2(i-2:i+2,j-2:j+2)
     &                         * filter_weights(-2:2,-2:2) )
     &                         * inv_weight_sum5
          endif  
#   endif
#   if defined LERAY_FILTER_9PTS | defined LERAY_FILTER_7PTS | defined LERAY_FILTER_5PTS | defined LERAY_FILTER_3PTS
          if(v_fwidth_array(i,j).eq.3) then 
            vrhs_smooth(i,j) = sum( work2d2(i-1:i+1,j-1:j+1)
     &                         * filter_weights(-1:1,-1:1) )
     &                         * inv_weight_sum3
          endif
#   endif
          if(v_fwidth_array(i,j).lt.3) then
            vrhs_smooth(i,j) = work2d2(i,j)
          endif
        enddo
      enddo
#   if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                   urhs_smooth(START_2D_ARRAY))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                   vrhs_smooth(START_2D_ARRAY))
#   endif 
