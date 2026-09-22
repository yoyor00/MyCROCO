
! ! gls_smooth.h begin
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          TRB_NEW(Istr-1,j,k)=TRB_NEW(Istr,j,k)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          TRB_NEW(Iend+1,j,k)=TRB_NEW(Iend,j,k)
        enddo
      endif
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          TRB_NEW(i,Jstr-1,k)=TRB_NEW(i,Jstr,k)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          TRB_NEW(i,Jend+1,k)=TRB_NEW(i,Jend,k)
        enddo
      endif
# ifndef EW_PERIODIC
      if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
        TRB_NEW(Istr-1,Jstr-1,k)=TRB_NEW(Istr,Jstr,k)
      endif
      if (WESTERN_EDGE.and.NORTHERN_EDGE) then
        TRB_NEW(Istr-1,Jend+1,k)=TRB_NEW(Istr,Jend,k)
      endif
      if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
        TRB_NEW(Iend+1,Jstr-1,k)=TRB_NEW(Iend,Jstr,k)
      endif
      if (EASTERN_EDGE.and.NORTHERN_EDGE) then
        TRB_NEW(Iend+1,Jend+1,k)=TRB_NEW(Iend,Jend,k)
      endif
# endif
#endif

         DO j=jstr-1,jend+1
            DO i=istr,iend+1
#ifndef OPENACC      
               FX(i,j  )=( TRB_NEW(i  ,j,k)
#else
               FX(i,j,k)=( TRB_NEW(i  ,j,k)
#endif               
     &                   -  TRB_NEW(i-1,j,k) )
#ifdef MASKING
     &                             *umask(i,j)
#endif
            ENDDO
         ENDDO
         DO j=jstr,jend+1
            DO i=istr-1,iend+1
#ifndef OPENACC      
               FE1(i,j,0)=( TRB_NEW(i,j  ,k)
#else
               FE1(i,j,k)=( TRB_NEW(i,j  ,k)
#endif
     &                    - TRB_NEW(i,j-1,k) )
#ifdef MASKING
     &                             *vmask(i,j)
#endif
            ENDDO
            DO i=istr,iend
#ifndef OPENACC      
              FE(i,j)=FE1(i,j,0)
     &                + smth_a*( FX(i+1,j)+FX(i  ,j-1)
     &                          -FX(i  ,j)-FX(i+1,j-1))
#else
              FE(i,j,k)=FE1(i,j,k)
     &                + smth_a*( FX(i+1,j,k)+FX(i  ,j-1,k)
     &                          -FX(i  ,j,k)-FX(i+1,j-1,k))
#endif              
            ENDDO
         ENDDO

         DO j=jstr,jend
            DO i=istr,iend+1
#ifndef OPENACC      
              FX(i,j)=FX(i,j  )
     &                + smth_a*( FE1(i,j+1,0)+FE1(i-1,j  ,0)
     &                          -FE1(i,j  ,0)-FE1(i-1,j+1,0))
#else
              FX(i,j,k)=FX(i,j,k)
     &                + smth_a*( FE1(i,j+1,k)+FE1(i-1,j  ,k)
     &                          -FE1(i,j  ,k)-FE1(i-1,j+1,k))
#endif              
            ENDDO
            DO i=istr,iend
               trb(i,j,k,nnew,ig)=TRB_NEW(i,j,k)
#ifndef OPENACC      
     &                        + smth_b*( FX(i+1,j)-FX(i,j)
     &                                  +FE(i,j+1)-FE(i,j) )
#else
     &                        + smth_b*( FX(i+1,j,k)-FX(i,j,k)
     &                                  +FE(i,j+1,k)-FE(i,j,k) )
#endif              
#ifdef MASKING
     &                                               *rmask(i,j)
#endif
            ENDDO
         ENDDO              !--> discard FX,FE,FE1
! ! gls_smooth.h end
