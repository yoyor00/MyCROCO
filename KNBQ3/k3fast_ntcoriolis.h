! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_ntcoriolis.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#  define Huw ntcoru
#  define Hvw ntcorv
#ifdef OPENACC
    stop 'to be completed OpenACC'
#endif
      do j=Jstr,Jend
        do i=Istr,Iend+1
          Huw(i,j,0)=0.5*qdmu_nbq(i,j,1)
          Huw(i,j,N)=0.5*qdmu_nbq(i,j,N)
        enddo
      enddo
      do j=Jstr,Jend+1
        do i=Istr,Iend
          Hvw(i,j,0)=0.5*qdmv_nbq(i,j,1)
          Hvw(i,j,N)=0.5*qdmv_nbq(i,j,N)
        enddo
      enddo
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend+1
            Huw(i,j,k)=0.5*(qdmu_nbq(i,j,k)+qdmu_nbq(i,j,k+1))
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=Istr,Iend
            Hvw(i,j,k)=0.5*(qdmv_nbq(i,j,k)+qdmv_nbq(i,j,k+1))
          enddo
        enddo
      enddo
      do k=0,N
        do j=Jstr,Jend
          do i=Istr,Iend
            ntcorw(i,j,k) = 0.5*e(i,j)*(
     &          cosa(i,j)*(Huw(i,j,k)+Huw(i+1,j,k))
     &        - sina(i,j)*(Hvw(i,j,k)+Hvw(i,j+1,k))
     &                                 )
          enddo
        enddo
      enddo
#  undef Huw
#  undef Hvw
!
      do k=1,N
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            UFx(i,j)=-0.5*e(i,j)*cosa(i,j)
     &                            *( qdmw_nbq(i,j,k-1)
     &                              +qdmw_nbq(i,j,k  )) 
            VFe(i,j)=+0.5*e(i,j)*sina(i,j)
     &                            *( qdmw_nbq(i,j,k-1)
     &                              +qdmw_nbq(i,j,k)  )
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            ntcoru(i,j,k)=0.5*(UFx(i,j)+UFx(i-1,j))
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            ntcorv(i,j,k)=0.5*(VFe(i,j)+VFe(i,j-1))
          enddo
        enddo
      enddo
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_ntcoriolis.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
