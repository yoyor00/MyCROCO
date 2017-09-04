 
#ifdef SOLVE3D

      cff1=weight(1,iif)
      cff2=weight(2,iif)

# ifndef NBQ_ZETAW
      if (FIRST_2D_STEP) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
          enddo
        enddo 
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
# endif

# ifdef NBQ_ZETAW
      if (iif==1) then
        do j=JstrR,JendR
          do i=IstrR,IendR
        !   DU_avg2(i,j)=cff2*DU_nbq(i,j)*on_u(i,j)   !DUon(i,j)
            DU_avg2(i,j)=cff2*DUon(i,j)
        !   DV_avg2(i,j)=cff2*DV_nbq(i,j)*om_v(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
          enddo
        enddo 
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
        !   DU_avg2(i,j)=DU_avg2(i,j)+cff2*DU_nbq(i,j)*on_u(i,j) !DUon(i,j)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j) 
        !   DV_avg2(i,j)=DV_avg2(i,j)+cff2*DV_nbq(i,j)*om_v(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
# endif

#endif  /* SOLVE3D */
