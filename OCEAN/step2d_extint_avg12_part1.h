 
#ifdef SOLVE3D

      cff1=weight(1,iif)
      cff2=weight(2,iif)

      if (FIRST_2D_STEP) then
        do j=JstrR,JendR
          do i=IstrR,IendR
# ifndef NBQ_ZETAW
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
# endif
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo 
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
# ifndef NBQ_ZETAW
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
# endif
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif

#endif  /* SOLVE3D */
