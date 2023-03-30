! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! m3fast_fbf.h
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
#  define UBOT UFx
#  define VBOT VFe
      if (mod(iif-1,inc_faststep).eq.0) then

        if (maxval(Zob).ne.0.) then
          do j=JstrV-1,Jend+1
            do i=IstrU-1,Iend+1
              UBOT(i,j)=2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
#   ifdef MVB
     &                  -u_mvb(i,j,knew2)
#   endif
              VBOT(i,j)=2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
#   ifdef MVB
     &                  -v_mvb(i,j,knew2)
#   endif
             enddo
          enddo
          do j=JstrV-1,Jend
            do i=IstrU-1,Iend
              cff1=z_r(i,j,1)-z_w(i,j,0)
              cff=vonKar/LOG(cff1/Zob(i,j))
              work(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
            enddo
          enddo
          do j=Jstr,Jend
            do i=IstrU,Iend
              cff=0.25*(VBOT(i  ,j)+VBOT(i  ,j+1)+
     &                  VBOT(i-1,j)+VBOT(i-1,j+1))
              bustr(i,j)=0.5*(work(i-1,j)+work(i,j))*UBOT(i,j)*
     &                   SQRT(UBOT(i,j)*UBOT(i,j)+cff*cff)
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              cff=0.25*(UBOT(i,j  )+UBOT(i+1,j  )+
     &                  UBOT(i,j-1)+UBOT(i+1,j-1))
              bvstr(i,j)=0.5*(work(i,j-1)+work(i,j))*VBOT(i,j)*
     &                   SQRT(VBOT(i,j)*VBOT(i,j)+cff*cff)
            enddo
          enddo
        elseif (rdrg2.gt.0.) then
          do j=JstrV-1,Jend+1
            do i=IstrU-1,Iend+1
              UBOT(i,j)=2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
#   ifdef MVB
     &                  -u_mvb(i,j,knew2)
#   endif
              VBOT(i,j)=2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
#   ifdef MVB
     &                  -v_mvb(i,j,knew2)
#   endif
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              cff=0.25*(VBOT(i  ,j)+VBOT(i  ,j+1)+
     &                  VBOT(i-1,j)+VBOT(i-1,j+1))
              bustr(i,j)=rdrg2*UBOT(i,j)*
     &                   SQRT(UBOT(i,j)*UBOT(i,j)+cff*cff)
            enddo
          enddo
          do j=Jstr,Jend
            do i=IstrU,Iend
              cff=0.25*(UBOT(i,j  )+UBOT(i+1,j  )+
     &                  UBOT(i,j-1)+UBOT(i+1,j-1))
              bvstr(i,j)=rdrg2*VBOT(i,j)*
     &                    SQRT(VBOT(i,j)*VBOT(i,j)+cff*cff)
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=IstrU,Iend
              bustr(i,j)=rdrg*(2.*qdmu_nbq(i,j,1)/(Hz(i,j,1)+Hz(i-1,j,1))
#   ifdef MVB
     &                  -u_mvb(i,j,knew2)
#   endif
     &                  )
            enddo
          enddo
          do j=JstrV,Jend
            do i=Istr,Iend
              bvstr(i,j)=rdrg*(2.*qdmv_nbq(i,j,1)/(Hz(i,j,1)+Hz(i,j-1,1))
#   ifdef MVB
     &                  -v_mvb(i,j,knew2)
#   endif
     &                  )
            enddo
          enddo
        endif
      endif
#  undef UBOT
#  undef VBOT
