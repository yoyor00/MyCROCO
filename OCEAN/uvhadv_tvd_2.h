        if (nnew.ne.3) then ! corrector     
        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cff = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
            rrp= (u(i+2,j,k,nstp)-u(i+1,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i+1,j)
#   endif             
            rr = (u(i+1,j,k,nstp)-u(i  ,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i  ,j)
#   endif
            rrm= (u(i  ,j,k,nstp)-u(i-1,j,k,nstp))
#   ifdef MASKING
     &         *rmask(i-1,j)
#   endif
            Fxl=( cff*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
     &       -abs(cff)*rr )*0.5
       vit=velux(i,j)/Hz(i,j,k) * pm(i,j)
       UFx(i,j)=Fxl+limiteur_h3(rrm,rr,rrp,vit,Fxl,UFx(i,j))       
          enddo
        enddo

        do j=JstrV-1,Jend
          do i=Istr,Iend
            cff = 0.5*(Hvom(i,j,k)+Hvom(i,j+1,k))
            rrp= (v(i,j+2,k,nstp)-v(i,j+1,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i,j+1)
#   endif
            rr = (v(i,j+1,k,nstp)-v(i,j,k,nstp)) 
#   ifdef MASKING
     &         *rmask(i,j  )
#   endif
            rrm= (v(i,j  ,k,nstp)-v(i,j-1,k,nstp))
#   ifdef MASKING
     &         *rmask(i,j-1)
#   endif                    
            Fxl=( cff*(v(i,j,k,nstp) + v(i,j+1,k,nstp))
     &        -abs(cff)*rr )*0.5
            vit=velve(i,j)/Hz(i,j,k)*pn(i,j)
            VFe(i,j)=Fxl+limiteur_h3(rrm,rr,rrp,vit,Fxl,VFe(i,j))       
          enddo
        enddo

        do j=JstrV,Jend+1
          do i=IstrU,Iend
            cff = 0.5*(Hvom(i,j,k)+Hvom(i-1,j,k))
            rrp= (u(i,j+1,k,nstp)-u(i,j  ,k,nstp)) 
#   ifdef MASKING
     &         *pmask(i,j+1)
#   endif
            rr = (u(i,j  ,k,nstp)-u(i,j-1,k,nstp)) 
#   ifdef MASKING
     &         *pmask(i,j  )
#   endif
            rrm= (u(i,j-1,k,nstp)-u(i,j-2,k,nstp))
#   ifdef MASKING
     &         *pmask(i,j-1)
#   endif
            Fxl=( cff*(u(i,j,k,nstp)+u(i,j-1,k,nstp))
     &        -abs(cff)*rr )*0.5
       vit=velue(i,j)/on_p(i,j)
     &    /( 0.25*( Hz(i  ,j,k)+Hz(i  ,j-1,k)
     &           +  Hz(i-1,j,k)+Hz(i-1,j-1,k) ) )     
       UFe(i,j)=Fxl+limiteur_h3(rrm,rr,rrp,vit,Fxl,UFe(i,j))       
          enddo
        enddo

        do j=JstrV,Jend
          do i=IstrU,Iend+1
            cff = 0.5*(Huon(i,j,k)+Huon(i,j-1,k))
            rrp=(v(i+1,j,k,nstp)-v(i  ,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i+1,j)
#   endif
            rr =(v(i  ,j,k,nstp)-v(i-1,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i  ,j)
#   endif
            rrm=(v(i-1,j,k,nstp)-v(i-2,j,k,nstp))
#   ifdef MASKING
     &         *pmask(i-1,j)
#   endif
            cff1=( cff*(v(i,j,k,nstp)+v(i-1,j,k,nstp))
     &       -abs(cff)*rr )*0.5
            vit=velvx(i,j)/om_p(i,j)
     &    /( 0.25*( Hz(i  ,j,k)+Hz(i  ,j-1,k)
     &           +  Hz(i-1,j,k)+Hz(i-1,j-1,k) ) )
            VFx(i,j)=Fxl+limiteur_h3(rrm,rr,rrp,vit,Fxl,VFx(i,j))       
          enddo
        enddo
        endif
