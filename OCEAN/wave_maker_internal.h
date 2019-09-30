!
!====================================================================
!                        Internal wave maker
!
! Mode-1 internal wave (e.g., Hall et al., JPO 2013)
!====================================================================
!
!--------------------------------------------------------------------
!  Configurations
!--------------------------------------------------------------------
!
        wp=12.4*3600           ! period
        wa=0.1                 ! amplitude
!
!--------------------------------------------------------------------
!  Initialisation
!--------------------------------------------------------------------
!
        ramp=tanh(dt/wp*float(iic-ntstart))
        wa=wa*ramp
        wf=2*pi/wp
!
!--------------------------------------------------------------------
!  Sea level zetabry
!--------------------------------------------------------------------
!
#   ifdef Z_FRC_BRY
        do j=JstrR,JendR
          zetabry_west(j)=0.
        enddo
#   endif /* Z_FRC_BRY */
!
!--------------------------------------------------------------------
!  XI velocity components ubry and ubarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrR,JendR
          h0=0.5*(h(0,j)+h(1,j))
          do k=1,N
            Zu=h0+0.5*(z_r(0,j,k)+z_r(1,j,k))
            ubry_west(j,k)=wa*cos(pi*Zu/h0)*sin(wf*time)
          enddo
        enddo
#   endif /* M3_FRC_BRY */

#   ifdef M2_FRC_BRY
        do j=JstrR,JendR
          cff4=0.
          cff5=0.
          do k=1,N
            cff4=cff4+ubry_west(j,k)*(Hz(0,j,k)+Hz(1,j,k))
            cff5=cff5+(Hz(0,j,k)+Hz(1,j,k))
          enddo
          ubarbry_west(j)=cff4/cff5
        enddo
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  ETA velocity components vbry and vbarbry
!--------------------------------------------------------------------
!
#   ifdef M3_FRC_BRY
        do j=JstrV,JendR
          do k=1,N
              vbry_west(j,k)=0.
          enddo
        enddo
#   endif /* M3_FRC_BRY */
#   ifdef M2_FRC_BRY
        do j=JstrV,JendR
          vbarbry_west(j)=0.
        enddo
#   endif /* M2_FRC_BRY */
!
!--------------------------------------------------------------------
!  Z velocity component wbry
!--------------------------------------------------------------------
!
#   ifdef W_FRC_BRY
        do j=JstrR,JendR
          h0=h(0,j)
          do k=1,N
            Zr=h0+z_r(0,j,k)
            wbry_west(j,k)=0. 
!           wbry_west(j,k)=wa*sin(pi*Zr/h0)*cos(wf*time)
!     &                      *sqrt(wf**2/(bvf(0,j,k)-wf**2))

          enddo
        enddo
#   endif /* W_FRC_BRY */
!
!--------------------------------------------------------------------
!  TRACERS tbry
!--------------------------------------------------------------------
!
#   ifdef T_FRC_BRY
        if (FIRST_TIME_STEP) then
          do k=1,N
            do j=JstrR,JendR
              do itrc=1,NT
                tbry_west(j,k,itrc)=t(0,j,k,1,itrc)
              enddo
            enddo
          enddo
        endif
#   endif
#   ifdef T_FRC_BRY0
        cff=2.e-4     ! thermal expansion coefficient °C-1
        do j=JstrR,JendR
          h0=h(0,j)
          do k=1,N
            Zr=h0+z_r(0,j,k)
            tbry_west(j,k,itemp)=24.+z_r(0,j,k)*0.06
     &                            -wa*sin(pi*Zr/h0)*sin(wf*time)
     &                               *sqrt(bvf(0,j,k))/(g*cff)
          enddo
        enddo
#   endif

