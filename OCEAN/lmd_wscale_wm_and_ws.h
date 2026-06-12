!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
! PART OF KPP2005 (Shchepetkin et al. 2005)
!
            if (Bfsfc.lt.0.) zscale=min(zscale,
     &                       my_hbl(i,j)*epssfc)
#ifdef MASKING
            zscale=zscale*rmask(i,j)
#endif
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
!
! Stable regime.
!
            if (zetahat.ge.0.) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5.*zetahat,
     &                                                    1.E-20)
              ws=wm
!
! Unstable regime: note that zetahat is always negative here, also
! negative are constants "zeta_m" and "zeta_s".
!
            else
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16.*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16.*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif

#ifdef LMD_LANGMUIR
!
! Enhanced turbulent velocity scale due to Langmuir turbulence
!
            cff1=max(eps,Langmuir(i,j))
            cff=sqrt(1+0.104/cff1**2+0.034/cff1**4)   ! Van Roekel et al. (2012)
            wm=wm*cff
            ws=ws*cff
#endif

