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
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Istr-1,j)=hwrk(Istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=J_EXT_RANGE
          hwrk(Iend+1,j)=hwrk(Iend,j)
        enddo
      endif
#endif
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jstr-1)=hwrk(i,Jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=I_EXT_RANGE
          hwrk(i,Jend+1)=hwrk(i,Jend)
        enddo
      endif
# ifndef EW_PERIODIC
      if (WESTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Istr-1,Jstr-1)=hwrk(Istr,Jstr)
      endif
      if (WESTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Istr-1,Jend+1)=hwrk(Istr,Jend)
      endif
      if (EASTERN_EDGE.and.SOUTHERN_EDGE) then
        hwrk(Iend+1,Jstr-1)=hwrk(Iend,Jstr)
      endif
      if (EASTERN_EDGE.and.NORTHERN_EDGE) then
        hwrk(Iend+1,Jend+1)=hwrk(Iend,Jend)
      endif
# endif
#endif

      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk(i,j)=0.25*(hwrk(i,j)  +hwrk(i-1,j)
     &                  +hwrk(i,j-1)+hwrk(i-1,j-1))
#ifdef MASKING
     &                   *pmask2(i,j)
#endif
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
#ifdef MASKING
          cff=0.25*(pmask2(i,j)   +pmask2(i+1,j)
     &             +pmask2(i,j+1) +pmask2(i+1,j+1))
#else
          cff=1.
#endif
          hwrk(i,j)=(1.-cff)*hwrk(i,j)+
     &              0.25*(wrk(i,j)  +wrk(i+1,j)
     &                   +wrk(i,j+1)+wrk(i+1,j+1))
#ifdef MASKING
          hwrk(i,j)=hwrk(i,j)*rmask(i,j)
#endif
        enddo
      enddo
