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
      integer imin,imax,ishft, jmin,jmax,jshft
#ifdef EW_PERIODIC
      if (NP_XI.eq.1) then                ! This means that if there
        imin=Istr-Npts                    ! is no partition in XI-
        imax=Iend+Npts                    ! direction, then periodic
      else                                ! margins are included into
        imin=Istr                         ! the message;
        imax=Iend                         ! otherwise strip them out.
      endif
#else
      if (ii.eq.0 .and. Istr.eq.1) then   ! Extra point on either
        imin=Istr-1                       ! side to accomodate ghost
      else                                ! points associated with
        imin=Istr                         ! PHYSICAL boundaries.
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
#endif
      ishft=imax-imin+1

#ifdef NS_PERIODIC
      if (NP_ETA.eq.1) then               ! This means that if there
        jmin=Jstr-Npts                    ! is no partition in ETA-
        jmax=Jend+Npts                    ! direction, then periodic
      else                                ! margins are included into
        jmin=Jstr                         ! the message;
        jmax=Jend                         ! otherwise strip them out.
      endif
#else
      if (jj.eq.0 .and. Jstr.eq.1) then   ! Extra point on either
        jmin=Jstr-1                       ! side to accomodate ghost
      else                                ! points associated with
        jmin=Jstr                         ! PHYSICAL boundaries.
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
#endif
      jshft=jmax-jmin+1
