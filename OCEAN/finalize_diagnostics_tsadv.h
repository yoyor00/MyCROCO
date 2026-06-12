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
! This block is inside  itrc,k,j,i loops

              cff1= Hz(i,j,k) / (pm(i,j)*pn(i,j))
     &              * 0.5*(t(i,j,k,nstp,itrc)+t(i,j,k,nnew,itrc))

              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff1
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff1
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff1
              TForc(i,j,k,itrc)=TForc(i,j,k,itrc)*cff1
              Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff1
              TVmix(i,j,k,itrc)=TVmix(i,j,k,itrc)*cff1
              THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff1
!              TVmixt(i,j,k,itrc)=TVmixt(i,j,k,itrc)*cff

