C
C
C***********************************************************************
C***********************************************************************
       subroutine RUNDEF
C
C Appele par main.f: lecture de runparam.i et initalisation des variables
C    du run - lecture de la grille verticale.
C
C ======================================================================
C
C -------------------------------------
C Definition de la configuration du run.
C -------------------------------------
C
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

       character*80 clpipo
C
       mflagrun=0
C
C Unit des fichiers output
C
       nundyn=20
       nundynfl=21
       nundyn1d=25
       nunbio=30
       nunbiofl=31
       nunbio1d=35
C
C Lecture de runparam.i
C
       open(90,file='../CONFIG_PHY/tke/inputs/runparam.i',status='old')
       read(90,'(a)') clpipo  ! reading "version" 
       read(90,'(a)') version  ! reading version number
       write(*,*) "version =",version
       read(90,'(a)') clpipo  ! reading "directory sauvegarde - adresse absolue:  cnmlog"
       read(90,'(a)') cnmlog  ! reading ./PHY/SORTIES
       read(90,'(a)') clpipo  ! reading "   directory local de travail: cdir"
       read(90,'(a)') cdir  ! reading "./PHY"
       read(90,'(a)') clpipo
       read(90,'(a)') ctit
       read(90,'(a)') clpipo
       read(90,'(a)') csbtit
       read(90,'(a)') clpipo
       read(90,*) dt
       read(90,'(a)') clpipo
       read(90,*) tdebut0
       read(90,'(a)') clpipo
       read(90,*) tfin0
       write(*,*) "tfin0 ", tfin0
       read(90,'(a)') clpipo
       read(90,*)  timesave0
       read(90,'(a)') clpipo
       read(90,*)  nwrite
       read(90,'(a)') clpipo
       write(*,*) clpipo
       read(90,*) mflagwrite
       read(90,'(a)') clpipo
       read(90,*)  nbrprono
       read(90,'(a)') clpipo
       if (nbrprono.gt.0) then
         do jtr=1,nbrprono
           read(90,'(a)') CNMTRA(jtr)
         enddo
       endif
       read(90,'(a)') clpipo
       read(90,*)  nbrdiag
       nbrbio=nbrprono+nbrdiag
       if (nbrbio.gt.jptract) then
         mflagrun=12
         return
       endif
       read(90,'(a)') clpipo
       if (nbrdiag.gt.0) then
         do jtr=nbrprono+1,nbrbio
           read(90,'(a)') CNMTRA(jtr)
         enddo
       endif
       read(90,'(a)') clpipo
       if (nbrbio.gt.0) then
         read(90,*) (MSTOCK(jtr),jtr=1,nbrbio)
       endif
       read(90,'(a)') clpipo
       if (nbrprono.gt.0) then
         read(90,*) (MBIODIF(jtr),jtr=1,nbrprono)
       endif
       read(90,'(a)') clpipo
       read(90,*)  mbiosed
       read(90,'(a)') clpipo
       read(90,*)  nbrflux
       if (nbrflux.gt.jpfluxt) then
         mflagrun=13
         return
       endif
       read(90,'(a)') clpipo
       if (nbrflux.gt.0) then
         do jf=1,nbrflux
           read(90,'(a)') CNMFLX(jf)
         enddo
       endif
       read(90,'(a)') clpipo
       if (nbrflux.gt.0) then
         read(90,*) (MFLUX(jf),jf=1,nbrflux)
       endif
       read(90,'(a)') clpipo
       read(90,*)  mtke
       read(90,'(a)') clpipo
       read(90,*)  mlength
       read(90,'(a)') clpipo
       read(90,*) mdifft
       read(90,'(a)') clpipo
       read(90,*)  mt
       read(90,'(a)') clpipo
       read(90,*)  ms
       read(90,'(a)') clpipo
       read(90,*)  mtket
       read(90,'(a)') clpipo
       read(90,*)  mtt
       read(90,'(a)') clpipo
       read(90,*)  mst
       read(90,'(a)') clpipo
       read(90,*) nzt
       if (nzt.gt.jpzt) then
         mflagrun=11
         return
       endif
       read(90,'(a)') clpipo
       read(90,'(a)')  cgrille
C
       close(90)

C Definition des directories de travail et des racines des fichiers out
       do jl=120,1,-1
         if (cdir(jl:jl).ne.' ') goto 1
       enddo
 1     ndir=jl
       do jl=ndir,1,-1
         if (cdir(jl:jl).eq.'/') goto 2
       enddo
       ndir0=1
       goto 3
 2     ndir0=jl+1
 3     continue
       do jl=50,1,-1
         if (cnmlog(jl:jl).ne.' ') goto 4
       enddo
 4     nmlog=jl
       croot=
     &    cnmlog(1:nmlog)//cdir(ndir0:ndir)
C MB     &    cnmlog(1:nmlog)//'/'//cdir(ndir0:ndir)//'_'

C MB       nroot=nmlog+1+(ndir-ndir0+1)+1
       nroot=nmlog+1+(ndir-ndir0+1)

C Definition de parametres temporels
       dts=dt*hourm
       timesave=timesave0/hourm
       timeyr=float(int(tdebut0/xyear))+1.
       timeday=tdebut0-(timeyr-1.)*xyear+1.
       npurday=int(daymn/dt+0.001)

C Verification de la coherence des entreees et des dimensions
       isum=mtke+mlength+mdifft+mt+ms
       isumb=0
       do j=1,nbrbio
         isumb=isumb+MSTOCK(j)
       enddo
       isum=isum+mtket+mtt+mts
       do j=1,nbrflux
         isumb=isumb+MFLUX(j)
       enddo

       tdebut=tdebut0*day*hourm
       tfin=tfin0*day*hourm
       ztime=timesave*hourm+0.01
       nsave=int(ztime/dt)
       if (isum.gt.0) then
         if (nsave.lt.1) mflagrun=2
       endif
       if (isumb.gt.0) then
         if (nsave.lt.1) mflagrun=2
       endif
       if (mflagrun.eq.2) return

       jstock=int(tfin0*day/timesave+0.001)
       write(*,*) jptemps
       write(*,*) jstock

       if (jstock.gt.jptemps) then
         mflagrun=4
         return
       endif

C Ecriture de la definition du RUN
       call WRUNPRINT

       return
       end
