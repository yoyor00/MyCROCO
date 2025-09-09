C     
C     
C***********************************************************************
C***********************************************************************
      subroutine BIODEF 
      Use comrunmod
      Use comdynmod
      USE mod_eco3m
      use mod_varphy_coupl
      
      Integer::ivar,ii
      character(200) :: clpipo
      real(8) :: vitsedfile

C     Lecture des parametres du modele biologique
      open(90,file=cdir(1:ndir)//'/inputs/bioparam.i',status='old')

C     Photo - phyto
      read(90,'(a)') clpipo
      read(90,*) teneurmin
      read(90,'(a)') clpipo
      read(90,*) xslopet
      read(90,'(a)') clpipo
      read(90,*) xzelight
      read(90,'(a)') clpipo
      read(90,*) flimlm0
      read(90,'(a)') clpipo
      read(90,*) chlaphea
      read(90,'(a)') clpipo
      read(90,*) cchlmin
      read(90,'(a)') clpipo
      read(90,*) cchlmax
      read(90,'(a)') clpipo
      read(90,*) xpurmx
      read(90,'(a)') clpipo
      read(90,*) ftimenut0
      read(90,'(a)') clpipo
      read(90,*) ftimenut1

C     Sedimentation
      read(90,'(a)') clpipo
      write(*,*) trim(clpipo)
      VITSEDfile = 0.0
      vitsed0 = 0.0
      read(90, *) vitsedfile
      write(*,*) 'vitsedfile =',vitsedfile
      do ivar = 1,nbvar
         if (VAR(ivar)%comp == 'mop') then 
            if (VAR(ivar)%scomp == 'MOPS') then 
                   vitsed0(ivar) = VITSEDfile
            elseif (VAR(ivar)%scomp == 'MOPL') then
                   vitsed0(ivar) = VITSEDfile * 25.d0  !100
            endif
!                   vitsed0(ivar) = VITSEDfile
         endif
      enddo

C     Passage de m/jour en SI (m/sec)
C MB correction bug      do jtr=1,nbrprono
      do jtr=1,nbvar
         VITSED(jtr)=VITSED0(jtr)/(day*hour)
         write(*,*) 'jtr,vitsed',jtr,vitsed(jtr)
      enddo

C     Rappels NO3 DOR - Initialisation     
      read(90,'(a)') clpipo
      read(90,'(a)') coptique
      read(90,'(a)') clpipo
      read(90,*) mrapbio
      read(90,'(a)') clpipo
      read(90,'(a)') crapbio
      read(90,'(a)') clpipo
      read(90,'(a)')  cbioinit
C     
      close(90)
C     
C     Clha: variabile pronostique ?
C     
      nchlapro=0
      do jt=1,nbrprono
         if (jt.eq.nicchl) nchlapro=1
      enddo
C     
C     Eta = uptake ammonium / uptake DOM labile par bacteries
C     
C     MB       eta=REDCN(nidom)/REDCN(nibac)-1.
C     
      return
      end
