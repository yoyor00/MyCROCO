PROGRAM my_layout
   !!---------------------------------------------------------------------
   !!
   !!                       PROGRAM MPP_OPTIMIZ_NC
   !!                     ***********************
   !!
   !!  PURPOSE :
   !!  ---------
   !!              This program is build to optimize the domain beakdown into
   !!              subdomain for mpp computing.
   !!              Once the grid size, and the land/sea mask is known, it looks
   !!              for all the possibilities within a range of setting parameters
   !!              and determine the optimal.
   !!
   !!              Optimization is done with respect to the maximum number of
   !!              sea processors and to the maximum numbers of procs (max_number_proc)
   !!                     
   !!
   !! history:
   !! --------
   !!       original  : 95-12 (Imbard M) for OPA8.1, CLIPPER
   !!       f90       : 03-06 (Molines JM), namelist as input
   !!                 : 05-05 (Molines JM), bathy in ncdf
   !!                 : 18-05 (Benshila R), adaptation for CROCO 
   !!----------------------------------------------------------------------
   !! * modules used
    USE netcdf

    IMPLICIT NONE

    INTEGER ::  max_number_proc=250   !: maximum number of proc. (Read from namelist)
       !
    INTEGER ::  &
         xi_rho  ,    & !: I-size of the model (namelist)
         eta_rho ,    & !: J-size of the model (namelist)
         Npts    =  2,    & !: number of ghost cells
         numnam  =  4       !: logical unit for the namelist
    INTEGER :: min_nb_proc_NP_XI,min_nb_proc_NP_ETA
    NAMELIST /namproc/ min_nb_proc_NP_XI,min_nb_proc_NP_ETA,max_number_proc
    INTEGER :: nb_proc_NP_XI,nb_proc_NP_ETA
    NAMELIST /namlayout/ nb_proc_NP_XI,nb_proc_NP_ETA

    INTEGER ::  jpnix ,jpnjx  
    REAL(kind=8) :: xlen,ylen

    !

    CHARACTER(LEN=80) :: cbathy, &       !: File name of the netcdf bathymetry (namelist)
        &                clvar           !: Variable name in netcdf for the bathy to be read
    CHARACTER(LEN=80) :: covdta, cdum
    NAMELIST /namfile/ cbathy, covdta
    
    INTEGER :: ji,jj,jni,jnj,jni2,jnj2
    INTEGER :: imoy,isurf,ivide
    INTEGER :: in
    INTEGER :: ipi,ipj
    INTEGER :: inf10,inf30,inf50,iptx,isw
    INTEGER :: iii,iij,iiii,iijj,iimoy,iinf10,iinf30,iinf50
    !
    INTEGER,DIMENSION(:,:),ALLOCATABLE     ::  ippdi, ippdj ,iidom, ijdom
    INTEGER ::  LLm, MMm, NP_XI,NP_ETA
    !
    REAL(KIND=4)                           ::  zmin,zmax,zper
    REAL(KIND=4)                           ::  zzmin,zzmax,zperx
    REAL(KIND=4),DIMENSION(:,:),ALLOCATABLE  ::  zmask ! xi_rho - eta_rho
    REAL(KIND=4),DIMENSION(:,:),ALLOCATABLE  ::  ztemp ! xi_rho - eta_rho


   ! CDF stuff
    INTEGER :: ncid, ivarid, dimid, istatus
    LOGICAL :: llbon=.FALSE.

    INTEGER :: jjc
    INTEGER :: chunk_size_X, margin_X, chunk_size_E, margin_E
    INTEGER :: Istrmpi, Iendmpi, Jstrmpi, Jendmpi, i_X, j_E
    INTEGER, DIMENSION(:), ALLOCATABLE :: nldi,nlei, nldj,nlej,icount
    INTEGER, DIMENSION(:), ALLOCATABLE :: nleiv, nldiv,nlejv,nldjv
    !
    ! 0. Initialisation
    ! -----------------
    OPEN(numnam,FILE='namelist')
!    REWIND(numnam)
!    READ(numnam,namspace)

    REWIND(numnam)
    READ(numnam,namfile)

    !REWIND(numnam)
    !READ(numnam,namparam)

    REWIND(numnam)
    READ(numnam,namproc)

    REWIND(numnam)
    READ(numnam,namlayout)

    
    jpnix = max_number_proc ; jpnjx= max_number_proc

    ALLOCATE (ippdi(jpnix,jpnjx), ippdj(jpnix,jpnjx) )
    ALLOCATE (iidom(jpnix,jpnjx), ijdom(jpnix,jpnjx) )
    ALLOCATE (nlei(max_number_proc), nldi(max_number_proc) )
    ALLOCATE (nlej(max_number_proc), nldj(max_number_proc) )
    ! empty processors
    ALLOCATE (nleiv(max_number_proc), nldiv(max_number_proc) )
    ALLOCATE (nlejv(max_number_proc), nldjv(max_number_proc) )
    ALLOCATE (ICOUNT(max_number_proc) ) 

    !
    ! * Read cdf mask file
    !
    clvar = 'mask_rho'  

    INQUIRE( FILE=cbathy, EXIST=llbon )
    IF( llbon ) THEN
       istatus=NF90_OPEN(cbathy,NF90_NOWRITE,ncid)
       istatus =NF90_OPEN(cbathy,NF90_NOWRITE,ncid)
       istatus = nf90_inq_dimid(ncid,'xi_rho',dimid)
       istatus = nf90_inquire_dimension(ncid,dimid,len=xi_rho)
       istatus = nf90_inq_dimid(ncid,'eta_rho',dimid)
       istatus = nf90_inquire_dimension(ncid,dimid,len=eta_rho) 
      ALLOCATE (zmask(0:xi_rho-1,0:eta_rho-1))
       ALLOCATE (ztemp(xi_rho,eta_rho))
       istatus = NF90_INQ_VARID(ncid,clvar,ivarid)
       istatus = NF90_GET_VAR(ncid,ivarid,ztemp)
       istatus = NF90_CLOSE(ncid)
    ELSE
        PRINT *,' File missing : ', TRIM(cbathy)
        STOP
    ENDIF

    !
    zmask(0:xi_rho-1,0:eta_rho-1)=ztemp
    DO jj=0,eta_rho-1
       DO ji=0,xi_rho-1
          zmask(ji,jj)=  MIN(REAL(1.,kind=4),MAX(REAL(0.,kind=4),zmask(ji,jj)))  ! Old vector coding rule ...
       END DO
    END DO
 
    PRINT *,'Number of pts     :', eta_Rho*eta_rho
    PRINT *,'Number of sea pts :', INT(SUM(zmask))
    PRINT *

    !
    !  0. Main loop on all possible combination of processors up to max_number_proc
    ! --------------------------------------------------------------------
    iii=1 ; iij=1
    iiii=xi_rho ; iijj=eta_rho
    iptx=0
    iimoy=0
    zzmin=0. ; zzmax=0.
    iinf10=0 ; iinf30=0 ; iinf50=0
    zperx=1.
    in=0
   
    LLm = xi_rho -2
    Mmm = eta_rho-2

    
    jni=nb_proc_NP_XI
    jnj=nb_proc_NP_ETA
        
    !  1. Global characteristics of the jni x jnj decomposition
    ! ---------------------------------------------------------
    !

    ! Limitation of the maximum number of PE's
    print *,'jni,jnj,max_number_proc : ',jni,jnj,max_number_proc
    IF(jni*jnj > max_number_proc) goto 1000
    !
    NP_XI=jni
    NP_ETA=jnj        
    chunk_size_X=(LLm+NP_XI-1)/NP_XI      
    chunk_size_E=(MMm+NP_ETA-1)/NP_ETA

    ! we requiere number of interior points > 3*Ghostcells : WHY +1 ???
    ! => to avoid too small domains at the boundaries
    IF (chunk_size_X < 2*Npts .OR. chunk_size_E < 2*Npts) go to 1000

    ipi=chunk_size_X+2*Npts  ! Interior + Ghost cells
    ipj=chunk_size_E+2*Npts   

!    IF(((ipi+2)/2 - (ipi+1)/2 )> 0)  go to 1000
!    IF(((ipj+2)/2 - (ipj+1)/2) > 0)  go to 1000

    ! Memory optimization ?
    isw=0
    IF(isw.EQ.1) go to 1000
    in=in+1
    !
    zper=(jni*jnj*ipi*ipj)/float(xi_rho*eta_rho)
    !
    ivide=0
    imoy=0
    zmin=1.e+20
    zmax=-1.e+20
    inf10=0
    inf30=0
    inf50=0
   
    !  2. Loop on the CPUS : Compute mpi stuff for each given decomposition
    ! -----------------------------------------------------------------------
    !
    xlen=xi_rho ; ylen= eta_rho
    DO jj=1,jnj
       DO ji=1,jni
    !DO jj=1,nb_proc_NP_ETA
    !   DO ji=1,nb_proc_NP_XI
          j_E=jj-1
          i_X=ji-1
          margin_X=(NP_XI*chunk_size_X-Llm)/2
          istrmpi=1+i_X*chunk_size_X-margin_X
          iendmpi=istrmpi+chunk_size_X-1
          istrmpi=MAX(istrmpi,1)
          iendmpi=MIN(iendmpi,LLm)
          ! 
          margin_E=(NP_ETA*chunk_size_E-MMm)/2
          jstrmpi=1+j_E*chunk_size_E-margin_E
          jendmpi=jstrmpi+chunk_size_E-1             
          jstrmpi=MAX(jstrmpi,1)   
          jendmpi=MIN(jendmpi,Mmm) 

          ! security chack, maybe useless by construction
          !if (margin_X >=chunk_size_X) go to 1000               
          !if (margin_E >=chunk_size_E) go to 1000

           xlen= min(xlen,real(iendmpi-istrmpi+1))
           ylen= min(ylen,real(jendmpi-jstrmpi+1))
           if(xlen<Npts ) go to 1000    
           if(ylen<Npts ) go to 1000    
          

          ! Check wet points over the entire domain to preserve the MPI communication stencil ???????
          isurf=0
          DO jnj2=Max(jstrmpi-Npts,1),Min(jendmpi+Npts,Mmm)
             DO  jni2=Max(istrmpi-Npts,1),Min(iendmpi+Npts,LLm)
                IF(zmask(jni2,jnj2).EQ.1.) isurf=isurf+1
             END DO
          END DO

          IF(isurf.EQ.0) THEN
             ivide=ivide+1
           ELSE
             imoy=imoy+isurf
            ENDIF
          zper=float(isurf)/float(ipi*ipj)   ! additional points for ghost cells
          IF(zmin.GT.zper.AND.isurf.NE.0) zmin=zper
          IF(zmax.LT.zper.AND.isurf.NE.0) zmax=zper
          IF(zper.LT.0.1.AND.isurf.NE.0) inf10=inf10+1
          IF(zper.LT.0.3.AND.isurf.NE.0) inf30=inf30+1
          IF(zper.LT.0.5.AND.isurf.NE.0) inf50=inf50+1
          !
          !
          ! 3. End of the loop on the CPUS, print
          ! ------------------------------------------------
          !
       END DO
    END DO
    zper=float((jni*jnj-ivide))*float(ipi*ipj)/float(xi_rho*eta_rho)

    ! 
    ! 4. Optimum search
    ! -------------------------
    !
    print*,'ivide,iptx : ',ivide,iptx
    IF(ivide.GT.iptx) THEN
       iii=jni
       iij=jnj
       iiii=ipi
       iijj=ipj
       iptx=ivide
       iimoy=imoy
       zzmin=zmin
       zzmax=zmax
       iinf10=inf10
       iinf30=inf30
       iinf50=inf50
       zperx=zper
    ELSE IF(ivide.EQ.iptx.AND.zperx.LT.zper) THEN
       iii=jni
       iij=jnj
       iiii=ipi
       iijj=ipj
       iimoy=imoy
       zzmin=zmin
       zzmax=zmax
       iinf10=inf10
       iinf30=inf30
       iinf50=inf50
       zperx=zper
    ENDIF
    !
    ! 5. End of loop on all possible decomposition
    ! --------------------------------------------
    !
  1000 continue

    !
    ! 6. loop on optimal cpus (iii x jjj) for plotting purposes
    ! ---------------------------------------------------------
    !
    jjc=0
    ivide=0
    imoy=0
    print *,'iij',iij
    print *,'iii',iii
    DO jj=1,iij
      DO ji=1,iii
          j_E=jj-1
          i_X=ji-1
          NP_XI=iii
          NP_ETA=iij
          chunk_size_X=(LLm+NP_XI-1)/NP_XI      
          margin_X=(NP_XI*chunk_size_X-Llm)/2
          istrmpi=1+i_X*chunk_size_X-margin_X
          iendmpi=istrmpi+chunk_size_X-1
          istrmpi=MAX(istrmpi,1)
          iendmpi=MIN(iendmpi,LLm)
          ! 
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
          margin_E=(NP_ETA*chunk_size_E-MMm)/2
          jstrmpi=1+j_E*chunk_size_E-margin_E
          jendmpi=jstrmpi+chunk_size_E-1             
          jstrmpi=MAX(jstrmpi,1)   
          jendmpi=MIN(jendmpi,Mmm) 
          
          ! Check wet points over the entire domain to preserve the MPI communication stencil
          isurf=0
          DO jnj2=Max(jstrmpi-Npts,1),Min(jendmpi+Npts,Mmm)
             DO  jni2=Max(istrmpi-Npts,1),Min(iendmpi+Npts,LLm)
                 IF(zmask(jni2,jnj2).EQ.1.) isurf=isurf+1
             END DO
          END DO

          IF (isurf.EQ.0) THEN
             ivide=ivide+1
             nldiv(ivide)=istrmpi
             nleiv(ivide)=iendmpi
             nldjv(ivide)=jstrmpi
             nlejv(ivide)=jendmpi
          ELSE
             imoy=imoy+isurf
             jjc=jjc+1
             icount(jjc)=isurf 
             nldi(jjc)=istrmpi
             nlei(jjc)=iendmpi
             nldj(jjc)=jstrmpi
             nlej(jjc)=jendmpi
          ENDIF
          if ( iendmpi ==-1) STOP !return
          !
          !
          ! End of the loop on the optimal CPUS
          ! -----------------------------------
          !
      END DO
    END DO


    !
    ! 8. Write optimum in a file
    ! --------------------------
    !
    !WRITE(cdum,'(a,1h-,i3.3,1hx,i3.3,1h_,i4.4)') TRIM(covdta),iii,iij,iii*iij-iptx
    print*,TRIM(covdta),iii,iij,iii*iij-iptx
    WRITE(cdum,'(a,1h-,i4.4,1hx,i4.4,1h_,i4.4)') TRIM(covdta),iii,iij,iii*iij-iptx
    OPEN (20,file=cdum)
    WRITE(20,'(a,i5)')'#',iii*iij -iptx
    DO jjc=1,iii*iij-iptx
       WRITE(20,'(a,i5)')'#',jjc
       WRITE(20,'(2i6)')nldi(jjc),nldj(jjc)
       WRITE(20,'(2i6)')nlei(jjc),nldj(jjc)
       WRITE(20,'(2i6)')nlei(jjc),nlej(jjc)
       WRITE(20,'(2i6)')nldi(jjc),nlej(jjc)
       WRITE(20,'(2i6)')nldi(jjc),nldj(jjc)
       WRITE(20,'(2i6)') 9999,9999
       ! Write warnings on standard entry if very few water points (arbitrary <20)
       ! We could test the ratio instead
       IF (icount(jjc).LT.20) THEN
          DO jj=nldj(jjc),nlej(jjc)
          ENDDO
  900           FORMAT(1x,i4,1x,9(10i1,1x))
       ENDIF
    END DO
    WRITE(20,'(a,i5)')'# vides:',iptx
    DO jjc=1,iptx
       WRITE(20,'(a,i5)')'# vide ',jjc
       WRITE(20,'(2i6)')nldiv(jjc),nldjv(jjc)
       WRITE(20,'(2i6)')nleiv(jjc),nldjv(jjc)
       WRITE(20,'(2i6)')nleiv(jjc),nlejv(jjc)
       WRITE(20,'(2i6)')nldiv(jjc),nlejv(jjc)
       WRITE(20,'(2i6)')nldiv(jjc),nldjv(jjc)
       WRITE(20,'(2i5)') 9999,9999
    END DO
    CLOSE(20)
    !
    STOP
    !
END PROGRAM my_layout
