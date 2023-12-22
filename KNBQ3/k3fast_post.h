! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_post.h (begin)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! !********************************
! !  Update total mass of water 
! !    volume Dnew 
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
# ifdef NBQ_MASS
     &               *rhobar_nbq(i,j,knew)
# endif
        enddo
      enddo
!$acc end kernels 
# if defined RVTK_DEBUG && defined KNBQ
C$OMP BARRIER
C$OMP MASTER
c       call check_tab3d(thetadiv_nbq,'step3d_fastthetadiv_nbq','rint')
# ifdef K3FAST_SEDLAYERS
      call check_tab3d_sedlay(rho_nbq,'step3d_fastrho_nbq_1',
     &  'r',N_sl+1,N,ondevice=.TRUE.)
# else
      call check_tab3d(rho_nbq,'step3d_fastrho_nbq_1','r',ondevice=.TRUE.)
# endif      
C$OMP END MASTER
# endif    
! !
! !********************************
! !  Update W KNBQ forcing for 
! !  internal mode (rw_nbq ~ rw_nbq_avg1)
! !  Note: here rw_nbq contains qdmw_nbq(m)
! !********************************
! !
# ifdef K3FAST_W
      if (LAST_FAST_STEP) then
!$acc kernels if(compute_on_device) default(present)
        do k=0,N 
          do j=Jstr,Jend              
            do i=Istr,Iend
              rw_nbq(i,j,k)=((qdmw_nbq(i,j,k)
     &                    -rw_nbq(i,j,k))
     &                       /dtfast
#  ifdef KNHINT_3M
     &                       /float(nsdtnbq)
#  endif
#   ifdef K3FAST_C3D_WSF
     &                     -rw_int_nbq(i,j,k)
#   endif     
     &                                  )/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo
!$acc end kernels 
!$acc update host( rw_nbq )         !! iif=last
      endif
# endif
! !
! !********************************
! ! Get filtered RHS terms
! ! and multiply by dx*dy to get 
! ! units of rho*Hz*dx*dy*ru
! !********************************
! !
      if (LAST_FAST_STEP) then
! !
! !********************************
! ! Store average fields AVG1 of rho 
! !    and rhobar
! !********************************
! !     
# ifdef NBQ_MASS
        do j=Jstr,Jend
          do i=Istr,Iend
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
          enddo
        enddo 
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rho_nbq_avg1(i,j,k)=1.d0
     &                            + ( rho_nbq(i,j,k)/Hzr(i,j,k)
     &                               +rho_grd(i,j,k) )  
            enddo
          enddo 
        enddo
# endif /* NBQ_MASS */
! !
! !********************************
! !  Compute average fields AVG2 of RHS NBQ forcing
! !  Note: here ru_nbq_avg2, ru_nbq_2d_old ... are working arrays
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)
# ifdef K3FAST_UV 
       do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend 
!             ru_nbq_avg1(i,j,k)= ru_nbq(i,j,k)
#  ifdef K3FAST_C3D_UVSF
              ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
#   ifndef K3FAST_COUPLING2D
     &                     -rubar(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))   ! CAUTION HERE
#   else
     &                     -ru_ext_nbq_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
#   endif     
#  endif     
#  ifdef K3FAST_C3D_UVFS
              ru_nbq_avg2(i,j,k)=
     &                       ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt  ! CAUTION: use define
#   ifdef K3FAST_C3D_UVSF
     &                   -ru_int_nbq(i,j,k)-(ru_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i-1,j,k))
#   endif
     &                       )*on_u(i,j)*om_u(i,j)
#  endif     
            enddo
          enddo 
        enddo  
        do k=1,N
          do j=JstrV,Jend
            do i=Istr,Iend             
    !          rv_nbq_avg1(i,j,k)= rv_nbq(i,j,k) 
#  ifdef K3FAST_C3D_UVSF
              rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
#   ifndef K3FAST_COUPLING2D
     &                     -rvbar(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
#   else
     &                     -rv_ext_nbq_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
#   endif
#  endif     
#  ifdef K3FAST_C3D_UVFS
              rv_nbq_avg2(i,j,k)=
     &                        ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt ! CAUTION: use define
#   ifdef K3FAST_C3D_UVSF
     &                   -rv_int_nbq(i,j,k)-(rv_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i,j-1,k))
#   endif
     &                        )*on_v(i,j)*om_v(i,j)
#  endif    
            enddo
          enddo 
        enddo
# endif /* K3FAST_UV */
# ifdef K3FAST_W
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
!             rw_nbq_avg1(i,j,k)= rw_nbq(i,j,k)
              rw_nbq_avg2(i,j,k)=
     &                          ((qdmw_nbq(i,j,k)
     &        -rw_nbq_avg2(i,j,k))/dt
#   ifdef K3FAST_C3D_WSF
     &                           -rw_int_nbq(i,j,k)
#   endif
     &                          )*on_r(i,j)*om_r(i,j)
            enddo
          enddo 
        enddo
# endif /* K3FAST_W */
!$acc end kernels
      endif ! LAST_FAST_STEP
! !
! !********************************
! ! Dismiss coupling of KNBQ, NBQ2EXT
! ! & NBQ2INT to debug
! !********************************
! !
# ifdef NBQ_NOCOUPLING 
      ru_nbq      =0.   ! 3D
      rv_nbq      =0.
      ru_nbq_avg2 =0.
      rv_nbq_avg2 =0.
#  ifdef K3FAST 
      rw_nbq      =0.
      rw_nbq_avg2 =0.
#  endif
#  ifdef NBQ_MASS
      rhobar_nbq  =1.
      rho_nbq     =0.
      rho_nbq_avg1=1.
#  endif
# endif /* NBQ_NOCOUPLING */
! !
! !********************************
! ! Exchange KNBQ coupling
! !********************************
! !
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
#  ifdef NBQ_MASS
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        rhobar_nbq_avg1(START_2D_ARRAY))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq_avg1(START_2D_ARRAY,1))
#  endif
! !
! !********************************
! ! ATTENTION exchange !!!!!!! FRANCIS
! !********************************
! !
#  ifdef K3FAST_C3D_UVFS
      if (LAST_FAST_STEP) then
       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                         ru_nbq(START_2D_ARRAY,1))
       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rv_nbq(START_2D_ARRAY,1))
       call exchange_u3d_tile (Istr,Iend,Jstr,Jend,  
     &                         ru_nbq_avg2(START_2D_ARRAY,1))
       call exchange_v3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rv_nbq_avg2(START_2D_ARRAY,1))
      endif
#  endif
! !
! !********************************
! ! ATTENTION exchange !!!!!!! FRANCIS
! !********************************
! !
#  ifdef K3FAST_C3D_WFS
      if (LAST_FAST_STEP) then
       call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rw_nbq(START_2D_ARRAY,0))
       call exchange_w3d_tile (Istr,Iend,Jstr,Jend,  
     &                         rw_nbq_avg2(START_2D_ARRAY,0))
      endif
#  endif 
# endif

# ifdef K3FAST_UV 
      if (LAST_FAST_STEP) then
!$acc update host( ru_int_nbq, ru_nbq_avg2      !! iif=last
!$acc&            ,rv_int_nbq, rv_nbq_avg2 
!$acc&            ,rw_nbq_avg2
!$acc&           ) !async(sync_ruv_nbq_avg2)
! ! ru_int_nbq => pre_step
! ! ru_nbq_avg2 => step3d_uv1
      endif
#endif
! !
! !********************************
! !  Depth-averaged velocity & 
! !  forcing from fast mode
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)
#  define Dstp DUon
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
#  ifndef MVB        
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
#  else
          Dstp(i,j)=zeta(i,j,kstp)+dh_mvb(i,j,kstp2)
#  endif          
        enddo
      enddo
!$acc end kernels
      cff=0.5*dtfast
    !  cff=2.*dtfast
      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)
! !
! !********************************
! ! Compute time averaged fields over 
! !  all short timesteps.
! !
! ! Reset/initialise arrays for averaged fields during the first
! ! barotropic time step; Accumulate averages after that. Include
! ! physical boundary points, but not periodic ghost points or
! ! computation  MPI computational margins.
! !********************************
! !
!$acc kernels if(compute_on_device) default(present)

! ! FA: TBT

# ifndef K3FAST_AVG_CLASSIC
      if (FIRST_2D_STEP) then  
        do j=JstrR,JendR
          do i=IstrR,IendR
            DU_avg1(i,j,nnew)=0.
            DV_avg1(i,j,nnew)=0.
            DU_avg2(i,j)=0.
            DV_avg2(i,j)=0. 
          enddo
        enddo
      endif
# endif
      do j=Jstr,Jend
        do i=IstrU,Iend
# ifdef K3FAST_ZETAW
          DUnew=DU_nbq(i,j) *2. 
# else
          DUnew=( (Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)
     &     +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                                *( rubar(i,j)
     &                                  +rufrc(i,j)
     &                          !   +ru_int2d_nbq(i,j)
     &                                 ))
#  ifdef MASKING
     &                                     *umask(i,j)
#  endif
# endif
# ifdef WET_DRY
     &                                     *umask_wet(i,j)
# endif
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) 
     &                             +cff1*on_u(i,j)*( DUnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#  ifdef WET_DRY
     &  *umask_wet(i,j)
#  endif
# endif
     &                                                   )
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(i,j)=DU_avg2(i,j)+cff2*on_u(i,j)*( DUnew 
#  ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i-1,j))*ust2d(i,j)
#   ifdef WET_DRY
     & *umask_wet(i,j)
#   endif
#  endif
     &                                                   )
# endif /*   K3FAST_AVG_CLASSIC */
        enddo
      enddo 
      do j=JstrV,Jend
        do i=Istr,Iend
# ifdef K3FAST_ZETAW
          DVnew=DV_nbq(i,j) *2.
# else
          DVnew=( (Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)
     &     +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                                *( rvbar(i,j)
     &                                  +rvfrc(i,j)
     &                          !   +rv_int2d_nbq(i,j)
     &                                 ))
#  ifdef MASKING
     &                                     *vmask(i,j)
#  endif
# endif
# ifdef WET_DRY
     &                                     *vmask_wet(i,j)
# endif
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) 
     &                              +cff1*om_v(i,j)*(DVnew
# ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#  ifdef WET_DRY
     & *vmask_wet(i,j)
#  endif
# endif
     &                                                   )
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(i,j)=DV_avg2(i,j)+cff2*om_v(i,j)*( DVnew
#  ifdef MRL_WCI
     &                +(Dnew(i,j)+Dnew(i,j-1))*vst2d(i,j)
#   ifdef WET_DRY
     & * vmask_wet(i,j)
#   endif
#  endif
     &                                                   )
# endif /*   K3FAST_AVG_CLASSIC */
        enddo
      enddo
!$acc end kernels
! !
! !********************************
! ! Apply point sources 
! !  for hydrostatic case
! !********************************
! !
# if defined PSOURCE && !defined K3FAST 
!$acc kernels if(compute_on_device) default(present)
      do is=1,Nsrc 
#  ifdef MPI
        i=Isrc_mpi(is,mynode)
        j=Jsrc_mpi(is,mynode)
#  else
        i=Isrc(is)
        j=Jsrc(is)
#  endif
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            ubar(i,j,knew)=2.*Qbar(is)/( on_u(i,j)
     &                       *(Dnew(i-1,j)+Dnew(i,j)) )
            DU_avg1(i,j,nnew)=Qbar(is)
          else
            vbar(i,j,knew)=2.*Qbar(is)/( om_v(i,j)
     &                       *(Dnew(i,j-1)+Dnew(i,j)) )
            DV_avg1(i,j,nnew)=Qbar(is)
          endif
        endif
      enddo
!$acc end kernels
# endif
! !
! !********************************
! !  Set 2D Momemtum nudging
! !********************************
! !
# if defined M2NUDGING && defined M2CLIMATOLOGY

#  ifdef ZONAL_NUDGING
      if (FIRST_TIME_STEP .or. mod(iic,10).eq.0) then
        if (FIRST_FAST_STEP) then
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubar(START_2D_ARRAY,knew),ubzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbar(START_2D_ARRAY,knew),vbzon)
        endif
      endif
      if (FIRST_TIME_STEP) then
        if (FIRST_FAST_STEP) then 
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   ubclm(START_2D_ARRAY),ubclmzon)
          call zonavg_2d(Istr,Iend,Jstr,Jend,
     &                   vbclm(START_2D_ARRAY),vbclmzon)
        endif
      endif
#  endif /* ZONAL_NUDGING */
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr,Jend
        do i=IstrU,Iend
#  ifdef ZONAL_NUDGING        
          DUnew = dtfast*M2nudgcof(i,j)*(ubclmzon(j)-ubzon(j))
#  else          
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
#  endif
#  ifdef MASKING
     &                 * umask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * umask_wet(i,j)
#  endif
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(i,j)     =DU_avg2(i,j)      +cff2*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
# endif /*   K3FAST_AVG_CLASSIC */    
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
#  if defined ZONAL_NUDGING  
          DVnew = dtfast*M2nudgcof(i,j)*(vbclmzon(j)-vbzon(j)) 
#  else
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))    
#  endif      
#  ifdef MASKING
     &                 * vmask(i,j)
#  endif 
#  ifdef WET_DRY
     &                 * vmask_wet(i,j)
#  endif
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(i,j)     =DV_avg2(i,j)      +cff2*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
# endif /*   K3FAST_AVG_CLASSIC */
        enddo
      enddo
!$acc end kernels      
# endif /* M2NUDGING */
! !
! !********************************
! ! Update BRY if necessary
! !********************************
! !
!# if defined M2_FRC_BRY || defined Z_FRC_BRY
!      if (FIRST_FAST_STEP) then
!!!$acc update device(                    
!!#   ifdef M2_FRC_BRY
!!!$acc&        ubarbry_west, ubarbry_east
!!#   endif
!!#  ifdef Z_FRC_BRY
!!!$acc&       ,zetabry_west, zetabry_east
!!#  endif
!!!$acc&              ) !! iif=1
!      endif
!# endif      
! !
! !********************************
! ! Set boundary conditions and 
! ! compute integral mass flux accross
! ! all open boundaries, if any.
! !********************************
! !
      M2bc_nbq_flag=.false.  
! !           skip wet/dry conditions
! !         & AGRIF
      call u2dbc_tile (Istr,Iend,Jstr,Jend, work) 
      call v2dbc_tile (Istr,Iend,Jstr,Jend, work)
# ifdef WET_DRY
!$acc kernels if(compute_on_device) default(present)
#  ifndef EW_COM_PERIODIC
      if (WESTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask_wet(Istr,j)
#   ifdef MRL_WCI
          ust2d(Istr,j)=ust2d(Istr,j)*umask_wet(Istr,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask_wet(Istr-1,j)
#   ifdef MRL_WCI
          vst2d(Istr-1,j)=vst2d(Istr-1,j)*vmask_wet(Istr-1,j)
#   endif
        END DO
      END IF
      if (EASTERN_EDGE) then
        DO j=Jstr,Jend
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask_wet(Iend+1,j)
#   ifdef MRL_WCI
          ust2d(Iend+1,j)=ust2d(Iend+1,j)*umask_wet(Iend+1,j)
#   endif
        END DO
        DO j=JstrV,Jend
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask_wet(Iend+1,j)
#   ifdef MRL_WCI
          vst2d(Iend+1,j)=vst2d(Iend+1,j)*vmask_wet(Iend+1,j)
#   endif
        END DO
      END IF
#  endif
#  ifndef NS_COM_PERIODIC
      if (SOUTHERN_EDGE) then
        DO i=IstrU,Iend
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask_wet(i,Jstr-1)
#   ifdef MRL_WCI
          ust2d(i,Jstr-1)=ust2d(i,Jstr-1)*umask_wet(i,Jstr-1)
#   endif
        END DO
        DO i=IstrU,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask_wet(i,Jstr)
#   ifdef MRL_WCI
          vst2d(i,Jstr)=vst2d(i,Jstr)*vmask_wet(i,Jstr)
#   endif
        END DO
      END IF
      if (NORTHERN_EDGE) then
        DO i=Istr,Iend
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask_wet(i,Jend+1)
#   ifdef MRL_WCI
          ust2d(i,Jend+1)=ust2d(i,Jend+1)*umask_wet(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask_wet(i,Jend+1)
#   ifdef MRL_WCI
          vst2d(i,Jend+1)=vst2d(i,Jend+1)*vmask_wet(i,Jend+1)
#   endif
        END DO
      END IF
#  endif
!$acc end kernels      
# endif
!
! zeta vill be recomputed via depth-integrated continuity equation
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined  MPI
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(START_2D_ARRAY,knew))
# endif
#ifdef RVTK_DEBUG
C$OMP BARRIER
C$OMP MASTER
       call check_tab2d(ubar(:,:,knew),'ubar step3d_fast','u'
     &    ,ondevice=.TRUE.)
       call check_tab2d(vbar(:,:,knew),'vbar step3d_fast','v'
     &    ,ondevice=.TRUE.)
#endif
         
# ifdef OBC_VOLCONS
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
# endif
! !
! !********************************
! ! Compute fast-time averaged barotropic 
! ! mass fluxes along physical
! ! boundaries.
! !********************************
! !

!$acc kernels if(compute_on_device) default(present)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=(h(Istr-1,j)+zeta(Istr-1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Istr-1,j,knew)
#  endif
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=(h(Iend+1,j)+zeta(Iend+1,j,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(Iend+1,j,knew)
#  endif
        enddo
      endif
# endif /* !EW_PERIODIC */

# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jstr-1)=(h(i,Jstr-1)+zeta(i,Jstr-1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jstr-1,knew)
#  endif
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=Istr-1,IendR
          Dnew(i,Jend+1)=(h(i,Jend+1)+zeta(i,Jend+1,knew))
#  ifdef NBQ_MASS
     &                          *rhobar_nbq(i,Jend+1,knew)
#  endif
        enddo
      endif
# endif /* !NS_PERIODIC */
!$acc end kernels

      cff1=0.5*weight(1,iif)
      cff2=0.5*weight(2,iif)

!$acc kernels if(compute_on_device) default(present)
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
# ifdef MRL_WCI
     &                                              +ust2d(IstrU-1,j)
# endif
     &                                               )*on_u(IstrU-1,j)
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(IstrU-1,j)=DU_avg2(IstrU-1,j)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */      
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Istr-1,j)+Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
# ifdef MRL_WCI
     &                                               +vst2d(Istr-1,j)
# endif
     &                                               )*om_v(Istr-1,j)
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(Istr-1,j)=DV_avg2(Istr-1,j)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */      
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=JstrR,JendR
          cff=(Dnew(Iend+1,j)+Dnew(Iend,j))*(ubar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                          +ust2d(Iend+1,j)
# endif
     &                                          )*on_u(Iend+1,j)
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(Iend+1,j)=DU_avg2(Iend+1,j)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */
        enddo
        do j=JstrV,Jend
          cff=(Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*(vbar(Iend+1,j,knew)
# ifdef MRL_WCI
     &                                              +vst2d(Iend+1,j)
# endif
     &                                              )*om_v(Iend+1,j)
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(Iend+1,j)=DV_avg2(Iend+1,j)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */     
        enddo
      endif
# endif /* !EW_PERIODIC */
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*(ubar(i,Jstr-1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jstr-1)
#  endif
     &                                              )*on_u(i,Jstr-1)
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(i,Jstr-1)=DU_avg2(i,Jstr-1)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */       
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
#  ifdef MRL_WCI
     &                                              +vst2d(i,JstrV-1)
#  endif
     &                                              )*om_v(i,JstrV-1)
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(i,JstrV-1)=DV_avg2(i,JstrV-1)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */     
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=IstrU,Iend
          cff=(Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*(ubar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                              +ust2d(i,Jend+1)
#  endif
     &                                              )*on_u(i,Jend+1)
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DU_avg2(i,Jend+1)=DU_avg2(i,Jend+1)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */       
        enddo
        do i=IstrR,IendR
          cff=(Dnew(i,Jend+1)+Dnew(i,Jend))*(vbar(i,Jend+1,knew)
#  ifdef MRL_WCI
     &                                          +vst2d(i,Jend+1)
#  endif
     &                                          )*om_v(i,Jend+1)
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)+cff1*cff
# ifndef K3FAST_AVG_CLASSIC
          DV_avg2(i,Jend+1)=DV_avg2(i,Jend+1)+cff2*cff
# endif /*   K3FAST_AVG_CLASSIC */     
        enddo
      endif
# endif /* !NS_PERIODIC */
!$acc end kernels
! !
! !********************************
! ! Update Free-slip
! !********************************
! !
# ifdef NBQ_FREESLIP
       if (FIRST_FAST_STEP) then
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=qdmw_nbq(i,j,0)
            enddo
           enddo
        elseif (LAST_FAST_STEP) then
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=qdmw0_nbq(i,j)+qdmw_nbq(i,j,0)
            enddo
           enddo
        else
          do j=JstrV-2,Jend+1
            do i=IstrU-2,Iend+1
              qdmw0_nbq(i,j)=(qdmw0_nbq(i,j)+qdmw_nbq(i,j,0))
     &               /float(ndtfast)
            enddo
           enddo
           
        endif
# endif
! !
! !********************************
! ! FRANCIS: EXCHANG to be tested
! !********************************
! !
#  if defined MRL_WCI && defined WET_DRY
      if (LAST_FAST_STEP) then  
       call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                         ust2d(START_2D_ARRAY))
       call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                         vst2d(START_2D_ARRAY))
      endif
#  endif
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! K3FAST_post.h (end)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !     
