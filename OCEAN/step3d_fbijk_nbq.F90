#include "cppdefs.h"
#if defined NBQ && defined NBQ_IJK
!
!======================================================================
!                      NBQ-Mode for NH-modeling
!                            Main Routine
!======================================================================
!
!> @note Main web documentation at:
! http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm
!
! DESCRIPTION: 
!
!> @brief SNBQ driver : Non-hydrostatic algorithm with the 
!                       Non-boussinesq solver.
!
!> @details NBQ time step. See SNBQ web pages :
!  http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/Restricted/NH-NBQ/Html_pages
!    Algorithme_NBQ.htm      --> SNBQ algorithm:               
!    Algebrique_SNBQ.htm     --> SNBQ algebric representation
!    Couplage_Numerique.htm  --> Numerical coupling
!    
!    Couplage_Split_NBQ.htm  --> Coupling Splitting
!
! REVISION HISTORY:
!
!> @authors
!> @date 2015 January
!> @todo
!
!======================================================================
!
      subroutine step3d_fbijk_nbq (Istr,Iend,Jstr,Jend, WORK)

      use module_nh 
      use module_nbq

      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"
# ifdef MPI
      include 'mpif.h'
# endif
      integer Istr,Iend,Jstr,Jend,i,j,k
      real    WORK(PRIVATE_2D_SCRATCH_ARRAY)
      real :: dum_s
      double precision :: a_m,b_m
# ifndef MPI
      integer mynode
      mynode=0
# endif
# undef DEBUG

!
!-------------------------------------------------------------------
!       Initialization of various test-cases
!-------------------------------------------------------------------
!       
        if (iif==1.and.iic==1) call initial_nh (3)
!
!-------------------------------------------------------------------
!  Get internal and external forcing terms for nbq equations:
!  ru+rubar (or lambda_ext+lambda_int)
!  dzdt*rhosurf
!-------------------------------------------------------------------
!
!     call ru_nbq(1)     ! MPI ! TBD
!
!------------------------------------------------------------------
!       Implicit part: system setup
!-------------------------------------------------------------------
!
# ifdef NBQ_IMP
      if (iif.eq.1.and.ifl_imp_nbq.eq.1) call implicitijk_nbq (1)
# endif
!
!*******************************************************************
!*******************************************************************
!              NBQ mode iteration (main loop)
!*******************************************************************
!*******************************************************************

      do iteration_nbq=1,iteration_nbq_max

!
!-------------------------------------------------------------------
!
!-------------------------------------------------------------------

        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh
              div_nbq(i,j,k)=-visc2_nbq*div_nbq(i,j,k)     &
                             +soundspeed2_nbq*rho_nbq(i,j,k)
            enddo
          enddo
        enddo

!
!-------------------------------------------------------------------
!      Horizontal Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!-------------------------------------------------------------------
!
!  XI- and ETA-Directions:    ETA-directions  ====> TBD
!
        do j=JstrU_nh,JendU_nh
          do i=IstrU_nh,IendU_nh

            k=1
            dum_s =       - (0.5*(Hzr_half_nbq(i-1,j,k)+Hzr_half_nbq(i,j,k)))* pm_u(i,j)        &
                            *(div_nbq(i,j,k) - div_nbq(i-1,j,k))                                &
                          + gdepth_u(i,j,k) * coefb_u(i  ,j,k) * ( div_nbq(i  ,j,k+1) -  div_nbq(i  ,j,k)   &
                                                                 + div_nbq(i-1,j,k+1) -  div_nbq(i-1,j,k) ) &
                          - gdepth_u(i,j,k) * coefa_u(i  ,j,k) * ( div_nbq(i  ,j,k-1) + div_nbq(i-1,j,k-1) )    

            rhssumu_nbq(i,j,k) = rhssumu_nbq(i,j,k) + dum_s
            qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k)                                                   & 
                            + dtnbq * ( dum_s +  rho0*(ruint_nbq(i,j,k)+ruext_nbq(i,j,k)))

            do k=2,N-1
               dum_s =       - (0.5*(Hzr_half_nbq(i-1,j,k)+Hzr_half_nbq(i,j,k))) * pm_u(i,j)    &
                               *(div_nbq(i,j,k) - div_nbq(i-1,j,k))                             &
                             + gdepth_u(i,j,k) * coefb_u(i  ,j,k) * ( div_nbq(i  ,j,k+1) - div_nbq(i  ,j,k)     &
                                                                    + div_nbq(i-1,j,k+1) - div_nbq(i-1,j,k) )   &
                             + gdepth_u(i,j,k) * coefa_u(i  ,j,k) * ( div_nbq(i  ,j,k  ) - div_nbq(i  ,j,k-1)   &  
                                                                    + div_nbq(i-1,j,k  ) - div_nbq(i-1,j,k-1) ) 

              rhssumu_nbq(i,j,k) = rhssumu_nbq(i,j,k) + dum_s
              qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k)                                                   & 
                              + dtnbq * ( dum_s +  rho0*(ruint_nbq(i,j,k)+ruext_nbq(i,j,k)))
            enddo

            k=N
            dum_s =       - (0.5*(Hzr_half_nbq(i-1,j,k)+Hzr_half_nbq(i,j,k))) * pm_u(i,j)        &
                            *(div_nbq(i,j,k) - div_nbq(i-1,j,k))                                 &
                          + gdepth_u(i,j,k  ) * coefa_u(i  ,j,k  ) * ( div_nbq(i  ,j,k  ) - div_nbq(i  ,j,k-1)   &  
                                                                     + div_nbq(i-1,j,k  ) - div_nbq(i-1,j,k-1) ) &
                          - gdepth_u(i,j,k+1) * coefb_u(i  ,j,k+1) * ( div_nbq(i  ,j,k  ) + div_nbq(i-1,j,k) )     
            rhssumu_nbq(i,j,k) = rhssumu_nbq(i,j,k) + dum_s
            qdmu_nbq(i,j,k) = qdmu_nbq(i,j,k)                                                   & 
                            + dtnbq * ( dum_s +  rho0*(ruint_nbq(i,j,k)+ruext_nbq(i,j,k)))

          enddo
        enddo
!
!  U-momentum open boundary conditions
!
# ifdef OBC_NBQ
!       call unbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)   ! TBD
!       call vnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)   ! TBD
# endif
          
!
!  Message passing: Send U (51) 
!
!       call parallele_nbq(51)   ! TBD
!
!  Message passing: Send V (52) 
!
!       call parallele_nbq(52)   ! TBD
!
!-------------------------------------------------------------------
!      Vertical Momentum equation: 
!         If explicit: (x,y,z) is dealt with here
!         If implicit: (x,y)   only
!-------------------------------------------------------------------
!
# ifndef NBQ_IMP
!
!  Z-Direction: Explicit
!
        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh                                                               
               dum_s =   float(mijk2lq_nh(i,j,k))                                          &
                      * (div_nbq(i,j,k) - div_nbq(i,j,k+1) * float(mijk2lq_nh(i,j,k+1)))
               rhssumw_nbq(i,j,k) = rhssumw_nbq(i,j,k) + dum_s                              
               qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)   &
                + dtnbq * ( dum_s + rho0 * rwint_nbq(i,j,k) )
             enddo             
           enddo
        enddo
        
# else

!       call parallele_nbq(151)  ! u only  ! TBD 
!       call parallele_nbq(152)  ! v only  ! TBD

        call implicitijk_nbq (2)    ! TBD

# endif
!
!      Vertical momentum open boundary conditions
!
# ifdef OBC_NBQ
!       call wnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!-------------------------------------------------------------------
!      Message passing 
!-------------------------------------------------------------------
!
!  Send
!       call parallele_nbq(53)    ! TBD
!
!  Receive
# ifndef NBQ_IMP
!       call parallele_nbq(151)   ! TBD 
!       call parallele_nbq(152)   ! TBD  
# endif 
!       call parallele_nbq(153)   ! TBD   
!
!-------------------------------------------------------------------
!      Message passing
!-------------------------------------------------------------------
!
!  Send
!       call parallele_nbq(7)     ! TBD

!  Receive
!       call parallele_nbq(17)    ! TBD

!
!-------------------------------------------------------------------
!      Acoustic wave emission
!-------------------------------------------------------------------
!
#  ifdef ACOUSTIC
!      call density_nbq(11)       ! TBD
#  endif
!
!-------------------------------------------------------------------
!      Mass equation: 
!-------------------------------------------------------------------
!
           do j=Jstr_nh,Jend_nh
              do i=Istr_nh,Iend_nh
                 k=1
#  ifndef NBQ_IMP
                 div_nbq(i,j,k) =                                       &
                        +  (- on_u(i,j)*pm(i,j)*pn(i,j)                 &
                          - coefb_u(i,j,k) * gdepth_u(i,j,k)            &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +              &
                                    Hzr_half_nbq(i,j,k) ) )  )          &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)       &
                              * mijk2lmom_nh(i,j,k,1)                   &

                         - coefa_u(i,j,k+1) * gdepth_u(i,j,k+1)         &
                          / (0.5*(Hzr_half_nbq(i-1,j,k+1)+              &
                                  Hzr_half_nbq(i,j,k+1)))               &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k+1)     & 
                              * mijk2lmom_nh(i,j,k+1,1)                 &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)                &
                          - coefb_u(i+1,j,k) * gdepth_u(i+1,j,k)        &
                          / (0.5*(Hzr_half_nbq(i,j,k)+                  &
                                  Hzr_half_nbq(i+1,j,k)))   )           &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)     &
                              * mijk2lmom_nh(i+1,j,k,1)                 &

                        - coefa_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)      &
                          / (0.5*(Hzr_half_nbq(i,j,k+1)+                &
                                  Hzr_half_nbq(i+1,j,k+1)))             &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k+1)   &
                              * mijk2lmom_nh(i+1,j,k+1,1)   
# endif            

                 divz_nbq(i,j,k) =                                      &
                        + 1. / Hzw_half_nbq(i,j,k)                      &
                             / Hzr_half_nbq(i,j,k) * qdmw_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,3)  
                                
                 div_nbq(i,j,k) = div_nbq(i,j,k) + divz_nbq(i,j,k)

                 rho_nbq(i,j,k) = rho_nbq(i,j,k)                        &
                      - dtnbq * div_nbq(i,j,k)
                               
                 do k=2,N-1
#  ifndef NBQ_IMP
                    div_nbq(i,j,k) =                                    &
                        +  ( - on_u(i,j)*pm(i,j)*pn(i,j)                &
                          - ( coefb_u(i,j,k) - coefa_u(i,j,k) )         &
                          * gdepth_u(i,j,k)                             &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +              &
                                    Hzr_half_nbq(i,j,k) ) )  )          &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)       &
                              * mijk2lmom_nh(i,j,k,1)                   &

                         - coefa_u(i,j,k+1) * gdepth_u(i,j,k+1)         &
                          / (0.5*(Hzr_half_nbq(i-1,j,k+1)+              &
                                  Hzr_half_nbq(i,j,k+1)))               &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k+1)     & 
                             * mijk2lmom_nh(i,j,k+1,1)                  &

                         + coefb_u(i,j,k-1) * gdepth_u(i,j,k-1)         &
                          / (0.5*(Hzr_half_nbq(i-1,j,k-1)+              &
                                  Hzr_half_nbq(i,j,k-1)))               &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k-1)     & 
                              * mijk2lmom_nh(i,j,k-1,1)                 &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)                &
                          - ( coefb_u(i+1,j,k) - coefa_u(i+1,j,k) )     &
                            * gdepth_u(i+1,j,k)                         &
                          / (0.5*(Hzr_half_nbq(i,j,k)+                  &
                                  Hzr_half_nbq(i+1,j,k)))   )           &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)     &
                              * mijk2lmom_nh(i+1,j,k,1)                 &

                        - coefa_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)      &
                          / (0.5*(Hzr_half_nbq(i,j,k+1)+                &
                                  Hzr_half_nbq(i+1,j,k+1)))             &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k+1)   &
                              * mijk2lmom_nh(i+1,j,k+1,1)               &

                        + coefb_u(i+1,j,k-1) * gdepth_u(i+1,j,k-1)      &
                          / (0.5*(Hzr_half_nbq(i,j,k-1)+                &
                                  Hzr_half_nbq(i+1,j,k-1)))             & 
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k-1)   &
                              * mijk2lmom_nh(i+1,j,k-1,1)   
# endif

                    divz_nbq(i,j,k) =                                   &
                        + 1. / Hzw_half_nbq(i,j,k)                      &
                             / Hzr_half_nbq(i,j,k) * qdmw_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,3)                   &

                        - 1. / Hzw_half_nbq(i,j,k-1)                    &
                             / Hzr_half_nbq(i,j,k) * qdmw_nbq(i,j,k-1)  &
                             * mijk2lmom_nh(i,j,k-1,3)  

                    div_nbq(i,j,k) = div_nbq(i,j,k) + divz_nbq(i,j,k)
                     rho_nbq(i,j,k) = rho_nbq(i,j,k)                     &
                      - dtnbq * div_nbq(i,j,k)

                 enddo

                 k=N

#  ifndef NBQ_IMP
                 div_nbq(i,j,k) = &
                        +  ( - on_u(i,j)*pm(i,j)*pn(i,j)                &
                          + ( coefa_u(i,j,k)   * gdepth_u(i,j,k)        &
                            - coefb_u(i,j,k+1) * gdepth_u(i,j,k+1) )    &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +              &
                                    Hzr_half_nbq(i,j,k) ) )  )          &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)       &
                              * mijk2lmom_nh(i,j,k,1)                   &

                         + coefb_u(i,j,k-1) * gdepth_u(i,j,k-1)         &
                          / (0.5*(Hzr_half_nbq(i-1,j,k-1)+              &
                                  Hzr_half_nbq(i,j,k-1)))               &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k-1)     & 
                              * mijk2lmom_nh(i,j,k-1,1)                 &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)                &
                           + ( coefa_u(i+1,j,k)   * gdepth_u(i+1,j,k)   &
                            - coefb_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)) &
                          / (0.5*(Hzr_half_nbq(i,j,k)+                  &
                                  Hzr_half_nbq(i+1,j,k)))   )           &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)     &
                              * mijk2lmom_nh(i+1,j,k,1)                 &

                        + coefb_u(i+1,j,k-1) * gdepth_u(i+1,j,k-1)      &
                          / (0.5*(Hzr_half_nbq(i,j,k-1)+                &
                                  Hzr_half_nbq(i+1,j,k-1)))             & 
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k-1)   &
                              * mijk2lmom_nh(i+1,j,k-1,1)  
#endif            
                 divz_nbq(i,j,k) =                                      &
                        + 1. / Hzw_half_nbq(i,j,k)                      &
                             / Hzr_half_nbq(i,j,k) * qdmw_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,3)                   &

                        - 1. / Hzw_half_nbq(i,j,k-1)                    &
                              / Hzr_half_nbq(i,j,k) * qdmw_nbq(i,j,k-1) &
                              * mijk2lmom_nh(i,j,k-1,3) 

                 div_nbq(i,j,k) = div_nbq(i,j,k) + divz_nbq(i,j,k)

                 rho_nbq(i,j,k) = rho_nbq(i,j,k)                        &
                      - dtnbq *  div_nbq(i,j,k)


              enddo
           enddo

!
!-------------------------------------------------------------------
!      Density open boundary conditions
!-------------------------------------------------------------------
!
# ifdef OBC_NBQ
!        call rnbq_bc_tile (Istr,Iend,Jstr,Jend, WORK)
# endif
!
!*******************************************************************
!*******************************************************************
      enddo    ! NBQ loop

!*******************************************************************
!*******************************************************************
!
!-------------------------------------------------------------------
!......Set NBQ/EXT coupling terms
!-------------------------------------------------------------------
!
      call ruijk_nbq(2)

#ifdef RVTK_DEBUG
!       call check_tab2d(rubar_nbq(:,:),'rubar_nbq step3d_nbq','u')
!       call check_tab2d(rvbar_nbq(:,:),'rvbar_nbq step3d_nbq','v')
!       call check_tab3d(rw_nbq_ext(:,:,0:N),'rw_nbq_ext step3d_nbq','r')
#endif      
      call densityijk_nbq(20)

   
      end subroutine step3d_fbijk_nbq

#else
      subroutine step3d_fbijk_nbq_empty
      end subroutine step3d_fbijk_nbq_empty
#endif
