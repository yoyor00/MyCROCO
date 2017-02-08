#include "cppdefs.h"
#if defined NBQ && defined NBQ_IMP && defined NBQ_IJK
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
!> @details Loops on the NBQ time step. See SNBQ web pages :
! REVISION HISTORY:
!
!> @authors
!> @date 2015 October
!> @todo
!
!======================================================================
!
      subroutine implicitijk_nbq ( icall )


!      use module_principal  
!      use module_parallele
       use module_nh 
       use module_nbq
!      use module_exp 
# ifdef TRACETXT
      use module_tracetxt_out
# endif
      implicit none

# include "param_F90.h"
# include "scalars_F90.h"
# include "ocean3d.h"
# include "grid.h"
# include "nbq.h"

      integer :: i,j,k
      integer :: icall,flag_i,ierr_i,ifl_ascii,nd_i
      integer :: ioff

      double precision :: t_i

      character *1 :: fact_i,trans_i
      integer :: info_i
      real :: rcond_i,ferr_i,berr_i


!.....Flag pour les sorties ASCII liees aux calculs implicites:
      ifl_ascii = 0


      if (icall.eq.0) then
!*******************************************************************************
! Initialisations
!*******************************************************************************

      pdv_nbq     = 1.
      puv_nbq     = 0.
      plv_nbq     = 0.
      
# ifdef NBQ_IMP_TRIDIAG
      pdvsave_nbq = 0.
      puvsave_nbq = 0.
      plvsave_nbq = 0.
#endif

      endif ! icall == 0

      if (icall.eq.1) then
!*******************************************************************************
!   Construction du systeme 
!*******************************************************************************
        ioff = neqmom_nh(1)+neqmom_nh(2)

        do j=Jstr_nh,Jend_nh
           do i=Istr_nh,Iend_nh 
              k = 1
              l_nbq  = ijk2lmom_nh(i,j,k  ,3) - ioff
              l2_nbq = ijk2lmom_nh(i,j,k+1,3) - ioff
              pdv_nbq(l_nbq ) = 1. + dtnbq * csvisc1_nbq           &
                                * (  1. / Hzr_half_nbq(i,j,k+1)    &
                                   + 1. / Hzr_half_nbq(i,j,k)   )  &     
                                        / Hzw_half_nbq(i,j,k) 
              puv_nbq(l_nbq) =    - dtnbq * csvisc1_nbq            &
                                        / Hzr_half_nbq(i,j,k+1)    &   
                                        / Hzw_half_nbq(i,j,k+1) 
              do k = 2, N-1
                l1_nbq = ijk2lmom_nh(i,j,k-1,3) - ioff - 1
                l_nbq  = ijk2lmom_nh(i,j,k  ,3) - ioff
                l2_nbq = ijk2lmom_nh(i,j,k+1,3) - ioff
                plv_nbq(l_nbq-1) =  -  dtnbq * csvisc1_nbq         &
                                        / Hzr_half_nbq(i,j,k)      &   
                                        / Hzw_half_nbq(i,j,k-1) 
                pdv_nbq(l_nbq ) =  1. + dtnbq * csvisc1_nbq        &
                               *  (  1. / Hzr_half_nbq(i,j,k+1)    &
                                   + 1. / Hzr_half_nbq(i,j,k)   )  &     
                                        / Hzw_half_nbq(i,j,k) 
                puv_nbq(l_nbq) =   -  dtnbq * csvisc1_nbq          &
                                        / Hzr_half_nbq(i,j,k+1)    &   
                                        / Hzw_half_nbq(i,j,k+1) 
              enddo

              k = N           
              l1_nbq = ijk2lmom_nh(i,j,k-1,3) - ioff - 1
              l_nbq  = ijk2lmom_nh(i,j,k  ,3) - ioff
              plv_nbq(l_nbq-1) = - dtnbq * csvisc1_nbq             &
                                        / Hzr_half_nbq(i,j,k)      &   
                                        / Hzw_half_nbq(i,j,k-1) 
              pdv_nbq(l_nbq ) = 1. +  dtnbq * csvisc1_nbq          &
                                        / Hzr_half_nbq(i,j,k)      &     
                                        / Hzw_half_nbq(i,j,k)        
           enddo
        enddo

        neqmimp_nbq = neqmom_nh(0) - ioff
       
# ifdef NBQ_IMP_LU
        if (neqmimp_nbq.ne.0) then
        call dgttrf (                          &
                    neqmimp_nbq                   &
                   ,plv_nbq    (1:neqmimp_nbq-1)  &
                   ,pdv_nbq    (1:neqmimp_nbq  )  &
                   ,puv_nbq    (1:neqmimp_nbq-1)  &
                   ,puv2_nbq   (1:neqmimp_nbq-2)  &
                   ,ipiv_nbq   (1:neqmimp_nbq  )  &
                   ,info_i )
        endif
# endif

      endif  ! icall == 1

      if (icall.eq.2) then
!*******************************************************************************
!   RHS 
!*******************************************************************************

!.......Produit MatxVect: CONTxQDM

          ioff = neqmom_nh(1)+neqmom_nh(2)
        
           do j=Jstr_nh,Jend_nh
              do i=Istr_nh,Iend_nh
                 k=1
                 div_nbq(i,j,k) =                                    &
                        +  (- on_u(i,j)*pm(i,j)*pn(i,j)              &
                          - coefb_u(i,j,k) * gdepth_u(i,j,k)         &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +           &
                                    Hzr_half_nbq(i,j,k) ) )  )       &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,1)                &

                         - coefa_u(i,j,k+1) * gdepth_u(i,j,k+1)      &
                          / (0.5*(Hzr_half_nbq(i-1,j,k+1)+           &
                                  Hzr_half_nbq(i,j,k+1)))            &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k+1)  & 
                              * mijk2lmom_nh(i,j,k+1,1)              &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)             &
                          - coefb_u(i+1,j,k) * gdepth_u(i+1,j,k)     &
                          / (0.5*(Hzr_half_nbq(i,j,k)+               &
                                  Hzr_half_nbq(i+1,j,k)))   )        &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)  &
                              * mijk2lmom_nh(i+1,j,k,1)              &

                        - coefa_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)   &
                          / (0.5*(Hzr_half_nbq(i,j,k+1)+             &
                                  Hzr_half_nbq(i+1,j,k+1)))          &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k+1)&
                              * mijk2lmom_nh(i+1,j,k+1,1)           

                 rhsijk_nbq(i,j,k) =  csvisc2_nbq * rho_nbq(i,j,k)   &
                                    - dtnbq       * div_nbq(i,j,k)  
                               
                 do k=2,N-1
                    div_nbq(i,j,k) =                                 &
                        +  ( - on_u(i,j)*pm(i,j)*pn(i,j)             &
                          - ( coefb_u(i,j,k) - coefa_u(i,j,k) )      &
                          * gdepth_u(i,j,k)                          &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +           &
                                    Hzr_half_nbq(i,j,k) ) )  )       &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,1)                &

                         - coefa_u(i,j,k+1) * gdepth_u(i,j,k+1)      &
                          / (0.5*(Hzr_half_nbq(i-1,j,k+1)+           &
                                  Hzr_half_nbq(i,j,k+1)))            &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k+1)  & 
                             * mijk2lmom_nh(i,j,k+1,1)               &

                         + coefb_u(i,j,k-1) * gdepth_u(i,j,k-1)      &
                          / (0.5*(Hzr_half_nbq(i-1,j,k-1)+           &
                                  Hzr_half_nbq(i,j,k-1)))            &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k-1)  & 
                              * mijk2lmom_nh(i,j,k-1,1)              &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)             &
                          - ( coefb_u(i+1,j,k) - coefa_u(i+1,j,k) )  &
                            * gdepth_u(i+1,j,k)                      &
                          / (0.5*(Hzr_half_nbq(i,j,k)+               &
                                  Hzr_half_nbq(i+1,j,k)))   )        &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)  &
                              * mijk2lmom_nh(i+1,j,k,1)              &

                        - coefa_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)   &
                          / (0.5*(Hzr_half_nbq(i,j,k+1)+             &
                                  Hzr_half_nbq(i+1,j,k+1)))          &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k+1)&
                              * mijk2lmom_nh(i+1,j,k+1,1)            &

                        + coefb_u(i+1,j,k-1) * gdepth_u(i+1,j,k-1)   &
                          / (0.5*(Hzr_half_nbq(i,j,k-1)+             &
                                  Hzr_half_nbq(i+1,j,k-1)))          & 
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k-1) &
                              * mijk2lmom_nh(i+1,j,k-1,1)   

                    rhsijk_nbq(i,j,k) =  csvisc2_nbq * rho_nbq(i,j,k) &
                                       - dtnbq       * div_nbq(i,j,k)  
                 enddo

                 k=N
                 div_nbq(i,j,k) = &
                        +  ( - on_u(i,j)*pm(i,j)*pn(i,j)             &
                          + ( coefa_u(i,j,k)   * gdepth_u(i,j,k)     &
                            - coefb_u(i,j,k+1) * gdepth_u(i,j,k+1) ) &
                          / (0.5*( Hzr_half_nbq(i-1,j,k) +           &
                                    Hzr_half_nbq(i,j,k) ) )  )       &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k)    &
                              * mijk2lmom_nh(i,j,k,1)                &

                         + coefb_u(i,j,k-1) * gdepth_u(i,j,k-1)      &
                          / (0.5*(Hzr_half_nbq(i-1,j,k-1)+           &
                                  Hzr_half_nbq(i,j,k-1)))            &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i,j,k-1)  & 
                              * mijk2lmom_nh(i,j,k-1,1)              &

                         + ( on_u(i+1,j)*pm(i,j)*pn(i,j)              &
                           + ( coefa_u(i+1,j,k)   * gdepth_u(i+1,j,k)   &
                            - coefb_u(i+1,j,k+1) * gdepth_u(i+1,j,k+1)) &
                          / (0.5*(Hzr_half_nbq(i,j,k)+               &
                                  Hzr_half_nbq(i+1,j,k)))   )        &
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k)  &
                              * mijk2lmom_nh(i+1,j,k,1)              &

                        + coefb_u(i+1,j,k-1) * gdepth_u(i+1,j,k-1)   &
                          / (0.5*(Hzr_half_nbq(i,j,k-1)+             &
                                  Hzr_half_nbq(i+1,j,k-1)))          & 
                          / Hzr_half_nbq(i,j,k) * qdmu_nbq(i+1,j,k-1)&
                              * mijk2lmom_nh(i+1,j,k-1,1)         
 
                 rhsijk_nbq(i,j,k) =  csvisc2_nbq * rho_nbq(i,j,k)   &
                                    - dtnbq       * div_nbq(i,j,k)     
              enddo
          enddo 


!.......RHS = MOMxRHS(cont) 
!.....Sauvegarde du second membre &
!     Ajouts des deux termes de l'equation QDM  du RHS: MOM * RHS(cont) + RHS(mom)
!     Equation de continuite:

        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh   
                l_nbq  = ijk2lmom_nh(i,j,k  ,3) - ioff                                          
                rhsimp_nbq (l_nbq) =                                 &
                        csvisc1_nbq * float(mijk2lq_nh(i,j,k))       &
                      * (rhsijk_nbq(i,j,k) - rhsijk_nbq(i,j,k+1)     &
                      * float(mijk2lq_nh(i,j,k+1)))                  &
                      + dtnbq * rho0 * rwint_nbq(i,j,k)              &
                      + qdmw_nbq(i,j,k)
             enddo
           enddo
        enddo

!*******************************************************************************
! Inversion du systeme implicite
!
!      matrice : tri-diagonale
!
!      rhs     : rhsimp_nbq
!
!      # de lignes: neqmimp_nbq
!
!
!
!*******************************************************************************

      nd_i = 1  ! Une seule inversion

# ifdef NBQ_IMP_TRIDIAG
      if (iif.eq.1.and.iteration_nbq.eq.1) then
!......Sauvegarde de la matrice tri-diagonale
!        necessaire car dgtsv modifie la matrice tri-diag...
       pdvsave_nbq(1:neqmimp_nbq) = pdv_nbq(1:neqmimp_nbq)
       plvsave_nbq(1:neqmimp_nbq) = plv_nbq(1:neqmimp_nbq)
       puvsave_nbq(1:neqmimp_nbq) = puv_nbq(1:neqmimp_nbq)
      else
!......Recopie de la matrice tri-diagonale...
       pdv_nbq(1:neqmimp_nbq) = pdvsave_nbq(1:neqmimp_nbq)
       plv_nbq(1:neqmimp_nbq) = plvsave_nbq(1:neqmimp_nbq)
       puv_nbq(1:neqmimp_nbq) = puvsave_nbq(1:neqmimp_nbq)
      endif

!..... Resolution du systeme tri-diagonal:
      if (neqmimp_nbq.ne.0) then
       call dgtsv  ( neqmimp_nbq                  &
                   ,nd_i                       &
                   ,plv_nbq    (1:neqmimp_nbq-1)  &
                   ,pdv_nbq    (1:neqmimp_nbq  )  &
                   ,puv_nbq    (1:neqmimp_nbq-1)  &
                   ,rhsimp_nbq (1:neqmimp_nbq)    &
                   ,neqmimp_nbq                   &
                   ,ierr_i )  
      endif
# endif

# ifdef NBQ_IMP_LU
      trans_i = 'N'
      if (neqmimp_nbq.ne.0) then
      call dgttrs (                            &
                    trans_i                    &
                   ,neqmimp_nbq                   &
                   ,nd_i                       &
                   ,plv_nbq    (1:neqmimp_nbq-1)  &
                   ,pdv_nbq    (1:neqmimp_nbq  )  &
                   ,puv_nbq    (1:neqmimp_nbq-1)  &
                   ,puv2_nbq   (1:neqmimp_nbq-2)  &
                   ,ipiv_nbq   (1:neqmimp_nbq  )  &
                   ,rhsimp_nbq (1:neqmimp_nbq  )  &
                   ,neqmimp_nbq                   &
                   ,info_i )
      endif
# endif
                   
!.....Recopie de la solution implicite pour w (update) et calcul rhs (sum)
      

        do j=Jstr_nh,Jend_nh
          do k=1,N
            do i=Istr_nh,Iend_nh      
               l_nbq  = ijk2lmom_nh(i,j,k  ,3) - ioff                                                                  
               rhssumw_nbq(i,j,k) = rhssumw_nbq(i,j,k)            &
                 + (rhsimp_nbq(l_nbq)-qdmw_nbq(i,j,k)) / dtnbq   &
                 - rho0 * rwint_nbq(i,j,k)    
               qdmw_nbq(i,j,k) = rhsimp_nbq(l_nbq)  
             enddo
           enddo
        enddo

      endif
      

       end subroutine implicitijk_nbq 
#else
      subroutine implicitijk_nbq_empty
      end subroutine implicitijk_nbq_empty
#endif
