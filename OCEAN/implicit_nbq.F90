#include "cppdefs.h"
#if defined NBQ && defined NBQ_IMP
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
      subroutine implicit_nbq ( icall )


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

      qdmimp_nbq  = 0. 
      pdv_nbq     = 0.
      puv_nbq     = 0.
      plv_nbq     = 0.
      pdvsave_nbq = 0.
      puvsave_nbq = 0.
      plvsave_nbq = 0.

      endif ! icall == 0

      if (icall.eq.1) then
!*******************************************************************************
!   Construction du systeme 
!*******************************************************************************
       flag_i = 1

       ptri_nbq = neqmimp_nbq
       call amub2_tri(                                             &
             neqmimp_nbq                                           &
            ,neqcimp_nbq                                           &
            ,flag_i                                                &
            ,mimpv_nbq(1:)                                         &  ! Taille!
            ,mimpj_nbq                                             &
            ,mimpi_nbq                                             &
            ,cimpv_nbq(1:)                                         &
            ,cimpj_nbq                                             &
            ,cimpi_nbq                                             &
                )
     

# ifdef NBQ_IMP_LU
      if (neqmimp_nbq.ne.0) then
      call dgttrf (                            &
                    ptri_nbq                   &
                   ,plv_nbq    (1:ptri_nbq-1)  &
                   ,pdv_nbq    (1:ptri_nbq  )  &
                   ,puv_nbq    (1:ptri_nbq-1)  &
                   ,puv2_nbq   (1:ptri_nbq-2)  &
                   ,ipiv_nbq   (1:ptri_nbq  )  &
                   ,info_i )
      endif
# endif

      endif  ! icall == 1

      if (icall.eq.2) then
!*******************************************************************************
!   RHS 
!*******************************************************************************

#ifdef IMP_SLOW
!.......QDM pour le RHS: terme CONTxQDM
!                        qdm_z must 0 here
        qdmimp_nbq(1:neqmom_nh(1)+neqmom_nh(2))=qdm_nbq_a(1:neqmom_nh(1)+neqmom_nh(2),vnnew_nbq) 

!.......Produit MatxVect: CONTxQDM
        call amux                (                                   &
               neqcont_nh                                            &
              ,qdmimp_nbq(1)                                         &
              ,rhsimp2_nbq(1)                                        &
              ,contv_nh(1)                                           &
              ,contj_nh(1)                                           &
              ,conti_nh(1)       )

        rhsimp2_nbq(1:neqcimp_nbq)=  csvisc2_nbq                    &   
                             * rhp_nbq_a(1:neqcimp_nbq,rnrhs_nbq)   &
                             - dtnbq*rhsimp2_nbq(1:neqcimp_nbq)  
#else

   qdmimp_nbq(1:neqmom_nh(1)+neqmom_nh(2))=qdm_nbq_a(1:neqmom_nh(1)+neqmom_nh(2),vnnew_nbq) 
   do l_nbq = 1, neqcont_nh 
     t_i = 0.0D+00
     do l1_nbq = conti_nh(l_nbq), conti_nh(l_nbq+1)-1
       t_i = t_i + contv_nh(l1_nbq) * qdmimp_nbq(contj_nh(l1_nbq))
     end do
     rhsimp2_nbq(l_nbq) = csvisc2_nbq * rhp_nbq_a(l_nbq,rnrhs_nbq)   &
                         -dtnbq*t_i
   end do
!        rhsimp2_nbq(1:neqcimp_nbq)=  csvisc2_nbq                    &   
!                             * rhp_nbq_a(1:neqcimp_nbq,rnrhs_nbq)   &
!                            -dtnbq*rhsimp2_nbq(1:neqcimp_nbq)  
#endif

#ifdef IMP_SLOW
!.......RHS = MOMxRHS(cont) 
        call amux(                                                  &
               neqmimp_nbq                                          &
              ,rhsimp2_nbq(1:neqcimp_nbq)                           & 
              ,rhsimp_nbq(1)                                        &
              ,mimpv_nbq(1)                                         &
              ,mimpj_nbq(1)                                         &
              ,mimpi_nbq(1)       )  

!.....Sauvegarde du second membre &
!     Ajouts des deux termes de l'equation QDM  du RHS: MOM * RHS(cont) + RHS(mom)
!     Equation de continuite:

       rhsimp_nbq(1:neqmimp_nbq) = rhsimp_nbq(1:neqmimp_nbq)             &
       + dtnbq * dqdmdt_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0))  &
       + qdm_nbq_a(neqmom_nh(1)+neqmom_nh(2)+1:neqmom_nh(0),vnrhs_nbq)   
#else
   do l_nbq = 1, neqmimp_nbq 
     t_i = 0.0D+00
     do l1_nbq = mimpi_nbq(l_nbq), mimpi_nbq(l_nbq+1)-1
       t_i = t_i + mimpv_nbq(l1_nbq) * rhsimp2_nbq(mimpj_nbq(l1_nbq))
     end do
     rhsimp_nbq(l_nbq) = t_i &
       + dtnbq * dqdmdt_nbq_a(neqmom_nh(1)+neqmom_nh(2)+l_nbq)  &
       + qdm_nbq_a(neqmom_nh(1)+neqmom_nh(2)+l_nbq,vnrhs_nbq)   
   end do
#endif
     
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
       pdvsave_nbq(1:ptri_nbq) = pdv_nbq(1:ptri_nbq)
       plvsave_nbq(1:ptri_nbq) = plv_nbq(1:ptri_nbq)
       puvsave_nbq(1:ptri_nbq) = puv_nbq(1:ptri_nbq)
      else
!......Recopie de la matrice tri-diagonale...
       pdv_nbq(1:ptri_nbq) = pdvsave_nbq(1:ptri_nbq)
       plv_nbq(1:ptri_nbq) = plvsave_nbq(1:ptri_nbq)
       puv_nbq(1:ptri_nbq) = puvsave_nbq(1:ptri_nbq)
      endif

!..... Resolution du systeme tri-diagonal:
      if (neqmimp_nbq.ne.0) then
       call dgtsv  ( ptri_nbq                  &
                   ,nd_i                       &
                   ,plv_nbq    (1:ptri_nbq-1)  &
                   ,pdv_nbq    (1:ptri_nbq  )  &
                   ,puv_nbq    (1:ptri_nbq-1)  &
                   ,rhsimp_nbq (1:ptri_nbq)    &
                   ,ptri_nbq                   &
                   ,ierr_i )  
      endif
# endif

# ifdef NBQ_IMP_LU
      trans_i = 'N'
      if (neqmimp_nbq.ne.0) then
      call dgttrs (                            &
                    trans_i                    &
                   ,ptri_nbq                   &
                   ,nd_i                       &
                   ,plv_nbq    (1:ptri_nbq-1)  &
                   ,pdv_nbq    (1:ptri_nbq  )  &
                   ,puv_nbq    (1:ptri_nbq-1)  &
                   ,puv2_nbq   (1:ptri_nbq-2)  &
                   ,ipiv_nbq   (1:ptri_nbq  )  &
                   ,rhsimp_nbq (1:ptri_nbq  )  &
                   ,ptri_nbq                   &
                   ,info_i )
      endif
# endif
                   
!.....Recopie de la solution implicite pour w (update) et calcul rhs (sum)
      
!     do l_nbq = neqw_nh(1)+1,neqw_nh(6) 
      do l_nbq = neqv_nh(7)+1,neqw_nh(7) 
         qdm_nbq_a(l_nbq,vnnew_nbq) = rhsimp_nbq(l_nbq-(neqmom_nh(1)+neqmom_nh(2)))  
         rhssum_nbq_a(l_nbq) = rhssum_nbq_a(l_nbq)                                    &
              + (qdm_nbq_a(l_nbq,vnnew_nbq)-qdm_nbq_a(l_nbq,vnrhs_nbq)) / dtnbq   &
              - dqdmdt_nbq_a(l_nbq)   
      enddo

      endif
      

       end subroutine implicit_nbq 
#else
      subroutine implicit_nbq_empty
      end subroutine implicit_nbq_empty
#endif
