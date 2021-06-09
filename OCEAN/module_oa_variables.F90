#include "cppdefs.h"
#ifdef ONLINE_ANALYSIS
!------------------------------------------------------------------------------
!                               NHOMS
!                Non Hydrostatic Ocean Modeling System      
!------------------------------------------------------------------------------
!
!> @note <a href="http://poc.obs-mip.fr/auclair/WOcean.fr/SNH/index_snh_home.htm"> Main web documentation </a>
!
! DESCRIPTION: 
!
!> @brief
!
!> @details 
!
!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Initial version
!! - B. Lemieux-Dudon
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!!  - BLXD nzvc_oa variables moved to the module_oa_interface with
!!    parametrized dimension maxtyp_oa.
!!  - Max # of convolution windows simultaneously openned pre-calculated 
!!    nmsimult_oa
!> @date 2015
!> @todo
!
!------------------------------------------------------------------------------

!************************************************************************
!.....module module_oa_variables
!************************************************************************

      module module_oa_variables

!      use module_interface_oa, only : maxtyp_oa                         ! Max. number of predef. configuration types
      
!.....parameters:

      ! Moved to module_oa_time 
      ! integer,parameter:: nmsimult_oa = 900000                          ! nombre maximum de fenetres ouvertes simultanement

!.....initialisation:
!STDALONE      integer :: ifl_init_oa = 0                               ! flag pour l'initialisation.

!.....variables:

      integer::                                               &
           nzv_oa                                                ! nombre de variables (assoc. aux configurations)
           !,nzvs_oa                                           & ! nombre total de points (taille de la structure spatiale 2d du vecteur d etat)
           !,nzvs3d_oa                                           ! nombre total de points (taille de la structure spatiale 3d du vecteur d etat)

! BLXD ltrec_oa and if_first_rec_oa are debugging variables
      integer,dimension(:),allocatable::                      & ! (nmv_oa)
           swt_d_oa                                           & ! caracteristiques spatiales de la variable (voir notebook_wavelet)
           ,swt_t_oa                                          & ! caracteristiques temporelles de la variable (voir notebook_wavelet)          
           ,tv_oa                                             & ! variable symphonie associee         
           ,cnb_oa                                            & ! call number: position of the call in the baroclinic / barotropic time step.
!           ,tvc_oa                                            & ! configuration de la variable associee
           ,updv_oa                                           & ! records 1st temporal OA in simulation 
           ,ltrec_oa                                          & ! ...
           ,if_first_rec_oa

      ! BLXD added in module_oa_variables croco2021 May 06
      !      removed from module variables declared in module_interface_oa
      !      with wrong size, ie
      !      integer, dimension(1:6000) :: des_oa  
      !
      ! Variable useful for NHOMS energy analysis if you wish to call OA at specific point in the code
      ! for specific variable that are updated at some irregular code place or time step (e.g., time splitting)
      ! depnds on iv_m so it should have nzv_oa size has cnb_oa 
      ! integer, dimension(1:6000) :: des_oa  

      integer,dimension(:),allocatable::                      &
          des_oa                                                ! (nmv_oa)

      integer,dimension(:),allocatable::                      & ! (nmv_oa+1)
          begvt_oa                                              ! structure temporelle du vecteur d etat

      !integer,dimension(:),allocatable::                      & ! (nmv_oa+1)
      !     begvs_oa                                           & ! structure 2d du vecteur d etat
      !     ,begvt_oa                                            ! structure temporelle du vecteur d etat
      !
      !integer,dimension(:),allocatable::                      & !(nmvs_oa)
      !     l2i_oa                                             & ! transformation l --> i
      !     ,l2j_oa                                              ! transformation l --> j
      !
      !integer,dimension(:),allocatable::                      & ! (nmvs_oa+1)
      !     begvs3d_oa                                         & ! structure 3d du vecteur d etat: debut de la colonne pour un point donne
      !     ,kmin3d_oa                                           ! structure 3d du vecteur d etat: premier niveau a considerer     
      !
      !integer,dimension(:,:,:),allocatable::                  & ! ( 0:meco+1,0:neco+1,nmv_oa)
      !     ij2l_oa                                              ! transformation (i,j) --> l (l=  structure 2d du vecteur d etat)

      !! #BLXD variables moved to the module_oa_interface with parametrized dimension, maxtyp_oa
      !integer:: nzvc_oa(200)
      !integer:: nzvc_oa(maxtyp_oa)  ! #BLXD

      integer                                                         &
           flag_nrj_oa                                                  ! flag pour le calcul de l'energie
 
      integer,dimension(:),allocatable::                              & !(nmv_oa)
           save_oa                                                    & ! sauvegarde finale de la variable dans un fichier
           ,tupd_oa                                                     ! pour variables communes

      ! BLXD parameter moved to the module_oa_interface
      !real :: pi_oa 

!     pour tester les sorties "online":

!      double precision ::                                              &
!       period_test_oa(10)                                             
      double precision,dimension(:,:,:),allocatable,target ::                                              &
       vardp_test_oa
!     ,vardp_test_oa(0:imax+1,0:jmax+1,0:kmax+1)

      integer :: ifl_test_oa
      logical :: isopycne_analysis

      !> Variables required to read OA namelists
      character(len=50), dimension(:), allocatable :: nzc_oa_names
      integer, dimension(:), allocatable           :: nzc_oa_vartyp

      !logical              :: isopycne_analysis
      !character(len=200)   :: directory_in, directory_out
      !character(len=1), parameter :: txtslash='/'
      !integer :: io_unit


      end module module_oa_variables

#else
      module module_oa_variables_empty
      end module
#endif
