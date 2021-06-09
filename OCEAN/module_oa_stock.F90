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
!> @brief Array of the structured type type_wf, which will store for each configuation and variable,
!! the statistical analysis (i.e., convolution product between the analysed variable and a given atome,
!! e.g., a daughter wavelet Psi(s,ta).
!
!> @details Array: 
!! - tpsi_oa stores the type of atome requested for each configuration,
!! - nzw_oa is the total size of the array wf_oa
!! Array wf_oa stores all the analysis requested in the simulation for:
!! - each configuration/variable, 
!! - each convolution window, 
!! - each scale of analysis,
!! - each horizontal grid point, and each vertical grid point.
!!
! REVISION HISTORY:
!
!> @authors
!! - Francis Auclair , Jochem Floor and Ivane Pairaud:
!!  - Initial version
!! - B. Lemieux-Dudon
!!  - Namelists (06/2014), Stand-alone version + optimization (01/2015)
!!  - doxygen comments
!! - Francis Auclair, B. Lemieux-Dudon, C. Nguyen
!!  - Croco-OnlineA module interface, 1st version, Spring 2020
!> @date 2015
!> @todo
!
!------------------------------------------------------------------------------

!************************************************************************
!.....module module_oa_stock
!************************************************************************

      module module_oa_stock


      use module_oa_type


      integer                                                         & 
           nzw_oa                                                        ! taille du vecteur resultat

      integer,dimension(:),allocatable::                              &  ! (nmv_oa)
           tpsi_oa                                                       ! type d atome utilise


      !> To store the online analysis
      type(type_wf),dimension(:),allocatable::wf_oa

      !> Conversion of the specific time of analysis lti_a to the corresponding wf_oa array entry l_a
      integer,dimension(:),allocatable,target::tallocated_oa


      end module module_oa_stock

#else
      module module_oa_stock
      end module
#endif
