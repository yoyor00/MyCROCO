! $Id: output_oa.F
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
!======================================================================
!
!> @brief Croco interface for Online Analysis (OA=OnlineAnalysis)
!! - Outputs the online spectral and wavelet analyses with XIOS
!!
!! @details written to support tile-threads and 
!!  possible dual OpenMP-MPI parallelization. 
!!  See Online_spectral_diags called by main.F and step.F
!!  More history/info in source module_interface_oa.F90
!!
!> @authors  
!! - B. Lemieux-Dudon
!!  - Croco Tile-thread compliant version (2021) 
!!    supporting tile-threads and possible dual OpenMP-MPI 
!!    parallelization of the horizontal domain
!!     => output_oa : OA updates global array looping on tiles 
!!        and finally send them to XIOS (out of the loop)
!!     => var2d_oa_out, var3d_oa_out : 
!!        2D using Croco working arrays work2d, work2d2
!!     => module_parameter_oa : allocating var3d_oa, var2d_oa 
!!        to the croco GLOBAL_3D_ARRAY, 2D_ARRAY 
!!  - based on a preliminary Croco-OA interface version : spring 2020.
!! - More history (authors, comments) in source module_interface_oa.F90
!> @todo BLXD TODO test var3d_oa_out with work.h as applied in
!!  var2d_oa_out
!
!  REFERENCE:
!  
!======================================================================
!
#include "cppdefs.h"
! BLXD includes :
! - cppdefs_dev.h
! - set_global_definitions.h
!   with GLOBAL_2D_ARRAY
!        PRIVATE_SCRATCH_2D_ARRAY... 
#ifdef ONLINE_ANALYSIS

      subroutine output_oa
!
!=======================================================================
!                                                                      !
! Output for OA Online analysis block
!                                                                      !
!=======================================================================
!
      use xios_module
      use module_interface_oa, only : tvar_oa
     & ,tgv3d_oa, tgv_oa, tgvnam_oa, nzvc_oa
     & ,verbose_oa, io_hist, directory_out_oa
      use module_oa_variables,  only : tv_oa     
     & ,tupd_oa  ! BLXD added for test
      use module_oa_periode,    only : begc_oa, tvc_oa, tc_oa, nzc_oa
     & ,swt_wfpf_oa

      implicit none
!
#include "param.h"
#include "scalars.h"
#include "work.h"
!
                integer :: ntrds,trd,range,
     &        my_first,my_last, tile, my_iif

      integer :: omp_get_num_threads, omp_get_thread_num

      character(len=22) :: vnamr, vnami 

      real :: workim(GLOBAL_2D_ARRAY,N)
      real :: workre(GLOBAL_2D_ARRAY,N)

      integer :: ic, iv, ivc

!      ilc=1+iic-ntstart   ! number of time step since restart

      oa_config_loop : do ic = 1, nzc_oa
      
        ivc = 0
        do iv = begc_oa(ic), begc_oa(ic+1) -1
            ivc = ivc + 1
        end do
        oa_stop1 : if ( ivc /= nzvc_oa( tc_oa(ic) ) ) then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR-iv-ic'
           stop
           !MPI_master_only write(stdout,*) 'Expected total # of var nzvc_oa( tc_oa(ic) ) ',nzvc_oa( tc_oa(ic) )
        end if oa_stop1
        !could enable to retrieve a name, requires to use tupd_oa
        ! # of variable per configuration ic
        ! do iv = 1, nzvc_oa( tc_oa(ic) )
        ivc = 0
        the_oa_var_loop : do iv = begc_oa(ic), begc_oa(ic+1) -1
           ivc = ivc + 1
        oa_stop2 : if ( 
     &     tvar_oa( tc_oa(ic), tvc_oa(ic), ivc ) /= tupd_oa(iv) ) 
     &     then
           MPI_master_only write(stdout,fmt='(a)') '...OA-ERR tupd_oa'
           stop
        end if oa_stop2
      
#ifdef OA_TRACES && ( MILES || IGW)
        call traces_oa1(verbose_oa)
#endif

! 3Dvar
        if_is_3dvar : if ( tgv3d_oa( tv_oa(iv) ) == 3 ) then

!
! XIOS Case only
!
#ifdef XIOS
          ntrds=omp_get_num_threads()
          trd=omp_get_thread_num()
          range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
          my_first=trd*range
          my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
          
          do tile=my_first,my_last
              ! BLXD TODO test with common space work3d_r and workr
              !call var3d_oa_out(tile,workr,workim,vnamr,vnami)
              call var3d_oa_out(tile,ic,iv,workre,workim,vnamr,vnami)
          end do
C$OMP BARRIER
          call xios_send_field( trim(vnamr), workre  )
          if ( swt_wfpf_oa(iv)==4 )
     &      call xios_send_field( trim(vnami), workim )
#else
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...Non XIOS outputs doesn''t work with OA'
#endif

! 2Dvar
        else if ( tgv3d_oa( tv_oa(iv) ) == 2 ) then if_is_3dvar

#ifdef XIOS
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...OA var2d outputs using croco work2d/work2d2'

          ntrds=omp_get_num_threads()
          trd=omp_get_thread_num()
          range=(NSUB_X*NSUB_E+ntrds-1)/ntrds
          my_first=trd*range
          my_last=min(my_first + range-1, NSUB_X*NSUB_E-1)
          
          do tile=my_first,my_last
              call var2d_oa_out(tile,ic,iv,work2d,work2d2,vnamr,vnami)
          end do
C$OMP BARRIER
          call xios_send_field( trim(vnamr), work2d  )
          if ( swt_wfpf_oa(iv)==4 )
     &      call xios_send_field( trim(vnami), work2d2 )
#else
           MPI_master_only write(stdout,fmt='(a)') 
     &      '...Non XIOS outputs doesn''t work with OA'
#endif


        else if_is_3dvar
           MPI_master_only write(stdout,fmt='(a)') '...ERR output_oa'
           stop
        endif if_is_3dvar

!
        end do the_oa_var_loop
      end do oa_config_loop
!   !
      contains

!      subroutine traces_oa1(verbose,iic)        
      subroutine traces_oa1(verbose)        

        implicit none
        integer, intent(in)    :: verbose
        if ( verbose >=6 ) then
            MPI_master_only write(stdout,fmt='(a13,(3(i3.3,1x)))')     
     &     'Cfg-#var-Var ',ic,ivc,iv
        end if
        return
      end subroutine traces_oa1



      end subroutine output_oa

#else /* ONLINE_ANALYSIS */

      subroutine output_oa_empty
      return
      end subroutine output_oa_empty

#endif /* ONLINE_ANALYSIS */

