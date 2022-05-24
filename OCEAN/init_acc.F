        subroutine initialisation_openacc
#ifdef JEANZAY
        USE openacc
        implicit none
#include "cppdefs.h"
#include "param.h"
 
        character(len=6) :: local_rank_env
        integer          :: local_rank_env_status
        integer          :: ngpus, local_rank, devicenum
 
        ! Initialisation d'OpenACC
        !$acc init
 
         ! Récupération du rang local du processus via la variable d'environnement
        ! positionnée par Slurm, l'utilisation de MPI_Comm_rank n'étant pas encore
        ! possible puisque cette routine est utilisée AVANT l'initialisation de MPI
!        call get_environment_variable(name="SLURM_LOCALID", &
        call get_environment_variable(name="SLURM_PROCID", 
     &      value=local_rank_env, status=local_rank_env_status)
 
        if (local_rank_env_status == 0) then
            read(local_rank_env, *) local_rank
                ngpus=acc_get_num_devices( acc_device_nvidia )
                !choose device to use by this rank
                devicenum = MOD( local_rank, ngpus )
            ! Définition du GPU à utiliser via OpenACC
            call acc_set_device_num( devicenum, acc_get_device_type())
            my_acc_device=devicenum
            compute_on_device=.TRUE.
!            write(100+local_rank,*) devicenum,ngpus
        else
            print *,    
     &   "Erreur : impossible de déterminer le rang local du processus"
            stop 1
        end if
#endif
        end subroutine initialisation_openacc
