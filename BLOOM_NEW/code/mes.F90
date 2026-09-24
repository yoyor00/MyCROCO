Module MES
#include "cppdefs.h"

    USE module_BIOLink
    USE comBIOLink
    USE COMBLOOM
    USE substance, only : name_var_r
    

    IMPLICIT NONE
    PRIVATE

    public :: get_total_MES

    contains
    
    function get_total_MES(i,j,k, c) result(TotMes)

        INTEGER :: i,j,k, i_sub, iv
        REAL(KIND=rsh), DIMENSION(nv_state):: c
        REAL(KIND=rsh):: TotMes
     
        INTEGER, PARAMETER   :: nb_part_subs = 3
        CHARACTER(len=100), parameter :: part_subs(nb_part_subs) = [CHARACTER(len=100):: "MUD", "GRAVEL","SAND"]
        INTEGER, DIMENSION(nb_part_subs) :: idx_sub
         

        DO i_sub=1,nb_part_subs
                idx_sub(i_sub) = findloc(name_var_r, TRIM(ADJUSTL(ADJUSTR(part_subs(i_sub)))), dim=1)
                !MPI_master_only WRITE(*,*) "l'index correcspondant à la variable ", TRIM(ADJUSTL(ADJUSTR(part_subs(i_sub)))), &
                !        " est " , idx_sub(i_sub)
        END DO

       
        TotMes = 0.0_rsh
        DO iv=1,nb_part_subs
            TotMes = TotMes + c(idx_sub(iv))*1000.0_rsh  !cmes_3dmgl passe en mg/l
        END DO
        !diag_3d_wat(irk_diag(id_spm_total),k,i,j)=cmes_3dmgl(k,i,j)


    end function get_total_MES

END MODULE MES




