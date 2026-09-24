
MODULE shared

#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"

    USE module_BIOLink
    USE comBIOLink
    USE comBIOLink_helping
    USE comBIOLink_physics
    USE COMBLOOM
    
    IMPLICIT NONE

    PUBLIC:: shared_parameters
    REAL(KIND=rsh) :: dtbiojour, temper,tempabs,sali,epn
    REAL(KIND=rsh) :: effetchaleur,effetturbidite
    REAL(KIND=rsh), PARAMETER :: ratio_mgO_to_mumolN_catabol=0.212, &
                                  ratio_mgO_to_mumolN_anabol=0.260
         
#if defined key_BLOOM_opt2
    LOGICAL, PARAMETER :: l_BLOOM_opt2 = .TRUE.
#else
    LOGICAL, PARAMETER :: l_BLOOM_opt2 = .FALSE.
#endif


#if defined key_BLOOM_insed
    LOGICAL :: l_BLOOM_insed = .TRUE.
#else
    LOGICAL :: l_BLOOM_insed = .FALSE.
#endif

    contains
        SUBROUTINE shared_parameters(i,j,k)

            INTEGER, intent(in)                 :: i,j,k

            dtbiojour=REAL(dtbio/86400.0_rlg,rsh)

            epn=thicklayerW_C(k,i,j)
            tempabs=TEMP_BIOLink(k,i,j)+273.15_rsh
            temper=max(0.0_rsh,TEMP_BIOLink(k,i,j))
            IF (temper > 30.0_rsh) THEN
                 write(*,*) '!!!!!tempwat>30 !!!!!!! (forcee a 30degres)',ijour_BIOLINK,imois_BIOLINK,iheure_BIOLINK,k,i,j,temper
                 temper=30.0_rsh
                 tempabs=temper+273.15_rsh
             ENDIF

             effetchaleur=exp(p_T_effect*temper)

             effetturbidite=1.0_rsh
             IF(l_ChlNratio_var) effetturbidite=EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct/ &
                 (sqrt(1+(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))

             IF(l_physadaptation) effetturbidite=EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct/ &
                 (sqrt(1+(EXTINCTION_RAD(k,i,j)/p_phyto_ChlN_ksmithextinct)**2))

            

        END SUBROUTINE shared_parameters



END MODULE shared

