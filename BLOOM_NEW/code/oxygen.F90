!&E---------------------------------------------------------------------
!&E                 ***  incellwat_bloom_oxygen  ***
!&E
!&E ** Purpose : Calcul de la variable oxygene
!&E
!&E       !  2010-02    (A.Chapelle) Original code
!&E
!&E      use from general modele : WIND_SPEED
!&E      use from BIOLink  :
!&E
!&E      use from basic bloom modele: c, epn,ibbenth,ibsurf, effetchaleur,ratio_mgO_to_mumolN_anabol
!&E               rationdiat,rationdino,rationnano ,rationpsnz,rationkarenia ,rationphaeocystiscolo,rationphaeocystiscell
!&E               effetlumierediat,effetlumieredino,effetlumierenano,..
!&E               ratio_mgO_to_mumolN_catabol,reminazdeteau,xnitrifeau,temper,tempabs
!&E      OUTPUT : diag_3d_wat(irk_diag(id_oxy_sat + dc
!&E---------------------------------------------------------------------


MODULE oxygen

#include "cppdefs.h"
#include "coupleur_define_BIOLink.h"

    USE module_BIOLink
    USE comBIOLink
    USE comBIOLink_helping
    USE COMBLOOM
    USE shared
    


    IMPLICIT NONE

    PUBLIC :: oxygen_dynamics

    contains
        SUBROUTINE oxygen_dynamics(i,j,k,c,dc,ibsurf)

            INTEGER, intent(in)                 :: i,j,k, ibsurf
            REAL(KIND=rsh),DIMENSION(nv_state), intent(in)  :: c
            REAL(KIND=rsh),DIMENSION(nv_state), intent(inout)  :: dc
            REAL(KIND=rsh) :: gst,o2sat,oair
            


            ! reoxygenation
            ! -------------
            gst=-173.4292_rsh+249.6329_rsh*100.0_rsh/tempabs+143.3483_rsh*log(tempabs/100.0_rsh)  &
                -21.8492_rsh*tempabs/100.0_rsh+sali*(-0.033096_rsh+0.014259_rsh*tempabs/100.0_rsh &
                -0.0017_rsh*(tempabs/100.0_rsh)**2) 
            o2sat=1.429_rsh*exp(gst)

            !Formule de Riley et Stephan (1988)
            oair=(0.64_rsh+0.0256_rsh*(WIND_SPEED(i,j)/0.447_rsh)**2)*(o2sat-c(iv_oxygen))/epn

            !Formule de Wanninkhof (1992)
            !schmidt=(1800.6_rsh-120.1_rsh*temper+3.7818_rsh*temper**2-0.047608_rsh*temper**3)*(1.0_rsh+0.00314_rsh*sali)
            !le facteur 0.31 cm/h est converti en 24*0.31/100=0.0744 m/jour
            !oair=0.0744_rsh*(WIND_SPEED(i,j)**2)/sqrt(schmidt/660.0_rsh)*(o2sat-c(iv_oxygen))/epn

            !Formule de Wanninkhof & McGillis (1999)
            !le facteur 0.0283 cm/h est converti en 24*0.0283/100=0.00679 m/jour
            !oair=0.00679_rsh*(WIND_SPEED(i,j)**3)/sqrt(schmidt/660.0_rsh)*(o2sat-c(iv_oxygen))/epn

            ! Evolution de l oxygene dissous
            ! ------------------------------
            dc(iv_oxygen) = dc(iv_oxygen) + oair*ibsurf

            ! pourcentage de saturation de l oxygene
            diag_3d_wat(irk_diag(id_oxy_sat),k,i,j)=100.0_rsh*c(iv_oxygen)/o2sat



        END SUBROUTINE oxygen_dynamics



END MODULE oxygen


















