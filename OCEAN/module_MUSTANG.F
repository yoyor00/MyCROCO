! $Id: module_MUSTANG  2018-08- bthouvenin 
#include "cppdefs.h"

!!! donner acces aux dimensions, variables, parametres dont a besoin MUSTANG
!!! + 
#if defined MUSTANG

      module module_MUSTANG

   !&E==========================================================================
   !&E                   ***    coupleur_MUSTANG dimensions hydro ***
   !&E==========================================================================



   ! si non MARS : parametres utilises pour dimensionner les tableaux spatiaux et autres, nombre de variables ..
   !               acces aux donnees fixes du modele hydro (non passe en argument pendant les calculs)
   !               les donnees normalement definis par le modele hydro mais qui n existent pas en fait
   !               doivent etre declares et definis dans comMUSTANG


   ! modele CROCO : ce fichier est remplace par module_MUSTANG.F
   !                qui permet de passer du format fixe (F) au format libre F90.
   ! CROCO : il y a acces grace au  module substance acces  aux donnees fixes
   !         par des include de  param.h , ncscrum.h 


      USE module_substance


      end module module_MUSTANG
#else
      module module_MUSTANG_empty
      end module module_MUSTANG_empty
      
   !
#endif /* MUSTANG */

