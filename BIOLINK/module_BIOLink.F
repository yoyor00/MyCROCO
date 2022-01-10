! $Id: module_BIOLink  2019-10- bthouvenin 
#include "cppdefs.h"

!!! donner acces aux dimensions, variables, parametres dont a besoin BIOLink
!!! + 
#if defined BIOLink

      module module_BIOLink

     !&E==========================================================================
     !&E                   ***    coupleur_BIOLink dimensions hydro ***
     !&E==========================================================================



     ! si non MARS : parametres utilises pour dimensionner les tableaux spatiaux et autres, nombre de variables ..
     !               acces aux donnees fixes du modele hydro (non passe en argument pendant les calculs)
     !               les donnees normalement definis par le modele hydro mais qui n existent pas en fait
     !               doivent etre declares et definis dans comBIOLink


     ! modele CROCO : ce fichier est remplace par module_BIOLink.F
     !                qui permet de passer du format fixe (F) au format libre F90.
     ! CROCO : il y a acces grace au  module substance acces  aux donnees fixes
     !         par des include de  param.h , ncscrum.h 
      USE module_substance

#include "ocean2d.h"
#include "ocean3d.h"
#if defined BLOOM && defined GLS_MIXING
#include "mixing.h"
#endif


      end module module_BIOLink
#else
      module module_BIOLink_empty
      end module module_BIOLink_empty
      
     !
#endif /* BIOLink */

