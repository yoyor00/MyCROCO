/*
!                       Single Column Example
!                       ====== ====== =======
! Seven sets up are encompassed:
! Choose ONE experiment by nml testcase_name: 
! KATO_PHILLIPS , WILLIS_DEARDORFF, DIURNAL_CYCLE, FORCED_EKBBL,
!FORCED_DBLEEK, FORCED_NONROTBBL,FORCED_OSCNONROTBBL  
*/
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define NEW_S_COORD
# define UV_COR
# define SOLVE3D
# undef  LMD_MIXING
# define GLS_MIXING
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define EW_PERIODIC
# define NS_PERIODIC


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
