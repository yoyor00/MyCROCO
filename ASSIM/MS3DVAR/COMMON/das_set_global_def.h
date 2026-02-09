/* This is "global_definitions.h": It contains a set of predetermined
 macro definitions which are inserted into the individual files by
 C-preprocessor. General user is strongly discouraged from attempts
 to modify anything below this line.
------------------------------------------------------------------ */

#ifdef MPI
# define GLOBAL_1D_XI -1:Lm+2+padd_X
# define GLOBAL_1D_ETA -1:Mm+2+padd_E
#else
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define GLOBAL_1D_XI -1:Lm+2+padd_X
#   define GLOBAL_1D_ETA -1:Mm+2+padd_E
#  else
#   define GLOBAL_1D_XI -1:Lm+2+padd_X
#   define GLOBAL_1D_ETA 0:Mm+1+padd_E
#  endif
# else
#  ifdef NS_PERIODIC
#   define GLOBAL_1D_XI 0:Lm+1+padd_X
#   define GLOBAL_1D_ETA -1:Mm+2+padd_E
#  else
#   define GLOBAL_1D_XI 0:Lm+1+padd_X
#   define GLOBAL_1D_ETA 0:Mm+1+padd_E
#  endif
# endif
#endif



