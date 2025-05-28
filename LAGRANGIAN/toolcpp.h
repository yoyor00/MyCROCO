#ifdef MPI
#define USE_MPI USE
#define IF_MPI IF
#define ELSE_MPI ELSE
#define ENDIF_MPI ENDIF
#define CALL_MPI CALL
#define STOP_MPI STOPavoir
#define MASTER (mynode == 0)
#else
#define MASTER 0
#define USE_MPI !USE
#define IF_MPI !IF
#define ELSE_MPI !ELSE
#define ENDIF_MPI !ENDIF
#define CALL_MPI !CALL
#define STOP_MPI !STOPavoir
#endif

#ifdef key_printdbg
#define PRINT_DBG PRINT
#else
#define PRINT_DBG !PRINT
#endif