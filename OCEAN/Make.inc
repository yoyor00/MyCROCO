# [2] Machine dependent definitions for SGI and SGI/CRAY Origin
#==== ======= ========= =========== === === === ======== ======
# CPP           Name of the C preprocessor
# CPPFLAGS      Flags to the C preprocessor
# CFT           Name of the fortran compiler to use
# FFLAGS        Flags to the fortran compiler
# LDR           Program to load the objects into an executable
# LDFLAGS       Flags to the loader
#
      CPP = $(CPP1)
 CPPFLAGS = $(CPPFLAGS1)
      CFT = $(CFT1)
      LDR = $(CFT)
   FFLAGS = $(FFLAGS1)
 FFLAGS90 = $(FFLAGS1)
  LDFLAGS =
     LCDF = $(LDFLAGS1)
