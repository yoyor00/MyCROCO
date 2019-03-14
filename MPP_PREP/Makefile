NETCDF_ROOT=$(NETCDF_DIR)
NETCDF_INC=$(NETCDF_ROOT)/include
NETCDF_LIB=$(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
#
F90=ifort
F90FLAGS= -O0 -C -debug -Wunused -traceback -w
LDFLAGS= 
#
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90))
EXEC=mpp_optimiz

$(EXEC): $(SRCS)
	$(F90) $(LDFLAGS) -o $@ $< -L$(NETCDF_LIB)

%.o: %.f90
	$(F90) $(F90FLAGS) -c  $< -I$(NETCDF_INC) 

clean:
	rm -f *.o *.mod $(EXEC)
