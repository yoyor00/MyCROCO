.SUFFIXES: .o .f90

FC = mpif90
FFLAGS = -c -O3 -frecord-marker=4 -fPIC -I/usr/include

.f90.o:
	$(FC) $(FFLAGS) $*.f90 -o $*.o

SRCS = mg_cst.f90 mg_mpi.f90 mg_tictoc.f90 mg_namelist.f90 mg_netcdf_out.f90 mg_grids.f90 mg_mpi_exchange.f90 mg_gather.f90 mg_horiz_grids.f90 mg_vert_grids.f90 mg_projection.f90 mg_relax.f90 mg_intergrids.f90 mg_solvers.f90 mg_autotune.f90 nhmg.f90

OBJS = $(SRCS:.f90=.o)

all: $(OBJS)
	$(FC) $(OBJS) $(LDFLAGS) -shared -o libnhmg.so.0
	ar cr libnhmg.a $(OBJS)

clean:
	rm -f *.o *.mod *.so
