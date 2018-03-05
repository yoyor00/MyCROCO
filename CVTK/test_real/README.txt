Scripts to run all test cases defined in cppdefs.h and make some plots :

- test_real.sh : main script calling run_real.sh and plot_real.sh
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            MPI/OPENMP
            MAX_PROCS
            Usage      : test_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR] [-p PARALLEL] [-m MAX_PROC]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all
 -p PARALLEL      : Type of parallelism (MPI or OPENMP), default : no
 -p MAX_PROC      : Max number of cpus available, default : 1
 -s MPI_RUN       : mpirun command, default : mpirun

- run_real.sh : compile and run a given test_case, called from test_real or alone
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            MPI/OPENMP
            MAX_PROCS
            Usage      : test_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR] [-p PARALLEL] [-m MAX_PROC]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all
 -p PARALLEL      : Type of parallelism (MPI or OPENMP), default : no
 -p MAX_PROC      : Max number of cpus available, default : 1
 -s MPI_RUN       : mpirun command, default : mpirun 

- plot_real.sh : compile and run a given test_case, called from test_real or alone
requirement : ghostscript and matlab
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            MPI/OPENMP
            MAX_PROCS
            Usage      : test_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all

 Method :
 - copy jobcomp, TEST_CASES and cppdefs.h from the Run repository unless there already there
 - copy .h .F files from OCEAN and Run
 - Compile and Run
 - for test cases with input .nc files, if not there, they are created (attempt)     
 - LOG files stored in LOG
 - final pdf file is merged.pdf