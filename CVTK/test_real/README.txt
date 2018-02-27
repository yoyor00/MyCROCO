Scripts to run all test cases defined in cppdefs.h and make some plots :

- test_real.sh : main script calling run_real.sh and plot_real.sh
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
=> change if needed

- run_real.sh : compile and run a given test_case, called from test_real or alone
input : arg1 : name of the test case
        arg2 : root repository for croco (default is the same thant CVTK) 

- plot_real.sh : compile and run a given test_case, called from test_real or alone
requirement : ghossctipt and matlab
input : arg1 : name of the test case
        arg2 : root repository for croco (default is the same thant CVTK) 

 Method :
 - copy jobcomp, TEST_CASES and cppdefs.h from the Run repository unless there already there
 - Compile and Run
 - for test cases with input .nc files, if not there, they are created (attempt)     
 - LOG files stored in LOG
 - final pdf files is mergeD.pdf