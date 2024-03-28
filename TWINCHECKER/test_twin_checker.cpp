#include <cstdlib>
#include <sys/wait.h>
#include "twin_checker.hpp"

twin::TwinAppChecker * gbl_checker = nullptr;

#define COUNT 100

void init_array(double * array, size_t count, double base)
{
    for (size_t i = 0 ; i < count ; i ++)
        array[i] = base + (double)i;
}

void sum_array(double * a, double * b, double * c, size_t count)
{
    for (size_t i = 0 ; i < count ; i ++)
        a[i] = a[i] + b[i] + c[i];
    twin_check_double_fixable_array(a, count, "a(:)", 4, __LINE__, __LINE__);
    twin_check_double_fixable_array(b, count, "b(:)", 4, __LINE__, __LINE__);
    twin_check_double_fixable_array(c, count, "c(:)", 4, __LINE__, __LINE__);
}

void sum_array_2(double * a, double * b, double * c, size_t count)
{
    for (size_t i = 0 ; i < count ; i ++)
        a[i] = a[i] + b[i] + c[i];
    twin_check_double_array(a, count, "a(:)", 4, __LINE__, __LINE__);
    twin_check_double_array(b, count, "b(:)", 4, __LINE__, __LINE__);
    twin_check_double_array(c, count, "c(:)", 4, __LINE__, __LINE__);
}

void check_arrays(double * a, double * b, double * c, size_t count)
{
    twin_check_double_array(a, count, "a(:)", 4, __LINE__, __LINE__);
    twin_check_double_array(b, count, "b(:)", 4, __LINE__, __LINE__);
    twin_check_double_array(c, count, "c(:)", 4, __LINE__, __LINE__);
}

void run_master()
{
    const size_t count = COUNT;
    double * array1 = (double*)malloc(sizeof(double) * count);
    double * array2 = (double*)malloc(sizeof(double) * count);
    double * array3 = (double*)malloc(sizeof(double) * count);

    //same
    init_array(array1, count, 1.0);
    init_array(array2, count, 100.0);
    init_array(array3, count, 10100.0);
    
    sum_array(array1, array2, array3, count);
    sum_array_2(array1, array2, array3, count);
    check_arrays(array1, array2, array3, count);
}

void run_slave()
{
    const size_t count = COUNT;
    double * array1 = (double*)malloc(sizeof(double) * count);
    double * array2 = (double*)malloc(sizeof(double) * count);
    double * array3 = (double*)malloc(sizeof(double) * count);

    //same
    init_array(array1, count, 1.0);
    init_array(array2, count, 100.0);
    init_array(array3, count, 10000.0);

    sum_array(array1, array2, array3, count);
    sum_array_2(array1, array2, array3, count);
    check_arrays(array1, array2, array3, count);
}

int main()
{
    setenv("TWIN_CHECKER_NAME", "croco-test", 1);
    setenv("TWIN_CHECKER_MEETPOINT","/tmp/twin_checker_croco_test", 1);
    mkfifo("/tmp/twin_checker_croco_test", 0666);

    pid_t child_pid = fork();
    if (child_pid == 0) {
        setenv("TWIN_CHECKER_MODE", "master", 1);
        run_master();
    } else {
        setenv("TWIN_CHECKER_MODE", "slave", 1);
        run_slave();
        int wstatus;
        waitpid(child_pid, &wstatus, WUNTRACED | WCONTINUED);
    }

    printf("Finished\n");

    return EXIT_SUCCESS;
}
