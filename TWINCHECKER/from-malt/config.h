#ifndef MALT_CONFIG_H
#define MALT_CONFIG_H

//versions
#define MALT_VERSION "1.2.2"
#define MALT_VERSION_NOTE "-dev"
#define MALT_JSON_FORMAT_VERSION "1.1"

//select one of the spinlock mode
#define MALT_PORTABILITY_SPINLOCK_PTHREAD
/* #undef MALT_PORTABILITY_SPINLOCK_DUMMY */

//select one of the mutex mode
#define MALT_PORTABILITY_MUTEX_PTHREAD
/* #undef MALT_PORTABILITY_MUTEX_DUMMY */

//select one of the OS mode
#define MALT_PORTABILITY_OS_UNIX

//select one of the compiler mode
#define MALT_PORTABILITY_COMPILER_GNU

//select one of the backtrace mode
/* #undef MALT_PORTABILITY_BACKTRACE_GLIBC */
/* #undef MALT_PORTABILITY_BACKTRACE_LIBUNWIND */

//Use code timing to help optimizing MALT code
/* #undef MALT_ENABLE_CODE_TIMING */

//if have libunwind
#define MALT_HAVE_LIBUNWIND

//arm RDTSC port selection
/* #undef HAVE_ARMV7A_CNTVCT */
/* #undef HAVE_ARMV7A_PMCCNTR */
/* #undef HAVE_ARMV8_CNTVCT_EL0 */
/* #undef HAVE_ARMV8_PMCCNTR_EL0 */

#endif //MALT_CONFIG_H
