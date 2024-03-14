#ifndef TWIN_CHECKER_HPP
#define TWIN_CHECKER_HPP

//stdc
#include <stdio.h>
//std C compat
#include <cstring>
#include <string>
#include <cerrno>
#include <cassert>
#include <cmath>
//std C++
#include <map>
#include <atomic>
#include <stdexcept>
#include <algorithm>
//Linux shm
#include <sys/mman.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
//Linux
#include <unistd.h>
#include <execinfo.h>
//from MALT
#include "from-malt/SymbolSolver.hpp"

#define TWIN_ARRAY_BATCH_SIZE 4096ul

namespace twin
{

enum TwinAppCheckerMode
{
    TWIN_MASTER,
    TWIN_SLAVE,
    TWIN_AUTO
};

struct TwinLocationInfos
{
    //members
    TwinLocationInfos(const char * equation, size_t equation_size, int64_t locationId, int sourceLine);
    //vars
    const char * equation{nullptr};
    size_t equation_size{0};
    int64_t locationId{-1};
    int sourceLine{-1};
};

template <class T>
struct TwinAppCheckerChannelValuesArray
{
    //members
    void set(const T * value, size_t count, bool isMaster);
    //vars
    T masterValues[TWIN_ARRAY_BATCH_SIZE];
    T slaveValues[TWIN_ARRAY_BATCH_SIZE];
};

template <class T>
struct TwinAppCheckerChannelValues
{
    //members
    void set(T value, bool isMaster);
    //vars
    volatile T masterValue{0.0};
    char __padd_3[64-sizeof(float)];
    volatile T slaveValue{0.0};
    char __padd_4[64-sizeof(float)];
};

struct TwinAppCheckerChannelBarrier
{
    //members
    void waitAll(bool isMaster, bool logging = false);
    //in
    std::atomic<int> masterIsIn{0};
    char __padd_1[64-sizeof(int)];
    std::atomic<int> slaveIsIn{0};
    char __padd_2[64-sizeof(int)];
};

struct TwinAppCheckerChannelMeta
{
    //in
    TwinAppCheckerChannelBarrier barrierIn;
    TwinAppCheckerChannelBarrier barrierOut;

    //dat
    TwinAppCheckerChannelValues<double> doubles;
    TwinAppCheckerChannelValues<float> floats;
    TwinAppCheckerChannelValues<int> integers;
    TwinAppCheckerChannelValues<bool> bools;

    //arrays
    TwinAppCheckerChannelValuesArray<double> doublesArray;
    TwinAppCheckerChannelValuesArray<float> floatsArray;
    TwinAppCheckerChannelValuesArray<int> integersArray;
    TwinAppCheckerChannelValuesArray<bool> boolsArray;
};

class TwinAppChecker
{
    public:
        TwinAppChecker(const std::string & meetingPointFifoFile, const std::string & name, TwinAppCheckerMode mode = TWIN_AUTO);
        ~TwinAppChecker(void);
        void registerSiteId(int64_t id, const std::string & sourceFile);
        void check(double & value, const TwinLocationInfos & infos, bool fixable = false);
        void check(float & value, const TwinLocationInfos & infos, bool fixable = false);
        void check(int & value, const TwinLocationInfos & infos, bool fixable = false);
        void check(bool & value, const TwinLocationInfos & infos, bool fixable = false);

        void checkArray(double * values, size_t count, const TwinLocationInfos & infos, bool fixable = false);
        void checkArray(float * values, size_t count, const TwinLocationInfos & infos, bool fixable = false);
        void checkArray(int * values, size_t count, const TwinLocationInfos & infos, bool fixable = false);
        void checkArray(bool * values, size_t count, const TwinLocationInfos & infos, bool fixable = false);

        template <class T> T checkGeneric(T value, const char * valueFormat, TwinAppCheckerChannelValues<T> & channel, const TwinLocationInfos & infos, bool reportError);
        template <class T> const T * checkGenericArrayBatch(const T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError);
        template <class T> void checkGenericArray(T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError, bool fixable);
        template <class T> void checkGenericArray(const T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError);
        void enableLoggin(void);
        void enableStopOnFirstError(void);
        void enableResolveStack(void);
    private:
        void meetup(void);
        void goodbye(void);
        void performBacktrace(void);
        template <class T> T checkGenericSingleValue(T masterValue, T slaveValue, T value, const char * valueFormat, const TwinLocationInfos & infos, bool reportError);
    private:
        bool logging{false};
        bool stopOnFirstError{false};
        bool resolveStack{false};
        std::string meetingPointFifoFile{};
        std::string name{};
        bool isMaster{false};
        size_t sharedSize{1024*1024};
        void * sharedMem{nullptr};
        int sharedFd{-1};
        TwinAppCheckerChannelMeta * channelMeta{nullptr};
        MALT::SymbolSolver solver;
        std::map<void*, bool> alreadySeen;
        std::map<int64_t, std::string> siteIdToFile;
};

}

extern "C" {

    twin::TwinAppChecker * twin_init(void);

    //////////////////////////////////////// BOOL

    void twin_check_bool(int value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_bool_fixable(int * value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_bool_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_bool_fixable_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);

    //////////////////////////////////////// FLOAT

    void twin_check_float(float value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_float_fixable(float * value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_float_array(float * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_float_fixable_array(float * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);

    //////////////////////////////////////// DOUBLE

    void twin_check_double(double value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_double_fixable(double * value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_double_array(double * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_double_fixable_array(double * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);

    //////////////////////////////////////// INT

    void twin_check_int(int value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_int_fixable(int * value, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_integer_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);
    void twin_check_integer_fixable_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line);

    //////////////////////////////////////// FILE NAMES

    void twin_register_site(int64_t id, const char * source_file, size_t source_file_size);

}

//template implementation
#include "twin_checker_impl.hpp"

#endif //TWIN_CHECKER_HPP
