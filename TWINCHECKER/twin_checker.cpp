//force to keep assert as used as error handling
#undef NDEBUG
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

template <class T>
void TwinAppCheckerChannelValues<T>::set(T value, bool isMaster)
{
    if (isMaster)
        this->masterValue = value;
    else
        this->slaveValue = value;
}

template <class T>
void TwinAppCheckerChannelValuesArray<T>::set(const T * value, size_t count, bool isMaster)
{
    assert(count < TWIN_ARRAY_BATCH_SIZE);

    if (isMaster)
        memcpy(this->masterValues, value, count * sizeof(T));
    else
        memcpy(this->slaveValues, value, count * sizeof(T));
}

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

void TwinAppCheckerChannelBarrier::waitAll(bool isMaster, bool logging)
{
    //inc enter
    const auto memoryOrder = std::memory_order_acq_rel;

    //enter
    if (logging) printf("------------- START BARRIER --------------\n");

    //case
    if (isMaster) {
        if (logging) printf("MASTER - enter\n");
        this->masterIsIn.store(1, memoryOrder);
        if (logging) printf("MASTER - wait\n");
        while (this->slaveIsIn.load() != 1){__builtin_ia32_pause();};
        if (logging) printf("MASTER - OK\n");
        this->slaveIsIn.store(0, memoryOrder);
    } else {
        if (logging) printf("SLAVE  - enter\n");
        this->slaveIsIn.store(1, memoryOrder);
        if (logging) printf("SLAVE  - wait\n");
        while (this->masterIsIn.load() != 1){__builtin_ia32_pause();};
        if (logging) printf("SLAVE  - OK\n");
        this->masterIsIn.store(0, memoryOrder);
    }

    if (logging) printf("------------- END BARRIER --------------\n");
}

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

TwinLocationInfos::TwinLocationInfos(const char * equation, size_t equation_size, int64_t locationId, int sourceLine)
{
    this->equation = equation;
    this->equation_size = equation_size;
    this->locationId = locationId;
    this->sourceLine = sourceLine;
}

TwinAppChecker::TwinAppChecker(const std::string & meetingPointFifoFile, const std::string & name, TwinAppCheckerMode mode)
{
    //init
    this->meetingPointFifoFile = meetingPointFifoFile;
    this->name = name;

    //mode
    char * env = getenv("TWIN_CHECKER_MODE");
    switch(mode)
    {
        case TWIN_MASTER:
            this->isMaster = true;
            break;
        case TWIN_SLAVE:
            this->isMaster = false;
            break;
        case TWIN_AUTO:
            assert(env != nullptr);
            assert(strcmp(env, "master") == 0 || strcmp(env, "slave") == 0);
            this->isMaster = (strcmp(env, "master") == 0);
            break;
        default:
            assert(false);
    }
    this->isMaster = isMaster;

    //solver
    this->solver.loadProcMap();

    //meet
    this->meetup();
}

TwinAppChecker::~TwinAppChecker(void)
{
    this->goodbye();
}

void TwinAppChecker::enableLoggin(void)
{
    this->logging = true;
}

void TwinAppChecker::enableStopOnFirstError(void)
{
    this->stopOnFirstError = true;
}

void TwinAppChecker::enableResolveStack(void)
{
    this->resolveStack = true;
}

void TwinAppChecker::registerSiteId(int64_t id, const std::string & sourceFile)
{
    //@todo optimize not to redo many times & report error if has different value when exist
    this->siteIdToFile[id] = sourceFile;
}

void TwinAppChecker::goodbye(void)
{
    //unmap
    int status = munmap(this->sharedMem, this->sharedSize);
    assert(status == 0);

    //unlink now we are connected
    if (this->isMaster) {
        int status = shm_unlink(this->name.c_str());
        assert(status == 0);
    }
}

void TwinAppChecker::meetup(void)
{
    //master allocate shared mem
    int shfd = -1;
    if (this->isMaster) {
        shfd = shm_open(this->name.c_str(), O_CREAT | O_RDWR, 0666);
        if (shfd < 0)
            throw std::runtime_error(std::string("Fail to create & open shm file : ") + this->name + " : " + std::string(std::strerror(errno)));
        /* configure the size of the shared memory object */
        int status = ftruncate(shfd, this->sharedSize);
        if (status < 0)
            throw std::runtime_error(std::string("Fail to ftruncate shm file : ") + this->name + " : " + std::string(std::strerror(errno)));
    }

    //open
    FILE * fp = NULL;
    if (this->isMaster)
        fp = fopen(this->meetingPointFifoFile.c_str(), "w");
    else
        fp = fopen(this->meetingPointFifoFile.c_str(), "r");

    //check
    if (fp == NULL)
        throw std::runtime_error(this->meetingPointFifoFile + " : " + std::string(std::strerror(errno)));

    //sync on fifo
    const char messageMaster[] = "webcome";
    char messageSlave[4096];
    size_t ioSize = 0;
    if (this->isMaster) {
        if (logging) printf("MASTER - write\n");
        ioSize = fwrite(messageMaster, 1, strlen(messageMaster) + 1, fp);
        fflush(fp);
        assert(ioSize == strlen(messageMaster) + 1);
        if (logging) printf("MASTER - done\n");
    } else {
        if (logging) printf("SLAVE  - read\n");
        ioSize = fread(messageSlave, 1, sizeof(messageSlave), fp);
        assert(ioSize == strlen(messageSlave) + 1);
        assert(strcmp(messageSlave, messageMaster) == 0);
        if (logging) printf("SLAVE  - read done\n");
    }

    //open mem slave
    if (this->isMaster == false) {
        if (logging) printf("SLAVE  - shmopen\n");
        shfd = shm_open(this->name.c_str(), O_RDWR, 0666);
        if (shfd < 0)
            throw std::runtime_error(std::string("Fail to open shm file : ") + this->name + " : " + std::string(std::strerror(errno)));
    }

    //mmap it
    this->sharedMem = mmap(NULL, this->sharedSize, PROT_READ|PROT_WRITE, MAP_SHARED, shfd, 0);
    if (this->sharedMem == MAP_FAILED)
        throw std::runtime_error(std::string("Fail to mmap shm file : ") + this->name + " : " + std::string(std::strerror(errno)));

    //keep
    this->sharedFd = shfd;
    this->channelMeta = static_cast<TwinAppCheckerChannelMeta*>(this->sharedMem);

    //close
    fclose(fp);
}

void TwinAppChecker::check(bool & value, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    bool masterValue = this->checkGeneric(value, "%B", this->channelMeta->bools, infos, reportError);
    if (masterValue != value)
        value = masterValue;
}

void TwinAppChecker::check(int & value, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    int masterValue = this->checkGeneric(value, "%d", this->channelMeta->integers, infos, reportError);
    if (masterValue != value)
        value = masterValue;
}

void TwinAppChecker::check(float & value, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    float masterValue = this->checkGeneric(value, "%.27g", this->channelMeta->floats, infos, reportError);
    if (value == NAN || value == -NAN) {
        fprintf(stderr, "Get NaN !\n");
        abort();
    }
    if (masterValue != value)
        value = masterValue;
}

void TwinAppChecker::check(double & value, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    double masterValue = this->checkGeneric(value, "%.27g", this->channelMeta->doubles, infos, reportError);
    if (value == NAN || value == -NAN) {
        fprintf(stderr, "Get NaN !\n");
        abort();
    }
    if (masterValue != value)
        value = masterValue;
}

void TwinAppChecker::checkArray(double * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%.27g", this->channelMeta->doublesArray, infos, reportError);
}

void TwinAppChecker::checkArray(float * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%.27g", this->channelMeta->floatsArray, infos, reportError);
}

void TwinAppChecker::checkArray(int * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%d", this->channelMeta->integersArray, infos, reportError);
}

void TwinAppChecker::checkArray(bool * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%B", this->channelMeta->boolsArray, infos, reportError);
}

void TwinAppChecker::performBacktrace(void)
{
    //stack
    void * frames[4096];
    int frameCount = 1;
    const int id_of_error_frame = 3;

    frameCount = backtrace(frames, 4096);
    for (int i = 3 ; i < frameCount - 3 ; i++) {
        const MALT::CallSite * site = this->solver.getCallSiteInfo(frames[i]);
        fprintf(stderr, "%s:%d (%s)\n", this->solver.getString(site->file).c_str(), site->line, this->solver.getString(site->function).c_str());
    }

    this->solver.solveNames();

    const MALT::CallSite * errorSite = this->solver.getCallSiteInfo(frames[id_of_error_frame]);
    FILE * fp = fopen(this->solver.getString(errorSite->file).c_str(), "r");
    if (fp != NULL)
    {
        fprintf(stderr, "-------------------------------------------------------------------\n");
        char lineContent[4096];
        for (int line = 0 ; line < errorSite->line ; line++)
            fgets(lineContent, sizeof(lineContent), fp);
        fprintf(stderr, "%s", lineContent);
    }
    fprintf(stderr, "-------------------------------------------------------------------\n\n");
    assert(errorSite->line != 0);
}

template <class T>
T TwinAppChecker::checkGenericSingleValue(T masterValue, T slaveValue, T value, const char * valueFormat, const TwinLocationInfos & infos, bool reportError)
{
    //both check
    T result = masterValue;
    if (masterValue != slaveValue) {
        char message[4096];
        char currentExe[4096];
        int status = readlink("/proc/self/exe", currentExe, sizeof(currentExe));
        assert(status > 0);
        currentExe[status] = '\0';
        T delta = slaveValue - masterValue;

        char valueStrMaster[256];
        char valueStrSlave[256];
        char valueDiff[256];
        memset(valueStrMaster, 0, 256);
        memset(valueStrSlave, 0, 256);
        snprintf(valueStrMaster, 256, valueFormat, masterValue);
        snprintf(valueStrSlave, 256, valueFormat, slaveValue);

        for (int i = 0 ; i < std::max(strlen(valueStrMaster), strlen(valueStrSlave)) ; i++) {
            if (valueStrMaster[i] == valueStrSlave[i])
                valueDiff[i] = ' ';
            else if (valueStrMaster[i] != ' ' && valueStrSlave[i] != ' ')
                valueDiff[i] = '^';
        }
        valueDiff[255] = '\0';


        char preparedFormatError[4096] = "";
        sprintf(preparedFormatError, "Not maching value in %%s :\n - master = %s\n - slave  = %s\n - err    = %%s\n - diff   = %s\n - local  = %s\n", valueFormat, valueFormat, valueFormat, valueFormat);
        sprintf(message, preparedFormatError, currentExe, masterValue, slaveValue, valueDiff, delta, value);

        //stack
        void * frames[4096];
        int frameCount = 1;
        const int id_of_error_frame = 3;
        void * frame_of_error = nullptr;

        //extract
        if (this->resolveStack) {
            frameCount = backtrace(frames, 4096);
            //char ** frameSymbols = backtrace_symbols(frames, frameCount);
            for (int i = 3 ; i < frameCount - 3 ; i++)
                this->solver.registerAddress(frames[i]);
            frame_of_error = frames[id_of_error_frame];
        } else {
            frame_of_error = (void*)infos.locationId;
        }

        if (this->alreadySeen[frame_of_error] == false)
        {
            //throw std::runtime_error(std::string(message));
            fprintf(stderr, "\n---------------------- Twin Checker Error -------------------------\n");
            fprintf(stderr, message);
            fprintf(stderr, "-------------------------------------------------------------------\n");
            char null_terminated_eq[4096];
            assert(infos.equation_size < 4096);
            strncpy(null_terminated_eq, infos.equation, std::min(4096ul, infos.equation_size));
            null_terminated_eq[infos.equation_size] = '\0';
            fprintf(stderr, "%s\n", null_terminated_eq);
            fprintf(stderr, "-------------------------------------------------------------------\n");
            fprintf(stderr, "At /home/svalat/Projects/minicroco/%s:%d\n", this->siteIdToFile[infos.locationId].c_str(), infos.sourceLine);
            fprintf(stderr, "-------------------------------------------------------------------\n");
            if (this->resolveStack)
                this->performBacktrace();
            if (this->stopOnFirstError)
                abort();
            this->alreadySeen[frame_of_error] = true;
        }
    }
}

template <class T>
T TwinAppChecker::checkGeneric(T value, const char * valueFormat, TwinAppCheckerChannelValues<T> & channel, const TwinLocationInfos & infos, bool reportError)
{
    //@todo: use atomics here
    T result;

    //master mark value
    channel.set(value, this->isMaster);

    //Barrier in
    this->channelMeta->barrierIn.waitAll(this->isMaster, this->logging);

    //display
    if (logging) {
        char preparedFormatLog[4096] = "";
        sprintf(preparedFormatLog, "%%s - value=%s\n", valueFormat);
        if (this->isMaster) {
            printf(preparedFormatLog, "MASTER", value);
        } else {
            printf(preparedFormatLog, "SLAVE", value);
        }
    }

    //both check
    result = this->checkGenericSingleValue(channel.masterValue, channel.slaveValue, value, valueFormat, infos, reportError);

    //barrier out
    this->channelMeta->barrierOut.waitAll(this->isMaster, this->logging);

    //ok
    return result;
}

template <class T>
void TwinAppChecker::checkGenericArray(T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError, bool fixable)
{
    for (size_t offset = 0 ; offset < count ; offset += TWIN_ARRAY_BATCH_SIZE) {
        size_t batch_size = std::min(count - offset, TWIN_ARRAY_BATCH_SIZE);
        const T * masterValue = this->checkGenericArrayBatch(value + offset, batch_size, valueFormat, channel, infos, reportError);
        if (fixable && masterValue != nullptr) {
            for (size_t i = 0 ; i < batch_size ; i++)
                value[offset + i] = masterValue[i];
        }
    }
}

template <class T>
void TwinAppChecker::checkGenericArray(const T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError)
{
    for (size_t offset = 0 ; offset < count ; offset += TWIN_ARRAY_BATCH_SIZE) {
        size_t batch_size = std::min(count - offset, TWIN_ARRAY_BATCH_SIZE);
        this->checkGenericArrayBatch(value + offset, batch_size, valueFormat, channel, infos, reportError);
    }
}

template <class T>
const T * TwinAppChecker::checkGenericArrayBatch(const T * value, size_t count, const char * valueFormat, TwinAppCheckerChannelValuesArray<T> & channel, const TwinLocationInfos & infos, bool reportError)
{
    //checj
    assert(value != nullptr);
    assert(count > 0);
    assert(valueFormat != nullptr);
    assert(count <= TWIN_ARRAY_BATCH_SIZE);

    //@todo: use atomics here
    T * result = nullptr;

    //master mark value
    channel.set(value, count, this->isMaster);

    //Barrier in
    this->channelMeta->barrierIn.waitAll(this->isMaster, this->logging);

    //display
    if (logging) {
        char preparedFormatLog[4096] = "";
        sprintf(preparedFormatLog, "%%s - value=%s\n", valueFormat);
        if (this->isMaster) {
            printf(preparedFormatLog, "MASTER", value);
        } else {
            printf(preparedFormatLog, "SLAVE", value);
        }
    }

    //make volatile
    const volatile T * masterValues = channel.masterValues;
    const volatile T * slaveValues = channel.slaveValues;

    //check all
    for (size_t i = 0 ; i < count ; i++) {
        T masterValue = this->checkGenericSingleValue(masterValues[i], slaveValues[i], value[i], valueFormat, infos, reportError);
        if (masterValue != value[i])
            result = channel.masterValues;
    }

    //barrier out
    this->channelMeta->barrierOut.waitAll(this->isMaster, this->logging);

    //ok
    return result;
}

}

extern "C" {

static twin::TwinAppChecker * gbl_twin_state = nullptr;

twin::TwinAppChecker * twin_init(void)
{
    //get env
    const char * meeting_point = getenv("TWIN_CHECKER_MEETPOINT");
    const char * name = getenv("TWIN_CHECKER_NAME");

    //build
    twin::TwinAppChecker * twin = new twin::TwinAppChecker(meeting_point, name, twin::TWIN_AUTO);
    //twin->enableLoggin();
    //twin->enableStopOnFirstError();
    //twin->enableResoveStack();

    //ok
    return twin;
}

void twin_check_bool(int value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(value, infos);
}

void twin_check_bool_fixable(int * value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(*value, infos);
}

void twin_check_float(float value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(value, infos);
}

void twin_check_float_fixable(float * value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(*value, infos);
}

void twin_check_double(double value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(value, infos);
}

void twin_check_double_fixable(double * value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(*value, infos);
}

void twin_check_int(int value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(value, infos);
}

void twin_check_int_fixable(int * value, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->check(*value, infos);
}

void twin_register_site(int64_t id, const char * source_file, size_t source_file_size)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    assert(source_file_size > 0 && source_file_size < 4096);
    char null_terminated[4096];
    strncpy(null_terminated, source_file, std::min(4096ul, source_file_size));
    null_terminated[source_file_size] = '\0';

    gbl_twin_state->registerSiteId(id, null_terminated);
}

void __attribute__((destructor)) twin_finish()
{
    if (gbl_twin_state != nullptr)
        delete gbl_twin_state;
}

}
