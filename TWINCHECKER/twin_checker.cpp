//force to keep assert as used as error handling
#undef NDEBUG
//stdc
#include "twin_checker.hpp"

#define TWIN_ARRAY_BATCH_SIZE 4096ul

namespace twin
{

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
        shm_unlink(this->name.c_str());
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
    if (masterValue != value && fixable)
        value = masterValue;
}

void TwinAppChecker::check(int & value, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    int masterValue = this->checkGeneric(value, "%d", this->channelMeta->integers, infos, reportError);
    if (masterValue != value && fixable)
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
    if (masterValue != value && fixable)
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
    if (masterValue != value && fixable)
        value = masterValue;
}

void TwinAppChecker::checkArray(double * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%.27g", this->channelMeta->doublesArray, infos, reportError, fixable);
}

void TwinAppChecker::checkArray(float * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%.27g", this->channelMeta->floatsArray, infos, reportError, fixable);
}

void TwinAppChecker::checkArray(int * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%d", this->channelMeta->integersArray, infos, reportError, fixable);
}

void TwinAppChecker::checkArray(bool * values, size_t count, const TwinLocationInfos & infos, bool fixable)
{
    bool reportError = (fixable || this->stopOnFirstError);
    this->checkGenericArray(values, count, "%B", this->channelMeta->boolsArray, infos, reportError, fixable);
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
    //twin->enableResolveStack();

    //ok
    return twin;
}

//////////////////////////////////////// BOOL

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
    gbl_twin_state->check(*value, infos, true);
}

void twin_check_bool_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, false);
}

void twin_check_bool_fixable_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, true);
}

//////////////////////////////////////// FLOAT

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
    gbl_twin_state->check(*value, infos, true);
}

void twin_check_float_array(float * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, false);
}

void twin_check_float_fixable_array(float * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, true);
}

//////////////////////////////////////// DOUBLE

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
    gbl_twin_state->check(*value, infos, true);
}

void twin_check_double_array(double * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, false);
}

void twin_check_double_fixable_array(double * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, true);
}

//////////////////////////////////////// INT

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
    gbl_twin_state->check(*value, infos, true);
}

void twin_check_integer_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, false);
}

void twin_check_integer_fixable_array(int * values, int count, const char * equation, size_t equation_size, int64_t location_id, int source_line)
{
    //lazy init
    if (gbl_twin_state == nullptr)
        gbl_twin_state = twin_init();

    //build infos
    twin::TwinLocationInfos infos(equation, equation_size, location_id, source_line);

    //check
    gbl_twin_state->checkArray(values, count, infos, true);
}

//////////////////////////////////////// FILE NAMES

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

//////////////////////////////////////// DESTRUCATOR

void __attribute__((destructor)) twin_finish()
{
    if (gbl_twin_state != nullptr)
        delete gbl_twin_state;
}

}
