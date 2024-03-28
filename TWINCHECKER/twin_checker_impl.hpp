#ifndef TWIN_CHECKER_IMPL_HPP
#define TWIN_CHECKER_IMPL_HPP

//header
#include "twin_checker.hpp"

namespace twin
{

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
    assert(count <= TWIN_ARRAY_BATCH_SIZE);

    if (isMaster)
        memcpy(this->masterValues, value, count * sizeof(T));
    else
        memcpy(this->slaveValues, value, count * sizeof(T));
}

template <class T>
T TwinAppChecker::checkGenericSingleValue(T masterValue, T slaveValue, T value, const char * valueFormat, const TwinLocationInfos & infos, bool reportError)
{
    //both check
    T result = masterValue;
    bool tooLargeGap = (fabs(masterValue - slaveValue) > 1e-07);
    if (masterValue != slaveValue || tooLargeGap) {
        char message[4096];
        static char currentExe[4096];
        static int status = readlink("/proc/self/exe", currentExe, sizeof(currentExe));
        assert(status > 0);
        currentExe[status] = '\0';
        T delta = slaveValue - masterValue;

        char valueStrMaster[256];
        char valueStrSlave[256];
        char valueDiff[256];
        memset(valueDiff, 0, 256);
        memset(valueStrMaster, 0, 256);
        memset(valueStrSlave, 0, 256);
        snprintf(valueStrMaster, 255, valueFormat, masterValue);
        snprintf(valueStrSlave, 255, valueFormat, slaveValue);

        for (int i = 0 ; i < std::min(strlen(valueStrMaster), strlen(valueStrSlave)) ; i++) {
            if (valueStrMaster[i] == valueStrSlave[i])
                valueDiff[i] = ' ';
            else if (valueStrMaster[i] != ' ' && valueStrSlave[i] != ' ' && valueStrMaster[i] != ' ')
                valueDiff[i] = '^';
        }
        valueDiff[255] = '\0';

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
            char preparedFormatError[4096] = "";
            sprintf(preparedFormatError, "Not maching value in %%s :\n - master = %s\n - slave  = %s\n - err    = %%s\n - diff   = %s\n - local  = %s\n", valueFormat, valueFormat, valueFormat, valueFormat);
            sprintf(message, preparedFormatError, currentExe, masterValue, slaveValue, valueDiff, delta, value);

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
            fprintf(stderr, "At %s:%d\n", this->siteIdToFile[infos.locationId].c_str(), infos.sourceLine);
            fprintf(stderr, "-------------------------------------------------------------------\n");
            if (this->resolveStack)
                this->performBacktrace();
            tooLargeGap = false;
            if (this->stopOnFirstError || tooLargeGap)
                abort();
            this->alreadySeen[frame_of_error] = true;
        }
    }

    return result;
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
        //fprintf(stderr, "Rand : %zu[%zu] (%zu)\n", offset, batch_size, count);
        const T * masterValue = this->checkGenericArrayBatch(value + offset, batch_size, valueFormat, channel, infos, reportError);
        /*if (masterValue != nullptr)
            fprintf(stderr, "!!!!!!!!!!!!!!!!!!!! ARRAY ERROR %p\n", masterValue);*/
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
        T retMasterValue = this->checkGenericSingleValue(masterValues[i], slaveValues[i], value[i], valueFormat, infos, reportError);
        if (retMasterValue != value[i])
            result = channel.masterValues;
    }

    //barrier out
    this->channelMeta->barrierOut.waitAll(this->isMaster, this->logging);

    //ok
    return result;
}

}

#endif //TWIN_CHECKER_IMPL_HPP
