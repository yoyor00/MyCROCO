
set(CMAKE_C_FLAGS "-O0")
set(CMAKE_CXX_FLAGS "-O0")

# build lib
include_directories(./from-malt)
include_directories(./from-malt/from-fftw)
include_directories(./from-malt/iniparser/src)
set(sources_from_malt ./from-malt/SymbolSolver.cpp
    ./from-malt/portability/CompilerGNU.cpp
    ./from-malt/portability/LinuxProcMapReader.cpp
    ./from-malt/portability/OSUnix.cpp
    ./from-malt/json/JsonState.cpp
    ./from-malt/json/ConvertToJson.cpp
    ./from-malt/common/Helpers.cpp
    ./from-malt/common/Options.cpp
    ./from-malt/common/Debug.cpp
    ./from-malt/common/CodeTiming.cpp
    ./from-malt/common/NoFreeAllocator.cpp
    ./from-malt/common/SimpleAllocator.cpp
    ./from-malt/common/FormattedMessage.cpp
    ./from-malt/iniparser/src/dictionary.c
    ./from-malt/iniparser/src/iniparser.c)
add_library(twin_checker SHARED twin_checker.cpp ${sources_from_malt})

add_executable(test_twin_checker test_twin_checker.cpp)
target_link_libraries(test_twin_checker twin_checker)