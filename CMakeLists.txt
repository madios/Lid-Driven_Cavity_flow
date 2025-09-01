cmake_minimum_required(VERSION 3.31)
project(QTTest LANGUAGES C CXX)

include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/release-1.12.1.zip
)

FetchContent_MakeAvailable(googletest)

set(CMAKE_C_COMPILER  "/opt/homebrew/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

include_directories(/opt/homebrew/opt/libomp/include)
link_directories(/opt/homebrew/opt/libomp/lib)
add_compile_options(-fopenmp)
set(OPENMP_LINK_TARGET omp)



find_package(Qt6 COMPONENTS
        Core
        Gui
        Widgets
        Charts
        REQUIRED
)
find_package(nlohmann_json CONFIG REQUIRED)

add_executable(QTTest main.cpp
        LidDrivenCavityFlow.cpp
        LidDrivenCavityFlow.h
        SimpleAlgorithme.cpp
        SimpleAlgorithme.h
        parameters.cpp
        parameters.h
        )
target_link_libraries(QTTest
    PRIVATE
        Qt::Core
        Qt::Gui
        Qt::Widgets
        Qt::Charts
        nlohmann_json::nlohmann_json
        ${OPENMP_LINK_TARGET}
)

enable_testing()

add_executable(run_tests
        Tests/test_Simple.cpp
        SimpleAlgorithme.cpp
        parameters.cpp
)

target_link_libraries(run_tests
        PRIVATE
        gtest_main
        gtest
        pthread
        Qt::Core
        Qt::Gui
        Qt::Widgets
        Qt::Charts
        nlohmann_json::nlohmann_json
        ${OPENMP_LINK_TARGET}
)

include(GoogleTest)
gtest_discover_tests(run_tests)