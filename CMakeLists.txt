cmake_minimum_required(VERSION 3.15)
project(CodeSynthesis)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(FetchContent)


set(CMAKE_CXX_STANDARD 11)

#SOURCE FILES ADDED HERE
set(CodeSynthesis_SOURCES
        src/utils/Utils.cpp
        src/utils/Utils.h
	src/Visitors.cc
	src/Visitors.h
	src/Permute.h
        src/CodeSynthesis.cc 
	src/CodeSynthesis.h)
set(CodSynthesis_t_SOURCES
        test/CodeSynthesisUnitTest.cpp
        test/unit_tests_main.cpp
       )

include_directories(src/)
enable_testing()
add_compile_options(-g)

#Add thread dependencies for gtest
find_package(Threads)
include_directories(${Threads_INCLUDE_DIR})



#ADD EXTERNAL IEGEN LIB
ExternalProject_Add(iegenlib_in
        GIT_REPOSITORY "https://github.com/BoiseState-AdaptLab/IEGenLib.git"
        GIT_TAG "6828a4940f75c734da5f67116c1a5271610901ac"
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/lib/installed
        BUILD_COMMAND ${MAKE})



# Add depeendence on GTEST

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.8.0
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()


# INCLUDE EXTERNAL
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/include)
link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/lib" )
link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/installed/lib" )
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/gtest/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/installed/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/include/iegenlib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/isl/include)

# BUILD CodeSynthesis_Driver
add_executable(CodeSynthesis_Driver ${CodeSynthesis_SOURCES} src/drivers/code_synthesis_main.cpp)


add_dependencies(CodeSynthesis_Driver iegenlib_in)
target_link_libraries(CodeSynthesis_Driver gtest iegenlib isl codegen omega )




#Set Target properties
set_target_properties(CodeSynthesis_Driver
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

#ADD KERNEL DRIVER
add_executable(CodeSynthesis_Kernels ${CodeSynthesis_SOURCES} src/drivers/conversion_kernels.cpp)
add_dependencies(CodeSynthesis_Kernels iegenlib_in)
target_link_libraries(CodeSynthesis_Kernels gtest iegenlib isl codegen omega )




#Set Target properties
set_target_properties(CodeSynthesis_Kernels
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

#ADD TEST
add_executable(CodeSynthesis_t 
	${CodeSynthesis_SOURCES} 
	${CodSynthesis_t_SOURCES} )
add_dependencies(CodeSynthesis_t iegenlib_in)
#Tell included gtest headers to NOT use tuple support.
add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
add_test(NAME unit_tests COMMAND CodeSynthesis_t)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

target_link_libraries(CodeSynthesis_t gtest iegenlib isl gmp codegen omega)
set_target_properties(CodeSynthesis_t
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
        )

add_subdirectory(scripts)
