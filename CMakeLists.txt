cmake_minimum_required(VERSION 3.17)
project(CodeSynthesis)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

set(CMAKE_CXX_STANDARD 11)

#SOURCE FILES ADDED HERE
set(CodeSynthesis_SOURCES src/smt/MinimalSatisfiablity.cpp src/smt/MinimalSatisfiablity.h)
set(CodSynthesis_t_SOURCES test/smt/MinimalUnitTest.cpp test/unit_tests_main.cpp)

include_directories(src/)
enable_testing()
add_compile_options(-g)

#Add thread dependencies for gtest
find_package(Threads)
include_directories(${Threads_INCLUDE_DIR})




#ADD EXTERNAL IEGEN LIB
ExternalProject_Add(iegenlib_in
        DOWNLOAD_COMMAND ""
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/lib/installed
        BUILD_COMMAND ${MAKE})

# INCLUDE EXTERNAL
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/include)
link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/lib" )
link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/installed/lib" )
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/gtest/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/IEGenLib/lib/installed/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/installed/include/iegenlib)

add_subdirectory(lib)

# BUILD CodeSynthesis_Driver
add_executable(CodeSynthesis_Driver ${CodeSynthesis_SOURCES} src/main.cpp)
add_dependencies(CodeSynthesis_Driver iegenlib_in)
target_link_libraries(CodeSynthesis_Driver gtest iegenlib isl)

#Set Target properties
set_target_properties(CodeSynthesis_Driver
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

#ADD TEST
add_executable(CodeSynthesis_t ${CodeSynthesis_SOURCES} ${CodSynthesis_t_SOURCES} )
add_dependencies(CodeSynthesis_t iegenlib_in)
#Tell included gtest headers to NOT use tuple support.
add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
add_test(NAME unit_tests COMMAND CodeSynthesis_t)

target_link_libraries(CodeSynthesis_t gtest iegenlib isl)
set_target_properties(CodeSynthesis_t
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
        )
