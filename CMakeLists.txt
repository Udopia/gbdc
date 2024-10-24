cmake_minimum_required(VERSION 3.15)
project(gbdc VERSION 1.0 LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CADICAL_DIR ${CMAKE_BINARY_DIR}/solvers/src/cadical_external/build)
set(CADICAL_LIB ${CADICAL_DIR}/libcadical.a)
if (EXISTS ${CADICAL_LIB})
    message(STATUS "Cadical already built")
    add_library(cadical STATIC IMPORTED)
    set_target_properties(cadical PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")
else()
    message(STATUS "Building Cadical")
    include(ExternalProject)
    ExternalProject_Add(cadical_external
        GIT_REPOSITORY https://github.com/arminbiere/cadical.git
        GIT_TAG master
        PREFIX solvers
        CONFIGURE_COMMAND ./configure "CXXFLAGS=-fPIC"
        BUILD_COMMAND make -j
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${CADICAL_LIB}"
        BUILD_IN_SOURCE 1
    )
endif()
    
add_library(solver STATIC IMPORTED)
add_dependencies(solver cadical_external)
set_target_properties(solver PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    execute_process(
        COMMAND brew --prefix libarchive
        OUTPUT_VARIABLE LIBARCHIVE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )
    set(LibArchive_INCLUDE_DIRS "${LIBARCHIVE_PREFIX}/include")
endif()

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIRS})
set(LIBS ${LIBS} md5 ${LibArchive_LIBRARIES})

include_directories(gbdc PUBLIC "${PROJECT_SOURCE_DIR}")

add_subdirectory("src")
add_subdirectory("test")

add_executable(gbdctool src/Main.cc)
target_link_libraries(gbdctool PUBLIC ${LIBS} solver util extract transform)
# target_include_directories(gbdctool PUBLIC "${PROJECT_SOURCE_DIR}")
set_property(TARGET gbdctool PROPERTY OUTPUT_NAME gbdc)

pybind11_add_module(gbdc src/gbdlib.cc)
target_link_libraries(gbdc PUBLIC ${LIBS} solver util extract transform)
install(TARGETS gbdc DESTINATION .)

add_test(NAME Test_StreamBuffer COMMAND "test/tests_streambuffer")
add_test(NAME Test_Feature_Extraction COMMAND "test/tests_feature_extraction")
add_test(NAME Test_StreamCompressor COMMAND "test/tests_streamcompressor")
add_test(NAME Test_GBDLib COMMAND "test/tests_gbdlib")