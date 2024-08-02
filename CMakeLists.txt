cmake_minimum_required(VERSION 3.10)

project(gbdc VERSION 1.0)

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

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND NOT DEFINED LA_PREFIX)
    execute_process(
        COMMAND brew --prefix libarchive
        OUTPUT_VARIABLE LIBARCHIVE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )
    message("Libarchive prefix: ${LIBARCHIVE_PREFIX}")
    set(LibArchive_INCLUDE_DIR "${LIBARCHIVE_PREFIX}/include")
    list(INSERT CMAKE_PREFIX_PATH 0 "${LIBARCHIVE_PREFIX}")
endif()

if(DEFINED LA_PREFIX)
    message(STATUS "Setting Libarchive include dir: ${LA_PREFIX}/include")
    set(LibArchive_INCLUDE_DIR "${LA_PREFIX}/include")
    list(INSERT CMAKE_PREFIX_PATH 0 "${LA_PREFIX}")
endif()

find_package(LibArchive REQUIRED)
if(LibArchive_FOUND)
    message(STATUS "LibArchive found:")
    message(STATUS "LibArchive_INCLUDE_DIRS: ${LibArchive_INCLUDE_DIRS}")
    message(STATUS "LibArchive_LIBRARIES: ${LibArchive_LIBRARIES}")
endif()

include_directories(${LibArchive_INCLUDE_DIRS})
set(LIBS ${LIBS} md5 ${LibArchive_LIBRARIES})

include_directories(gbdc PUBLIC "${PROJECT_SOURCE_DIR}")

add_subdirectory("src")
add_subdirectory("lib/md5")

add_executable(gbdc src/Main.cc)
add_dependencies(gbdc solver)
target_link_libraries(gbdc PUBLIC ${LIBS} solver $<TARGET_OBJECTS:gates> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:transform> $<TARGET_OBJECTS:extract>)

target_include_directories(gbdc PUBLIC "${PROJECT_SOURCE_DIR}")

add_test(NAME Test_StreamBuffer COMMAND "src/test/tests_streambuffer")
add_test(NAME Test_Feature_Extraction COMMAND "src/test/tests_feature_extraction")
add_test(NAME Test_StreamCompressor COMMAND "src/test/tests_streamcompressor")
add_test(NAME Test_ThreadPool COMMAND "src/test/tests_threadpool")
add_test(NAME Test_GBDLib COMMAND "src/test/tests_gbdlib")