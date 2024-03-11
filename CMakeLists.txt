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
        GIT_REPOSITORY git@github.com:arminbiere/cadical.git 
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

if(APPLE)
    set(LibArchive_INCLUDE_DIR "/usr/local/opt/libarchive/include")
endif()

find_package(LibArchive REQUIRED)
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
add_test(NAME Test_CNFBaseFeatures COMMAND "src/test/tests_cnfbasefeatures")
