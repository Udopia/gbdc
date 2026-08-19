cmake_minimum_required(VERSION 3.15)
project(gbdc_gate VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Reach back to the gbdc root for shared headers and sources
set(GBDC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
include_directories("${GBDC_ROOT}")

# --- CaDiCaL ---
set(CADICAL_ARCH_SUFFIX "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    string(REPLACE ";" "_" CADICAL_ARCH_SUFFIX "macos-${CMAKE_OSX_ARCHITECTURES}")
endif()

set(CADICAL_PREFIX "${CMAKE_BINARY_DIR}/solvers-${CADICAL_ARCH_SUFFIX}")
set(CADICAL_DIR ${CADICAL_PREFIX}/src/cadical_external/build)
set(CADICAL_LIB ${CADICAL_DIR}/libcadical.a)

set(CADICAL_CFLAGS "-fPIC")
set(CADICAL_CXXFLAGS "-fPIC")

if(APPLE)
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CADICAL_CFLAGS "${CADICAL_CFLAGS} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
        set(CADICAL_CXXFLAGS "${CADICAL_CXXFLAGS} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    endif()
    if(CMAKE_OSX_ARCHITECTURES)
        foreach(CADICAL_ARCH IN LISTS CMAKE_OSX_ARCHITECTURES)
            set(CADICAL_CFLAGS "${CADICAL_CFLAGS} -arch ${CADICAL_ARCH}")
            set(CADICAL_CXXFLAGS "${CADICAL_CXXFLAGS} -arch ${CADICAL_ARCH}")
        endforeach()
    endif()
endif()

if(EXISTS ${CADICAL_LIB})
    message(STATUS "CaDiCaL found at ${CADICAL_LIB}")
    add_library(cadical STATIC IMPORTED)
    set_target_properties(cadical PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")
else()
    message(STATUS "Building CaDiCaL")
    include(ExternalProject)
    ExternalProject_Add(cadical_external
        GIT_REPOSITORY https://github.com/arminbiere/cadical.git
        GIT_TAG master
        PREFIX ${CADICAL_PREFIX}
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env "CXXFLAGS=${CADICAL_CXXFLAGS}" "CFLAGS=${CADICAL_CFLAGS}" ./configure
        BUILD_COMMAND make -j
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${CADICAL_LIB}"
        BUILD_IN_SOURCE 1
    )
endif()

add_library(solver STATIC IMPORTED)
if(NOT EXISTS ${CADICAL_LIB})
    add_dependencies(solver cadical_external)
endif()
set_target_properties(solver PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIRS})

# Gate extractor + the util .cc files it needs
add_library(gate_extract OBJECT
    CNFGateFeatures.cc
    "${GBDC_ROOT}/src/util/CaptureDistribution.cc"
)
set_property(TARGET gate_extract PROPERTY POSITION_INDEPENDENT_CODE ON)
if(NOT EXISTS ${CADICAL_LIB})
    add_dependencies(gate_extract cadical_external)
endif()

add_executable(gbdc-gate main.cc)
target_link_libraries(gbdc-gate PRIVATE gate_extract solver ${LibArchive_LIBRARIES})
