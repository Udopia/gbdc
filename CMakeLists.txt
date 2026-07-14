cmake_minimum_required(VERSION 3.15)
project(gbdc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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

if (EXISTS ${CADICAL_LIB})
    message(STATUS "Cadical found at ${CADICAL_LIB}")
    add_library(cadical STATIC IMPORTED)
    set_target_properties(cadical PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")
else()
    message(STATUS "Building Cadical")
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
if (NOT EXISTS ${CADICAL_LIB})
    add_dependencies(solver cadical_external)
endif()

set_target_properties(solver PROPERTIES IMPORTED_LOCATION "${CADICAL_LIB}")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    execute_process(
        COMMAND brew --prefix libarchive
        OUTPUT_VARIABLE LIBARCHIVE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )
    set(LibArchive_INCLUDE_DIR "${LIBARCHIVE_PREFIX}/include")
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

# Per-tool commands: thin `gbd-<tool>` executables realised as symlinks to the single `gbdc`
# binary, which dispatches on its invocation name (argv[0]). See src/Main.cc.
set(GBD_TOOL_NAMES
    gbd-extract-base gbd-extract-gate gbd-extract-wcnf gbd-extract-opb
    gbd-checksani gbd-isohash gbd-isohash2 gbd-identify
    gbd-cnf2kis gbd-cnf2bip gbd-sanitize gbd-normalize)

foreach(toolname IN LISTS GBD_TOOL_NAMES)
    add_custom_command(TARGET gbdctool POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                "$<TARGET_FILE_NAME:gbdctool>" "$<TARGET_FILE_DIR:gbdctool>/${toolname}"
        COMMENT "Creating symlink ${toolname} -> gbdc")
endforeach()

# Install the binary and per-tool symlinks. Prefer the wheel scripts dir (scikit-build-core)
# so the tools land on PATH; fall back to a conventional bin directory otherwise.
if(DEFINED SKBUILD_SCRIPTS_DIR)
    set(GBD_TOOL_BINDIR "${SKBUILD_SCRIPTS_DIR}")
else()
    set(GBD_TOOL_BINDIR "bin")
endif()

install(TARGETS gbdctool RUNTIME DESTINATION "${GBD_TOOL_BINDIR}")

# For local/source installs, also expose the per-tool `gbd-<tool>` commands as symlinks to the
# single `gbdc` binary. Wheel builds (scikit-build-core, SKBUILD_SCRIPTS_DIR defined) ship only the
# single `gbdc` binary to keep the distribution small; gbd invokes tools as `gbdc <tool>` there.
if(NOT DEFINED SKBUILD_SCRIPTS_DIR)
    install(CODE "
        set(_bindir \"${GBD_TOOL_BINDIR}\")
        if(NOT IS_ABSOLUTE \"\${_bindir}\")
            set(_bindir \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/\${_bindir}\")
        endif()
        foreach(toolname ${GBD_TOOL_NAMES})
            execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink gbdc \"\${_bindir}/\${toolname}\")
        endforeach()
    ")
endif()

pybind11_add_module(gbdc src/gbdlib.cc)
target_link_libraries(gbdc PUBLIC ${LIBS} solver util extract transform)
install(TARGETS gbdc DESTINATION .)

add_test(NAME Test_StreamBuffer COMMAND "test/tests_streambuffer")
# add_test(NAME Test_Feature_Extraction COMMAND "test/tests_feature_extraction")
add_test(NAME Test_StreamCompressor COMMAND "test/tests_streamcompressor")
add_test(NAME Test_GBDLib COMMAND "test/tests_gbdlib")
add_test(NAME Test_IsoHash2 COMMAND "test/tests_isohash2")