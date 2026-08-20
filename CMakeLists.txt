cmake_minimum_required(VERSION 3.15)
project(gbdc VERSION 1.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Prefer Homebrew's keg-only libarchive for local/source builds; CI wheel builds
    # supply libarchive via CMAKE_PREFIX_PATH, so a missing brew formula is not fatal.
    execute_process(
        COMMAND brew --prefix libarchive
        OUTPUT_VARIABLE LIBARCHIVE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE LIBARCHIVE_BREW_RESULT
        ERROR_QUIET
    )
    if(LIBARCHIVE_BREW_RESULT EQUAL 0 AND EXISTS "${LIBARCHIVE_PREFIX}/include/archive.h")
        set(LibArchive_INCLUDE_DIR "${LIBARCHIVE_PREFIX}/include")
    endif()
endif()

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIRS})
set(LIBS ${LIBS} md5 ${LibArchive_LIBRARIES})

include_directories(gbdc PUBLIC "${PROJECT_SOURCE_DIR}")

add_subdirectory("src")
add_subdirectory("test")

add_executable(gbdctool src/Main.cc)
target_link_libraries(gbdctool PUBLIC ${LIBS} util extract transform)
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
target_link_libraries(gbdc PUBLIC ${LIBS} util extract transform)
install(TARGETS gbdc DESTINATION .)

add_test(NAME Test_StreamBuffer COMMAND "test/tests_streambuffer")
# add_test(NAME Test_Feature_Extraction COMMAND "test/tests_feature_extraction")
add_test(NAME Test_StreamCompressor COMMAND "test/tests_streamcompressor")
add_test(NAME Test_GBDLib COMMAND "test/tests_gbdlib")
add_test(NAME Test_IsoHash2 COMMAND "test/tests_isohash2")