cmake_minimum_required(VERSION 3.10)

project(cnftools VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIRS})
set(LIBS ${LIBS} md5 ${LibArchive_LIBRARIES})

add_subdirectory("lib/md5")

add_executable(cnftools Main.cc)

target_link_libraries(cnftools PUBLIC ${LIBS})

target_include_directories(cnftools PUBLIC "${PROJECT_SOURCE_DIR}")