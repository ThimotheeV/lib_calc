cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_CXX_STANDARD 17)

project(simulateur VERSION 0.0.1 LANGUAGES CXX)

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#add path search for #include <>
include_directories ("${PROJECT_SOURCE_DIR}/source")
include_directories(BEFORE SYSTEM "${PROJECT_SOURCE_DIR}")

#include subprojet in main cmake projet
add_subdirectory(source)
add_subdirectory(unit_test)

add_executable(libstat main.cpp)
target_link_libraries(libstat libsource)