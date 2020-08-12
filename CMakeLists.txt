cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_CXX_STANDARD 17)

project(simulateur VERSION 0.0.1 LANGUAGES CXX)

include(CTest)
include(lib/Catch/Catch.cmake)

enable_testing()

option(CMAKE_FUNC_TEST OFF)
# pour pouvoir lui donner des options de compilation (e.g. dans sources/CMakeLists.txt : target_compile_options() )
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

#add path search for #include <> in CXX files for the compiler
include_directories ("${PROJECT_SOURCE_DIR}/sources")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Catch")
include_directories ("${PROJECT_SOURCE_DIR}/lib/random")
include_directories(BEFORE SYSTEM "${PROJECT_SOURCE_DIR}")

#include subprojet in main cmake projet, so that cmake can find other cmake files.
add_subdirectory(sources)

if (CMAKE_FUNC_TEST)
    add_subdirectory(func_tests)
else()
    add_subdirectory(unit_tests)
endif()

# fait le .o
add_executable(libstat main.cpp)
# lie les .o/.a ...
target_link_libraries(libstat libsources)