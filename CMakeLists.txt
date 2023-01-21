cmake_minimum_required(VERSION 3.10)
project(AVBMC)

# Set compiler
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Compilation flags used
set(CMAKE_CXX_FLAGS_RELEASE "-ggdb3 -O3 -o -DNDEBUG -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -ggdb3 -g -O0 -std=c+11")


# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set default build configuration
if( NOT CMAKE_BUILD_TYPE )
   SET(CMAKE_BUILD_TYPE Release)
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/boundary
    ${PROJECT_SOURCE_DIR}/src/constraint
    ${PROJECT_SOURCE_DIR}/src/forcefield
    ${PROJECT_SOURCE_DIR}/src/moves
    ${PROJECT_SOURCE_DIR}/src/rng
    ${PROJECT_SOURCE_DIR}/src/sampler
)

# Specify which functions to compile
file(GLOB SOURCES
    src/*.cpp
    src/boundary/*.cpp
    src/constraint/*.cpp
    src/forcefield/*.cpp
    src/moves/*.cpp
    src/rng/*.cpp
    src/sampler/*.cpp
)

set(PROJECT_EXECUTABLE avbmc)


# Add executable file VMaChine
add_executable(${PROJECT_EXECUTABLE} ${SOURCES})
