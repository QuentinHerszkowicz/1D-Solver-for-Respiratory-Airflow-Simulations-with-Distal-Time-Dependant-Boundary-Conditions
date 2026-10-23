cmake_minimum_required(VERSION 3.5)

project(OneDSolver C CXX)

# Linear solver used for the finite element system.
#   skyline : built-in serial solver, no external dependency (default, used in the paper)
#   csparse : serial sparse solver, sources bundled in Code/Source/sparse/csparse
#   superlu : parallel sparse solver, requires an external SuperLU_MT installation
set(sparseSolverType "skyline" CACHE STRING "Linear solver: skyline, csparse or superlu")
set_property(CACHE sparseSolverType PROPERTY STRINGS skyline superlu csparse)

if(sparseSolverType STREQUAL "skyline")
  add_definitions("-DUSE_SKYLINE")
elseif(sparseSolverType STREQUAL "superlu")
  add_definitions("-DUSE_SUPERLU")
  set(SUPERLU_DIR " " CACHE PATH "SuperLU_MT installation folder")
elseif(sparseSolverType STREQUAL "csparse")
  add_definitions("-DUSE_CSPARSE")
else()
  message(FATAL_ERROR "Unknown sparseSolverType '${sparseSolverType}' (use skyline, csparse or superlu)")
endif()

# Source files
if(sparseSolverType STREQUAL "skyline")
  file(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx")
  file(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h")
elseif(sparseSolverType STREQUAL "superlu")
  file(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.cxx"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.cxx")
  file(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.h"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.h")
elseif(sparseSolverType STREQUAL "csparse")
  file(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.cxx"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.cxx")
  file(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.h"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.h")
endif()

# Compiler flags (same as the build used for the published simulations)
set(CMAKE_CXX_FLAGS "-g -m64 -O3 -std=c++0x -fPIC")

# Executable goes to <build>/bin/OneDSolver
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(Code/Source)
