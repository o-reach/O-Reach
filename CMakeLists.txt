cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

add_definitions()

# tweak compiler flags
#CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
#if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
#  add_definitions(-funroll-loops)
#endif()
#CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
#if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
#  add_definitions(-fno-stack-limit)
#endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-Wextra COMPILER_SUPPORTS_WEXTRA)
if(COMPILER_SUPPORTS_WEXTRA)
  add_definitions(-Wextra)
endif()
CHECK_CXX_COMPILER_FLAG(-pedantic COMPILER_SUPPORTS_PEDANTIC)
if(COMPILER_SUPPORTS_PEDANTIC)
  add_definitions(-pedantic)
endif()
#CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
#if(COMPILER_SUPPORTS_MARCH_NATIVE)
#  add_definitions(-march=native)
#endif()
#CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
#if(COMPILER_SUPPORTS_FPERMISSIVE)
#  add_definitions(-fpermissive)
#endif()

# check dependencies
# find_package(MPI REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/flat_hash_map)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/tools)


set(LIB_COMMON_SOURCE_FILES
	extern/KaHIP/lib/algorithms/strongly_connected_components.cpp
  extern/KaHIP/lib/io/graph_io.cpp
  extern/KaHIP/lib/tools/random_functions.cpp
  extern/KaHIP/lib/tools/macros_assertions.h
  extern/KaHIP/lib/tools/macros_common.h
  extern/KaHIP/lib/definitions.h
  extern/KaHIP/lib/data_structure/union_find.h
  app/query_generator.h
  lib/algorithms/topological_sort.cpp
  extern/argtable3-3.0.3/argtable3.c
  extern/flat_hash_map/flat_hash_map.hpp
  extern/flat_hash_map/bytell_hash_map.hpp
  extern/murmur_hash/MurmurHash3.cpp
  )
add_library(libcommon OBJECT ${LIB_COMMON_SOURCE_FILES})


set(LIBREACHABILITY_SOURCE_FILES OBJECT 
  lib/algorithms/bfs.cpp
  lib/reachability/oreach.cpp
  lib/reachability/matrix/full_reach.h
)
add_library(libreachability ${LIBREACHABILITY_SOURCE_FILES})



## reachability
add_executable(reachability app/reachability.cpp $<TARGET_OBJECTS:libcommon> $<TARGET_OBJECTS:libreachability> )
target_link_libraries(reachability ${OpenMP_CXX_LIBRARIES})
install(TARGETS reachability DESTINATION bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address,undefined,leak -g -Wall")
# target_link_options(reachability PUBLIC -g)
# target_compile_options(reachability PUBLIC -g)
#target_compile_definitions(reachability PUBLIC NUM_SUPPORTS=16)
#target_compile_definitions(libreachability PUBLIC NUM_SUPPORTS=16)
#target_compile_definitions(reachability PUBLIC MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS=75)
#target_compile_definitions(libreachability PUBLIC MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS=75)
#add_compile_definitions(RDEBUG)
#add_compile_definitions(NUM_SUPPORTS=16)
#add_compile_definitions(MAX_SUPPORT_CANDIDATE_FACTOR_LEVELS=75)
# add_compile_definitions(CHECK_EDGE)
# add_compile_definitions(FULL_REACH)
