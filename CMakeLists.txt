cmake_minimum_required (VERSION 3.1)
project (ipws)

set(LATTICE_RING "int" CACHE STRING "Lattice ring")
set(DIMENSION 4 CACHE STRING "Dimension")
set(R_NUMERATOR 1 CACHE STRING "Index numerator")
set(R_DENOMINATOR 1 CACHE STRING "Index denominator")

set(SUFFIX ${DIMENSION}-${R_NUMERATOR}-${R_DENOMINATOR})

option(LINK_STATICALLY "Link statically" OFF)
if(LINK_STATICALLY)
  link_libraries(-static) # valgrind reports problems with this!?
endif()

# Signed integer overflow check for the types 'int' and 'long'.
set(FLAGS "${FLAGS} -fsanitize=signed-integer-overflow -fsanitize-undefined-trap-on-error")
set(FLAGS "${FLAGS} -pedantic -Wall -Wextra -Wconversion -Wno-maybe-uninitialized")

set(FLAGS_DEBUG "-O0 -g") # -fsanitize=address
set(FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS "${FLAGS}")
set(CMAKE_CXX_FLAGS "${FLAGS}")

set(CMAKE_C_FLAGS_DEBUG "${FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_RELEASE "${FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${FLAGS_RELEASE}")

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_HASH)
git_local_changes(GIT_CHANGES)
if(${GIT_CHANGES} STREQUAL CLEAN)
  add_definitions("-DGIT_REVISION=\"${GIT_HASH}\"")
else()
  add_definitions("-DGIT_REVISION=\"${GIT_HASH}(dirty)\"")
endif()

include_directories("${PROJECT_SOURCE_DIR}/src/include")

configure_file(
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

find_package(Boost 1.46 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_library(palp src/palp_glue.cpp src/palp/Polynf.c src/palp/Coord.c
  src/palp/Rat.c src/palp/Vertex.c)
target_compile_options(palp PRIVATE -Wno-conversion -Wno-unused-parameter
  -Wno-misleading-indentation -Wno-unused-but-set-variable -Wno-parentheses
  -Wno-int-in-bool-context)

add_executable(ipws-${SUFFIX} src/ipws.cpp src/buffered_reader.cpp
  src/buffered_writer.cpp src/history.cpp src/point.cpp src/settings.cpp
  src/weight_system.cpp src/weight_system_builder.cpp
  src/weight_system_pair.cpp src/polytope.cpp)
target_link_libraries(ipws-${SUFFIX} palp)

add_executable(add-half-${SUFFIX} src/add-half.cpp src/buffered_reader.cpp
  src/buffered_writer.cpp src/history.cpp src/point.cpp src/settings.cpp
  src/weight_system.cpp src/weight_system_builder.cpp
  src/weight_system_pair.cpp src/polytope.cpp)
target_link_libraries(add-half-${SUFFIX} palp)
