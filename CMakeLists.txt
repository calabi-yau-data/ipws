cmake_minimum_required (VERSION 3.1)
project (ipws)

option(LINK_STATICALLY "Link statically" OFF)
set(LATTICE_RING "int" CACHE STRING "Lattice ring")
set(DIMENSION 4 CACHE STRING "Dimension")
set(R_NUMERATOR 1 CACHE STRING "Index numerator")
set(R_DENOMINATOR 1 CACHE STRING "Index denominator")

# Signed integer overflow check for the types 'int' and 'long'.
set(FLAGS "${FLAGS} -fsanitize=signed-integer-overflow -fsanitize-undefined-trap-on-error")
set(FLAGS "${FLAGS} -pedantic -Wall -Wextra -Wconversion")

set(FLAGS_DEBUG "-O0 -g")
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
else(${GIT_CHANGES} STREQUAL CLEAN)
  add_definitions("-DGIT_REVISION=\"${GIT_HASH}(dirty)\"")
endif(${GIT_CHANGES} STREQUAL CLEAN)

add_library(ip_check ip_check.c palp/Polynf.c palp/Coord.c palp/Rat.c palp/Vertex.c)
target_compile_options(ip_check PRIVATE -Wno-conversion -Wno-misleading-indentation -Wno-unused-but-set-variable -Wno-unused-parameter)

configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

add_executable(ipws ipws.cpp file.cpp history.cpp point.cpp settings.cpp weight_system.cpp weight_system_builder.cpp weight_system_pair.cpp)
target_include_directories(ipws PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(ipws ip_check)

if(LINK_STATICALLY)
  target_link_libraries(ipws -static) # requires: glibc-static libstdc++-static
endif(LINK_STATICALLY)
