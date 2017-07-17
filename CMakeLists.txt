cmake_minimum_required (VERSION 3.1)
project (cws)

option(LINK_STATICALLY "Link statically" OFF)

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

add_library(wsip wsip.c palp/Polynf.c palp/Coord.c palp/Rat.c palp/Vertex.c)
target_compile_options(wsip PRIVATE -Wno-conversion -Wno-misleading-indentation -Wno-unused-but-set-variable -Wno-unused-parameter)

add_executable(cws cws.cpp file.cpp point.cpp settings.cpp weight_system.cpp weight_system_builder.cpp weight_system_pair.cpp)
target_include_directories(cws PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(cws wsip)

if(LINK_STATICALLY)
  target_link_libraries(cws -static) # requires: glibc-static libstdc++-static
endif(LINK_STATICALLY)
