cmake_minimum_required (VERSION 2.6)
project (cws)

# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)

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

set(CMAKE_CXX_STANDARD_REQUIRED 14)
set(CMAKE_CXX_STANDARD 14)

add_library(wsip wsip.c palp/Polynf.c palp/Coord.c palp/Rat.c palp/Vertex.c)
target_compile_options(wsip PRIVATE -Wno-conversion -Wno-misleading-indentation -Wno-unused-but-set-variable -Wno-unused-parameter)

add_executable(cws cws.cpp point.cpp weight_system.cpp weight_system_builder.cpp weight_system_pair.cpp)
target_include_directories(cws PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(cws wsip)
