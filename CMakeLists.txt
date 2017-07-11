cmake_minimum_required (VERSION 2.6)
project (PALP)

# set(CMAKE_C_COMPILER clang)
# set(CMAKE_CXX_COMPILER clang++)

# Signed integer overflow check for the types 'int' and 'long'.
set(FLAGS "${FLAGS} -fsanitize=signed-integer-overflow -fsanitize-undefined-trap-on-error")

set(FLAGS_DEBUG "-O0 -g -pedantic -Wall -Wextra -Wconversion")
set(FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS "${FLAGS}")
set(CMAKE_CXX_FLAGS "${FLAGS}")

set(CMAKE_C_FLAGS_DEBUG "${FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_DEBUG "${FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_RELEASE "${FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${FLAGS_RELEASE}")

set(CMAKE_CXX_STANDARD_REQUIRED 14)
set(CMAKE_CXX_STANDARD 14)

add_library(common Coord.c Polynf.c Rat.c Vertex.c)
target_compile_options(common PRIVATE -Wno-misleading-indentation)

add_executable(poly poly.c LG.c)
target_link_libraries(poly common)

add_executable(class class.c Subpoly.c Subadd.c Subdb.c)
target_link_libraries(class common)
target_compile_definitions(class PUBLIC _FILE_OFFSET_BITS=64 _LARGEFILE_SOURCE)

add_executable(cws cws.c LG.c)
target_link_libraries(cws common)

add_executable(nef nef.c E_Poly.c Nefpart.c LG.c)
target_link_libraries(nef common)

add_executable(mori mori.c MoriCone.c SingularInput.c LG.c)
target_link_libraries(mori common)

add_executable(cws0 cws0.cpp weight_system_store_set.cpp)
target_link_libraries(cws0 common)

add_executable(cws1 cws1.cpp weight_system.cpp weight_system_builder.cpp weight_system_pair.cpp)
target_include_directories(cws1 PRIVATE "${PROJECT_SOURCE_DIR}/include")

find_package(JNI)
if (JNI_FOUND)
  add_library(WsIp SHARED WsIp.c Polynf.c Coord.c Rat.c Vertex.c)
  target_include_directories(WsIp PRIVATE ${JNI_INCLUDE_DIRS})
  set_target_properties(WsIp PROPERTIES C_VISIBILITY_PRESET hidden)
endif (JNI_FOUND)
