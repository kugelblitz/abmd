cmake_minimum_required(VERSION 3.8)
project(dde)

set(CMAKE_CXX_STANDARD 11)

include_directories(src test .)

set(SOURCE_FILES test/main.c src/poly.c test/plot.c src/coeffs.c src/abm.c
        src/queue.c src/api.c src/rk.c src/rk.h)
add_executable(dde ${SOURCE_FILES})

target_compile_definitions(dde PUBLIC USE_LONG_DOUBLE)
