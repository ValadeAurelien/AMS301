cmake_minimum_required(VERSION 2.8.12)


project(AMS301-Projet1)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CFLAGS -Wall -O3)

add_executable(seq seq.c)

add_executable(parallel justParallel.c)

target_link_libraries(seq
  m
)

target_link_libraries(parallel
  m
)

