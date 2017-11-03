cmake_minimum_required(VERSION 3.8.2)


project(AMS301-Projet1)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CFLAGS -Wall -O3)

add_executable(seq
  seq.c
  )

target_link_libraries(seq
  m
  )

add_executable(para
  justParallel.c
  )

target_link_libraries(para
  m
  )
