cmake_minimum_required(VERSION 3.15)
project(LRmethod C)

set(CMAKE_C_STANDARD 99)

add_executable(LRmethod main_01_03.c sim_01_03.c evc_01_03.c task_01_03.h)
target_link_libraries(LRmethod m)
add_executable(test testgen_s.c)
target_link_libraries(test m)