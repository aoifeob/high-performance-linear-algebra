cmake_minimum_required(VERSION 3.26)
project(high_performance_linear_algebra)

set(CMAKE_CXX_STANDARD 14)

add_executable(high_performance_linear_algebra
        matrix-multiplication/matrix_mul.c
        matrix-multiplication/testing-openmp/simple_parallel_sum.c)