cmake_minimum_required(VERSION 3.26)
project(high_performance_linear_algebra)

set(CMAKE_CXX_STANDARD 14)

add_executable(high_performance_linear_algebra
        matrix-multiplication/atlas_blocked_matrix_mul.c
        matrix-multiplication/blas_matrix_mul.c
        matrix-multiplication/simple_matrix_mul.c)
