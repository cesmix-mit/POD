################################################################################                   
#                    Proper Orthogonal Descriptors (POD)
#                           CESMIX-MIT Project  
# 
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
################################################################################
# mkdir exec
# cd exec
# cmake -D CMAKE_C_COMPILER=gcc-10 -D CMAKE_CXX_COMPILER=g++-10 ../cmake 
# cmake --build .
################################################################################

cmake_minimum_required(VERSION 3.12)

project(pod LANGUAGES CXX)

get_filename_component(POD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

set(POD_SOURCE_DIR     ${POD_DIR}/src)
set(POD_MAIN_DIR       ${POD_DIR}/src)
file(GLOB MAIN_SOURCES ${POD_MAIN_DIR}/main.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

add_executable(pod ${MAIN_SOURCES})
target_compile_options(pod PRIVATE -std=c++11 -Wall -ffast-math -O3 -DNDEBUG) 
target_link_libraries(pod PRIVATE  ${lapackblas_libraries})        


