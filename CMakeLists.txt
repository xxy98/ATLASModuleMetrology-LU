cmake_minimum_required(VERSION 3.23)
project(Ana)

find_path(EIGEN_DIR NAMES signature_of_eigen3_matrix_library
        //PATHS C:/Users/Xiangyu/Documents/CppPkg/eigen-3.4.0/
        PATHS ./external/eigen-3.4.0/
        PATH_SUFFIXES eigen3 eigen
        )
if(NOT EXISTS ${EIGEN_DIR})
    message(FATAL_ERROR "Please check that the set Eigen directory is valid!")
endif()
include_directories(${EIGEN_DIR})

include_directories(C:/Users/Xiangyu/Documents/CppPkg/date-3.0.1/include)
include_directories(C:/Users/Xiangyu/Documents/CppPkg/date-3.0.1/)

add_subdirectory(external/date-3.0.1)
include_directories(external/date-3.0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(Ana main.cpp structures.h functions.h mathematics.h measurements.h)
target_link_libraries(Ana date)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
