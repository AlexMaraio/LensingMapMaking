CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(thingy CXX)

# Set the flags for all builds to include aggressive optimisation settings and increased warnings
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -flto -Wall -Wextra -fPIC")

# Find OpenMP, MPFR, and Boost as required modules for the code
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}  "/home/amaraio/Documents/PhD/Codes/LensingMapMaking/lib/cpp")
FIND_PACKAGE(MPFR 3.0 REQUIRED)

FIND_PACKAGE(OpenMP)

FIND_PACKAGE(Boost REQUIRED)

# Compile the code as a shared library and as an executable
ADD_LIBRARY(thingy SHARED alm_to_cl.cpp ComputeAlms.h)
ADD_EXECUTABLE(thingy2  alm_to_cl.cpp )

TARGET_LINK_LIBRARIES(thingy2  ${MPFR_LIBRARIES} OpenMP::OpenMP_CXX)
TARGET_COMPILE_OPTIONS(thingy2 PRIVATE -std=c++17 -fopenmp -lpthread)
