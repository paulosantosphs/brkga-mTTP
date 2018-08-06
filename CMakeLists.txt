cmake_minimum_required(VERSION 3.9)
project(brkga_mTTP)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wextra -Wall -Weffc++ -ansi -pedantic -Woverloaded-virtual -Wcast-align -Wpointer-arith -std=c++11")

add_executable(brkga_mTTP main.cpp mTTPInstance.cpp mTTPInstance.h mTTPDecoder.cpp mTTPDecoder.h brkgaAPI/Population.cpp mTTPChromosome.cpp mTTPChromosome.h)