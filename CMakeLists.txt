
cmake_minimum_required(VERSION 2.8.9)


# Defining the project name
project (EllipticPDEwithOpenMP-master)
										#check
set(CMAKE_CXX_FLAGS "-Wall -Winline -Wshadow -Wextra -DNDEBUG -pedantic -ffast-math -march=native -mavx -funit-at-a-time -fstrength-reduce -fexpensive-optimizations -fschedule-insns2 -funroll-all-loops -O3 -fopenmp")

add_definitions(-std=c++11)

# Bringing in the include directories
include_directories(include)


#Adding the sources using the set command
set(SOURCES ./src/PDE.cc ./test/main.cc)  

#Setting the executable file
add_executable(rbgs ${SOURCES})

