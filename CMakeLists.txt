
cmake_minimum_required(VERSION 2.8.9)


# Defining the project name
project (ex02)

set(CMAKE_CXX_FLAGS "-Wall -Winline -Wshadow -pedantic -fopenmp -march=native -ffast-math -mavx")

add_definitions(-std=c++11)

# Bringing in the include directories
include_directories(include)


#Adding the sources using the set command
set(SOURCES src/PDE.cc test/main.cc)  

#Setting the executable file
add_executable(rbgs ${SOURCES})

