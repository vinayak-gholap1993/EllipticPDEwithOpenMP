#!/bin/bash
clear
make clean
#cmake CMakeLists.txt
make
set OMP_NUM_THREADS=32
./rbgs 2049 2049 500 
echo 
