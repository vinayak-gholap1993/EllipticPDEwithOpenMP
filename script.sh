#!/bin/bash
clear
make clean
#cmake CMakeLists.txt
make
set OMP_NUM_THREADS = 16
./rbgs 2049 2049 500
 
