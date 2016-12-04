#!/bin/bash
clear
make clean
#cmake CMakeLists.txt
make
for i in 32 33 1024 1025 2048 2049
do
echo nx = $i ny = $i
#echo nx = $i ny = $i >> "plot.txt"
for threads in 1 2 4 8 16 32
do
echo
echo Num. of Threads $threads >> "plot.txt"
set OMP_NUM_THREADS=$threads
./rbgs $i $i 500 >> "plot.txt"
echo Run with number of threads $threads finished
echo 
done
done
 
