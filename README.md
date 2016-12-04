# Simulation and Scientific Computing 1 
#Assignment 2
#Group: Mayank Patel, Vinayak Gholap

To run the code follow the steps. 
1. Open the terminal in the code directory. 
2. Then to compile all the files type "cmake ." (without the inverted commas). This will generate the Makefile for the code and set up the environment.
3. Type "make clean" to remove any initial output files. 
4. Type "make". Then type the command to set the OpenMP threads "set OMP_NUM_THREADS= 32".
5. We now will have an executable file "rbgs". 
6. Finally run the executable with inputs as nx, ny, c. "./rbgs nx ny c" 
7. Alternatively you can do above all step 3,4,5,6 using "script.sh". Open the script, you can change the paramters nx, ny and c in sciprt and also num of threads. Run the script "./script.sh" to compile all the .h and .cpp files. Appropriate flags are given in CMakefile as required. 
8. An output file named solution.txt will be generated. This file contains x y and data value at this point in grid.
9. A Performance report is present named "Performance plots for Elliptic PDE Solver using Open MP.pdf".

