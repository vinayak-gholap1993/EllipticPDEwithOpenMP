#!/bin/bash
#PBS -l nodes=2:ppn=32
#PBS -l walltime=00:15:00
#PBS -q siwir
#PBS -M up28itig@fau.de -m abe
#PBS -N test

set OMP_NUM_THREADS=16

cd ./script.sh
