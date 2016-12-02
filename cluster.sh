#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l walltime=00:05:00
#PBS -q siwir
#PBS -M up28itig@fau.de -m abe
#PBS -N test

./script.sh
