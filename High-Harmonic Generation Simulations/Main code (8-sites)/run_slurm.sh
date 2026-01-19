#!/bin/bash
#SBATCH -J GBae_HHG         # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -p g5
#SBATCH -w n042
#SBATCH -N 1              # total number of node
#SBATCH -n 1              # total number of node

## HPC ENVIRONMENT 
. /etc/profile.d/TMI.sh
##

mpicc -g -Wall -o hhg_cal hhg_main.c hhg_mpi.c hhg_mpi.h rk4_kxky.c rk4_kxky.h -llapacke -llapack -lgsl -lgslcblas -lm -O3
mpiexec -n $SLURM_NTASKS ./hhg_cal
