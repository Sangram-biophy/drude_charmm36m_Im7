#!/bin/bash -l
#$ -l h_rt=120:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N em_1CEI #job name
#$ -e em_1CEI.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe mpi_28_tasks_per_node 28

module purge
module load openmpi/3.1.4
module load cuda/9.2
module load gromacs/2018.3

mpirun -np 28 gmx_mpi mdrun -deffnm step4.0_minimization -s step4.0_minimization.tpr
