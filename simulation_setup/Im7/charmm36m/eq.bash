#!/bin/bash -l
#$ -l h_rt=120:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N eq_1CEI #job name
#$ -e eq_1CEI.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe mpi_28_tasks_per_node 56

module purge
module load openmpi/3.1.4
module load cuda/9.2
module load gromacs/2018.3

mpirun -np 56 gmx_mpi mdrun -deffnm step4.2_equilibration -s step4.2_equilibration.tpr
