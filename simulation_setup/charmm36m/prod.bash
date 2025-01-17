#!/bin/bash -l
#$ -l h_rt=120:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N prod1charmm #job name
#$ -e prod1charmm.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe mpi_28_tasks_per_node 84

module purge
module load openmpi/3.1.4
module load cuda/9.2
module load gromacs/2018.3

mpirun -np 84 gmx_mpi mdrun -deffnm step5_production -s step5_production_extend_till_2.5us.tpr -cpi step5_production.cpt
