#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -P mcprot
#$ -j y
#$ -m e
#$ -N prodCBD1_C36
#$ -pe omp 16
#$ -l gpu_type=L40S
#$ -l gpus=1

module purge
module load openmpi/3.1.4_gnu-8.1
module load cuda/11.2
module load gromacs/2021.5

export OMP_NUM_THREADS=16
#gmx grompp -f ../bench.mdp -c equil.gro -r equil.gro -p ../system4x.top -n ../system4x.ndx -o bench.tpr -maxwarn 1
gmx mdrun -pin on -v -nt 16 -ntomp 16 -ntmpi 1 -deffnm step5_production -s step5_production.tpr -maxh 48

