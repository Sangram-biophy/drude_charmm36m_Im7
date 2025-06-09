#!/bin/csh
#$ -l h_rt=24:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N CBD1drude #job name
#$ -e CBD1drude.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe omp 16
#$ -l gpus=1
#$ -l gpu_type=L40S

source ~/.bashrc
module load cuda
module load miniconda
conda activate openmm-8.2

cd /projectnb/cui-buchem/sangram/CBD1_simulation/drude2019/set1

set init = step3_charmm2omm
set equi_prefix = step4.1_equilibration
set prod_prefix = step5_production
set prod_step   = step5
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str"

setenv OPENMM_CPU_THREADS 16

# Equilibration = step 4.1
python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -irst step4_equilibration.rst -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcd > ${equi_prefix}.out
