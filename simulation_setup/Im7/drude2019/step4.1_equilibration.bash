#!/bin/csh
#$ -l h_rt=24:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N Im7drude #job name
#$ -e Im7drude.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe omp 28
#$ -l gpus=1
#$ -l gpu_c=3.5

module load miniconda
conda activate openmm_plumed

cd /projectnb/cui-buchem/sangram/Im7_simulation/drude2019/1CEI_Im7_modelled_Cterminus/openmm

set init = step3_charmm2omm
set equi_prefix = step4.1_equilibration
set prod_prefix = step5_production
set prod_step   = step5
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str"

setenv OPENMM_CPU_THREADS 28

# Equilibration = step 4.1
python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -irst step4_equilibration.rst -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcd > ${equi_prefix}.out
