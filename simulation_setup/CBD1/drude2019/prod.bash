#!/bin/csh
#$ -l h_rt=48:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N s1n1 #job name
#$ -m e
#$ -j y #merge console out and error in same file
#$ -pe omp 16
#$ -l gpus=1
#$ -l gpu_type=L40S

module load miniconda
conda activate openmm_plumed

#cd /projectnb/cui-buchem/sangram/Im7_simulation/drude2019/1CEI_Im7_modelled_Cterminus/openmm

set init = step3_charmm2omm
set equi_prefix = step4.1_equilibration    #previous state similar to cpt file in gromacs
set prod_prefix = step5_production
set prod_step   = step5_production_1
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str"

setenv OPENMM_CPU_THREADS 16
set FILE=${prod_step}.out     
if ( -f ${FILE} ) then
   rm -f $FILE #rmf is the alias I am using for rm -f
endif

# Equilibration = step 4.1
python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -irst ${equi_prefix}.rst -orst ${prod_step}.rst -odcd ${prod_step}.dcd > ${prod_step}.out

