#!/bin/csh
#$ -l h_rt=48:00:00 #hard time limit for job
#$ -P cui-buchem #project name
#$ -N Im7drude #job name
#$ -e Im7drude.err #error file
#$ -j y #merge console out and error in same file
#$ -V
#$ -pe omp 28
#$ -l gpus=1
#$ -l gpu_memory=32G
#$ -l gpu_c=6.0

module load miniconda
conda activate openmm_plumed

cd /projectnb/cui-buchem/sangram/Im7_simulation/drude2019/1CEI_Im7_modelled_Cterminus/openmm

set init = step3_charmm2omm
set equi_prefix = step4.1_equilibration
set prod_prefix = step5_production
set prod_step   = step5
#set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -b ${init}.str"


# Equilibration = step 4.1
#python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -irst step4_equilibration.rst -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcd > ${equi_prefix}.out

# Production
# The OpenMM check point file (.chk) cannot be used in a different machine environment.
# So please make sure if you are using the same GPU and CUDA version of machine while doing additional
# production steps with the check point file.
set file = round.txt
set cntmax = 22

#Check current round
if ( -f ${file} ) then
    set cnt = (`cat ${file}`)
    echo "round file exists, reading in current round to do: ${cnt}"
    @ pcnt = ${cnt} - 1
    set pstep = ${prod_step}_production_${pcnt}
else
    echo "round file does NOT exist, round set to 1"
    set cnt = 1
    set pstep = ${equi_prefix}
endif

setenv OPENMM_CPU_THREADS 28

#Set up and run current round
set istep = ${prod_step}_production_${cnt}
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst"
python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out

#Do NOT resubmit if the output file for this round doesn't exist
if ( -f ${istep}.rst ) then
    echo "Output file exists: it is possible to proceed"
else
    echo "Output file does NOT exist, please check this run"
    echo "Terminating daisy-chain for now"
    exit
endif

#Increase the round
#Update the round file 
@ cnt += 1
touch ${file}
rm -f ${file}
echo ${cnt} > ${file}

#Resubmit self if round is not greater than $cntmax
if ( ${cnt} <= ${cntmax} ) then
    echo "Resubmitting to do round ${cnt}"
    qsub step5_production.bash
else
    echo "Max round reached, terminating daisy-chain for now"
endif

