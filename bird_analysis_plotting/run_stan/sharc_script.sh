#!/bin/bash
# Cores
#$ -pe smp 15
# Memory (per core)
#$ -l rmem=3G
# Output
#$ -o runInfo/script_oe_$JOB_ID.txt
# Time
#$  -l h_rt=48:00:00 
#$ -j y
# Email 
#$ -m bea
#$ -M jacob.socolar@gmail.com

# chains and threads
nchains=1
nthreads=15

# Specify threads again
export OMP_NUM_THREADS=15

# run
module load apps/R/3.6.3/gcc-8.2.0
Rscript code/run_breakpoint_obsvr_model_cmdstanr.R $nchains $nthreads
