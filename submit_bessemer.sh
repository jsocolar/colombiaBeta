#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=3G
#SBATCH --mail-user=s.c.mills@sheffield.ac.uk
#SBATCH --output=output_model_run1.txt

# run (currently in staging- update when it's been moved)
module use /usr/local/modulefiles/staging/eb/all
module load R/4.0.0-foss-2020a

Rscript code/run_full_model.R