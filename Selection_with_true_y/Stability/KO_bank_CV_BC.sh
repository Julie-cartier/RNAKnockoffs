#!/bin/bash

#SBATCH --job-name=KO_bank_CV                                           # Job name
#SBATCH --output=/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_bank_CV_sets_BC/%j.log    # Output log
#SBATCH --error=/cluster/CBIO/data1/jcartier/Knockoffs/KO_Matrices_Bank/KO_bank_CV_sets_BC/%j.err     # Error log
#SBATCH --mem 20000                                                            # Job memory request
#SBATCH -p cbio-cpu                                                            # Name of the partition to use
#SBATCH --cpus-per-task=10                                                      # CPU cores per process (default 1, typically 4 or 5 - do not use more to let space for other people)

export PATH=~/anaconda3/bin:$PATH

source activate r-environment

Rscript KO_bank_CV_BC.R --job_id $SLURM_JOB_ID
