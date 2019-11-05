#!/bin/bash -l

#SBATCH -J smys_np
#SBATCH -o smys_np-%j.output
#SBATCH -e smys_np-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r SMYS_Naples_pre2012


