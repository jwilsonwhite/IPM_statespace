#!/bin/bash -l

#SBATCH -J smys_vb
#SBATCH -o smys_vb-%j.output
#SBATCH -e smys_vb-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r SMYS_Vandenberg_pre2007


