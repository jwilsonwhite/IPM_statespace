#!/bin/bash -l

#SBATCH -J smys_am
#SBATCH -o smys_am-%j.output
#SBATCH -e smys_am-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r SMYS_Andrew_Molera_pre2007


