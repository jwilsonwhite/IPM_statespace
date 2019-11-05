#!/bin/bash -l

#SBATCH -J smys_100ss
#SBATCH -o smys_100ss-%j.output
#SBATCH -e smys_100ss-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r SMYS_F02_100ss_12July


