#!/bin/bash -l

#SBATCH -J smys_100bss
#SBATCH -o smys_100bss-%j.output
#SBATCH -e smys_100bss-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r SMYS_F02_100ss_12bJuly


