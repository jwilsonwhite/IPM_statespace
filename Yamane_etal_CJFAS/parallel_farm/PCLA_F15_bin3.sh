#!/bin/bash -l

#SBATCH -J pcla_f15
#SBATCH -o pcla_f15-%j.output
#SBATCH -e pcla_f15-%j.output

module load matlab
hostname
matlab -nodisplay -nosplash -r PCLA_F15_bin3


