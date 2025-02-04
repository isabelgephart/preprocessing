#!/bin/bash

#SBATCH --job-name=preprocessing 
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=20
#SBATCH --account=pi-mdrosenberg
#SBATCH --partition=caslake
##SBATCH --mem=200GB
#SBATCH --output=%j.out       # output log file
#SBATCH --error=%j.err        # error file
#SBATCH --mail-type=ALL  # Email notification options: ALL, BEGIN, END, FAIL, ALL, NONE
#SBATCH --mail-user=igephart@rcc.uchicago.edu 


# Load all required modules below. As an example, we load cuda/10.1
module load python
# module unload afni
module load afni
module load fsl
# Launch your run
#jobnum=$1
python /project/mdrosenberg/IG/sleep_networks/code/python/run_preprocessing.py --sub=$sub --ses=$ses --dir=$dir