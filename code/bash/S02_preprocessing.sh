#!/bin/bash

#SBATCH --job-name=preproc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=pi-mdrosenberg
#SBATCH --time=2:00:00
#SBATCH --partition=caslake
##SBATCH --mem=200GB
#SBATCH --output=%j.out       # output log file
#SBATCH --error=%j.err        # error file
#SBATCH --mail-type=FAIL  # Email notification options: ALL, BEGIN, END, FAIL, ALL, NONE
#SBATCH --mail-user=igephart@rcc.uchicago.edu 

#set a random sleep delay so that the jobs do not submit at the same time
sleep $((RANDOM % 600))

# Load all required modules below. As an example, we load cuda/10.1
module load python
# module unload afni
module load afni
module load fsl
# Launch your run
#jobnum=$1
python /project/mdrosenberg/IG/ongoing_thought/code/python/S02_Preprocessing.py --sub=$sub --dir=$dir --task=$task --fd=$fd_thresh --low_pass=$low_pass --high_pass=$high_pass