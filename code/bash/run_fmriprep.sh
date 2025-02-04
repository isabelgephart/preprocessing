#!/bin/bash
#SBATCH --job-name=fmriprep
#SBATCH --account=pi-mdrosenberg
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=20
#SBATCH --mem=20GB
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=ALL  # Email notification options: ALL, BEGIN, END, FAIL, ALL, NONE
#SBATCH --mail-user=igephart@rcc.uchicago.edu 

module load singularity

FMRIPREP_VER=latest
subject=$subj
session=$ses
PROJ_DIR=$dir
DATA_DIR=$data

export SINGULARITYENV_FS_LICENSE=/project2/mdrosenberg/ZZ/paranoia/license.txt
export FS_LICENSE=/project2/mdrosenberg/ZZ/paranoia/license.txt

BIDS_DIR=${DATA_DIR}
OUTPUT_DIR="${PROJ_DIR}/fmriprep/" 
WORK_DIR="${PROJ_DIR}/work/" 
FS_DIR="${OUTPUT_DIR}/freesurfer"

#write job info to top of error and output files
echo "JOB INFO: " >> $SLURM_JOB_ID.out
echo "Subject: $subject" >> $SLURM_JOB_ID.out
echo "Session: $session" >> $SLURM_JOB_ID.out
echo "PROJ_DIR: $PROJ_DIR" >> $SLURM_JOB_ID.out
echo "BIDS_DIR: $BIDS_DIR" >> $SLURM_JOB_ID.out
echo "OUTPUT_DIR: $OUTPUT_DIR" >> $SLURM_JOB_ID.out
echo "----------------------" >> $SLURM_JOB_ID.out

echo "JOB INFO: " >> $SLURM_JOB_ID.err
echo "Subject: $subject" >> $SLURM_JOB_ID.err
echo "Session: $session" >> $SLURM_JOB_ID.err
echo "----------------------" >> $SLURM_JOB_ID.err

singularity run \
--cleanenv \
-B $DATA_DIR \
-B $FS_LICENSE:/project2/mdrosenberg/ZZ/paranoia/license.txt \
-B $OUTPUT_DIR \
-B $WORK_DIR \
-e /home/igephart/fmriprep.simg \
${BIDS_DIR} ${OUTPUT_DIR} participant \
-w ${WORK_DIR} \
--participant-label ${subject} \
--output-spaces anat fsnative MNI152NLin2009cAsym:res-2 fsaverage:den-10k fsLR \
--skip_bids_validation \
--n_cpus 3 \
--omp-nthreads 20 \
--ignore {fieldmaps,slicetiming}
--fs-subjects-dir ${FS_DIR}
