#!/bin/bash
#SBATCH --job-name=fmriprep
#SBATCH --account=pi-mdrosenberg
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=igephart@rcc.uchicago.edu

# ---------------------------
# Load modules
# ---------------------------
module load singularity

# ---------------------------
# Define variables
# ---------------------------
FMRIPREP_VER=latest
subject=${subj}
PROJ_DIR=${dir}
DATA_DIR=${data}

export FS_LICENSE=/project2/mdrosenberg/ZZ/paranoia/license.txt
export SINGULARITYENV_FS_LICENSE=$FS_LICENSE

BIDS_DIR=${DATA_DIR}
OUTPUT_DIR="${PROJ_DIR}/derivatives/fmriprep"
WORK_DIR="${PROJ_DIR}/derivatives/work"
FS_DIR="${OUTPUT_DIR}/freesurfer"

# Create necessary directories if missing
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${FS_DIR}"

# ---------------------------
# Log job info
# ---------------------------
{
  echo "===================================="
  echo "JOB INFO:"
  echo "SLURM Job ID: ${SLURM_JOB_ID}"
  echo "Subject: ${subject}"
  echo "Project Directory: ${PROJ_DIR}"
  echo "BIDS Directory: ${BIDS_DIR}"
  echo "Output Directory: ${OUTPUT_DIR}"
  echo "Work Directory: ${WORK_DIR}"
  echo "fMRIPrep Version: ${FMRIPREP_VER}"
  echo "Framewise Displacement Cutoff (FD_THRESH): ${FD_THRESH} mm"
  echo "===================================="
} >> "${SLURM_JOB_ID}.out"

singularity run \
  --cleanenv \
  --bind /proc:/proc \
  --bind ${DATA_DIR}:${DATA_DIR} \
  --bind ${OUTPUT_DIR}:${OUTPUT_DIR} \
  --bind ${WORK_DIR}:${WORK_DIR} \
  --bind ${FS_LICENSE}:${FS_LICENSE} \
  /home/igephart/fmriprep.simg \
  ${BIDS_DIR} ${OUTPUT_DIR} participant \
  --participant-label ${subject} \
  --fs-subjects-dir ${OUTPUT_DIR}/freesurfer \
  --fs-license-file ${FS_LICENSE} \
  -w ${WORK_DIR} \
  --output-spaces anat MNI152NLin2009cAsym:res-2 \
  --skip_bids_validation \
  --n_cpus 8 \
  --omp-nthreads 8 \
  --ignore fieldmaps 

# ---------------------------
# Completion message
# ---------------------------
if [ $? -eq 0 ]; then
  echo "fMRIPrep completed successfully for subject ${subject}" >> "${SLURM_JOB_ID}.out"
else
  echo "ERROR: fMRIPrep failed for subject ${subject}" >> "${SLURM_JOB_ID}.err"
fi
