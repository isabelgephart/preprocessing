# preprocessing

This repo contains the scripts I generally use to preprocess fMRI data. I do this in a few steps:
1. Run fmriprep
2. Run additional preprocessing steps (ultimately generating functional connectivity matrices) in afni and fsl
3. Run some quick sanity checks on the FC matrices that might show something is up with the data (but you should still do visual QC!)

This code is all writen to be run on midway, but may be useful to do on a local computer. Everything is organized in steps (S0#) and the general workflow for each step is: 
- A jupyter notebook generates an SBATCH script to submit jobs for each subject/scan 
- A .sh script submits the jobs 
- In some cases a python script is called within the jobs

My general recommendation is that you copy all of /code into your own folder and it should all work as long as you set the [roject name correctly at the top of each notebook. Note that you will need to have the atlases accessible as well. 

