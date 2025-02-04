import numpy as np
import pandas as pd
import os
import os.path as osp
import argparse

#import sleep_networks.proj_info as info

def read_command_line():
    parser = argparse.ArgumentParser(description='Run preprocessing on given subject and session.')
    parser.add_argument('-s','--sub', type=str, help="subject number", required=True, dest='sub')
    parser.add_argument('-S','--ses', type=str, help="session number", required=False, dest='ses')
    
    return parser.parse_args()

opts = read_command_line()

#Subjects:
# sub = 'sub-01'
# ses = 'ses-01'
sub = opts.sub
ses = opts.ses

print('test: ', sub)
print('test: ', ses)

print(make error)

#Setup:
#TODO: eventaully should be loaded in from info
#PRJ_DIR = '/project/mdrosenberg/IG/sleep_networks'
# OUTPUT_DIR= osp.join(PRJ_DIR, 'fmriprep')
# PREPROC_DIR= osp.join(PRJ_DIR, 'preprocessed')
PRJ_DIR = '/project/mdrosenberg/IG/Flanker_onesub'
OUTPUT_DIR= osp.join(PRJ_DIR, 'derivatives',sub)
PREPROC_DIR= osp.join(PRJ_DIR, 'preprocessed',sub)
fd_thresh=0.5 #TODO: check w prepreg
TR = 2

#set up directories
if not osp.exists(PREPROC_DIR):
    os.makedirs(PREPROC_DIR)

#check if regressor file has already been created, if not make summary files
if not osp.exists(osp.join(PREPROC_DIR,sub+'_task-flanker_run-1_confounds_regressors.1D')):
    print('Creating confounds file for: ',sub,ses)
    conf_file = osp.join(OUTPUT_DIR, sub, 'func',sub+'_task-flanker_run-1_desc-confounds_timeseries.tsv')
    cf = pd.read_csv(conf_file,delimiter='\t')

    # time course of framewise displacement (FD)
    fd = cf['framewise_displacement']

    # proportion of frames where FD >= fd_thres
    outliers = np.sum(fd >= fd_thresh) / (len(fd)-np.sum(pd.isna(fd)))

    print('   meanFD: '+str(np.round(np.nanmean(fd),2))+', %FD>='+str(fd_thresh)+': '+str(np.round(outliers*100, 2)))

    #TODO: figure out if this is missing parameters? Ask Monica
    # Select variables to regress out column names of *_desc-confounds_timeseries.tsv
    variables=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
               'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1',
               'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1',
               'trans_x_derivative1_power2', 'trans_y_derivative1_power2',
               'trans_z_derivative1_power2', 'rot_x_derivative1_power2',
               'rot_y_derivative1_power2', 'rot_z_derivative1_power2',
               'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2',
               'rot_y_power2', 'rot_z_power2', 'global_signal', 'csf', 'white_matter']
                # 24 motion parameters, global signal, CSF, white matter

    #TODO: fix this filepath to osp
    t = pd.read_csv(OUTPUT_DIR+'/'+sub+'/func/'+sub+'_task-flanker_run-1_desc-confounds_timeseries.tsv',
                        delimiter='\t')

        # ------ OUTLIERS ------ #
    # head motion outliers (FD, DVARS) & non-steady-state outliers
    # default fMRIprep outlier detection: exceeding 0.5 mm FD or 1.5 standardized DVARS
    outlier_columns = [i for i in t.columns if i.startswith('motion_outlier') or i.startswith('non_steady_state_outlier')]

       
    # time points to be included (1), to be excluded (0)
    outliers = np.zeros((len(t),))+np.nan
    outliers[np.where(np.sum(np.array(t[outlier_columns]),1)>0)[0]] = 0
    outliers[np.where(np.sum(np.array(t[outlier_columns]),1)==0)[0]] = 1

       
    # save outliers file
    outliers_file = osp.join(PREPROC_DIR, sub+'_task-flanker_run-1_outliers.1D')
    pd.DataFrame(outliers).to_csv(outliers_file)
    print('File saved: ', outliers_file)
       
    # ------ REGRESSORS ------ #
    # duplicate values in second row if the first row is NaN
    for varb in variables:
        if pd.isna(t[varb][0]):
            t[varb][0] = t[varb][1]
    regressors = t[variables]

    # save regressors file
    reg_1D_file = osp.join(PREPROC_DIR, sub+'_task-flanker_run-1_confounds_regressors.1D')
    reg_csv_file = osp.join(PREPROC_DIR, sub+'_task-flanker_run-1_confounds_regressors.csv')
    regressors.to_csv(reg_1D_file, index=False, header=False)
    regressors.to_csv(reg_csv_file)
    print('File saved: ', reg_1D_file)
    print('File saved: ', reg_csv_file)
else:
    print('Outlier and regressors files found.')

# -------------------------------------------------------------------------------------------- #
print('STARTING PREPROCESSING')
print('Running preprocessing on ', sub)

# From fmriprep: preprocessed, MNI-registered EPI
prep_img = osp.join(OUTPUT_DIR,sub,'func',sub+'_task-flanker_run-1_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz')
print('Image: ', prep_img)
# brain mask
prep_mask_img = osp.join(OUTPUT_DIR,sub,'func',sub+'_task-flanker_run-1_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz')
print('Mask: ', prep_mask_img)

# ------------ #
print('1. Mask fMRIprep output')
out_file = '1_'+sub+'_task-flanker_fmriprep_output.nii.gz'
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('Masking already complete: ', out_file)
else:
    run_masking = (f'3dcalc -a {prep_img} -b {prep_mask_img} -expr "(a*b)" -prefix {PREPROC_DIR}/{out_file}')
    os.system(run_masking)
    print('File written: ', out_file)
 
 # ------------ #
print('2. Intensity normalization')
out_file = '2_'+sub+'_task-flanker_intensity_normalization.nii.gz'
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('Intensity norm already complete: ', out_file)
else:
    run_norm = (f'fslmaths {PREPROC_DIR}/1_{sub}_task-flanker_fmriprep_output.nii.gz -inm 10000 {PREPROC_DIR}/{out_file}')
    os.system(run_norm)
    print('File written: ', out_file)


 # ------------ #
print('3. Create low-frequency signal regressors for high-pass filtering')
# 128 s cutoff
out_file = (f'{sub}_task-flanker_highpass_regressors.1D')
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('High-pass filtering complete: ', out_file)
else:
    run_tmp_hpass = (f'1dBport -input {PREPROC_DIR}/2_{sub}_task-flanker_intensity_normalization.nii.gz -TR {TR} -band 0 0.0078 -nozero > {PREPROC_DIR}/tmp_rm.hpass.1D')
    os.system(run_tmp_hpass)
    run_hpass = (f'1d_tool.py -infile {PREPROC_DIR}/tmp_rm.hpass.1D -write {PREPROC_DIR}/{sub}_task-flanker_highpass_regressors.1D')
    os.system(run_hpass)
    delete_tmp = (f'rm {PREPROC_DIR}/tmp_rm.hpass.1D')
    os.system(delete_tmp)
    print('File written: ', out_file)

# ------------ #
print('4. Create regressors')
# creating *_xmat.1D matrix that combines low-frequency signals, linear trend, fMRIprep output confound matrix, and an intercept
out_file = (f'{sub}_task-flanker_xmat.1D')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regressors already found: ',out_file)
else:
    os.system('3dDeconvolve -input '+PREPROC_DIR+'/2_'+sub+'_task-flanker_intensity_normalization.nii.gz '+
              '-ortvec '+PREPROC_DIR+'/'+sub+'_task-flanker_highpass_regressors.1D highpass '+
              '-ortvec '+PREPROC_DIR+'/'+sub+'_task-flanker_run-1_confounds_regressors.1D confounds '+
              '-polort 1 '+
              '-fout -tout -x1D '+PREPROC_DIR+'/'+sub+'_task-flanker_xmat.1D '+
              '-fitts fitts -errts errts -x1D_stop -bucket stats')
    print('Regressor files written: ', out_file)

# ------------ #
print('5. Regressing out nuissance variables')
# inputs
#    regressors: *_xmat.1D
#    censor frames: *_outliers.1D (1: included frames, 0: excluded frames)
#    (censored frames are treated as ZEROS: this should be converted to NaNs during analyses)
# outputs
#    3_*_nuissance_regressed.nii.gz
out_file = (f'3_{sub}_task-flanker_nuissance_regressed.nii.gz')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regression already complete: ', out_file)
else:
    os.system('3dTproject -polort 0 -input '+PREPROC_DIR+'/2_'+sub+'_task-flanker_intensity_normalization.nii.gz '+
              '-ort '+PREPROC_DIR+'/'+sub+'_task-flanker_xmat.1D '+
              '-censor '+PREPROC_DIR+'/'+sub+'_task-flanker_run-1_outliers.1D '+
              '-cenmode ZERO '+
              '-prefix '+PREPROC_DIR+'/3_'+sub+'_task-flanker_nuissance_regressed.nii.gz')
    print('Regression complete: ', out_file)



















