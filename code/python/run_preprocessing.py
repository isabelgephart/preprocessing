import numpy as np
import pandas as pd
import os
import os.path as osp
import argparse

def read_command_line():
    parser = argparse.ArgumentParser(description='Run preprocessing on given subject and session for given task and run.')
    parser.add_argument('-s','--sub', type=str, help="subject number", required=True, dest='sub')
    parser.add_argument('-S','--ses', type=str, help="session number", required=False, dest='ses')
    parser.add_argument('-d','--dir', type=str, help="project directory", required=True, dest='dir')
    parser.add_argument('-t','--task', type=str, help="task", required=True, dest='task')
    parser.add_argument('-r','--run', type=str, help="run", required=False, dest='run')
    
    return parser.parse_args()

opts = read_command_line()

#set variables:
sub = opts.sub
ses = opts.ses
PRJ_DIR = opts.dir
task = opts.task
run = opts.run

#filename prefix
if ses == '' and run =='':
    prefix = (f'{sub}_task-{task}_')
elif run == '':
    prefix = (f'{sub}_{ses}_task-{task}_')
elif ses =='':
    prefix = (f'{sub}_task-{task}_run-{run}_')
else:
    prefix = (f'{sub}_{ses}_task-{task}_run-{run}_')

#Setup:
OUTPUT_DIR= osp.join(PRJ_DIR,'derivatives','fmriprep', sub, ses)
PREPROC_DIR= osp.join(PRJ_DIR,'derivatives','preprocessed',sub, ses)

fd_thresh=0.5 #TODO: check w prepreg
TR = 2


#set up directories
if not osp.exists(PREPROC_DIR):
    os.makedirs(PREPROC_DIR)

# -------------------------------------------------------------------------------------------- #
#check if regressor file has already been created, if not make summary files
if not osp.exists(osp.join(PREPROC_DIR,prefix+'confounds_regressors.1D')):
    print('Creating confounds file for: ',sub,ses)
    conf_file = osp.join(OUTPUT_DIR, 'func',prefix+'desc-confounds_timeseries.tsv')
    cf = pd.read_csv(conf_file,delimiter='\t')

    # time course of framewise displacement (FD)
    fd = cf['framewise_displacement']

    # proportion of frames where FD >= fd_thres
    outliers = np.sum(fd >= fd_thresh) / (len(fd)-np.sum(pd.isna(fd)))

    print('   meanFD: '+str(np.round(np.nanmean(fd),2))+', %FD>='+str(fd_thresh)+': '+str(np.round(outliers*100, 2)))

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

    t = pd.read_csv(osp.join(OUTPUT_DIR, 'func', prefix+'desc-confounds_timeseries.tsv'),delimiter='\t')

        # ------ OUTLIERS ------ #
    # head motion outliers (FD, DVARS) & non-steady-state outliers
    # default fMRIprep outlier detection: exceeding 0.5 mm FD or 1.5 standardized DVARS
    outlier_columns = [i for i in t.columns if i.startswith('motion_outlier') or i.startswith('non_steady_state_outlier')]

       
    # time points to be included (1), to be excluded (0)
    outliers = np.zeros((len(t),))+np.nan
    outliers[np.where(np.sum(np.array(t[outlier_columns]),1)>0)[0]] = 0
    outliers[np.where(np.sum(np.array(t[outlier_columns]),1)==0)[0]] = 1

       
    # save outliers file
    outliers_file = osp.join(PREPROC_DIR, prefix+'outliers.1D')
    pd.DataFrame(outliers).to_csv(outliers_file)

    #check that code was successful
    if osp.exists(outliers_file):
        print('File saved: ', outliers_file)
    else:
        raise ValueError('****ERROR creating: ', outliers_file)
       
    # ------ REGRESSORS ------ #
    # duplicate values in second row if the first row is NaN
    for varb in variables:
        if pd.isna(t[varb][0]):
            t[varb][0] = t[varb][1]
    regressors = t[variables]

    # save regressors file
    reg_1D_file = osp.join(PREPROC_DIR, prefix+'confounds_regressors.1D')
    reg_csv_file = osp.join(PREPROC_DIR, prefix+'confounds_regressors.csv')
    regressors.to_csv(reg_1D_file, index=False, header=False)
    regressors.to_csv(reg_csv_file)


    #check that code was successful
    if osp.exists(reg_1D_file):
        print('File saved: ', reg_1D_file)
    else:
        raise ValueError('****ERROR creating: ', reg_1D_file)


    if osp.exists(reg_csv_file):
        print('File saved: ', reg_csv_file)
    else:
        raise ValueError('****ERROR creating: ', reg_csv_file)

else:
    print('Outlier and regressors files found.')

# -------------------------------------------------------------------------------------------- #
print('STARTING PREPROCESSING')
print('Running preprocessing on ', sub, ses)

# From fmriprep: preprocessed, MNI-registered EPI
prep_img = osp.join(OUTPUT_DIR,'func',prefix+'space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz')
print('Image: ', prep_img)
# brain mask
prep_mask_img = osp.join(OUTPUT_DIR,'func',prefix+'space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz')
print('Mask: ', prep_mask_img)

# -------------------------------------------------------------------------------------------- #
print('1. Mask fMRIprep output')
out_file = '1_'+prefix+'fmriprep_output.nii.gz'
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('Masking already complete: ', out_file)
else:
    run_masking = (f'3dcalc -a {prep_img} -b {prep_mask_img} -expr "(a*b)" -prefix {PREPROC_DIR}/{out_file}')
    os.system(run_masking)

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('File saved: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)
   
 
# -------------------------------------------------------------------------------------------- #
print('2. Intensity normalization')
out_file = '2_'+prefix+'intensity_normalization.nii.gz'
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('Intensity norm already complete: ', out_file)
else:
    run_norm = (f'fslmaths {PREPROC_DIR}/1_{prefix}fmriprep_output.nii.gz -inm 10000 {PREPROC_DIR}/{out_file}')
    os.system(run_norm)

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('File saved: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)
   


# -------------------------------------------------------------------------------------------- #
print('3. Create low-frequency signal regressors for high-pass filtering')
# 128 s cutoff
out_file = (f'{prefix}highpass_regressors.1D')
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('High-pass filtering complete: ', out_file)
else:
    run_tmp_hpass = (f'1dBport -input {PREPROC_DIR}/2_{prefix}intensity_normalization.nii.gz -TR {TR} -band 0 0.0078 -nozero > {PREPROC_DIR}/tmp_rm.hpass.1D')
    os.system(run_tmp_hpass)
    run_hpass = (f'1d_tool.py -infile {PREPROC_DIR}/tmp_rm.hpass.1D -write {PREPROC_DIR}/{prefix}highpass_regressors.1D')
    os.system(run_hpass)
    delete_tmp = (f'rm {PREPROC_DIR}/tmp_rm.hpass.1D')
    os.system(delete_tmp)

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('File saved: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)
   

# -------------------------------------------------------------------------------------------- #
print('4. Create regressors')
# creating *_xmat.1D matrix that combines low-frequency signals, linear trend, fMRIprep output confound matrix, and an intercept
out_file = (f'{prefix}xmat.1D')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regressors already found: ',out_file)
else:
    os.system('3dDeconvolve -input '+PREPROC_DIR+'/2_'+prefix+'intensity_normalization.nii.gz '+
              '-ortvec '+PREPROC_DIR+'/'+prefix+'highpass_regressors.1D highpass '+
              '-ortvec '+PREPROC_DIR+'/'+prefix+'confounds_regressors.1D confounds '+
              '-polort 1 '+
              '-fout -tout -x1D '+PREPROC_DIR+'/'+prefix+'xmat.1D '+
              '-fitts fitts -errts errts -x1D_stop -bucket stats')

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('Regressor files written: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)
   

# -------------------------------------------------------------------------------------------- #
print('5. Regressing out nuissance variables')
# inputs
#    regressors: *_xmat.1D
#    censor frames: *_outliers.1D (1: included frames, 0: excluded frames)
#    (censored frames are treated as ZEROS: this should be converted to NaNs during analyses)
# outputs
#    3_*_nuissance_regressed.nii.gz
out_file = (f'3_{prefix}nuissance_regressed.nii.gz')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regression already complete: ', out_file)
else:
    os.system('3dTproject -polort 0 -input '+PREPROC_DIR+'/2_'+prefix+'intensity_normalization.nii.gz '+
              '-ort '+PREPROC_DIR+'/'+prefix+'xmat.1D '+
              '-censor '+PREPROC_DIR+'/'+prefix+'outliers.1D '+
              '-cenmode ZERO '+
              '-prefix '+PREPROC_DIR+'/3_'+prefix+'nuissance_regressed.nii.gz')
    

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('Regression complete: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)


# -------------------------------------------------------------------------------------------- #
print('6. Resample Shen atlas')
out_file = 'shen_2mm_2mm_2_75mm_268_parcellation_resampled+tlrc'
check_file = 'shen_2mm_2mm_2_75mm_268_parcellation_resampled+tlrc.BRIK' #this step writes multiple files, this checks for one of them
if not osp.exists(osp.join(PREPROC_DIR, check_file)):
    cmd = '3dresample -master ' + PREPROC_DIR+'/3_'+prefix+'nuissance_regressed.nii.gz' + ' -prefix ' + PREPROC_DIR + '/' + out_file + ' -inset /project/mdrosenberg/IG/preprocessing/atlases/shen_2mm_268_parcellation.nii.gz' + ' -rmode NN'
    os.system(cmd)
    print(cmd)

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, check_file)):
        print('Resampling complete: ', out_file, '*')
    else:
        raise ValueError('****ERROR creating: ', out_file)
else:
    print('Resampling already complete: ', out_file, '*')


# -------------------------------------------------------------------------------------------- #
print('7. Parcellate the brain and generate FC')
out_file = (f'4_{prefix}LPI')
check_file = '4_' + prefix + 'LPI_000.netcc'
if osp.exists(osp.join(PREPROC_DIR, check_file)):
    print('FC calculation already complete: ', out_file, '*')
else:
    cmd = '3dNetCorr -prefix '+ PREPROC_DIR + '/' + out_file +'\
    -mask ' + prep_mask_img + ' \
    -inset '+ PREPROC_DIR + '/3_' + prefix + 'nuissance_regressed.nii.gz \
    -in_rois '+ PREPROC_DIR + '/shen_2mm_2mm_2_75mm_268_parcellation_resampled+tlrc \
    -fish_z \
    -ts_out \
    -push_thru_many_zeros'
    os.system(cmd)


    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, check_file)):
        print('FC calculation complete: ', out_file, '*')
    else:
        raise ValueError('****ERROR creating: ', out_file)


print('-------------------------------- ALL PREPROCESSING COMPLETED! --------------------------------')
















