import numpy as np
import pandas as pd
import os
import os.path as osp
import argparse
import json

def read_command_line():
    parser = argparse.ArgumentParser(description='Run preprocessing on given subject and session for given task and run.')
    parser.add_argument('-s','--sub', type=str, help="subject number", required=True, dest='sub')
    parser.add_argument('-S','--ses', type=str, help="session number", required=False, default='', dest='ses')
    parser.add_argument('-d','--dir', type=str, help="project directory", required=True, dest='dir')
    parser.add_argument('-t','--task', type=str, help="task", required=True, dest='task')
    parser.add_argument('-fd','--fd_thresh', type=float, help="framewise displacement threshold", required=False, default=0.25, dest='fd_thresh')
    parser.add_argument('-lp','--low_pass', type=float, help="lower limit of bandpass", required=False, default=0.01, dest='low_pass')
    parser.add_argument('-hp','--high_pass', type=float, help="upper limit of bandpass", required=False, default=0.1, dest='high_pass')
    
    return parser.parse_args()

opts = read_command_line()

#set variables:
sub = opts.sub
ses = opts.ses
PRJ_DIR = opts.dir
task = opts.task
fd_thresh = opts.fd_thresh
low_pass = opts.low_pass
high_pass = opts.high_pass

if ses == '':
    #filename prefix
    prefix = (f'{sub}_task-{task}_')
    
    #Setup:
    OUTPUT_DIR= osp.join(PRJ_DIR,'derivatives','fmriprep', sub)
    PREPROC_DIR= osp.join(PRJ_DIR,'derivatives','preprocessed',sub)
else:
    #filename prefix
    prefix = (f'{sub}_{ses}_task-{task}_')
    
    #Setup:
    OUTPUT_DIR= osp.join(PRJ_DIR,'derivatives','fmriprep', sub, ses)
    PREPROC_DIR= osp.join(PRJ_DIR,'derivatives','preprocessed',sub, ses)

json_file = f"{sub}_{'ses-' + ses + '_' if ses else ''}task-{task}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.json"
with open(osp.join(OUTPUT_DIR, 'func', json_file), 'r') as file:
    data = json.load(file)
    TR = data['RepetitionTime']

print('TR: ', TR)
print('Bandpass: ', low_pass, '-', high_pass)

if TR < 0.1:
    raise ValueError('****ERROR with TR: TR=', TR)
elif TR > 5:
    raise ValueError('****ERROR with TR: TR=', TR)

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
    prop_outliers = np.sum(fd >= fd_thresh) / (len(fd)-np.sum(pd.isna(fd)))

    print('   meanFD: '+str(np.round(np.nanmean(fd),2))+', %FD>='+str(fd_thresh)+': '+str(np.round(prop_outliers*100, 2)))

    # Select variables to regress out column names of *_desc-confounds_timeseries.tsv
    variables=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
               'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1',
               'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1',
               'trans_x_derivative1_power2', 'trans_y_derivative1_power2',
               'trans_z_derivative1_power2', 'rot_x_derivative1_power2',
               'rot_y_derivative1_power2', 'rot_z_derivative1_power2',
               'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2',
               'rot_y_power2', 'rot_z_power2', 'global_signal', 'csf', 'white_matter']
                # 27 parameters =  [6 motion] * derivatives, power, & derivatives^2, GS, CSF, WM

    t = pd.read_csv(osp.join(OUTPUT_DIR, 'func', prefix+'desc-confounds_timeseries.tsv'),delimiter='\t')

        # ------ OUTLIERS ------ #
    # head motion outliers (FD, DVARS) & non-steady-state outliers
    # default fMRIprep outlier detection: exceeding 0.5 mm FD or 1.5 standardized DVARS

    outlier_columns = [c for c in t.columns if c.startswith('motion_outlier') or c.startswith('non_steady_state_outlier')]

    # fMRIprep-defined outliers 
    # time points to be included (1), to be excluded (0) - initialize as Nans
    outliers = np.zeros((len(t),))+np.nan
    if len(outlier_columns) > 0:
        outliers[np.where(np.sum(np.array(t[outlier_columns]),1)>0)[0]] = 0
        outliers[np.where(np.sum(np.array(t[outlier_columns]),1)==0)[0]] = 1

    # Here i want to do more stringent excusion than the fmriprep default: exclude frames where FD ≥ 0.25 mm
    fd_bad = np.where(t['framewise_displacement'] >= fd_thresh)[0]
    outliers[fd_bad] = 0  # mark additional volumes as outliers

    # remove the first 3 TRs to match trimmed functional data
    outliers = outliers[3:]
       
    # save outliers file
    outliers_file = osp.join(PREPROC_DIR, prefix+'outliers.1D')
    pd.DataFrame(outliers).to_csv(outliers_file, index=False, header=False, sep='\t')

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
    
    # remove the first 3 TRs to match trimmed functional data
    regressors = regressors.iloc[3:, :].reset_index(drop=True)
    
    # save regressors file
    reg_1D_file = osp.join(PREPROC_DIR, prefix+'confounds_regressors.1D')
    reg_csv_file = osp.join(PREPROC_DIR, prefix+'confounds_regressors.csv')
    regressors.to_csv(reg_1D_file, index=False, header=False)
    regressors.to_csv(reg_csv_file, index=False)


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
print('0. Remove first 3 TRs')
trimmed_img = osp.join(PREPROC_DIR, '0_' + prefix + 'trimmed.nii.gz')
if osp.exists(trimmed_img):
    print('Already removed first 3 TRs: ', trimmed_img)
else:
    cmd = f'3dTcat -prefix {trimmed_img} {OUTPUT_DIR}/func/{prefix}space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz[3..$]'
    os.system(cmd)
    if osp.exists(trimmed_img):
        print('First 3 TRs removed: ', trimmed_img)
    else:
        raise ValueError('****ERROR removing first 3 TRs')

# Update prep_img variable to point to the trimmed version
prep_img = trimmed_img

# -------------------------------------------------------------------------------------------- #
print('1. Mask fMRIprep output')
out_file = '1_'+prefix+'masked.nii.gz'
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
print('2. Despiking')
despike_out_file = osp.join(PREPROC_DIR, '2_' + prefix + 'despiked.nii.gz')
if not osp.exists(despike_out_file):
    cmd = f'3dDespike -prefix {despike_out_file} {PREPROC_DIR}/1_{prefix}masked.nii.gz'
    os.system(cmd)
    if osp.exists(despike_out_file):
        print('Despiking complete: ', despike_out_file)
    else:
        raise ValueError('****ERROR during despiking')
else:
    print('Despiking already complete: ', despike_out_file)

# -------------------------------------------------------------------------------------------- #
#can be added back in if this step is required. Note that you need to update the numbers of the files below to make sure this gets used 

# print('#. Intensity normalization')
# out_file = '2_'+prefix+'intensity_normalization.nii.gz'
# if osp.exists(osp.join(PREPROC_DIR,out_file)):
#     print('Intensity norm already complete: ', out_file)
# else:
#     run_norm = (f'fslmaths {PREPROC_DIR}/1_{prefix}fmriprep_output.nii.gz -inm 10000 {PREPROC_DIR}/{out_file}')
#     os.system(run_norm)

#     #check that code was successful
#     if osp.exists(osp.join(PREPROC_DIR, out_file)):
#         print('File saved: ', out_file)
#     else:
#         raise ValueError('****ERROR creating: ', out_file)
   

# -------------------------------------------------------------------------------------------- #
print('3. Create signal regressors for bandpass filtering')

out_file = (f'{prefix}bandpass_regressors.1D')
if osp.exists(osp.join(PREPROC_DIR,out_file)):
    print('Bandpass filtering complete: ', out_file)
else:
    run_tmp_bpass = (f'1dBport -input {PREPROC_DIR}/2_{prefix}despiked.nii.gz -TR {TR} -band {low_pass} {high_pass} -nozero > {PREPROC_DIR}/tmp_rm.bpass.1D')
    os.system(run_tmp_bpass)
    run_bpass = (f'1d_tool.py -infile {PREPROC_DIR}/tmp_rm.bpass.1D -write {PREPROC_DIR}/{prefix}bandpass_regressors.1D')
    os.system(run_bpass)
    delete_tmp = (f'rm {PREPROC_DIR}/tmp_rm.bpass.1D')
    os.system(delete_tmp)

    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, out_file)):
        print('File saved: ', out_file)
    else:
        raise ValueError('****ERROR creating: ', out_file)
   

# -------------------------------------------------------------------------------------------- #
print('4. Create regressors')
# creating *_xmat.1D matrix that combines selected frequency signals, linear trend, fMRIprep output confound matrix, and an intercept
out_file = (f'{prefix}xmat.1D')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regressors already found: ',out_file)
else:
    os.system('3dDeconvolve -input '+PREPROC_DIR+'/2_'+prefix+'despiked.nii.gz '+
              '-ortvec '+PREPROC_DIR+'/'+prefix+'bandpass_regressors.1D highpass '+
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
print('5. Regressing out nuisance variables')
# inputs
#    regressors: *_xmat.1D
#    censor frames: *_outliers.1D (1: included frames, 0: excluded frames)
#    (censored frames are treated as ZEROS: this should be converted to NaNs during analyses)
# outputs
#    3_*nuisance.nii.gz
out_file = (f'3_{prefix}nuisance_regressed.nii.gz')
if osp.exists(osp.join(PREPROC_DIR, out_file)):
    print('Regression already complete: ', out_file)
else:
    os.system('3dTproject -polort 0 -input '+PREPROC_DIR+'/2_'+prefix+'despiked.nii.gz '+
              '-ort '+PREPROC_DIR+'/'+prefix+'xmat.1D '+
              '-censor '+PREPROC_DIR+'/'+prefix+'outliers.1D '+
              '-cenmode KILL '+
              '-prefix '+PREPROC_DIR+'/3_'+prefix+'nuisance_regressed.nii.gz')
    

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
    cmd = '3dresample -master ' + PREPROC_DIR+'/3_'+prefix+'nuisance_regressed.nii.gz' + ' -prefix ' + PREPROC_DIR + '/' + out_file + ' -inset /project/mdrosenberg/IG/preprocessing/atlases/shen_2mm_268_parcellation.nii.gz' + ' -rmode NN'
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
    -inset '+ PREPROC_DIR + '/3_' + prefix + 'nuisance_regressed.nii.gz \
    -in_rois '+ PREPROC_DIR + '/shen_2mm_2mm_2_75mm_268_parcellation_resampled+tlrc \
    -fish_z \
    -ts_out \
    '
    os.system(cmd)


    #check that code was successful
    if osp.exists(osp.join(PREPROC_DIR, check_file)):
        print('FC calculation complete: ', out_file, '*')
    else:
        raise ValueError('****ERROR creating: ', out_file)


print('-------------------------------- ALL PREPROCESSING COMPLETED! --------------------------------')
















