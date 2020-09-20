### This script executes all preprocessing steps.
### In order to run this script the following data has to be provided:
### the T1-weighted raw structural image (highres_head.nii) and
### the T2-weighted raw functional image (func_raw.nii.gz)
### The first part includes motion correction, spacial filtering and
### temporal filtering.
### The second part handles the WM and CSF regression.

# equation needed for later calculations
calc () {
    bc -l <<< "$@"
}

### PARAMETERS
# number of volumes
TR=2.
T=140
T_half=70
fwhm=6.
highpass_f=-1
lowpass_f=0.09


######## preprocessing steps (motion correction and filtering) #######

### directories and file names
# data directory
data_path=/home/bachmann/data/test/ #YourDataDirectory
# path of structural image (T1)
T1_name=${data_path}highres_head
# raw functional image
func_data_raw=${data_path}func_raw.nii.gz
# path of the truncated functional data file
func_data_trunc=${data_path}func_data_trunc.nii.gz
# functional brain data
func_brain=${data_path}func_brain
# initial functional data in float
prefiltered_func_data=${data_path}prefiltered_func_data_trunc
# example volume taking at T/2
example_func=${data_path}example_func
# functional data after motion correction
prefiltered_func_data_mcf=${data_path}prefiltered_func_data_mcf
# mean image of motion correction functional brain 
mean_func=${data_path}func_mean
# mask created from mean image
mask=${data_path}mask
# functional data after bet
prefiltered_func_data_bet=${data_path}prefiltered_func_data_bet
# thresholded functional image
prefiltered_func_data_thresh=${data_path}prefiltered_func_data_thresh
# functional image after spatial smoothing with susan
prefiltered_func_data_smooth=${data_path}prefiltered_func_data_smooth
# functional image after normalizing whole brain mode to 10000
prefiltered_func_data_intnorm=${data_path}prefiltered_func_data_intnorm
# functional image after temporal bandpass filter
prefiltered_func_data_tempfilt=${data_path}prefiltered_func_data_tempfilt
# temporary header for hd image
tmpHeader=${data_path}tmpHeader
# final preprocessed image
filtered_func_data=${data_path}filtered_func_data

### preprocessing steps (motion correction and filtering)
# removing first 10 images
fslroi $func_data_raw $func_data_trunc 0 63 0 63 0 32 10 149

# removing skull of functional image
bet $func_data_trunc $func_brain -F

fslmaths $func_brain $prefiltered_func_data -odt float
# example volume chosen
fslroi $prefiltered_func_data $example_func $T_half 1 
# motion correction
mcflirt -in $prefiltered_func_data -out $prefiltered_func_data_mcf -mats -plots  -rmsrel -rmsabs 

# creating mean image of motion corrected functional image
fslmaths $prefiltered_func_data_mcf -Tmean $mean_func
# creating mask of mean functional image
bet2 $mean_func $mask -f 0.3 -n -m
immv ${mask}_mask $mask
fslmaths $prefiltered_func_data_mcf -mas $mask $prefiltered_func_data_bet
k=$(fslstats $prefiltered_func_data_bet -p 2 -p 98)
n=${k:8}
m=$(calc $n*0.1) 
# creating thresholded mask
fslmaths $prefiltered_func_data_bet -thr $m -Tmin -bin mask -odt char
# estimating median on masked data
cmd=$(fslstats $prefiltered_func_data_mcf -k $mask -p 50)
median=${cmd} 
# creating mask functional image
fslmaths $mask -dilF $mask
fslmaths $prefiltered_func_data_mcf -mas $mask $prefiltered_func_data_thresh
# creating new mean functional image
fslmaths $prefiltered_func_data_thresh -Tmean $mean_func
# spatial smoothing
median_part=$(calc $median*0.75) #0.75*float(median) 
fwhm_real=$(calc $fwhm/2.355)
susan $prefiltered_func_data_thresh ${median_part} $fwhm_real 3 1 1 $mean_func ${median_part} $prefiltered_func_data_smooth
# applying mask to smoothed image
fslmaths $prefiltered_func_data_smooth -mas $mask $prefiltered_func_data_smooth
# normalization of smoothed image 
median_rat=$(calc 10000./$median)
fslmaths $prefiltered_func_data_smooth -mul $median_rat $prefiltered_func_data_intnorm 
# temporal filtering
hp_freq_pre=$(calc $highpass_f*2.*$TR)
lp_freq_pre=$(calc $lowpass_f*2.*$TR)
hp_freq=$(calc 1./$hp_freq_pre)
lp_freq=$(calc 1./$lp_freq_pre)
fslmaths $prefiltered_func_data_intnorm -bptf $hp_freq $lp_freq $prefiltered_func_data_tempfilt
# creating HD image
fslmaths $prefiltered_func_data_tempfilt $filtered_func_data
fslhd -x $filtered_func_data | sed 's/ dt = .*/ dt = '0.7'/g' > $tmpHeader
fslcreatehd $tmpHeader $filtered_func_data
# creating mean functional image
fslmaths $filtered_func_data -Tmean $mean_func


##########  White matter and CSF regression ############

### directories and file names
# results paths segmentation (gray and white matter and liquor)
dataDIR=/home/bachmann/data/test/
# path of structural image
head=${dataDIR}/highres_head.nii
# structural image includin GM WM and liquor
brain_high=${dataDIR}/highres.nii

# filtered functional image
func_image=${dataDIR}filtered_func_data.nii
# mean functional image
mean_func=${dataDIR}/mean_func.nii
# filtered and denoised functional image
func_image_denoised=${dataDIR}filtered_func_data_denoised.nii

# GM, WM and CSF in structural space derived from segmentation
graymat=${dataDIR}/1_0.0001_060.nii
whitemat=${dataDIR}/2_0.0001_060.nii
liqmat=${dataDIR}/3_0.0001_060.nii
# threshholding white mattter and CSF
whitemat_th=${dataDIR}/2_0.0001_060th.nii
liqmat_th=${dataDIR}/3_0.0001_060th.nii
# whole brain, CSF and WM mask in functional space
brain_mask=${dataDIR}/brain_mask.nii
WM_mask=${dataDIR}/WM_mask.nii
CSF_mask=${dataDIR}/CSF_mask.nii
# structural white matter and CSF in functional space
white_func2=${dataDIR}/wm_func2.nii
CSF_func2=${dataDIR}/CSF_func2.nii
# mean signal of extracted white matter and CSF signal and both
mean_WM=${dataDIR}/mean_WM.txt
mean_CSF=${dataDIR}/mean_CSF.txt
nuisance_regressors=${dataDIR}/nuisance_regressors.txt


### script

# segmentation of the structural image into GM, WM, liquor and rest
# input for segmentation made by SPM in Matlab having the following format
# 'general_datapath' 'SPM_path' 'Script_path' 'folder_name' 'bias regularisation' 'bias fwhm'
fslchfiletype NIFTI ${head}
$HOME/matlab-R2013b/bin/matlab -nosplash -r "segment_input '/home/bachmann/data/' '/home/bachmann/spm12/' '/home/bachmann/Schizo/scripts_in_schoen/segmentation' 'test' 0.0001 60"

# creating structural image of the brain (including gray, white matter and liquor structures) 
fslmaths ${graymat} -add ${whitemat} -add ${liqmat} ${brain_high}
# thresholding WM and CSF structural image
$FSLDIR/bin/fslmaths $whitemat -thr 0.99 $whitemat_th
$FSLDIR/bin/fslmaths $liqmat -thr 0.99 $liqmat_th
# creating mean functional image
fslmaths ${func_image} -Tmean ${mean_func}
## go to data directory
cd ${dataDIR}
ls
## register high resolution brain image to mean functional image (outcome are images with brain_low...Warped.nii.gz and brain_low_0GenericAffine.mat
/usr/lib/ants/antsRegistrationSyNQuick.sh -d 3 -t sr -m ${brain_high}.gz  -f ${mean_func}.gz -o brain_low
## checking results via sliced overlaying images
slices  ${mean_func}.gz brain_lowWarped.nii.gz -o brain_low_slices
## register WM and CSF structural image to functional space applying same transformation as before
/usr/lib/ants/antsApplyTransforms -d 3 -i ${whitemat_th}.gz -r ${mean_func}.gz  -o ${white_func2} -t brain_low0GenericAffine.mat
/usr/lib/ants/antsApplyTransforms -d 3 -i ${liqmat_th}.gz -r ${mean_func}.gz  -o ${CSF_func2} -t brain_low0GenericAffine.mat
## creating brain mask of functional brain
fslmaths ${mean_func} -thr 300 -bin $brain_mask
## multiply the thresholded segmentations with the global mask to create final WM, CSF masks
$FSLDIR/bin/fslmaths $white_func2  -thrP 99.99 -bin -mul $brain_mask  $WM_mask
$FSLDIR/bin/fslmaths $CSF_func2 -thrP 99.0 -bin -mul $brain_mask $CSF_mask
## get mean WM and CSF time-course
$FSLDIR/bin/fslmeants -i $func_image -o $mean_WM -m $WM_mask
$FSLDIR/bin/fslmeants -i $func_image -o $mean_CSF -m $CSF_mask
## combine the WM and CSF text files 
paste $mean_WM $mean_CSF > $nuisance_regressors
## run nuisance regression
$FSLDIR/bin/fsl_regfilt -i $func_image -o $func_image_denoised -d $nuisance_regressors -f "1,2"



