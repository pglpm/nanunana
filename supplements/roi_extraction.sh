#!/bin/bash

# This script registers the functional images first to structural
# and then to standard space. In standard space ROIs are extracted
# according to the Oxford atlas.

#### paths and file names
# data path
dataDIR=/home/bachmann/data/test/
# script path
resultDIR=/home/bachmann/Schizo/scripts_in_schoen/roi_extraction
# functional image path
func_image=${dataDIR}filtered_func_data_denoised.nii
# mean functional image
mean_func=${dataDIR}/mean_func.nii
# thresholded functional image
mean_func_th=${dataDIR}/mean_func_th.nii
# functional image in structural space
func_in_struc=${dataDIR}/func_in_struc
# transformation matrix which transforms from functional to structural space
func_in_struc_mat=${dataDIR}/func_in_struc.mat
# functional image in standard space
func_in_stand=${dataDIR}/func_in_stand
# structural image including GM, WM and liquor
brain_high=${dataDIR}/highres.nii.gz
# path of the MNI152 standard image
gmpath=/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz

# going to data directory
cd ${dataDIR}
# separating all functional images
fslsplit $func_image -t
# thresholding meam image
fslmaths $mean_func -thr 350 $mean_func_th

# register (thresholded) mean functional image to structural space
flirt -in $mean_func_th -ref ${brain_high} -out func_in_teno -omat $func_in_struc_mat -bins 256 -cost corratio -searchrx -140 140 -searchry -140 140 -searchrz -140 140 -dof 12  -interp trilinear
flirt -in $mean_func -ref ${brain_high} -applyxfm -init ${func_in_struc_mat} -out $func_in_struc

# registration from structural space to standard space
/usr/lib/ants/antsRegistrationSyNQuick.sh -d 3 -t s -m ${brain_high}  -f ${gmpath} -o ${func_in_stand}
/usr/lib/ants/antsApplyTransforms -d 3 -i $func_in_struc.nii.gz -r ${gmpath}  -o ${func_in_stand2} -t ${func_in_stand}1Warp.nii.gz -t ${func_in_stand}0GenericAffine.mat


# applying transformation from functional to structural space to all volumes
for i in $(seq 0 1 9); do  flirt -in $dataDIR/vol000${i}.nii.gz -ref ${brain_high} -applyxfm -init ${func_in_struc_mat} -out $dataDIR/vol000${i}_struc.nii.gz ; done
for i in $(seq 10 1 99); do  flirt -in $dataDIR/vol00${i}.nii.gz -ref ${brain_high} -applyxfm -init ${func_in_struc_mat} -out $dataDIR/vol00${i}_struc.nii.gz ; done
for i in $(seq 100 1 139); do  flirt -in $dataDIR/vol0${i}.nii.gz -ref ${brain_high} -applyxfm -init ${func_in_struc_mat} -out $dataDIR/vol0${i}_struc.nii.gz ; done

# applying transformation from structural to standard space to all volumes
for i in $(seq 0 1 9); do /usr/lib/ants/antsApplyTransforms -d 3 -i $dataDIR/vol000${i}_struc.nii.gz -r ${gmpath}  -o  $dataDIR/vol000${i}_stand.nii.gz -t ${func_in_stand}1Warp.nii.gz -t ${func_in_stand}0GenericAffine.mat ; done
for i in $(seq 10 1 99); do /usr/lib/ants/antsApplyTransforms -d 3 -i $dataDIR/vol00${i}_struc.nii.gz -r ${gmpath}  -o  $dataDIR/vol00${i}_stand.nii.gz -t ${func_in_stand}1Warp.nii.gz -t ${func_in_stand}0GenericAffine.mat ; done
for i in $(seq 100 1 139); do /usr/lib/ants/antsApplyTransforms -d 3 -i $dataDIR/vol0${i}_struc.nii.gz -r ${gmpath}  -o  $dataDIR/vol0${i}_stand.nii.gz -t ${func_in_stand}1Warp.nii.gz -t ${func_in_stand}0GenericAffine.mat ; done



cd $resultDIR
echo $resultDIR

# running python script that extract the mean signal of the areas
python extract_areas.py ${dataDIR}


# deleting all volumes

for i in $(seq 0 1 9); do rm $dataDIR/vol000${i}.nii.gz; done
for i in $(seq 0 1 9); do rm $dataDIR/vol000${i}_struc.nii.gz; done
for i in $(seq 0 1 9); do rm $dataDIR/vol000${i}_stand.nii.gz; done

for i in $(seq 10 1 99); do rm $dataDIR/vol00${i}.nii.gz; done
for i in $(seq 10 1 99); do rm $dataDIR/vol00${i}_struc.nii.gz; done
for i in $(seq 10 1 99); do rm $dataDIR//vol00${i}_stand.nii.gz; done

for i in $(seq 100 1 139); do rm $dataDIR/vol0${i}.nii.gz; done
for i in $(seq 100 1 139); do rm $dataDIR/vol0${i}_struc.nii.gz; done
for i in $(seq 100 1 139); do rm $dataDIR/vol0${i}_stand.nii.gz; done

