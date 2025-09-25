#!/bin/bash
set -euo pipefail

# fMRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.


# Paths – adjust as needed
denoised_dir="/home/ubuntu/volume/Antinomics/subjects_fsl_dir/denoised"
fmri_dir="/home/ubuntu/volume/Antinomics/raws/fMRI"
sessions=("s1" "s2")     # list all sessions here


# Function to preprocess a single subject
preprocess_subject() {
    local subject=$1
    cd $fmri_dir/$subject
    echo "[0] Reorient to MNI"
    fslreorient2std "$denoised_dir/${subject}_denoised.nii" "t1_reoriented"
    fslreorient2std "fmri_dir/${subject}_s1.nii" "t2_s1_reoriented" # fix this one
    fslreorient2std "${subject}_s2.nii" "t2_s2_reoriented"

    echo "[1] FAST segmentation"
    fast -B -o . "t1_reoriented"

    echo "[2] Brain extraction of T1"
    bet t1_reoriented "t1_brain"

    echo "[3] Creating mean functional for registration"
    fslmaths "t2_s1_reoriented" -Tmean "meanfunc"

    echo "[4] Brain extraction of functional"
    bet "meanfunc" "meanfunc_brain" -F -m

    echo "[4] Registration T1 -> mean functional linear 6 and 12 dof"
    flirt -in "t1_brain" -ref "meanfunc_brain" -dof 6 -omat "t1_to_func_6dof.mat"
    flirt -in "t1_brain" -ref "meanfunc_brain" -dof 12 -init "t1_to_func_6dof.mat" -omat "t1_to_func_12dof.mat"
    
    echo "[5] Registration T1 -> mean functional nonlinear"
    fnirt --in="t1_brain" --aff="t1_to_func_12dof.mat" --ref="meanfunc_brain" --iout="t1_in_func_fnirt" --cout="func_warpcoef"

    echo "[6] Unwarping masks"
    applywarp --in="pve_0" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="csf_in_func" \
            --interp=trilinear

    applywarp --in="pve_1" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="gm_in_func" \
            --interp=trilinear
    
    applywarp --in="pve_2" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="wm_in_func" \
            --interp=trilinear
    
    echo "[7] Slice time correction with slicetimer"
    slicetimer -i "t2_s1_reoriented" -o "t2_s1_st"
    slicetimer -i "t2_s2_reoriented" -o "t2_s2_st"

    echo "[8] Motion correction with mcflirt"
    mcflirt -in "t2_s1_st" -out "t2_s1_mc" -plots
    mcflirt -in "t2_s2_st" -out "t2_s2_mc" -plots

    echo "[9] Intensity normalization"
    fslmaths "t2_s1_mc" -ing 1000 "t2_s1_norm"
    fslmaths "t2_s2_mc" -ing 1000 "t2_s2_norm"

    echo "=== Finished preprocessing subject: $subject ==="
}


# Loop over subjects
for t1_file in "$denoised_dir/*_denoised.nii"; do
    subject=$(basename "$t1_file" _denoised.nii)
    if [ "$subject" = "asjt" ]; then
        preprocess_subject "$subject"
    fi
done

## for visual check
## fsleyes mean_func_brain.nii.gz ${reg_dir}/${subject}_t1_in_func_fnirt.nii.gz &
## fsleyes wm_in_func.nii.gz gm_in_func.nii.gz csf_in_func.nii.gz t1_in_func_fnirt.nii.gz &