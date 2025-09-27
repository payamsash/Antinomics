#!/bin/bash
set -e

# fMRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.

denoised_dir="/home/ubuntu/volume/Antinomics/raws/t1_denoised"
fmri_dir="/home/ubuntu/volume/Antinomics/raws/fMRI"
subjects_fsl_dir="/home/ubuntu/volume/Antinomics/subjects_fsl_dir"
sessions=("s1" "s2")

preprocess_subject() {
    local subject=$1
    mkdir -p "$subjects_fsl_dir/$subject"
    cd "$subjects_fsl_dir/$subject"

    echo "[0] Reorient to MNI"
    fslreorient2std "$denoised_dir/${subject}_denoised.nii" "t1_reoriented"
    fslreorient2std "$fmri_dir/s1/${subject}.nii" "t2_s1_reoriented"
    fslreorient2std "$fmri_dir/s2/${subject}.nii" "t2_s2_reoriented"

    echo "[1] Brain extraction of T1"
    bet t1_reoriented t1_brain -m

    echo "[2] FAST segmentation"
    fast -B -o t1 t1_brain

    echo "[3] Creating mean functional for registration"
    fslmaths "t2_s1_reoriented" -Tmean "meanfunc"

    echo "[4] Brain extraction of functional"
    bet "meanfunc" "meanfunc_brain" -F -m

    echo "[5] Registration T1 -> mean functional linear 6 and 12 dof"
    flirt -in "t1_brain" -ref "meanfunc_brain" -dof 6 -omat "t1_to_func_6dof.mat"
    flirt -in "t1_brain" -ref "meanfunc_brain" -dof 12 -init "t1_to_func_6dof.mat" -omat "t1_to_func_12dof.mat"
    
    echo "[6] Registration T1 -> mean functional nonlinear"
    fnirt --in="t1_brain" --aff="t1_to_func_12dof.mat" --ref="meanfunc_brain" --iout="t1_in_func_fnirt" --cout="func_warpcoef"

    echo "[7] Unwarping masks"
    applywarp --in="t1_pve_0" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="csf_in_func" \
            --interp=trilinear

    applywarp --in="t1_pve_1" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="gm_in_func" \
            --interp=trilinear
    
    applywarp --in="t1_pve_2" \
            --ref="meanfunc_brain" \
            --warp="func_warpcoef" \
            --out="wm_in_func" \
            --interp=trilinear

    applywarp --in="t1_brain_mask" \
                --ref="meanfunc_brain" \
                --warp="func_warpcoef" \
                --out="brain_mask_in_func" \
                --interp=nn
    
    echo "[8] Slice time correction with slicetimer"
    slicetimer -i "t2_s1_reoriented" -o "t2_s1_st" --odd -r 2.5
    slicetimer -i "t2_s2_reoriented" -o "t2_s2_st" --odd -r 2.5

    echo "[9] Motion correction with mcflirt"
    mcflirt -in "t2_s1_st" -out "t2_s1_mc" -plots
    mcflirt -in "t2_s2_st" -out "t2_s2_mc" -plots

    echo "[10] Intensity normalization"
    fslmaths "t2_s1_mc" -ing 1000 "t2_s1_norm"
    fslmaths "t2_s2_mc" -ing 1000 "t2_s2_norm"

    echo "[11] Cleaning directory"
    rm t2_s1_st.nii.gz t2_s2_st.nii.gz t2_s1_mc.nii.gz t2_s2_mc.nii.gz
    echo "=== Finished preprocessing subject: $subject ==="

}

# Loop over subjects
for t1_file in "$denoised_dir"/*_denoised.nii; do
    subject=$(basename "$t1_file" _denoised.nii)
    if [ "$subject" != "asjt" ]; then
        preprocess_subject "$subject"
    fi
done

## for visual check
## fsleyes meanfunc_brain t1_in_func_fnirt.nii.gz &
## fsleyes wm_in_func gm_in_func csf_in_func t2_s1_norm &