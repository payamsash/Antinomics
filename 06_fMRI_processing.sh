#!/bin/bash
set -euo pipefail

# Diffusion MRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.


# Paths – adjust as needed
t1_dir="/home/ubuntu/volume/Antinomics/subjects_fsl_dir"
fmri_dir="/home/ubuntu/volume/Antinomics/raws/fMRI"
sessions=("s1" "s2")     # list all sessions here

mkdir -p "$t1_dir/02_fast" "$t1_dir/03_reg" "$t1_dir/04_funcprep"
t1_fast_dir="$t1_dir/02_fast"
reg_dir="$t1_dir/03_reg"
funcprep_dir="$t1_dir/04_funcprep"


# Function to preprocess a single subject
preprocess_subject() {
    local subject=$1
    echo "=== Preprocessing subject: $subject ==="

    # ------------------ 0. Reorient to MNI ------------------
    fslreorient2std asjt_denoised.nii reoriented_asjt_denoised.nii.gz # fix here

    # ------------------ 1. FAST segmentation ------------------
    t1_file="$t1_dir/01_denoised/${subject}_denoised.nii"
    echo "[1] FAST segmentation"
    fast -B -o "$t1_fast_dir/${subject}" "$t1_file"

    # ------------------ 2. Registration reference ------------------
    first_session="${sessions[0]}"
    bold_ref="$fmri_dir/$first_session/${subject}.nii"
    mean_func="$reg_dir/${subject}_${first_session}_meanfunc.nii.gz"
    echo "[2] Creating mean functional for registration"
    fslmaths "$bold_ref" -Tmean "$mean_func"

    # Brain extraction
    bet "$mean_func" "${mean_func%.nii.gz}_brain.nii.gz" -F -m
    mean_func_brain="${mean_func%.nii.gz}_brain.nii.gz"

    # ------------------ 3. Registration T1 -> mean functional ------------------
    echo "[3] FLIRT 6 dof"
    flirt -in "$t1_brain" -ref "$mean_func_brain" -dof 6 \
      -omat "$reg_dir/${subject}_6dof.mat"

    ## flirt -in t1_brain.nii.gz -ref mean_func_brain.nii.gz -dof 6 -omat t1_to_func_6dof.mat
    ## flirt -in t1_brain.nii.gz -ref mean_func_brain.nii.gz -dof 12 -init t1_to_func_6dof.mat -omat t1_to_func_12dof.mat
    # fnirt --in="asjt_denoised_reoriented.nii.gz" \
    #   --aff="t1_to_func_12dof.mat" \
    #   --ref="mean_func_brain.nii.gz" \
    #   --iout="t1_in_func_fnirt.nii.gz" \
    #   --cout="func_warpcoef.nii.gz"

    ## fsleyes mean_func_brain.nii.gz ${reg_dir}/${subject}_t1_in_func_fnirt.nii.gz &
    ## fsleyes wm_in_func.nii.gz gm_in_func.nii.gz csf_in_func.nii.gz t1_in_func_fnirt.nii.gz &
    applywarp --in="pve_0.nii.gz" \
            --ref="mean_func_brain.nii.gz" \
            --warp="func_warpcoef.nii.gz" \
            --out="csf_in_func.nii.gz" \
            --interp=trilinear

    applywarp --in="pve_1.nii.gz" \
            --ref="mean_func_brain.nii.gz" \
            --warp="func_warpcoef.nii.gz" \
            --out="gm_in_func.nii.gz" \
            --interp=trilinear
    
    applywarp --in="pve_2.nii.gz" \
            --ref="mean_func_brain.nii.gz" \
            --warp="func_warpcoef.nii.gz" \
            --out="wm_in_func.nii.gz" \
            --interp=trilinear
    
    # ------------------ 4. Loop over sessions ------------------
    for session in "${sessions[@]}"; do
        echo "[4] Processing session: $session"
        bold="$fmri_dir/$session/${subject}.nii"
        outdir="$funcprep_dir/${session}_${subject}"
        mkdir -p "$outdir"

        # ------------------ 4b. Slice timing correction ------------------
        echo "[5] Slice time correction with slicetimer"
        slicetimer -i "asjt_s1_reoriented.nii.gz" -o "asjt_st"

        # ------------------ 4c. Motion correction ------------------
        echo "[6] Motion correction with mcflirt"
        mcflirt -in "asjt_st" -out "asjt_mc" -plots

        # ------------------ 4d. Intensity normalization ------------------
        echo "[7] Intensity normalization"
        fslmaths "asjt_mc" -ing 1000 "asjt_norm"
    done
    echo "=== Finished preprocessing subject: $subject ==="
}


# Loop over subjects
for t1_file in "$t1_dir"/01_denoised/*_denoised.nii; do
    subject=$(basename "$t1_file" _denoised.nii)
    if [ "$subject" = "asjt" ]; then
        preprocess_subject "$subject"
    fi
done
