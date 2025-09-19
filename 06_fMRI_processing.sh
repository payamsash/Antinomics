#!/bin/bash
set -euo pipefail

# -------------------------------------------------------------------------
# Paths – adjust as needed
# -------------------------------------------------------------------------
t1_dir="/home/ubuntu/volume/Antinomics/subjects_fsl_dir"
fmri_dir="/home/ubuntu/volume/Antinomics/raws/fMRI"
sessions=("s1" "s2")     # list all sessions here

mkdir -p "$t1_dir/02_fast" "$t1_dir/03_reg" "$t1_dir/04_funcprep"
t1_fast_dir="$t1_dir/02_fast"
reg_dir="$t1_dir/03_reg"
funcprep_dir="$t1_dir/04_funcprep"

# -------------------------------------------------------------------------
# Function to preprocess a single subject
# -------------------------------------------------------------------------
preprocess_subject() {
    local subject=$1
    echo "=== Preprocessing subject: $subject ==="

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
    flirt -in "$t1_file" -ref "$mean_func_brain" -dof 6 \
          -omat "$reg_dir/${subject}_6dof.mat"

    echo "[3] FLIRT 12 dof"
    flirt -in "$t1_file" -ref "$mean_func_brain" -dof 12 \
          -omat "$reg_dir/${subject}_12dof.mat"

    echo "[3] FNIRT nonlinear warp"
    fnirt --in="$t1_file" \
          --aff="$reg_dir/${subject}_12dof.mat" \
          --ref="$mean_func_brain" \
          --iout="$reg_dir/${subject}_fnirt.nii.gz" \
          --cout="$reg_dir/${subject}_warpcoef.nii.gz"

    # ------------------ 4. Loop over sessions ------------------
    for session in "${sessions[@]}"; do
        echo "[4] Processing session: $session"
        bold="$fmri_dir/$session/${subject}.nii"
        outdir="$funcprep_dir/${session}_${subject}"
        mkdir -p "$outdir"

        # ------------------ 4a. Skip applytopup if acqparams.txt is missing ------------------
        if [[ -f "$funcprep_dir/acqparams.txt" && -f "$funcprep_dir/topup_results_field.nii.gz" ]]; then
            echo "[4a] BOLD unwarping with applytopup"
            applytopup --imain="$bold" \
                       --datain="$funcprep_dir/acqparams.txt" \
                       --inindex=1 \
                       --topup="$funcprep_dir/topup_results" \
                       --out="$outdir/${subject}_unwarped" \
                       --method=jac
        else
            echo "[4a] Skipping applytopup, copying original BOLD"
            cp "$bold" "$outdir/${subject}_unwarped.nii.gz"
        fi

        # ------------------ 4b. Slice timing correction ------------------
        slicetimer -i "$outdir/${subject}_unwarped.nii.gz" \
                   -o "$outdir/${subject}_st"

        # ------------------ 4c. Motion correction ------------------
        mcflirt -in "$outdir/${subject}_st" \
                -out "$outdir/${subject}_mc" -plots -mats

        # ------------------ 4d. Intensity normalization ------------------
        fslmaths "$outdir/${subject}_mc" -ing 1000 \
                 "$outdir/${subject}_norm"
    done

    echo "=== Finished preprocessing subject: $subject ==="
}

# Loop over subjects in 01_denoised folder
for t1_file in "$t1_dir"/01_denoised/*_denoised.nii; do
    subject=$(basename "$t1_file" _denoised.nii)
    if [ "$subject" = "asjt" ]; then
        preprocess_subject "$subject"
    fi
done
