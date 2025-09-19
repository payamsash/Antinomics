#!/bin/bash

# Diffusion MRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.

set -euo pipefail
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/Antinomics/subjects_fs_dir
export PATH=/home/ubuntu/ants-2.6.2/ants-2.6.2/bin:$PATH
export LUT_DIR=/usr/local/mrtrix3/share/mrtrix3/labelconvert
export ANTSPATH=/home/ubuntu/ants-2.6.2/ants-2.6.2/bin

processDTI () {
    local subject_id=$1
    echo "Processing subject: $subject_id"
    echo "dMRI processing started at $(date '+%Y-%m-%d %H:%M:%S')"
    start_time=$(date +%s)

    ## set paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    mkdir "$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    raw_dwi="$ANTINOMICS_DIR/raws/dMRI/${subject_id}.rec"
    raw_anat="$ANTINOMICS_DIR/raws/sMRI_T1/${subject_id}.nii"

    ## Conversion
    dcm2niix -f raw_dwi -o $subject_dwi_dir $raw_dwi
    cd $subject_dwi_dir
    mrconvert raw_dwi.nii raw_dwi.mif -fslgrad raw_dwi.bvec raw_dwi.bval
    echo -e "\e[32mConverted to nifti and mif formats successfuly!"

    ### Preprocessing
    dwidenoise raw_dwi.mif dwi_den.mif -noise noise.mif
    mrcalc raw_dwi.mif dwi_den.mif -subtract residual.mif
    mrdegibbs dwi_den.mif dwi_gibb.mif
    dwifslpreproc dwi_gibb.mif dwi_preproc.mif -pe_dir ap -rpe_none  # roger says its okay
    dwibiascorrect ants dwi_preproc.mif dwi_unbiased.mif -bias bias.mif 
    echo -e "\e[32mPreprocessing is done successfuly!"

    ### Constrained Spherical Deconvolution
    ## lets see how dwi2mask works, if bad: 1. I'll compute the mask from raw_anat then register it to the diffusion space
    ## lets see how dwi2mask works, if bad: 2. I'll provide the template image and corresponding mask from T2 data.
    dwi2mask dwi_preproc.mif mask.mif
    dwi2response tournier dwi_unbiased.mif wm_response.txt -voxels voxels.mif
    dwi2fod csd dwi_unbiased.mif -mask mask.mif wm_response.txt wmfod.mif
    mtnormalise wmfod.mif wmfod_norm.mif -mask mask.mif
    echo -e "\e[32mCSD is done successfuly!"

    ### Registration to anatomical image
    mrconvert $raw_anat raw_anat.mif
    5ttgen hsvs $subject_fs_dir 5tt_nocoreg.mif
    dwiextract dwi_unbiased.mif mean_b0.mif -bzero
    mrconvert mean_b0.mif mean_b0.nii.gz
    mrconvert 5tt_nocoreg.mif 5tt_nocoreg.nii.gz
    
    fslroi 5tt_nocoreg.nii.gz 5tt_vol0.nii.gz 0 1
    flirt -in mean_b0.nii.gz \
            -ref 5tt_vol0.nii.gz \
            -interp nearestneighbour \
            -dof 6 \
            -omat diff2struct_fsl.mat
    
    ### Generate GM–WM interface for ACT 
    transformconvert diff2struct_fsl.mat mean_b0.nii.gz 5tt_nocoreg.nii.gz flirt_import diff2struct_mrtrix.txt
    mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif
    5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif

    ### Tractography
    tckgen -act 5tt_coreg.mif \
                -backtrack \
                -seed_gmwmi gmwmSeed_coreg.mif \
                -select 10000000 \
                wmfod_norm.mif \
                tracks_10M.tck
    tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck
    tcksift2 -act 5tt_coreg.mif \
                -out_mu sift_mu.txt \
                -out_coeffs sift_coeffs.txt \
                tracks_10M.tck \
                wmfod_norm.mif \
                sift_1M.txt

    ## Create a Connectome for Atlases (pending from here)
    labelconvert $subject_fs_dir/mri/aparc+aseg.mgz \
                    $FREESURFER_HOME/FreeSurferColorLUT.txt \
                    $LUT_DIR/fs_default.txt \
                    aparc_parcels.mif
    labelconvert $subject_fs_dir/mri/aparc.a2009s+aseg.mgz \
                    $FREESURFER_HOME/FreeSurferColorLUT.txt \
                    $LUT_DIR/fs_a2009s.txt \
                    aparc2009s_parcels.mif

    parcels=("aparc" "aparc.a2009s" "schaefer")
    for parc in "${parcels[@]}"; do
        if [[ $parc == "aparc" || $parc == "aparc.a2009s"]]; then
            labelconvert $subject_fs_dir/mri/${parc}+aseg.mgz \
                    $FREESURFER_HOME/FreeSurferColorLUT.txt \
                    $LUT_DIR/fs_default.txt \
                    ${parc}_parcels.mif
            
            tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in \
                        sift_1M.txt tracks_10M.tck ${parc}_parcels.mif ${parc}_parcels_connectome.csv \
                        -out_assignment ${parc}_parcels_assignments.txt
            label2mesh ${parc}_parcels.mif ${parc}_parcels_mesh.obj
            meshfilter ${parc}_parcels_mesh.obj smooth ${parc}_parcels_mesh_smoothed.obj
            connectome2tck tracks_10M.tck ${parc}_parcels_assignments.txt ${parc}_parcels_edge_exemplar.tck \
                            -files single -exemplars ${parc}_parcels.mif
        fi

        if [[ $parc == "schaefer" ]]; then
            for n in 400 800 1000; do
                
                    mrconvert $subject_fs_dir/schaefer/${n}Parcels_7Networks.mgz sch_${n}Parcels_7Networks.mif

                    tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in \
                                        sift_1M.txt tracks_10M.tck sch_${n}Parcels_7Networks.mif \
                                        sch_${n}Parcels_7Networks_connectome.csv \
                                        -out_assignment sch_${n}Parcels_7Networks_assignments.txt
                    label2mesh sch_${n}Parcels_7Networks.mif sch_${n}Parcels_7Networks_mesh.obj
                    meshfilter sch_${n}Parcels_7Networks_mesh.obj smooth sch_${n}Parcels_7Networks_mesh_smoothed.obj
                    connectome2tck tracks_10M.tck sch_${n}Parcels_7Networks_assignments.txt \
                                    sch_${n}Parcels_7Networks_edge_exemplar.tck -files single \
                                    -exemplars sch_${n}Parcels_7Networks.mif
                done
            done
        fi
    done
mrtrix_cleanup $subject_dwi_dir
}

for t1_file in /home/ubuntu/volume/Antinomics/raws/sMRI_T1/*.nii; do
    subject_id=$(basename "$t1_file" .nii)
    processDTI "$subject_id"
done
