#!/bin/bash

# Diffusion MRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.

set -e
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/Antinomics/subjects_fs_dir
export PATH=/home/ubuntu/ants-2.6.2/ants-2.6.2/bin:$PATH
export LUT_DIR=/usr/local/mrtrix3/share/mrtrix3/labelconvert
export ANTSPATH=/home/ubuntu/ants-2.6.2/ants-2.6.2/bin
source $FREESURFER_HOME/SetUpFreeSurfer.sh


create_tractography () {
    local subject_id=$1
    echo "Processing subject: $subject_id"

    ### set paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    raw_anat="$ANTINOMICS_DIR/raws/sMRI_T1/${subject_id}.nii"
    tian_atlas_dir="/home/ubuntu/volume/Tian_atlas"
    mni_hcp_dir="/home/ubuntu/volume/MNI_HCP"
    lut_dir="/home/ubuntu/volume/Schaefer_atlas"
    conn_dir=$subject_dwi_dir/conn
    mrtrix_home="/usr/local/mrtrix3"
    mkdir -p $conn_dir

    cd $subject_dwi_dir

    ### get the b0 image
    mrconvert $raw_anat raw_anat.mif
    dwiextract dwi_unbiased.mif mean_b0.mif -bzero
    mrconvert mean_b0.mif mean_b0.nii.gz

    ## create 5 tissue type and remove pons from lesion
    5ttgen hsvs "$subject_fs_dir" 5tt_nocoreg.mif -thalami nuclei -hippocampi aseg
    for i in 0 1 2 3 4; do
        mrconvert 5tt_nocoreg.mif -coord 3 $i vol${i}.mif
    done
    mrcalc vol4.mif 0 -mul vol4_zero.mif
    mrcat -axis 3 vol0.mif vol1.mif vol2.mif vol3.mif vol4_zero.mif 5tt_fixed.mif
    rm vol0.mif vol1.mif vol2.mif vol3.mif vol4.mif vol4_zero.mif

    ## Registration to anatomical image
    mrconvert 5tt_fixed.mif 5tt_fixed.nii.gz
    fslroi 5tt_fixed.nii.gz 5tt_vol0.nii.gz 0 1
    flirt -in mean_b0.nii.gz \
            -ref 5tt_vol0.nii.gz \
            -interp nearestneighbour \
            -dof 6 \
            -omat diff2struct_fsl.mat

    ## Generate GM–WM interface for global tractography 
    transformconvert diff2struct_fsl.mat mean_b0.nii.gz 5tt_fixed.nii.gz flirt_import diff2struct_mrtrix.txt
    mrtransform 5tt_fixed.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif
    5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif
    rm 5tt_fixed.nii.gz 5tt_fixed.mif 5tt_vol0.nii.gz

    rm -f 5tt_fixed.mif gmwmSeed_bin.mif subcortical_gmwmi.mif \
        Tian_subcortical.mif Tian_subcortical_dwi_resampled.mif hull_exclude.mif hull_exclude.nii.gz \
        cortical_ribbon_reg.mif hull_exclude_reg.mif tracks_subcortical_10M.tck sift_subcortical_mu.txt \
        sift_subcortical_coeffs.txt sift_subcortical_1M.txt

    ## Register Tian atlas to subject T1 and DWI space
    mkdir $subject_dwi_dir/t1_temp/
    mkdir $subject_dwi_dir/tian/
    mrconvert raw_anat.mif $subject_dwi_dir/t1_temp/raw_anat.nii.gz
    cd $subject_dwi_dir/t1_temp
    fslreorient2std raw_anat.nii.gz raw_anat_std.nii.gz
    bet raw_anat_std.nii.gz sub-T1_brain.nii.gz
    
    antsRegistrationSyNQuick.sh -d 3 \
                                -f $mni_hcp_dir/tpl-MNI152NLin6Asym_res-01_desc-brain_T1w.nii.gz \
                                -m sub-T1_brain.nii.gz \
                                -o sub2MNI_ \
                                -x $mni_hcp_dir/tpl-MNI152NLin6Asym_res-01_desc-brain_mask.nii.gz

    for scale in S1 S2 S3 S4; do
        antsApplyTransforms -d 3 \
                            -i $tian_atlas_dir/Tian_Subcortex_${scale}_3T_1mm.nii.gz \
                            -r sub-T1_brain.nii.gz \
                            -o $subject_dwi_dir/tian/tian_${scale}_T1.nii.gz \
                            -t "[sub2MNI_0GenericAffine.mat,1]" \
                            -t sub2MNI_1InverseWarp.nii.gz \
                            -n NearestNeighbor
    
        mrconvert $subject_dwi_dir/tian/tian_${scale}_T1.nii.gz tian_${scale}_T1.mif
        mrtransform tian_${scale}_T1.mif -linear $subject_dwi_dir/diff2struct_mrtrix.txt -inverse -interp nearest $subject_dwi_dir/tian/tian_${scale}_dwi.mif
    done
    
    cd $subject_dwi_dir
    rm -r $subject_dwi_dir/t1_temp

    ## Generate GM–WM interface for subcortical tractography 
    mrcalc gmwmSeed_coreg.mif 0.5 -gt gmwmSeed_bin.mif
    mrgrid $subject_dwi_dir/tian/tian_S3_dwi.mif \
            regrid \
            -template gmwmSeed_bin.mif \
            -interp nearest \
            tian_S3_dwi_resampled.mif

    mrcalc gmwmSeed_bin.mif tian_S3_dwi_resampled.mif -mul subcortical_gmwmi.mif

    ## Create exclusion masks (cortical ribbon exclusion)
    mrconvert 5tt_coreg.mif -coord 3 0 ribbon_core.mif
    mrconvert ribbon_core.mif -axes 0,1,2 cortical_ribbon.mif

    ## Create exclusion masks (convex hull)
    mrconvert tian_S3_dwi_resampled.mif tian_S3_dwi_resampled.nii.gz
    fslmaths tian_S3_dwi_resampled.nii.gz -thr 0.5 -bin Tian_sub_bin.nii.gz
    fslmaths Tian_sub_bin.nii.gz -kernel sphere 3 -dilM Tian_sub_dil1.nii.gz
    fslmaths Tian_sub_dil1.nii.gz -kernel sphere 3 -dilM Tian_sub_dil2.nii.gz
    fslmaths Tian_sub_dil2.nii.gz -kernel sphere 2 -ero Tian_sub_hull.nii.gz
    mrconvert Tian_sub_hull.nii.gz Tian_sub_hull.mif
    mrconvert mask.mif mask.nii.gz
    mrgrid Tian_sub_hull.nii.gz regrid -template mask.mif -interp nearest Tian_sub_hull_resampled.mif
    mrconvert Tian_sub_hull_resampled.mif Tian_sub_hull_resampled.nii.gz
    fslmaths mask.nii.gz -sub Tian_sub_hull_resampled.nii.gz -bin hull_exclude.nii.gz
    mrconvert hull_exclude.nii.gz hull_exclude.mif

    ## create exclusion mask (cortical ribbon + convex hull)
    mrgrid cortical_ribbon.mif regrid -template mask.mif -interp nearest cortical_ribbon_reg.mif
    mrgrid hull_exclude.mif regrid -template mask.mif -interp nearest hull_exclude_reg.mif

    rm \
        gmwmSeed_bin.mif \
        ribbon_core.mif \
        Tian_sub_bin.nii.gz \
        Tian_sub_dil1.nii.gz \
        Tian_sub_dil2.nii.gz \
        Tian_sub_hull.mif \
        Tian_sub_hull.nii.gz \
        Tian_sub_hull_resampled.mif \
        Tian_sub_hull_resampled.nii.gz \
        tian_S3_dwi_resampled.nii.gz \
        cortical_ribbon.mif \
        hull_exclude.nii.gz \
        mask.nii.gz

    ## subcortical tractography
    voxsize=$(mrinfo wmfod.mif -spacing | awk '{print $1}')
    step=$(echo "$voxsize * 0.25" | bc -l)
    tckgen -algorithm iFoD2 \
            -act 5tt_coreg.mif \
            -backtrack \
            -angle 45 \
            -select 10000000 \
            -step $step \
            -exclude cortical_ribbon_reg.mif \
            -exclude hull_exclude_reg.mif \
            -seed_image subcortical_gmwmi.mif \
            wmfod.mif \
            tracks_subcortical_10M.tck

    tcksift2 -act 5tt_coreg.mif \
            -out_mu sift_subcortical_mu.txt \
            -out_coeffs sift_subcortical_coeffs.txt \
            tracks_subcortical_10M.tck \
            wmfod.mif \
            sift_subcortical_1M.txt

    ## global Tractography
    tckgen -algorithm iFoD2 \
            -act 5tt_coreg.mif \
            -backtrack \
            -angle 45 \
            -select 10000000 \
            -seed_gmwmi gmwmSeed_coreg.mif \
            wmfod.mif \
            tracks_10M.tck

    tcksift2 -act 5tt_coreg.mif \
            -out_mu sift_mu.txt \
            -out_coeffs sift_coeffs.txt \
            tracks_10M.tck \
            wmfod.mif \
            sift_1M.txt

    ################
    cd $subject_dwi_dir
    transformcalc diff2struct_mrtrix.txt invert struct2diff_mrtrix.txt

    #### Schaefer atlas ###
    mkdir -p $subject_dwi_dir/atlases
    cd $subject_dwi_dir/atlases
    n_network=7
    for n_roi in 100 200 400 800 1000; do
        mrconvert --datatype uint32 \
                    $subject_fs_dir/mri/Schaefer2018_${n_roi}_7Networks.mgz \
                    Schaefer2018_${n_roi}_${n_network}Networks_T1.mif

        mrtransform Schaefer2018_${n_roi}_${n_network}Networks_T1.mif \
                    -linear $subject_dwi_dir/struct2diff_mrtrix.txt \
                    -interp nearest \
                    Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif
        mrcalc Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif 0.5 -add -floor Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif
        labelconvert Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif \
                    $lut_dir/Schaefer2018_${n_roi}Parcels_${n_network}Networks_order_LUT.txt \
                    $lut_dir/Schaefer2018_${n_roi}Parcels_${n_network}Networks_order.txt \
                    Schaefer2018_${n_roi}_${n_network}Networks_DWI_converted.mif
        
        rm  Schaefer2018_${n_roi}_${n_network}Networks_T1.mif \
            Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif \
            Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif
        
        tck2connectome $subject_dwi_dir/tracks_10M.tck \
                    Schaefer2018_${n_roi}_${n_network}Networks_DWI_converted.mif \
                    $conn_dir/sch_${n_roi}_conn.csv \
                    -tck_weights_in $subject_dwi_dir/sift_1M.txt \
                    -out_assignment $conn_dir/sch_${n_roi}_assign.txt \
                    -symmetric
    done
    
    #### Glasser atlas ###
    mrconvert --datatype uint32 \
                $subject_fs_dir/mri/HCPMMP1.mgz \
                HCPMMP1_T1.mif
    mrtransform HCPMMP1_T1.mif \
                -linear $subject_dwi_dir/struct2diff_mrtrix.txt \
                -interp nearest \
                HCPMMP1_DWI.mif
    # mrcalc HCPMMP1_DWI.mif 0.5 -add -floor HCPMMP1_DWI_int.mif
    labelconvert HCPMMP1_DWI.mif \
                $mrtrix_home/share/mrtrix3/labelconvert/hcpmmp1_original.txt \
                $mrtrix_home/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt \
                HCPMMP1_DWI_converted.mif
    rm HCPMMP1_DWI.mif HCPMMP1_T1.mif
    tck2connectome $subject_dwi_dir/tracks_10M.tck \
                    HCPMMP1_DWI_converted.mif \
                    $conn_dir/glasser_conn.csv \
                    -tck_weights_in $subject_dwi_dir/sift_1M.txt \
                    -out_assignment $conn_dir/glasser_assign.txt \
                    -symmetric \
                    -zero_diagonal

    #### Tian atlas ###
    for scale in S1 S2 S3 S4; do
        tck2connectome $subject_dwi_dir/tracks_subcortical_10M.tck \
                        $subject_dwi_dir/tian/tian_${scale}_dwi.mif \
                        $conn_dir/tian_${scale}_conn.csv \
                        -tck_weights_in $subject_dwi_dir/sift_subcortical_1M.txt \
                        -out_assignment $conn_dir/tian_${scale}_assign.txt \
                        -symmetric \
                        -zero_diagonal
    done

}


for t1_file in /home/ubuntu/volume/Antinomics/raws/sMRI_T1/*.nii; do
    subject_id=$(basename "$t1_file" .nii)
    create_tractography "$subject_id"
done