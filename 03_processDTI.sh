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

preprocessDTI () {
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
    dwi2mask dwi_preproc.mif mask.mif
    mrtrix_cleanup $subject_dwi_dir
}

compute_fixel_fixel_connectome () {
    ### set Paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    AFD_DIR="/home/ubuntu/volume/Antinomics/AFD"
    mkdir -p $AFD_DIR
    mkdir -p "$AFD_DIR/unbiased"
    cd "$ANTINOMICS_DIR/subjects_mrtrix_dir"
    
    ### compute average response
    for_each * : ln -sr IN/dwi_unbiased.mif ../AFD/unbiased/IN.mif
    for_each * : ln -sr IN/mask.mif ../AFD/masks/IN.mif
    dwinormalise group ../AFD/unbiased/ ../AFD/masks/ ../AFD/normalized/ ../AFD/fa_template.mif ../AFD/fa_template_wm_mask.mif
    for_each ../AFD/normalized/* : ln -sr IN PRE/dwi_normalised.mif
    for_each * : dwi2response tournier IN/dwi_normalised.mif IN/response.txt
    responsemean */response.txt ../AFD/group_average_response.txt
    rm -r ../AFD/normalized
    rm -r ../AFD/masks

    ### upsample the dwi data
    for_each * : mrgrid IN/dwi_normalised.mif regrid -vox 1.25 IN/dwi_upsampled.mif
    for_each * : dwi2mask IN/dwi_upsampled.mif IN/mask_upsampled.mif

    ### Fibre Orientation Distribution estimation
    for_each * : dwiextract IN/dwi_upsampled.mif - \| dwi2fod msmt_csd - ../AFD/group_average_response.txt IN/wmfod.mif -mask IN/mask_upsampled.mif

    ### Generate a study-specific unbiased FOD template
    mkdir -p ../diffusion_template/fod_input
    mkdir ../diffusion_template/mask_input
    for_each * : ln -sr IN/wmfod.mif ../diffusion_template/fod_input/PRE.mif
    for_each * : ln -sr IN/mask_upsampled.mif ../diffusion_template/mask_input/PRE.mif
    population_template ../diffusion_template/fod_input -mask_dir ../diffusion_template/mask_input ../diffusion_template/wmfod_template.mif -voxel_size 1.25

    ### Register all subject FOD images to the FOD template
    for_each * : mrregister IN/wmfod.mif -mask1 IN/mask_upsampled.mif ../diffusion_template/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif

    ### Compute the template mask
    for_each * : mrtransform IN/mask_upsampled.mif -warp IN/subject2template_warp.mif -interp nearest -datatype bit IN/mask_in_template_space.mif
    mrmath */mask_in_template_space.mif min ../diffusion_template/template_mask.mif -datatype bit

    ### Compute a white matter template analysis fixel mask
    fod2fixel -mask ../diffusion_template/template_mask.mif -fmls_peak_value 0.10 ../diffusion_template/wmfod_template.mif ../diffusion_template/fixel_mask
    for_each * : mrtransform IN/wmfod.mif -warp IN/subject2template_warp.mif -reorient_fod no IN/fod_in_template_space_NOT_REORIENTED.mif
    for_each * : fod2fixel -mask ../diffusion_template/template_mask.mif IN/fod_in_template_space_NOT_REORIENTED.mif IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif
    for_each * : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif IN/fixel_in_template_space
    for_each * : rm -r IN/fixel_in_template_space_NOT_REORIENTED

    ### Assign subject fixels to template fixels
    for_each * : fixelcorrespondence IN/fixel_in_template_space/fd.mif ../diffusion_template/fixel_mask ../diffusion_template/fd PRE.mif

    ### Compute the fibre cross-section (FC) metric
    for_each * : warp2metric IN/subject2template_warp.mif -fc ../diffusion_template/fixel_mask ../diffusion_template/fc IN.mif
    mkdir ../diffusion_template/log_fc
    cp ../diffusion_template/fc/index.mif ../diffusion_template/fc/directions.mif ../diffusion_template/log_fc
    for_each * : mrcalc ../diffusion_template/fc/IN.mif -log ../diffusion_template/log_fc/IN.mif

    ### Compute a combined measure of fibre density and cross-section (FDC)
    mkdir ../diffusion_template/fdc
    cp ../diffusion_template/fc/index.mif ../diffusion_template/fdc
    cp ../diffusion_template/fc/directions.mif ../diffusion_template/fdc
    for_each * : mrcalc ../diffusion_template/fd/IN.mif ../diffusion_template/fc/IN.mif -mult ../diffusion_template/fdc/IN.mif

    ### Perform whole-brain fibre tractography on the FOD template
    cd ../diffusion_template
    tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif -mask template_mask.mif -select 10000000 -cutoff 0.10 tracks_10_million.tck
    tcksift tracks_10_million.tck wmfod_template.mif tracks_1_million_sift.tck -term_number 1000000

    ### Generate fixel-fixel connectivity matrix
    fixelconnectivity fixel_mask/ tracks_1_million_sift.tck matrix/
    fixelfilter fd smooth fd_smooth -matrix matrix/
    fixelfilter log_fc smooth log_fc_smooth -matrix matrix/
    fixelfilter fdc smooth fdc_smooth -matrix matrix/
}

create_tractography () {

    ### Constrained Spherical Deconvolution
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
}








processDTI_4 () {
    cd ../diffusion_template 
    ### Perform statistical analysis of FD, FC, and FDC
    fixelcfestats fd_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_fd/
    fixelcfestats log_fc_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_log_fc/
    fixelcfestats fdc_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_fdc/
}






    ### Constrained Spherical Deconvolution
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


# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/st_fibre_density_cross-section.html#introduction
# FBA pinpoints where microstructure changes are happening; connectome analysis shows whether those local changes propagate into altered network connectivity.