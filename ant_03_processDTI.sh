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

######################## PREPROCESSING ########################

preprocessDTI () {
    local subject_id=$1
    echo "Processing subject: $subject_id"
    echo "dMRI processing started at $(date '+%Y-%m-%d %H:%M:%S')"
    start_time=$(date +%s)

    ## set paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    mkdir -p "$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    raw_dwi="$ANTINOMICS_DIR/raws/dMRI/${subject_id}.rec"
    raw_anat="$ANTINOMICS_DIR/raws/sMRI_T1/${subject_id}.nii"

    ## Conversion
    dcm2niix -f raw_dwi -o $subject_dwi_dir $raw_dwi
    cd $subject_dwi_dir
    mrconvert raw_dwi.nii raw_dwi.mif -fslgrad raw_dwi.bvec raw_dwi.bval

    ### Preprocessing
    dwidenoise raw_dwi.mif dwi_den.mif -noise noise.mif
    mrcalc raw_dwi.mif dwi_den.mif -subtract residual.mif
    mrdegibbs dwi_den.mif dwi_gibb.mif
    dwifslpreproc dwi_gibb.mif dwi_preproc.mif -pe_dir ap -rpe_none  # roger says its okay
    dwibiascorrect ants dwi_preproc.mif dwi_unbiased.mif -bias bias.mif 
    echo -e "\e[32mPreprocessing is done successfuly!"
    dwi2mask dwi_preproc.mif mask.mif

    rm noise.mif residual.mif dwi_den.mif dwi_gibb.mif dwi_preproc.mif bias.mif
    mrtrix_cleanup $subject_dwi_dir
}

######################## TRACTOGRAPHY ########################

create_tractography () {
    local subject_id=$1
    echo "Processing subject: $subject_id"

    ### set paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    raw_anat="$ANTINOMICS_DIR/raws/sMRI_T1/${subject_id}.nii"
    tian_atlas_dir="/home/ubuntu/volume/tools/Tian_atlas"
    mni_ref="/home/ubuntu/volume/MNI_HCP/MNI152NLin6Asym_1mm.nii.gz"
    mni_ref_mask="/home/ubuntu/volume/MNI_HCP/MNI152NLin6Asym_1mm_brain.nii.gz"

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

    ## Register Tian atlas to subject T1 and DWI space
    mkdir $subject_dwi_dir/t1_temp/
    mkdir $subject_dwi_dir/tian/
    mrconvert raw_anat.mif t1_temp/raw_anat.nii.gz
    cd $subject_dwi_dir/t1_temp
    fslreorient2std raw_anat.nii.gz raw_anat_std.nii.gz
    bet raw_anat_std.nii.gz sub-T1_brain.nii.gz
    
    ## this part is wrong ... maybe switch to ants
    flirt -in sub-T1_brain.nii.gz \
            -ref $mni_ref \
            -out sub2MNI_lin \
            -omat sub2MNI_lin.mat

    fnirt --in=raw_anat_std.nii.gz \
        --aff=sub2MNI_lin.mat \
        --ref=$mni_ref \
        --refmask=$mni_ref_mask \
        --cout=sub2MNI_warp

    ##

    invwarp -w sub2MNI_warp \
            -r raw_anat_std.nii.gz \
            -o MNI2sub_warp

    for scale in S1 S2 S3 S4; do
        atlas_file="$tian_atlas_dir/Tian_Subcortex_${scale}_3T_1mm.nii.gz"
        out_t1=$subject_dwi_dir/tian/tian_${scale}_T1.nii.gz
        
        # Apply warp: MNI -> T1
        out_dwi=$subject_dwi_dir/tian/tian_${scale}_dwi.mif 
        applywarp -i "$atlas_file" -r "raw_anat_std.nii.gz" -o "$out_t1" -w "MNI2sub_warp" --interp=nn

        # Convert to MRtrix and transform to DWI space
        mrconvert "$out_t1" "${out_t1%.nii.gz}.mif"
        mrtransform "${out_t1%.nii.gz}.mif" -linear $subject_dwi_dir/diff2struct_mrtrix.txt -inverse -interp nearest "$out_dwi"
        
    done
        
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
    rm ribbon_core.mif

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

    rm Tian_sub_bin.nii.gz \
        Tian_sub_dil1.nii.gz \
        Tian_sub_dil2.nii.gz \
        Tian_sub_hull.mif \
        Tian_sub_hull.nii.gz \
        Tian_sub_hull_resampled.mif \
        Tian_sub_hull_resampled.nii.gz \
        tian_S3_dwi_resampled.nii.gz \
        cortical_ribbon.mif \
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
}

######################## CONNECTOME ########################

create_connectome () {
    local subject_id=$1
    echo "Processing subject: $subject_id"

    ### define paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    lut_dir="/home/ubuntu/volume/Schaefer_atlas"
    conn_dir=$subject_dwi_dir/conn
    mrtrix_home="/usr/local/mrtrix3"
    mkdir -p $conn_dir
    
    # now convert it and bring it to diffusion space
    cd $subject_dwi_dir
    transformcalc diff2struct_mrtrix.txt invert struct2diff_mrtrix.txt

    #### Schaefer atlas ###
    mkdir -p $subject_dwi_dir/atlases
    cd $subject_dwi_dir/atlases
    n_network=7
    for n_roi in 400 800 1000; do
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
                    -symmetric \
                    -zero_diagonal \
                    -scale_invnodevol
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
                    -zero_diagonal \
                    -scale_invnodevol

    #### Tian atlas ###
    for scale in S1 S2 S3 S4; do
        tck2connectome $subject_dwi_dir/tracks_subcortical_10M.tck \
                        $subject_dwi_dir/tian/tian_${scale}_dwi.mif \
                        $conn_dir/tian_${scale}_conn.csv \
                        -tck_weights_in $subject_dwi_dir/sift_subcortical_1M.txt \
                        -out_assignment $conn_dir/tian_${scale}_assign.txt \
                        -symmetric \
                        -zero_diagonal \
                        -scale_invnodevol
    done
}

######################## FIXEL-FIXEL ########################

compute_fixel_fixel_connectome () {
    ### set Paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    AFD_DIR="/home/ubuntu/volume/Antinomics/AFD"
    mkdir -p $AFD_DIR
    mkdir -p "$AFD_DIR/unbiased"
    mkdir -p "$AFD_DIR/masks"
    cd "$ANTINOMICS_DIR/subjects_mrtrix_dir"
    
    ### compute average response
    for_each * : ln -sr IN/dwi_unbiased.mif ../AFD/unbiased/IN.mif
    for_each * : ln -sr IN/mask.mif ../AFD/masks/IN.mif
    dwinormalise group ../AFD/unbiased/ ../AFD/masks/ ../AFD/normalized/ ../AFD/fa_template.mif ../AFD/fa_template_wm_mask.mif
    for_each ../AFD/normalized/* : ln -sr IN PRE/dwi_normalised.mif
    for_each * : dwi2response tournier IN/dwi_normalised.mif IN/response.txt
    responsemean */response.txt ../AFD/group_average_response.txt

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
    tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif -mask template_mask.mif -select 20000000 -cutoff 0.10 tracks_20_million.tck
    tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number 2000000

    ### Generate fixel-fixel connectivity matrix
    fixelconnectivity fixel_mask/ tracks_2_million_sift.tck matrix/
    fixelfilter fd smooth fd_smooth -matrix matrix/
    fixelfilter log_fc smooth log_fc_smooth -matrix matrix/
    fixelfilter fdc smooth fdc_smooth -matrix matrix/
}

run_stats_on_fixels () {
    ## stats
    nohup bash -c 'fixelcfestats fd_smooth/ subjects.txt design_matrix.txt contrast_CT.txt matrix/ stats_fd_CT/ && \
    fixelcfestats log_fc_smooth/ subjects.txt design_matrix.txt contrast_CT.txt matrix/ stats_log_fc_CT/ && \
    fixelcfestats fdc_smooth/ subjects.txt design_matrix.txt contrast_CT.txt matrix/ stats_fdc_CT/ && \
    fixelcfestats fd_smooth/ subjects.txt design_matrix.txt contrast_TC.txt matrix/ stats_fd_TC/ && \
    fixelcfestats log_fc_smooth/ subjects.txt design_matrix.txt contrast_TC.txt matrix/ stats_log_fc_TC/ && \
    fixelcfestats fdc_smooth/ subjects.txt design_matrix.txt contrast_TC.txt matrix/ stats_fdc_TC/' \
    > fixelcfestats_all.out 2>&1 &
}

######################## TOOLS ########################

create_exempler_files () {
    local subject_id=$1
    echo "Processing subject: $subject_id"
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"

    connectome2tck tracks_10M.tck \
                    sch_400_7n_assign.txt \
                    sch_${n}Parcels_7Networks_edge_exemplar.tck \
                    -files single \
                    -exemplars Schaefer2018_400_7Networks_DWI_converted.mif

    label2mesh Schaefer2018_400_7Networks_DWI_converted.mif mesh_400.obj
    meshfilter mesh_400.obj smooth mesh_400_smooth.obj
}


## extract some metrics
    fod2fixel wmfod_norm.mif fixel_dir/ -afd fd.mif -mask mask.mif
    tcksample -stat_tck mean -weight sift_coeffs.txt tracks_20M.tck fd.mif mean_fd_weighted.txt


for t1_file in /home/ubuntu/volume/Antinomics/raws/sMRI_T1/*.nii; do
    subject_id=$(basename "$t1_file" .nii)
    if [ "$subject_id" = "asjt" ]; then
        create_tractography "$subject_id"   
    fi
done


antsRegistrationSyNQuick.sh \
  -d 3 \
  -f /Users/payamsadeghishabestari/Antinomics/data/tpl-MNI152NLin6Asym_res-01_desc-brain_T1w.nii.gz \
  -m sub-T1_brain.nii.gz \
  -o sub2MNI_ \
  -x /Users/payamsadeghishabestari/Antinomics/data/MNI152NLin6Asym_1mm_brain_mask.nii.gz

################
antsApplyTransforms -d 3 \
  -r sub-T1_brain.nii.gz \
  -o sub2MNI_CompositeInvWarp.nii.gz \
  -t "[sub2MNI_0GenericAffine.mat,1]" \
  -t sub2MNI_1InverseWarp.nii.gz \
  --compose

for scale in S1 S2 S3 S4; do
  antsApplyTransforms -d 3 \
    -i /Users/payamsadeghishabestari/Antinomics/data/Tian2020MSA/3T/Subcortex-Only/Tian_Subcortex_${scale}_3T_1mm.nii.gz \
    -r sub-T1_brain.nii.gz \
    -o Tian_${scale}_inSubj.nii.gz \
    -t sub2MNI_CompositeInvWarp.nii.gz \
    -n NearestNeighbor
done

