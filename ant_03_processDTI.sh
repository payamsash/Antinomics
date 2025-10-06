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


## subject level
create_tractography () {
    local subject_id=$1
    echo "Processing subject: $subject_id"

    ### define paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    raw_anat="$ANTINOMICS_DIR/raws/sMRI_T1/${subject_id}.nii"
    atlas_dir="/home/ubuntu/volume/Tian_atlas"
    mni_ref="/home/ubuntu/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
    mni_ref_mask="/home/ubuntu/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz"

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

    ## Register Tian S3 atlas to subject T1 space
    mkdir t1_temp/
    mrconvert raw_anat.mif t1_temp/raw_anat.nii.gz
    cd t1_temp
    fslreorient2std raw_anat.nii.gz raw_anat_std.nii.gz
    bet raw_anat_std.nii.gz sub-T1_brain.nii.gz
    flirt -in sub-T1_brain.nii.gz \
            -ref $mni_ref \
            -out sub2MNI_lin \
            -omat sub2MNI_lin.mat

    fnirt --in=raw_anat_std.nii.gz \
        --aff=sub2MNI_lin.mat \
        --ref=$mni_ref \
        --refmask=$mni_ref_mask \
        --cout=sub2MNI_warp

    invwarp -w sub2MNI_warp \
            -r raw_anat_std.nii.gz \
            -o MNI2sub_warp

    applywarp -i "$atlas_dir/Tian_Subcortex_S3_3T_1mm.nii.gz" \
                -r "raw_anat_std.nii.gz" \
                -o "tian_subspace.nii.gz" \
                -w "MNI2sub_warp" \
                --interp=nn
    
    ## Register Tian S3 atlas to subject DWI space
    mrconvert tian_subspace.nii.gz ../Tian_subcortical.mif
    cd $subject_dwi_dir
    mrtransform Tian_subcortical.mif \
                -linear diff2struct_mrtrix.txt \
                -inverse \
                -interp nearest \
                Tian_subcortical_dwi.mif

    ## Generate GM–WM interface for subcortical tractography 
    mrcalc gmwmSeed_coreg.mif 0.5 -gt gmwmSeed_bin.mif
    mrgrid Tian_subcortical_dwi.mif regrid -template gmwmSeed_bin.mif -interp nearest Tian_subcortical_dwi_resampled.mif
    mrcalc gmwmSeed_bin.mif Tian_subcortical_dwi_resampled.mif -mul subcortical_gmwmi.mif

    ## Create exclusion masks (cortical ribbon exclusion)
    mrconvert 5tt_coreg.mif -coord 3 0 ribbon_core.mif
    mrconvert ribbon_core.mif -axes 0,1,2 cortical_ribbon.mif
    

    ## Create exclusion masks (convex hull)
    mrconvert Tian_subcortical_dwi_resampled.mif Tian_subcortical_dwi_resampled.nii.gz
    fslmaths Tian_subcortical_dwi_resampled.nii.gz -thr 0.5 -bin Tian_sub_bin.nii.gz
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

    ## cleaning
    rm -r t1_temp
    rm 5tt_fixed.nii.gz 5tt_vol0.nii.gz
    rm Tian_sub_bin.nii.gz \
        Tian_sub_dil1.nii.gz \
        Tian_sub_dil2.nii.gz \
        Tian_sub_hull.mif \
        Tian_sub_hull.nii.gz \
        Tian_sub_hull_resampled.mif \
        Tian_sub_hull_resampled.nii.gz \
        Tian_subcortical_dwi.mif \
        Tian_subcortical_dwi_resampled.nii.gz
    rm bias.mif cortical_ribbon.mif \
        dwi_den.mif dwi_gibb.mif noise.mif mask.nii.gz \
        mean_b0.mif residual.mif ribbon_core.mif
}


create_connectome () {
    local subject_id=$1
    echo "Processing subject: $subject_id"

    ### define paths
    ANTINOMICS_DIR="/home/ubuntu/volume/Antinomics"
    subject_fs_dir="$SUBJECTS_DIR/$subject_id"
    subject_dwi_dir="$ANTINOMICS_DIR/subjects_mrtrix_dir/$subject_id"
    lut_dir="/home/ubuntu/volume/Schaefer_atlas"


    # now convert it and bring it to diffusion space
    cd $subject_dwi_dir
    transformcalc diff2struct_mrtrix.txt invert struct2diff_mrtrix.txt

    for n_roi in 100 200 400 800 1000:
            mrconvert $subject_fs_dir/mri/Schaefer2018_${n_roi}_7Networks.mgz \
                        Schaefer2018_${n_roi}_${n_network}Networks_T1.mif

            mrtransform \
                        Schaefer2018_${n_roi}_${n_network}Networks_T1.mif \
                        -linear struct2diff_mrtrix.txt \
                        -interp nearest \
                        Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif

            mrcalc Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif 0.5 -add -floor Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif

            labelconvert \
                Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif \
                $lut_dir/Schaefer2018_${n_roi}Parcels_${n_network}Networks_order_LUT.txt \
                $lut_dir/Schaefer2018_${n_roi}Parcels_${n_network}Networks_order.txt \
                Schaefer2018_${n_roi}_${n_network}Networks_DWI_int_converted.mif

            rm Schaefer2018_${n_roi}_${n_network}Networks_DWI.mif \
                Schaefer2018_${n_roi}_${n_network}Networks_DWI_int.mif

            tck2connectome tracks_10M.tck \
                            Schaefer2018_400_7Networks_DWI_converted.mif \
                            sch_400_conn.csv \
                            -tck_weights_in sift_1M.txt \
                            -out_assignment sch_400_7n_assign.txt \
                            -symmetric \
                            -zero_diagonal \
                            -scale_invnodevol
            
            ## only for 1 subject
            connectome2tck tracks_10M.tck \
                            sch_400_7n_assign.txt \
                            sch_${n}Parcels_7Networks_edge_exemplar.tck \
                            -files single \
                            -exemplars Schaefer2018_400_7Networks_DWI_converted.mif


            label2mesh Schaefer2018_400_7Networks_DWI_converted.mif mesh_400.obj
            meshfilter mesh_400.obj smooth mesh_400_smooth.obj

}



# processDTI_4 () {

#     ## Create a Connectome for Atlases (pending from here)
#     labelconvert $subject_fs_dir/mri/aparc+aseg.mgz \
#                     $FREESURFER_HOME/FreeSurferColorLUT.txt \
#                     $LUT_DIR/fs_default.txt \
#                     aparc_parcels.mif
#     labelconvert $subject_fs_dir/mri/aparc.a2009s+aseg.mgz \
#                     $FREESURFER_HOME/FreeSurferColorLUT.txt \
#                     $LUT_DIR/fs_a2009s.txt \
#                     aparc2009s_parcels.mif

#     parcels=("aparc" "aparc.a2009s" "schaefer")
#     for parc in "${parcels[@]}"; do
#         if [[ $parc == "aparc" || $parc == "aparc.a2009s"]]; then
#             labelconvert $subject_fs_dir/mri/${parc}+aseg.mgz \
#                     $FREESURFER_HOME/FreeSurferColorLUT.txt \
#                     $LUT_DIR/fs_default.txt \
#                     ${parc}_parcels.mif
            
#             tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in \
#                         sift_1M.txt tracks_10M.tck ${parc}_parcels.mif ${parc}_parcels_connectome.csv \
#                         -out_assignment ${parc}_parcels_assignments.txt
#             label2mesh ${parc}_parcels.mif ${parc}_parcels_mesh.obj
#             meshfilter ${parc}_parcels_mesh.obj smooth ${parc}_parcels_mesh_smoothed.obj
#             connectome2tck tracks_10M.tck ${parc}_parcels_assignments.txt ${parc}_parcels_edge_exemplar.tck \
#                             -files single -exemplars ${parc}_parcels.mif
#         fi

#         if [[ $parc == "schaefer" ]]; then
#             for n in 400 800 1000; do
                
#                     mrconvert $subject_fs_dir/schaefer/${n}Parcels_7Networks.mgz sch_${n}Parcels_7Networks.mif

#                     tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in \
#                                         sift_1M.txt tracks_10M.tck sch_${n}Parcels_7Networks.mif \
#                                         sch_${n}Parcels_7Networks_connectome.csv \
#                                         -out_assignment sch_${n}Parcels_7Networks_assignments.txt
#                     label2mesh sch_${n}Parcels_7Networks.mif sch_${n}Parcels_7Networks_mesh.obj
#                     meshfilter sch_${n}Parcels_7Networks_mesh.obj smooth sch_${n}Parcels_7Networks_mesh_smoothed.obj
#                     connectome2tck tracks_10M.tck sch_${n}Parcels_7Networks_assignments.txt \
#                                     sch_${n}Parcels_7Networks_edge_exemplar.tck -files single \
#                                     -exemplars sch_${n}Parcels_7Networks.mif
#                 done
#             done
#         fi
#     done
# mrtrix_cleanup $subject_dwi_dir
# }

## extract some metrics
    fod2fixel wmfod_norm.mif fixel_dir/ -afd fd.mif -mask mask.mif
    tcksample -stat_tck mean -weight sift_coeffs.txt tracks_20M.tck fd.mif mean_fd_weighted.txt


for t1_file in /home/ubuntu/volume/Antinomics/raws/sMRI_T1/*.nii; do
    subject_id=$(basename "$t1_file" .nii)
    if [ "$subject_id" = "asjt" ]; then
        create_tractography "$subject_id"   
    fi
done



'''
Goal:
We want to reconstruct white matter tracts that connect subcortical nuclei (thalamus, basal ganglia, hippocampus, amygdala, etc.) while avoiding cortical regions or other unwanted areas.
1️⃣ Start with an atlas of subcortical regions
The Tian atlas provides a detailed map of subcortical nuclei in standard (MNI) space.
These nuclei are the “regions of interest” we want to seed for tractography.
We choose a scale (e.g., S3) that balances anatomical detail and seed robustness.
2️⃣ Bring the atlas into the subject’s native anatomical space
Every subject’s brain is slightly different in size, shape, and orientation.
We warp the atlas to the subject’s T1-weighted MRI so that the labels correctly match their subcortical anatomy.
Nonlinear registration is ideal here because it allows small nuclei to line up accurately with the subject’s brain.
3️⃣ Bring the labels into diffusion space
Tractography is performed on diffusion-weighted images (DWI).
We need to map the subcortical labels from T1 space into DWI space, so the seeds are in the correct location relative to the diffusion data.
4️⃣ Restrict seeds to the GM–WM interface
Streamlines are generated from the interface between gray matter and white matter, not deep inside gray matter or in CSF.
By multiplying the subcortical labels with a GM–WM interface mask, we ensure that seeds start from just inside the white matter connected to the subcortical nuclei.
5️⃣ Optional exclusion masks
Cortical ribbon mask: prevents streamlines from entering cortex.
Convex hull mask around subcortical structures: prevents streamlines from leaving the subcortical area too far.
These help constrain the tractogram to subcortical–subcortical connections, avoiding contamination from cortical or CSF regions.

When you seed from the GM–WM interface near subcortical regions (e.g. thalamus, striatum, hippocampus), streamlines can sometimes:
Leak into the cortical ribbon (i.e. grey matter sheet on the surface).
Shoot into CSF or outside the brain if the mask has imperfections.
Create spurious loops across sulci/gyri because the seeds are near boundaries.
👉 To prevent this, you add exclusion masks: regions where streamlines are not allowed to go.
If a streamline enters an exclusion mask, it gets terminated.



6️⃣ Generate the tractogram
Using the subcortical GM–WM seeds, we run probabilistic tractography with the ACT framework (anatomically constrained tractography).
ACT ensures that streamlines follow white matter pathways, terminate in GM, and avoid non-brain tissues.
The result is a tractogram connecting subcortical nuclei, which can be used for connectivity analyses.
'''

## remove 17 network
## add 100, 200 schaefer
## make proper foldering
## run aparc2aseg for all and move it to subsegment
## run connectome and proper foldering for it
## tian atlas
## glasser atlas
## connectome2tck plots

## mrview Schaefer2018_400_7Networks_DWI_converted.mif -connectome.init Schaefer2018_400_7Networks_DWI_converted.mif -connectome.load sch_400_conn.csv
## So for instance: You could load a structural connectome file as your connectome matrix and show only those edges where the connection density is above a certain threshold, but then set the colour of each edge based on a different matrix file that contains functional connectivity values.