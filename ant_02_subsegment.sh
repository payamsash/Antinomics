#!/bin/bash

# Written by Payam S. Shabestari, Zurich, 09.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.

## set Paths
# set -euo pipefail
export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export SUBJECTS_DIR=/home/ubuntu/volume/Antinomics/subjects_fs_dir
export LD_LIBRARY_PATH=$FREESURFER_HOME/MCRv97/runtime/glnxa64:$FREESURFER_HOME/MCRv97/bin/glnxa64:$FREESURFER_HOME/MCRv97/sys/os/glnxa64:$FREESURFER_HOME/MCRv97/extern/bin/glnxa64
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=42
source $FREESURFER_HOME/SetUpFreeSurfer.sh

export PATH=/usr/lib/mrtrix3/bin:$PATH
export PATH=/home/ubuntu/fsl/bin:$PATH
export PATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin:$PATH
export ANTSPATH=/home/ubuntu/data/src_codes/ants-2.5.4/bin

subsegment () {
    local subject_id=$1
    echo "Processing subject: $subject_id"


    ## set variables
    sch_gcs_dir="$FREESURFER_HOME/gcs"
    n_threads=30
    mask_options=("Tight" "Loose")
    hemis=("lh" "rh")

    ## check if all subjects are segmented by recon-all
    echo -e "\e[32mChecking for missing recon-all.done files..."
    for dir in "$SUBJECTS_DIR"/*/; do
        if [ -d "$dir" ]; then
            if [ ! -f "$dir/scripts/recon-all.done" ]; then
            echo "Missing: $(basename "$dir")"
            fi
        fi
    done

    ## the real part
    log_file="$SUBJECTS_DIR/$subject_id/scripts/sub_segmentation.log"
    echo "$(date): Sub segmentation for subject $subject_id has been started." >> "$log_file"

    # hippocampal subfields and nuclei of the amygdala
    echo -e "\e[32mSegmentation of hippocampal subfields and nuclei of the amygdala!"
    T2_scan=/home/ubuntu/volume/Antinomics/raws/sMRI_T2/${subject_id}.nii
    segmentHA_T2.sh $subject_id $T2_scan "T2" 1
    echo "$(date): Hippocampal segmentation completed for subject $subject_id" >> "$log_file"

    # brainstem
    echo -e "\e[32mSegmentation of Brainstem Substructures!"
    segmentBS.sh $subject_id
    echo "$(date): Brainstem segmentation completed for subject $subject_id" >> "$log_file"

    # ascending arousal network
    echo -e "\e[32mSegmentations of brainstem nuclei that are part of the Ascending Arousal Network!"
    # sudo chmod +x $FREESURFER_HOME/bin/segmentNuclei
    segmentAAN.sh $subject_id
    echo "$(date): AAN segmentation completed for subject $subject_id" >> "$log_file"

    # thalamic nuclei
    echo -e "\e[32mSegmentation of the thalamic nuclei"
    segmentThalamicNuclei.sh $subject_id

    # hypothalamus
    echo -e "\e[32mSegmentation of the hypothalamus and its associated subunits"
    mri_segment_hypothalamic_subunits --s $subject_id --threads $n_threads
    echo "$(date): Hypothalamus segmentation completed for subject $subject_id" >> "$log_file"

    # striatum
    echo -e "\e[32mStriatal parcellation!"
    for mask_option in "${mask_options[@]}"; do
        mri_vol2vol --mov $SUBJECTS_DIR/$subject_id/mri/norm.mgz \
                    --s $subject_id \
                    --targ $SUBJECTS_DIR/MNI152/choi_atlas/17_network_${mask_option}_mask.nii.gz \
                    --m3z $SUBJECTS_DIR/MNI152/mri/transforms/talairach.m3z \
                    --noDefM3zPath \
                    --o $SUBJECTS_DIR/$subject_id/mri/striatum_17_network_${mask_option}_mask.nii.gz \
                    --inv-morph \
                    --interp nearest
        echo "$(date): Striatal parcellation with $mask_option mask completed for subject $subject_id" >> "$log_file"
    done

    # cerebellum
    echo -e "\e[32mCerebellum parcellation!"
    for mask_option in "${mask_options[@]}"; do
        mri_vol2vol --mov $SUBJECTS_DIR/$subject_id/mri/norm.mgz \
                    --s $subject_id \
                    --targ $SUBJECTS_DIR/MNI152/buckner_atlas/17_network_${mask_option}_mask.nii.gz \
                    --m3z $SUBJECTS_DIR/MNI152/mri/transforms/talairach.m3z \
                    --noDefM3zPath \
                    --o $SUBJECTS_DIR/$subject_id/mri/cerebellum_17_network_${mask_option}_mask.nii.gz \
                    --inv-morph \
                    --interp nearest
        echo "$(date): Cerebellum parcellation with $mask_option mask completed for subject $subject_id" >> "$log_file"
    done

    # schaefer atlas
    echo -e "\e[32mSchaefer2018 parcellation in individual surface space!"
    for n in 100 200 400 800 1000; do
        for hemi in "${hemis[@]}"; do
            mris_ca_label -l $SUBJECTS_DIR/$subject_id/label/${hemi}.cortex.label \
                                $subject_id \
                                ${hemi} \
                                $SUBJECTS_DIR/$subject_id/surf/${hemi}.sphere.reg \
                                $sch_gcs_dir/${hemi}.Schaefer2018_${n}Parcels_7Networks.gcs \
                                $SUBJECTS_DIR/$subject_id/label/${hemi}.Schaefer2018_${n}Parcels_7Networks_order.annot
            echo "$(date): Schaefer parcellation with $n labels in $hemi completed for subject $subject_id" >> "$log_file"
        done

        # bring schaefer from surface to volume
        mri_aparc2aseg --s $subject_id \
                        --annot Schaefer2018_${n}Parcels_7Networks_order \
                        --o $SUBJECTS_DIR/$subject_id/mri/Schaefer2018_${n}_7Networks.mgz

        echo "$(date): Schaefer parcellation with $n labels brought to volume space for subject $subject_id" >> "$log_file"
    done
    
    ## Glasser atlas
    mri_surf2surf \
                    --srcsubject fsaverage \
                    --trgsubject $subject_id \
                    --hemi lh \
                    --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.HCPMMP1.annot \
                    --tval $SUBJECTS_DIR/$subject_id/label/lh.HCPMMP1.annot

    mri_surf2surf \
                    --srcsubject fsaverage \
                    --trgsubject $subject_id \
                    --hemi rh \
                    --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.HCPMMP1.annot \
                    --tval $SUBJECTS_DIR/$subject_id/label/rh.HCPMMP1.annot

    mri_aparc2aseg \
                --s $subject_id \
                --annot HCPMMP1 \
                --o $SUBJECTS_DIR/$subject_id/mri/HCPMMP1.mgz
    



}

## main part
for subj_dir in "$SUBJECTS_DIR"/*; do
    if [ -d "$subj_dir" ]; then
        subject_id=$(basename "$subj_dir")
        if [ "$subject_id" != "fsaverage" ]; then
            subsegment "$subject_id"
        fi
    fi
done







#map the annotation files of the HCP MM1.0 atlas from fsaverage to the subject (first left, then right hemisphere)
      mri_surf2surf --srcsubject fsaverage --trgsubject $subjectID --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.hcpmmp1.annot --tval $subjectDIR/label/lh.hcpmmp1.annot
      mri_surf2surf --srcsubject fsaverage --trgsubject $subjectID --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.hcpmmp1.annot --tval $subjectDIR/label/rh.hcpmmp1.annot
      
      #Map the HCP MMP 1.0 annotations onto the volumetric image and add (Freesurfer-specific) subcortical segmentation.
      #Convert the resulting file to .mif format (use datatype uint32 --> liked best by mrtrix)
      mri_aparc2aseg --old-ribbon --s $subjectID --annot hcpmmp1 --o $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1.mgz"
      mrconvert --datatype uint32 $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1.mgz" $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1.mif"
      
      #Replace the random integers of the hcpmmp1.mif file with integers that start at 1 and increase by 1
      labelconvert $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1.mif" /home/julia/mrtrix3/share/mrtrix3/labelconvert/hcpmmp1_original.txt /home/julia/mrtrix3/share/mrtrix3/labelconvert/hcpmmp1_ordered.txt $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_nocoreg.mif"

      Register the ordered atlas-based volumetric parcellation to diffusion space
      first calculate transformation
     
      flirt -in $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_dir-APPA_meanB0brain.nii.gz" -ref $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/anat/"$subjectID"_acq-mprage_T1wBrain.nii.gz" -dof 6 -omat $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_diff2struct_fsl.mat"
      transformconvert -force $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_diff2struct_fsl.mat" $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_dir-APPA_meanB0brain.nii.gz" $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/anat/"$subjectID"_acq-mprage_T1wBrain.nii.gz" flirt_import $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_diff2struct_mrtrix.txt"
      mrtransform -force $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_nocoreg.mif" -linear $derivatives_path"/dwi_oct2020/"$sub"/"$ses"/dwi/"$subjectID"_diff2struct_mrtrix.txt" -inverse -datatype uint32 $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_coreg.mif"
      mrstats $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_nocoreg.mif"
      mrstats $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_coreg.mif"

      #connectomics
      connectome_filename=$derivatives_path"/dwi/"$sub"/"$ses"/dwi/"$subjectID"_connectome_Glasser.csv"
      tck2connectome $derivatives_path"/dwi/"$sub"/"$ses"/dwi/"$subjectID"_iFOD2.tck" $derivatives_path"/dwi/"$sub"/"$ses"/anat/"$subjectID"_hcpmmp1_parcels_coreg.mif" $connectome_filename --assignment_radial_search 3.2 -zero_diagonal --force