#!/bin/bash

# T1 + T2 anatomical MRI Analysis Pipeline
# Written by Payam S. Shabestari, Zurich, 01.2025
# email: payam.sadeghishabestari@uzh.ch
# This script is written mainly for Antinomics project. However It could be used for other purposes.


## echo with color
## chmod u+x sMRI_processing.sh
## ./sMRI_processing.sh
## echo sth >> log.txt
## copy these files to work with linux : license, buckner atlas, choi atlas, tracula config file

set -e
display_usage() {
	echo "$(basename $0) [raw_t1] [raw_t2] [raw_dwi] [subject_id] [saving_dir] [recon_all]"
	echo "This script uses Freesurfer for cortical and subcortical segmentation
            as well as extracting probabilistic white matter tracts:
			1) The structural T1 image (.rec / .nii format);
			2) The structural T2 image (.rec / .nii format);
			3) The diffusion image (.rec / .nii format);
			4) The subject ID number;
			5) Path to a directory to save anatomical data of the subject. If not provided, the default directory will be used;
            6) If false, the cortical segmentation is already done. By default, recon-all function from FS will be run." 
	}

if [[ "$1" == "--h" || $# -lt 5 ]]; then
	display_usage
	exit 1
fi

raw_t1=$1
raw_t2=$2
raw_dwi=$3
subject_id=$4
saving_dir=$5
recon_all=${6:-false}

## set Paths
spath="../subjects/${subject_id}/sMRI"
gcs_dir="./gcs"
saving_dir="${5:-$spath}"
export FREESURFER_HOME=/Applications/freesurfer/dev
export SUBJECTS_DIR=$saving_dir
source $FREESURFER_HOME/SetUpFreeSurfer.sh

## cortical and subcortical segmentation
if [[ "$recon_all" == "true" ]]; then
	recon-all -s $subject_id -i $raw_t1
	recon-all -all -subjid $subject_id  
fi

## subsegmentation (hippocampus + amygdala, brainstem)
segmentHA_T2.sh $subject_id raw_t2 T1_T2 1 
segmentBS.sh $subject_id # -> fixed

## second round (AAN, hypothalamus, thalamic nuclei, cerebellum)
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8
segmentAAN.sh $subject_id # fix later (not found) -> in the development built -> working with 2019 -> looks like fixed -> so local payam


###
mri_segment_hypothalamic_subunits $subject_id # fix later (illegal hardware instruction) -> probably linux will fix.





#####
segmentThalamicNuclei_DTI.sh -s $subject_id # fix later (not found) -> in the development built -> changed to matlab R2014b -> cant

'''
Unrecognized function or variable 'atlasInitialisationFiles'.

Error in TS_fnc_thalamus_seg_gem_joint (line 777)

Error in TS_fnc_jointSegmentationWrapper (line 66)


looks like its under development so not yet ...
'''

mri_vol2vol --mov 0002/mri/norm.mgz --o 0002/mri/norm.dwispace.mgz --lta 0002/dmri/xfms/anatorig2diff.bbr.lta  --no-resample --targ 0002/dmri/dtifit_FA.nii.gz
mri_vol2vol --mov 0002/mri/aseg.mgz --o 0002/mri/aseg.dwispace.mgz --lta 0002/dmri/xfms/anatorig2diff.bbr.lta  --no-resample --targ 0002/dmri/dtifit_FA.nii.gz
mri_segment_thalamic_nuclei_dti_cnn --t1 0002/mri/norm.dwispace.mgz --aseg 0002/mri/aseg.dwispace.mgz --fa 0002/dmri/dtifit_FA.nii.gz --v1 0002/dmri/dtifit_V1.nii.gz --o 0002/mri/thalamic_dti.nii.gz --vol 0002/tables/thalamic_dti_volumes.csv

'''
so now error is again illegal hardware instruction which is probably due to TF. -> probably linux will fix

'''




#### cerebellum parcelation

# Step 1: Upsample Buckner atlas from 2mm to 1mm resolution (matching MNI152 1mm template)
mri_vol2vol --mov Buckner_atlas.nii.gz \
            --targ MNI152/mri/norm.mgz \
            --regheader \
            --o Buckner_atlas1mm.nii.gz \
            --no-save-reg \
            --interp nearest

# Step 2: Warp the upsampled atlas to FreeSurfer's nonlinear volumetric space
mri_vol2vol_used --mov Buckner_atlas1mm.nii.gz \
                 --s MNI152_FS \
                 --targ $FREESURFER_HOME/average/mni305.cor.mgz \
                 --m3z talairach.m3z \
                 --o Buckner_atlas_freesurfer_internal_space.nii.gz \
                 --interp nearest

# Step 3: Warp the atlas from FreeSurfer internal space to subject’s native space
mri_vol2vol --mov $SUBJECTS_DIR/SUBJECT_FS/mri/norm.mgz \
            --s SUBJECT_FS \
            --targ Buckner_atlas_freesurfer_internal_space.nii.gz \
            --m3ztalairach.m3z \
            --o Buckner_atlas_subject.nii.gz \
            --interp nearest \
            --inv-morph

# freeview -v ${SUBJECTS_DIR}/${SUBJECT_ID}/mri/orig.mgz \
#            ${OUTPUT_DIR}/Buckner_atlas_subject.nii.gz:colormap=lut


#### Striatal Parcellation
# same as cerebellum

# freeview -v FSL_MNI152_FreeSurferConformed_1mm.nii.gz Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz:colormap=lut:lut=Choi2012_7Networks_ColorLUT.txt

# freeview -v FSL_MNI152_FreeSurferConformed_1mm.nii.gz Choi2012_7Networks_MNI152_FreeSurferConformed1mm_TightMask.nii.gz:colormap=lut:lut=Choi2012_7Networks_ColorLUT.txt Choi2012_7NetworksConfidence_MNI152_FreeSurferConformed1mm_TightMask.nii.gz:colormap=heat:heatscale=0,0.5,1



## schafer atlas (we might need it for fMRI or dMRI)
hemis=("lh" "rh")
for hemi in "${hemis[@]}"; do
	for n in 100 200 300 400; do
		for net_option in 7 17; do
			mris_ca_label -l $SUBJECTS_DIR/$subject_id/label/${hemi}.cortex.label \
			$subject_id ${hemi} $SUBJECTS_DIR/$subject_id/surf/${hemi}.sphere.reg \
			$gcs_dir/${hemi}.Schaefer2018_${n}Parcels_${net_option}Networks.gcs \
			$SUBJECTS_DIR/$subject_id/label/${hemi}.Schaefer2018_${n}Parcels_${net_option}Networks_order.annot
		done
	done
done

## extracting probabilistic white matter tracts



## bem watershed

mne watershed_bem -s 0002 -d /Applications/freesurfer/dev/subjects

##### tables

## extract DK atlas information
measures=("area" "volume" "thickness")
parcels=("aparc" "aparc.a2009s")
for hemi in "${hemis[@]}"; do
    for meas in "${measures[@]}"; do
        for parc in "${parcels[@]}"; do
            # Construct the table file name
            if [[ $parc == "aparc" ]]; then
                tablefile="tables/${parc}_${meas}_${hemi}.txt"
            else
                tablefile="tables/${parc}_${meas}_${hemi}.txt"
            fi
            aparcstats2table --subjects "$subject_id" --hemi "$hemi" --meas "$meas" --parc="$parc" --tablefile "$tablefile"
        done
    done
done

# extract segmentation information
asegstats2table --subjects $subject_id --common-segs --meas volume --stats=aseg.stats --table=tables/segstats.txt
asegstats2table --subjects $subject_id --statsfile=brainstem.v13.stats --tablefile=tables/brainstem.txt

for hemi in "${hemispheres[@]}"; do
	statsfile="thalamic-nuclei.${hemi}.v13.T1.stats"
	tablefile="tables/thalamic-nuclei_${hemi}.txt"
	asegstats2table --subjects $SUB_ID --statsfile="$statsfile" --tablefile="$tablefile"
	
	statsfile="hipposubfields.${hemi}.T2.v22.T2.stats"
	tablefile="tables/hipposubfields_${hemi}.txt"
	asegstats2table --subjects $SUB_ID --statsfile=hipposubfields.lh.T2.v22.T2.stats --tablefile=tables/hipposubfields_lh.txt 
	
	statsfile="amygdalar-nuclei.${hemi}.T2.v22.T2.stats"
	tablefile="tables/amygdalar_${hemi}.txt"
	asegstats2table --subjects $SUB_ID --statsfile=amygdalar-nuclei.lh.T2.v22.T2.stats --tablefile=tables/amygdalar_lh.txt 
done










## for report 

## AAN
freeview -v $SUBJECTS_DIR/0002/mri/T1.mgz -v  $SUBJECTS_DIR/0002/mri/arousalNetworkLabels.v10.mgz:colormap=lut:lut=$FREESURFER_HOME/average/AAN/atlas/freeview.lut.txt

