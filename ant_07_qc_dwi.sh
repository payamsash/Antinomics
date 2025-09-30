## check tissue type segmentation
mrview raw_anat.mif -overlay.load 5tt_fixed.mif -size 1000,800 -mode 2

## Check coregistration
mrview mean_b0.mif -overlay.load gmwmSeed_coreg.mif -size 1000,800 -mode 2

## check tian atlas coregistration in T1 space
mrview raw_anat.mif -overlay.load tian_subspace.nii.gz -size 1000,800 -mode 2

## check tian atlas coregistration in DWI space
mrview mean_b0.mif -overlay.load Tian_subcortical_dwi_resampled.mif -size 1000,800 -mode 2

## check subcortical seeding
mrview mean_b0.mif -overlay.load subcortical_gmwmi_raw.mif -size 1000,800 -mode 2

## exclusion = brainmask - hull
mrview mean_b0.mif -overlay.load hull_exclude.nii.gz -size 1000,800 -mode 2

## sum of masks
mrview mean_b0.mif -overlay.load tmp_excl1.mif -size 1000,800 -mode 2