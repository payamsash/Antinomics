mrview mean_b0.mif \
        -tractography.geometry points \
        -mode 2 \
        -size 1000,800 \
        -noannotations \
        -comments false \
        -voxelinfo false \
        -colourbar false \
        -tractography.load tracks_subcortical_10M.tck \
        -capture.prefix raw_diffusion_subcortical_tck \
        -capture.grab -exit


transformconvert diff2struct_fsl.mat mean_b0.nii.gz raw_anat.mif flirt_import diff2struct_mrtrix.txt


tcktransform tracks_10M.tck \
            tracks_10M_T1.tck \
            -matrix diff2struct_mrtrix.txt