#!/bin/bash
set -e

# Base directory containing all subjects
base_dir="/Volumes/G_USZ_ORL$/Research/ANTINOMICS/payam/subjects_mrtrix_dir"

# Loop through each subject folder
for subject_dir in "$base_dir"/*/; do
  subject=$(basename "$subject_dir")
  # if [[ "$subject" != "vfav" ]]; then
  #   continue
  # fi
  subject_mrtrix_dir="$base_dir/$subject"

  echo ">>> Processing subject: $subject"

  mkdir -p "$subject_mrtrix_dir/report"
  echo ">>> Generating QC snapshots..."

  capture_scene () {
      local base=$1
      local main=$2
      local overlay=$3
      local mode=$4

      for plane in 0 1 2; do
          mrview "$main" \
          -overlay.load "$overlay" \
          -mode "$mode" \
          -plane "$plane" \
          -size 1000,800 \
          -noannotations \
          -comments false \
          -voxelinfo false \
          -colourbar false \
          -capture.folder "$subject_mrtrix_dir/report" \
          -capture.prefix "${base}_plane${plane}" \
          -capture.grab -exit
      done
  }

  cd "$subject_mrtrix_dir"

  # for i in 0 1 2 3 4; do
  #     mrconvert 5tt_nocoreg.mif -coord 3 $i vol${i}.mif -force
  # done
  # mrcalc vol4.mif 0 -mul vol4_zero.mif -force

  # # # ----------------------
  # # # Run each QC step
  # # # ----------------------
  # capture_scene "gm_ribbon" raw_anat.mif vol0.mif 4
  # capture_scene "subcortical_gm" raw_anat.mif vol1.mif 4
  # capture_scene "wm" raw_anat.mif vol2.mif 4
  # capture_scene "csf" raw_anat.mif vol3.mif 4
  # rm vol0.mif vol1.mif vol2.mif vol3.mif vol4.mif vol4_zero.mif

  capture_scene "gmwmi_coreg" mean_b0.nii.gz gmwmSeed_coreg.mif 4
  # capture_scene "tian_t1" raw_anat.mif Tian_subcortical.mif 4
  # capture_scene "tian_dwi" mean_b0.nii.gz Tian_subcortical_dwi_resampled.mif 4
  # capture_scene "subcortical_seed" mean_b0.nii.gz subcortical_gmwmi.mif 4
  # capture_scene "exclusion_hull" mean_b0.nii.gz hull_exclude_reg.mif 4
  # capture_scene "exclusion_ribbon" mean_b0.nii.gz cortical_ribbon_reg.mif 4

  

  echo ">>> Writing HTML report..."

  cat > report/qc_report.html <<EOF

<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>Tractography QC Report - $subject</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    h2 { margin-top: 40px; }
    .imgrow { display: flex; gap: 10px; }
    .imgrow img { width: 32%; border: 1px solid #ccc; border-radius: 6px; }
  </style>
</head>
<body>

<h1>Tractography QC Report - $subject</h1>

EOF

  # append all report sections (reuse same HTML body as before)
  cat <<'HTML' >> report/qc_report.html

<h2>1. Tissue Type Segmentation (all tissues)</h2>

<h2>1a. GM Ribbon</h2>
<div class="imgrow">
  <img src="gm_ribbon_plane00000.png">
  <img src="gm_ribbon_plane10000.png">
  <img src="gm_ribbon_plane20000.png">
</div>

<h2>1b. Subcortical GM</h2>
<div class="imgrow">
  <img src="subcortical_gm_plane00000.png">
  <img src="subcortical_gm_plane10000.png">
  <img src="subcortical_gm_plane20000.png">
</div>

<h2>1c. White Matter</h2>
<div class="imgrow">
  <img src="wm_plane00000.png">
  <img src="wm_plane10000.png">
  <img src="wm_plane20000.png">
</div>

<h2>1d. CSF</h2>
<div class="imgrow">
  <img src="csf_plane00000.png">
  <img src="csf_plane10000.png">
  <img src="csf_plane20000.png">
</div>

<h2>2. Coregistration (GMWMI)</h2>
<div class="imgrow">
  <img src="gmwmi_coreg_plane00000.png">
  <img src="gmwmi_coreg_plane10000.png">
  <img src="gmwmi_coreg_plane20000.png">
</div>

<h2>3. Tian Atlas in T1 Space</h2>
<div class="imgrow">
  <img src="tian_t1_plane00000.png">
  <img src="tian_t1_plane10000.png">
  <img src="tian_t1_plane20000.png">
</div>

<h2>4. Tian Atlas in DWI Space</h2>
<div class="imgrow">
  <img src="tian_dwi_plane00000.png">
  <img src="tian_dwi_plane10000.png">
  <img src="tian_dwi_plane20000.png">
</div>

<h2>5. Subcortical Seeding</h2>
<div class="imgrow">
  <img src="subcortical_seed_plane00000.png">
  <img src="subcortical_seed_plane10000.png">
  <img src="subcortical_seed_plane20000.png">
</div>

<h2>6. Exclusion Mask (Hull)</h2>
<div class="imgrow">
  <img src="exclusion_hull_plane00000.png">
  <img src="exclusion_hull_plane10000.png">
  <img src="exclusion_hull_plane20000.png">
</div>

<h2>7. Exclusion Mask (Ribbon)</h2>
<div class="imgrow">
  <img src="exclusion_ribbon_plane00000.png">
  <img src="exclusion_ribbon_plane10000.png">
  <img src="exclusion_ribbon_plane20000.png">
</div>

</body>
</html>
HTML

  echo ">>> QC report ready at: $subject_mrtrix_dir/report/qc_report.html"
  echo "-------------------------------------------------------------"
done