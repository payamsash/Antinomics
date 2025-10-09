# Antinomics

## Subcortical Connectivity in Tinnitus

This repository provides the **analysis pipeline** and documentation for our study:

> **“Multimodal Characterisation of Subcortical Networks in Tinnitus”**

The project investigates **structural** and **functional connectivity** in tinnitus by comparing the **cortical and subcortical organization** of tinnitus patients and healthy controls.

We combine **structural MRI (sMRI)**, **diffusion MRI (dMRI)**, and **resting-state fMRI** to characterise both the brain’s structural wiring and functional networks.

The pipeline supports preprocessing of structural, functional, and diffusion MRI data by integrating established tools such as *FreeSurfer*, *FSL*, *ANTs*, and *MRtrix*. It also facilitates atlas handling—projecting atlases like Schaefer, Glasser, and Tian subcortical atlases into each subject’s native space—and provides a complete workflow for cortical and subcortical tractography, enabling the construction of structural connectomes across both cortical and subcortical regions.

---

## Pipeline Overview

The processing scripts are numbered to indicate their execution order. Each script name follows the pattern **`ant_{i}_*`**, where `{i}` is the step number.

| Step | Script      | Description |
|------|------------|------------|
| 0    | `ant_00_*` | Creates the project directory structure and organizes raw **sMRI**, **fMRI**, and **dMRI** data acquired on the Philips MRI scanner at University Hospital Zurich. |
| 1    | `ant_01_*` | Runs **FreeSurfer `recon-all`** in parallel on the structural MRI data to produce the `SUBJECTS_DIR` for FreeSurfer. |
| 2    | `ant_02_*` | Uses **FreeSurfer Bayesian segmentation tools** to subsegment the **hippocampus**, **amygdala**, **thalamic nuclei** **ascending arousal network (AAN)**, and register the **Schaefer** and **Glasser** atlas. to both surface and volume. |
| 3    | `ant_03_*` | Performs **tractography** at both population and individual levels to generate **global** and **subcortical structural connectivity (SC)** matrices. |
| 4    | `ant_04_*` | Jointly segments **thalamic nuclei** using combined **T1** and **dMRI** data. |
| 5    | `ant_05_*` | Applies **Pierrick Coupé’s non-local means algorithm** to denoise structural MRI data prior to downstream fMRI analysis. |
| 6    | `ant_06_*` | Handles **fMRI preprocessing** and denoising for both sessions. |

---

<p align="center">
  <img src="images/picture_1.svg" alt="Pipeline illustration" width="800">
</p>

*For questions or collaboration inquiries, please contact payam.*
