% fMRI denoising pipeline
% Written by Payam S. Shabestari, Zurich, 01.2025
% email: payam.sadeghishabestari@uzh.ch
% This script is written mainly for Antinomics project. However It could be used for other purposes.

%% --- Paths
fmri_dir = "/Volumes/Extreme_SSD/payam_data/antinomics_data/subjects_fsl_dir";
smri_dir = "/Volumes/Extreme_SSD/payam_data/antinomics_data/t1_denoised";
cwd      = pwd;  % folder for Antinomics.mat

nsub = 37;
nses = 2;
ncond = 1;  % only "rest"
condition_name = 'rest';

subjects = dir(fmri_dir);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name},'.'));
if numel(subjects) ~= nsub
    error('Expected %d subjects, found %d', nsub, numel(subjects));
end

clear batch;
batch.Setup.functionals = repmat({{}}, [nsub, 1]);
structurals = cell(1, nsub);
gm_masks    = cell(1,nsub);
wm_masks    = cell(1,nsub);
csf_masks   = cell(1,nsub);
brain_masks = cell(1,nsub);
motion_files = cell(nsub, nses);

for i = 1:nsub
    sid = subjects(i).name;
    
    % Functional sessions
    batch.Setup.functionals{i}{1} = fullfile(fmri_dir,sid,'t2_s1_norm.nii.gz');
    batch.Setup.functionals{i}{2} = fullfile(fmri_dir,sid,'t2_s2_norm.nii.gz');
    
    % Structural
    structurals{i} = fullfile(smri_dir,[sid '_denoised.nii']);
    
    % Masks
    gm_masks{i}    = fullfile(fmri_dir,sid,'gm_in_func.nii.gz');
    wm_masks{i}    = fullfile(fmri_dir,sid,'wm_in_func.nii.gz');
    csf_masks{i}   = fullfile(fmri_dir,sid,'csf_in_func.nii.gz');
    brain_masks{i} = fullfile(fmri_dir,sid,'brain_mask_in_func.nii.gz');
    
    % Motion files
    motion_files{i,1} = fullfile(fmri_dir,sid,'t2_s1_mc.par');
    motion_files{i,2} = fullfile(fmri_dir,sid,'t2_s2_mc.par');
end


batch.filename = fullfile(cwd,'Antinomics.mat');
batch.Setup.nsubjects       = nsub;
batch.Setup.RT              = 2.5;
batch.Setup.acquisitiontype = 1;
batch.Setup.structurals     = structurals;

batch.Setup.masks.Grey  = gm_masks;
batch.Setup.masks.White = wm_masks;
batch.Setup.masks.CSF   = csf_masks;
batch.Setup.masks.Brain = brain_masks;

batch.Setup.done = 1;

%% Denoising
batch.Denoising.confounds.names      = { ...
    'White Matter', ...
    'CSF', ...
    'Grey Matter', ...   % whole-brain signal
    'motion' ...
    };
batch.Denoising.confounds.dimensions = {1, 1, 1, 6};
batch.Denoising.confounds.deriv      = {1, 1, 1, 1};  
batch.Denoising.filter = [0.01 0.15];
batch.Denoising.detrending = 1;   
batch.Denoising.despiking   = 0;
batch.Denoising.done        = 1;
batch.Denoising.overwrite   = 'Yes';

%% Run all analyses
conn_batch(batch);


