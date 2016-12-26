%% This script automatically reads in the SPM.mat file with a certain image mask
% and display the results in the Graphics window. 

matlabbatch{1}.spm.stats.results.spmmat = cellstr([SPMDir 'SPM.mat']);
matlabbatch{1}.spm.stats.results.conspec(1).titlestr = current_mask;
matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec(1).thresh = treshold;
matlabbatch{1}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{1}.spm.stats.results.conspec(1).mask.contrasts = cellstr([MaskDir current_mask]);
matlabbatch{1}.spm.stats.results.conspec.mask.thresh = {};
matlabbatch{1}.spm.stats.results.conspec.mask.mtype = 0;
matlabbatch{1}.spm.stats.results.units = 1; % Data type: Volumetric, 1; Scalp-Time, 2; ...
matlabbatch{1}.spm.stats.results.print = false; % true or false

%% Run job file
spm_jobman('run', matlabbatch);