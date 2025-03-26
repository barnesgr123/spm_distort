
function coreg_fname=coreg_subject_distort(subject, fname,surf_name,coreg_fname,HeadModel,datadirectory);
%function coreg_fname=coreg_subject_distort(subject, fname,surf_name,coreg_fname,HeadModel,datadirectory);
%% GRB. 26/03/2025. basic copy and coreg to fit in with distort directory structure

%% subject- one item from structure - create_anon_subject_structure 
%% fname  - MEG dataset
%% surf_name - full path to cortical surface mesh
%% coreg_fname - full path to coregistered copy of fname
%% HeadModel - SPM Headmodel choice e.g. 'Single Shell'

mri_dir=fullfile(datadirectory,'mri'); %% mri directory
subject_id=subject.subj_id;
spmname = fullfile(datadirectory,'meg',subject_id,fname); %% meg dataset



%% now load in dataset and coregister
D     = spm_eeg_load(spmname);

Mtest=gifti(surf_name); %% just to check that surf_name exists

if exist(coreg_fname),
    fprintf('\n Deleting existing file')
    delete(coreg_fname)
end;
% Coregister to mesh
clear jobs
matlabbatch={};
% Copy datafile
matlabbatch{1}.spm.meeg.other.copy.D = {spmname};
matlabbatch{1}.spm.meeg.other.copy.outfile = coreg_fname;
spm_jobman('run', matlabbatch);


spm_get_defaults('cmdline',1); %% need this bit to keep things running in unix, not sure why
spm('CmdLine')
spm_jobman('initcfg');

matlabbatch={};
% Coregister dataset to reconstruction mesh
matlabbatch{1}.spm.meeg.source.headmodel.D = {coreg_fname};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {fullfile(mri_dir, subject.subj_id, [subject.headcast_t1 ',1'])};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {surf_name};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = subject.nas;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = subject.lpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = subject.rpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = HeadModel;
spm_jobman('run', matlabbatch);

D1=spm_eeg_load(coreg_fname);
M2=gifti(D1.inv{D1.val}.mesh.tess_ctx);
%% sometimes if selected cortex not found/loaded then uses default
if length(M2.vertices)~=length(Mtest.vertices),
    error('Mesh not correctly coregistered')
end;