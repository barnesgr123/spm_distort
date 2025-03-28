%% DEMO SCRIPT TO RUN LOPEZ ET AL. 2025 distortion paper
%% GRB March 2025

%% reproduces first trajectory from figure 3 in paper

clear all;
close all;
PLOTSTUFF=0;

addpath('D:\spm'); %% spm directory
addpath('D:\spm_distort'); %% directory with distort code from github
rootdir='C:\Users\gbarnes\Documents\jimmyupload\'; %% where the downloaded data sits

spm('defaults','eeg');
spm_jobman('initcfg');

%% set up the directory structure
datadirectory=fullfile(rootdir,'data');
mri_dir=fullfile(datadirectory,'mri');

outputdirectory=fullfile(rootdir,'output')
mkdir(outputdirectory);

%% where the PCA template for the diffeomorphic code sits
PCAtemplate=fullfile(spm('Dir'),'tpm','shp','Template_0.nii');

%% load in information on the 8 subjects
subjects=create_anon_subject_structure();


%%%% FILES OF AVERAGED MEG DATA AND Time-Frequency WINDOWS
allfname=strvcat('mpinstr_rcinstr_Tafdf.mat','mpdots_rcinstr_Tafdf.mat','mprcresp_Tafdf.mat');
allwoi=[ 0 500; -2500 -2000;-200 300]; %% approx where evoked responses sit in these files
foi=[5 90]; %% frequency window in Hz


%%%%%%%%%%%%%%%%%%%%% DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%%
subind=1; %% subject data to look at (1-8)
fileind=1; %% file to look at (1-3)


Npoints=17; %% points on trajectory (17 in paper)
RandSeeds=[1:8]; %% seed which sets the trajectory (runs from 1 to 8 in paper)
RunClean='Yes'; %% if 'No' already have all surfaces worked out

%% distortion parameters
Zmax=3; %% amount of distortion per component (as Z score)
DistortIndices=[8:100]; %% which PCAs to distort (as in paper)
HeadModel='Single Shell'; %% head model to use
loc_surface='white'; %% cortical surface to use white (in paper) or pial


%% inversion parameters for MEG Data
woi=allwoi(fileind,:); %% time window in ms
patch_size=0.6; % default smoothness of reconstruction
n_temp_modes=16; %% number of temporal modes to use
invmethods={'EBB','IID','GS'}; %% SPM inversion algorithms
Nblocks=1; pctest=0; %% not testing cross validation
%% for cross validation: Nblocks=10; pctest=10; for example


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%

allF_vals=[];
allR2=[];
allcrosserr=[];
allfiles=[];


%% loop through the list of random seeds- each seed defines a trajectory
RandSeed=RandSeeds(1);

for subind=1:8,
    %%  load in average evoked dataset
    fname = fullfile(datadirectory,'meg',subjects(subind).subj_id,deblank(allfname(fileind,:)));
    D     = spm_eeg_load(fname);

    %% Decide on the cortical surface to use for this subject (freesurfer/laMEG output)
    subj_surf_dir=fullfile(datadirectory, 'surf',sprintf('%s-synth',...
        subjects(subind).subj_id),'surf');
    %%%% FROM COREG SUBJ ERF
    fprintf('\n Using original cortex %s \n',sprintf('%s.ds.gii', loc_surface))

    tmp_loc_surface=gifti(fullfile(subj_surf_dir, sprintf('%s.ds.gii', loc_surface)));
    %% make new surface without the surface normals (just faces and vertices)..
    ds_loc_surface=[]; ds_loc_surface.vertices=tmp_loc_surface.vertices;
    ds_loc_surface.faces=tmp_loc_surface.faces; ds_loc_surface.mat=tmp_loc_surface.mat;
    fprintf('Removing orienation info from original surface (if it exists)');
    surf_name=sprintf('ori_norm_%s_loc.gii', loc_surface);
    dosstr=sprintf('copy %s %s',fullfile(subj_surf_dir, sprintf('%s.ds.gii', loc_surface)),fullfile(subj_surf_dir(1:end-5),sprintf('%s.ds.gii', loc_surface)))
    dos(dosstr);
    dosstr=sprintf('copy %s %s',fullfile(subj_surf_dir, sprintf('%s.ds.gii', 'pial')),fullfile(subj_surf_dir(1:end-5),sprintf('%s.ds.gii', 'pial')))
    dos(dosstr);
    dosstr=sprintf('cd %s', subj_surf_dir);
    eval(dosstr)
    dosstr=sprintf('delete *.*');
    eval(dosstr)
    dosstr=sprintf('move %s %s',fullfile(subj_surf_dir(1:end-5),sprintf('%s.ds.gii', loc_surface)),fullfile(subj_surf_dir, sprintf('%s.ds.gii', loc_surface)))
    dos(dosstr)
    dosstr=sprintf('move %s %s',fullfile(subj_surf_dir(1:end-5),sprintf('%s.ds.gii', 'pial')),fullfile(subj_surf_dir, sprintf('%s.ds.gii', 'pial')))
    dos(dosstr)
end; % for subind
