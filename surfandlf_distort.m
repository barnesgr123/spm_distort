
function [aM,aL]=surfandlf_distort(spmfilename,RandSeed,HeadModel,Npoints,Zmax,DistortIndices,RunClean);
%function [aM,aL]=surfandlf_distort(spmfilename,RandSeed,HeadModel,Npoints,Zmax,DistortIndices,RunClean);
%% spmfilename - full path to dataset
%% random seed to use to generate trajectory
%% HeadModel to use eg 'Single Shell'
%% Npoints - number of points on trajectory (defaults to 17 as in paper)
%% Zmax - maximum deviation (+Zmax to -Zmax over Npoints) as a z score 
%% DistortIndices - which surfaces will be distorted - eg [8:100] in paper.
%% RunClean - if 'Yes' all volumes and lead fields recomputed 


if nargin<7,
    RunClean='Yes';
end;

if isempty(Npoints)
    Npoints=17;
end;


spm('defaults','eeg') % or 'eeg' or 'pet'
spm_jobman('initcfg');
spm_get_defaults('cmdline',1);
spm('Cmdline')
PLOTSTUFF=0; %% no graphics 
HeadModel= replace(HeadModel,'_',' ');

if isstr(RandSeed), %% batch call can be string in unix
    RandSeed=str2num(RandSeed);
end;



Zlimit=Zmax+3+Zmax/6; %% arbitrary limit to stop components going too far from normal range
fprintf('\n Using a threshold of %3.2f for Zmax=%3.2f\n ',Zlimit,Zmax);

LFRedo=RunClean;
WriteClean=RunClean;
D=spm_eeg_load(spmfilename);


PCAtemplate=fullfile(spm('Dir'),'tpm','shp','Template_0.nii');

Zrange=[-Zmax,Zmax];
fprintf('\n Running surf and lf batch SHP no coreg for %s',spmfilename)

matlabbatch=[];
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.D = {spmfilename};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.val = {D.val};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.PCAtemplate = {PCAtemplate};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.TemplateRedo = {'No'};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.DistortIndices = {DistortIndices};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Npoints = {Npoints}; %% only 4 points on trajectory
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zrange = {Zrange}; %% move from minimal to maximal distortion
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zlimit = {Zlimit};
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.Zcentre = {'sub'}; %% keep it subject specific
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.RandSeed = {RandSeed}; %% random seed defines trajectory
matlabbatch{1}.spm.meeg.source.eeg_shp_distort.WriteClean = {WriteClean}; %% force delete of any existing seed directory
[aM,b]=spm_jobman('run', matlabbatch);

%% now compute distances of brians to brain.
M0=gifti(aM{1}.brain); %% get original cortex
for i=1:Npoints,
    M1=gifti(aM{1}.brians{i});
    d0=M1.vertices-M0.vertices;
    d1(i)=mean(sqrt(dot(d0',d0'))); %% mean vertex-vertex distance to M0
end;% for i

if PLOTSTUFF,
    figure;
    plot(linspace(Zrange(1),Zrange(end),Npoints),d1)
    ylabel('Distance to true cortex (mm)')
    xlabel('Trajectory (in Z scores)')
end;
%% now make lead fields for each surface (looks in directory seedXXX as specified by RandomSeed above)

LFsubdir=sprintf('seed%03d',RandSeed);
matlabbatch=[];
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.D = {spmfilename};
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.val = 1;
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFheadmodel = HeadModel;
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFRedo = LFRedo;
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.LFsubdir = LFsubdir;
matlabbatch{1}.spm.meeg.source.eeg_shp_gainmat.WriteClean = WriteClean;
[aL,b]=spm_jobman('run', matlabbatch);

%% now quantify how much the lead fields for the distorted surfaces differ from the original
%% first LF is the original
G=load(aL{1}.gainmat{1}); %% original lead fields
L0=G.G;
diffL=[];
for i=1:Npoints,
    G=load(aL{1}.gainmat{i+1});
    L1=G.G;
    diffL(i)=mean(dot(L1-L0,L1-L0));
end; % for i

if PLOTSTUFF,
    figure;
    plot(linspace(Zrange(1),Zrange(end),Npoints),diffL)
    ylabel('Distance to true LF (arb units)')
    xlabel('Trajectory (in Z scores)')
end;






