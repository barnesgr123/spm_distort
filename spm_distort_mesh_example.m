%% DEMO SCRIPT TO RUN LOPEZ ET AL. REVISION December 2025 distortion paper
%% GRB December 2025

%% Run code for subject 1

clear all;
close all;

DISTORTANAT=0; %% distort anatomies


addpath('D:\spm'); %% spm directory
addpath('D:\spm_distort'); %% directory with distort code from github
rootdir='C:\Users\gbarnes\Documents\jimmytest\'; %% where the downloaded data sits

spm('defaults','eeg');
spm_jobman('initcfg');

[v,r]=spm('ver'); %% Code is currently in the 2025 development version of spm (r=0.0); but not the 2025 'release' version
%% next year (2026) it will be in both versions.
if  (str2num(v(4:end))>=26) | (str2num(r)==0),
    fprintf('SPM version OK\n')
else
    error('Need to use current (>March 2025) development version of SPM : https://github.com/spm/spm')
end;


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
%% inversion choices
allwoi=[ 0 500; -2500 -2000;-200 300]; %% approx where evoked responses sit in these files
invmethods={'EBB','IID','GS'}; %% SPM inversion algorithms
foi=[5 90]; %% frequency window in Hz
patch_size=0.6; %% extent of source (arb units) on cortical surface
n_temp_modes=16; %% number of temporal modes
Nblocks=1; pctest=0; %% no cross validation
swap_subj=[];
%% for cross validation: Nblocks=10; pctest=10; for example



%%%%%%%%%%%%%%%%%%%%% DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%%
subind=1; %% individual subject data to look at (1-8)
fileind=1; %% file to look at

if length(subind)*length(fileind)>1,
    error('Only set up to deal with one subject and one dataset in this example code')
    % NB although simple to modify if needed
end;

Npoints=17; %% points on trajectory (17 in paper)
RandSeeds=[1:8]; %% seed which sets the trajectory (runs from 1 to 8 in paper)
%% RandSeeds=[3] to run just one representative seed (faster)
RunClean='Yes'; %% if 'No' already have all surfaces worked out

%% distortion parameters
Zmax=3; %% amount of distortion per component (as Z score)
DistortIndices=[8 100]; %% which PCAs to distort (as in paper)
HeadModel='Single Shell'; %% head model to use
loc_surface='white'; %% cortical surface to use white (in paper) or pial




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%

allFvals=[];
allR2=[];
allcrosserr=[];
allfiles=[];
gainmatfiles={};


if DISTORTANAT,
    %% loop through the list of random seeds- each seed defines a trajectory
    for rs=1:length(RandSeeds),
        RandSeed=RandSeeds(rs);

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
        full_surf_name=fullfile(subj_surf_dir, surf_name);
        save(gifti(ds_loc_surface),full_surf_name);

        %% now make a unique file name based on datafile and surface
        subj_outputdirectory=fullfile(outputdirectory,subjects(subind).subj_id);
        mkdir(subj_outputdirectory);
        [a1,b1,c1]=spm_fileparts(fname)
        coreg_fname=fullfile(subj_outputdirectory, sprintf('loc_%s_%s_norm.mat',loc_surface,b1))
        %% make copy of original fname and coregister this file to the chosen surface
        fprintf('\n Running initial coreg')
        coreg_fname=coreg_subject_distort(subjects(subind), deblank(allfname(fileind,:)),full_surf_name,coreg_fname,HeadModel,datadirectory);
        fprintf('\n Now running distortion')
        %% now make some distorted brians (distorted brain)
        %% for each brain and brian compute a lead-field
        [aM,aL]=surfandlf_distort(coreg_fname,RandSeed,HeadModel,Npoints,Zmax,DistortIndices,RunClean);

        %% now compute distances of brians to brain.
        M0=gifti(aM{1}.brain); %% get original cortex
        for i=1:Npoints,
            M1=gifti(aM{1}.brians{i});
            d0=M1.vertices-M0.vertices;
            d1(i)=mean(sqrt(dot(d0',d0'))); %% mean vertex-vertex distance to M0
        end;% for i


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


        %%&%%%%%%%%%%%%%%%%% Now run inversion models %%%%%%%%%%%
        %%% trajectory includes the new brians and the original brain

        %% get a list of the lead-field files
        gainmatfiles(rs,1:Npoints+1)=aL{1}.gainmat;
        spmfilename=aL{1}.Dc{1};
    end; % for rs
    save([rootdir filesep 'distortwkspace.mat']); % keep useful info for this subject
else
    load([rootdir filesep 'distortwkspace.mat']); % load info on distorted surfaces
end; % if DISTORTANAT


INVDISTORT=0;  %% run inversions on the distorted anatomies


%% make a copy (MAY NOT BE NECESSARY)
[a1,b1,c1]=spm_fileparts(spmfilename);
newspmfilename=[a1 filesep b1 '_copy' '.mat']
Dorig=spm_eeg_load(spmfilename);
Dorig.copy(newspmfilename)
Dnew=spm_eeg_load(newspmfilename); %% new spmfilename is a copy of Dorig


idstr=sprintf('w%d-%d_f%d-%d_Nb%d_pc%d_pt%d_Nt%d%s',round(allwoi(fileind,1)),round(allwoi(fileind,2)),...
    round(foi(1)),round(foi(2)),Nblocks,pctest,patch_size*10,n_temp_modes,swap_subj);
fprintf('\n id is %s\n',idstr)

if INVDISTORT, %% run inversions on these surfaces
    everyresults=[];
    for rs=1:length(RandSeeds),
        RandSeed=RandSeeds(rs);
        gmfiles=strvcat(gainmatfiles{rs,1:Npoints+1});

        spatialmodesname=fullfile(a1, sprintf('Grpmod_%s_seed%03dNb%d.mat',idstr,RandSeed,Npoints));
        [F_vals,R2,VE,crosserr,allresultsfiles]=invert_group_matched(Dnew,gmfiles,spatialmodesname,invmethods,allwoi(fileind,:),foi,idstr,Nblocks,pctest,patch_size,n_temp_modes)
        allFvals(rs,:,:)=F_vals;
        allR2(rs,:,:)=R2;
        allVE(rs,:,:)=VE;
        allcrosserr(rs,:,:)=crosserr;
        everyresults=strvcat(everyresults,allresultsfiles);
    end;
    save([rootdir filesep 'inv_distortwkspace.mat']); % keep useful info for this subject
else %% DONT REDO INVERSION JUST LOAD
    load([rootdir filesep 'inv_distortwkspace.mat']); % load info on distorted surfaces
end; % INVDISTORT



markerstr={'d','o','s'}
colstr={'b','r','g'}
mspairs=strvcat('bd','ro','gs');
Lw=4; %% line width
Ms=10; %% marker size

figure; %% part 1 of figure for slides
for im1=1:numel(invmethods),
    [shortF,shortdist]=simplify_dist_metric(allFvals(:,:,im1),[0 d1]);
    hold on;
    plotSeed=8; %% shows local maximum
    h01=plot(shortdist,squeeze(shortF(plotSeed,:)),colstr{im1});
    [dum,ind]=max(squeeze(shortF(plotSeed,:)));
    h1=plot(shortdist(ind),shortF(plotSeed,ind),mspairs(im1,:))
    set(h01,'LineWidth',Lw);
    set(h1(1),'LineWidth',Lw);
    set(h1(1),'MarkerSize',Ms);
end;

xlabel('Distortion (mm)')
ylabel('Free Energy')
set(gca,'Fontsize',18);


Nseed=length(RandSeeds)

figure; hold on;
for im1=1:numel(invmethods),
    plotSeed=[1:Nseed]; %% show all seeds
    [shortF,shortdist]=simplify_dist_metric(allFvals(:,:,im1),[0 d1]);
    [dum,ind]=max(squeeze(shortF(plotSeed,:))');
    for f=1:length(ind)

        h1=plot(shortdist(ind(f)),shortF(plotSeed(f),ind(f))',mspairs(im1,:))

        set(h1,'LineWidth',Lw);
        set(h1,'MarkerSize',Ms);
    end;
end;

hold on;
msumdist=zeros(numel(invmethods),1);
for im1=1:numel(invmethods),

    [shortF,shortdist]=simplify_dist_metric(allFvals(:,:,im1),[0 d1]);
    [dum,ind]=max(squeeze(shortF(plotSeed,:))');
    msumdist(im1)=0;
    for f=1:length(ind),
        msumdist(im1)=msumdist(im1)+shortdist(ind(f));
    end;
    msumdist(im1)=msumdist(im1)/length(ind);
    ax=axis;
    h=plot(msumdist(im1),ax(3),mspairs(im1,:))
    set(h,'Markersize',Ms*2)
    set(h,'LineWidth',Lw)
end;

xlabel('Distortion (mm)')
ylabel('Free Energy')
set(gca,'Fontsize',18);
axis('tight')



rs=8;
RandSeed=RandSeeds(rs);
invmethod_example={'EBB'}
gmfiles=strvcat(gainmatfiles{rs,1:Npoints+1});
spatialmodesname=fullfile(a1, sprintf('Grpmod_%s_seed%03dNb%d.mat',idstr,RandSeed,Npoints));
[F_vals,R2,VE,crosserr,allresultsfiles,Mall]=invert_group_matched(Dnew,gmfiles,spatialmodesname,invmethod_example,allwoi(fileind,:),foi,idstr,Nblocks,pctest,patch_size,n_temp_modes)
allFvals(rs,:,:)=F_vals;
allR2(rs,:,:)=R2;
allVE(rs,:,:)=VE;
allcrosserr(rs,:,:)=crosserr;
everyresults=strvcat(everyresults,allresultsfiles);


%% PLOT SOME REAL DATA AND INVERSIONS
figure;
megind=setdiff(Dnew.indchantype('MEG'),Dnew.badchannels)
data=squeeze(D(megind,:,1));
plot(D.time,data)
set(gca,'Fontsize',18);
xlabel('Time (s)')
ylabel('Field (fT)')

ind0=1;
indmaxdist=Npoints+1;
J0=squeeze(Mall(ind0,1,:,:))*data;
Jmax=squeeze(Mall(indmaxdist,1,:,:))*data;

rmsJ0=std(J0');
rmsJmax=std(Jmax');
rmsJtot=rmsJ0+rmsJmax;
[dum,maxsurfind]=max(rmsJtot);

M0=Dnew.inv{1}.mesh.tess_mni;
M0.vertices=M0.vert;
M0.faces=M0.face;
M0i=spm_mesh_inflate(gifti(M0));
figure;
h=trisurf(M0i.faces,M0i.vertices(:,1),M0i.vertices(:,2),M0i.vertices(:,3),rmsJ0)
set(h,'EdgeColor','none');
F_valsrel=F_vals-min(F_vals);
colorbar;
title(sprintf('F rel=%3.2f',F_valsrel(ind0)))
xlabel('x'); ylabel('y');zlabel('z');
set(gca,'Fontsize',18)
caxis([0 0.065])
figure;
h=trisurf(M0i.faces,M0i.vertices(:,1),M0i.vertices(:,2),M0i.vertices(:,3),rmsJmax)
set(h,'EdgeColor','none')
title(sprintf('F rel=%3.2f',F_valsrel(indmaxdist)))
xlabel('x'); ylabel('y');zlabel('z');
set(gca,'Fontsize',18)
caxis([0 0.065])
colorbar;

figure;
h=plot(D.time,J0(maxsurfind,:),D.time,Jmax(maxsurfind,:))
set(h,'Linewidth',3)
legend('True','Distorted')
title(sprintf('At MNI %3.2f %3.2f %3.2f',M0.vertices(maxsurfind,1),M0.vertices(maxsurfind,2),M0.vertices(maxsurfind,3)));
set(gca,'FontSize',18)


