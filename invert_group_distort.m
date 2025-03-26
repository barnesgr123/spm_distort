function [F_vals,R2,VE,crosserr,allresultsfiles]=invert_group_distort(Dnew,gainmatfiles,spatialmodesname,invmethods,woi,foi,idstr,Nblocks,pctest,patch_size,n_temp_modes)
%%function [F_vals,R2,VE,crosserr,allresultsfiles]=invert_group_distort(Dnew,gainmatfiles,spatialmodesname,invmethods,woi,foi,idstr,Nblocks,pctest,patch_size,n_temp_modes)

% Setup spatial modes for cross validation and to ensure same modes used
% across the group of inversions (not biased to first or last etc)

val=Dnew.val;
[spatialmodesname,Nmodes,pctest]=spm_eeg_inv_prep_modes_xval(Dnew.fullfile, [], spatialmodesname, Nblocks, pctest,gainmatfiles);


% now invert these data using this headmodel

Ngainmat=size(gainmatfiles,1);

F_vals=zeros(1,numel(invmethods));
R2=zeros(1,numel(invmethods));
allresultsfiles=[];

for f=1:Ngainmat,
    
    D=spm_eeg_load(Dnew.fullfile);
    % load in a dataset but swap in different gainmat (lead field) files
    gainmatfile=deblank(gainmatfiles(f,:));
    [gpath,gname,gext]=spm_fileparts(gainmatfile);
    copyfile(gainmatfile,[D.path filesep gname gext]); % have to put gainmat in spmfile directory
    D.inv{val}.gainmat=[gname gext]; %%  just change name of gainmat in file
    D.save;

    % now run the inversion
    matlabbatch={};
    batch_idx=1;    

    resultsfile=[gpath filesep 'Inv_' idstr gname '.mat']
    allresultsfiles=strvcat(allresultsfiles,resultsfile);
    for i1=1:numel(invmethods),
      %%  Source reconstruction
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.D = {Dnew.fullfile};
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.val = val;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.whatconditions.all = 1;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.invtype = invmethods{i1}; %;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.woi = woi;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.foi = foi;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.hanning = 0;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = patch_size; %% NB A fiddle here- need to properly quantify
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = Nmodes;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.umodes = {spatialmodesname};
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = n_temp_modes;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.modality = {'All'};
        matlabbatch{batch_idx}.spm.meeg.source.invertiter.crossval = [pctest Nblocks];
        spm_jobman('run', matlabbatch);

     %   Get metrics for inversion
        Drecon=spm_eeg_load(Dnew.fullfile);
        F_vals(f,i1)=Drecon.inv{1}.inverse.F;
        R2(f,i1)=Drecon.inv{1}.inverse.R2;
        VE(f,i1)=Drecon.inv{1}.inverse.VE;
        crosserr(f,i1)=mean(Drecon.inv{1}.inverse.crosserr(:));
    
    end; % for invmethods
    
    fprintf('\n Done %d of %d',f,Ngainmat);
    save(resultsfile,'invmethods','R2','VE','crosserr','F_vals','woi','foi','Ngainmat','patch_size','Nblocks','pctest','spatialmodesname','gname','n_temp_modes');
    fprintf('\n Deleting Gainmat')
    delete([D.path filesep gname gext]); %% remove copied gainmat file to save disk space
    D.inv{val}.gainmat=''; %%  just change name of gainmat in spmfile
    D.save;
    
end; % for gainmat


