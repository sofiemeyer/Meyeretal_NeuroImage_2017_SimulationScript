%% Script for simulating and reconstructing hippocampal activity
%% Results presented in Meyer et al., 2017 NeuroImage
%% Created by Sofie Meyer 12/2016

clear all

%% Set parameters
nRuns           = 1;                    % number of simulation iterations
snrVector       = [0 -5 -10 -15 -20];   % signal-to-noise raio, X log10 (RMSsource/RMSnoise) sensor level
coregerrVector  = [0 1 2 3];            % co-registration error: standard deviation of error (mm) to each of three fiducial locations in x,y,z dimensions
simMesh         = 1;                    % 1=hippocampal simulations, 2=cortical
algorithm       = 'EBB';                % 'EBB' for Empirical Bayes Beamformer; 'IID' for Minimum Norm; 'GS' (Greedy Search) for Multiple Sparse Priors

allmeshes=char('combinedmesh.gii','corticalmesh.gii');

%% Create simulation locations and MSP priors
load('verts_mni'); % mesh vertices in MNI space
nCortverts = 10595;
nHippverts = 162;

for i = 1:3
    rng(i,'twister'); % initialize
    
    % generate 100 random cortical priors
    cort100i = randi([1 nCortverts],100,1);
    cort100v = verts_mni(cort100i,:);
    
    % generate 10 random hippocampal priors
    rng(i,'twister');
    hipp10i = randi([nCortverts nCortverts+nHippverts],10,1);
    hipp10v = verts_mni(hipp10i,:);
    
    % save vertex values as simulation locations
    cortex_locs = cort100v(1:90,:);
    hipp_locs   = hipp10v;
    name = ['simverts' num2str(i)];
    save(name,'cortex_locs', 'hipp_locs');
    
    allhipplocs((i-1)*10+1:10*i,1:3) = hipp_locs;
    allcortlocs((i-1)*10+1:10*i,1:3) = cortex_locs(1:10,:);
    
    % save indices as priorlist for nested model
    Ip = [cort100i(1:90); hipp10i]';
    name = ['nestedpriors' num2str(i)];
    save(name,'Ip');
    
    % save indices as priorlist for nonnested model
    Ip = cort100i';
    name = ['nonnestedpriors' num2str(i)];
    save(name,'Ip');
end
name = 'allsimlocs';
save(name,'allcortlocs', 'allhipplocs');

%% Simulate data
rawfilename = 'bemgfbinfspmeeg_fingertest_LUCIA_20120223_06.mat';

load('C:\Users\smeyer\Dropbox\work\SimulationPaper\Priors+Simlocs\allsimlocs.mat')
if simMesh == 1
    simverts = allhipplocs;
elseif simMesh == 2
    simverts = allcortlocs;
end

for SNR = snrVector
    for run = 1:nRuns
        
        % simulate data:
        matlabbatch=[];
        matlabbatch{1}.spm.meeg.source.simulate.D = {rawfilename};
        matlabbatch{1}.spm.meeg.source.simulate.val = 1;
        matlabbatch{1}.spm.meeg.source.simulate.prefix = ['run_' num2str(run) '_snr_' num2str(SNR)]; % prefix for datafile
        matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = [0 300];               % time window, ms
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = 20;              % simulation frequency,  Hz
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = 20;                 % dipole strength, nAm
        matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = [simverts(run,1) simverts(run,2) simverts(run,3)];   % simulation location, mm
        matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;                                 % signal-to-noise raio, X log10 (RMSsource/RMSnoise) sensor level
        [a,~]=spm_jobman('run', matlabbatch);
        simfile=cell2mat(a{1}.D);
        
        % save single trial
        Dfull       = spm_eeg_load(simfile);
        mfilename   = ['st' num2str(run) '_snr' num2str(SNR) '.mat']; % single trial
        D           = clone(Dfull, mfilename, [Dfull.nchannels, Dfull.nsamples, 1], 2);
        D(:, :, :)  = Dfull(:,:,1); % copy the data, save only one trial
        D           = conditions(D, ':', Dfull.conditions(1));
        save(mfilename, 'D'); % prefix is ss, single trial
    end
end


%% Invert data
for SNR = snrVector
    for coregerror = coregerrVector
        for run = 1:nRuns
            
            filename = ['st' num2str(run) '_snr' num2str(SNR) '.mat']; % prefix st for single trial
            
            for anatModel = 1:2 %combined, cortical
                
                % Co-register and add co-registration error
                matlabbatch = [];
                matlabbatch{1}.spm.meeg.source.headmodel.D = {filename};
                matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
                matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {'msMQ0484_orig.img'};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {deblank(allmeshes(anatModel,:))};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
                matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = [9.9898 142.5147 -8.1787];
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = [-55.7659 49.4636 -26.2089];
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = [88.7153 62.4787 -29.1394];
                matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
                matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg = 'EEG BEM';
                matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg = 'Single Shell';
                matlabbatch{2}.spm.meeg.source.coregshift.D(1) = cfg_dep('MEG helmet head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
                matlabbatch{2}.spm.meeg.source.coregshift.val = 1;
                matlabbatch{2}.spm.meeg.source.coregshift.meanshift = [0 0 0];
                matlabbatch{2}.spm.meeg.source.coregshift.sdshift = [0 0 0];
                matlabbatch{2}.spm.meeg.source.coregshift.meanangshift = [0 0 0];
                matlabbatch{2}.spm.meeg.source.coregshift.sdangshift = [0 0 0];
                matlabbatch{2}.spm.meeg.source.coregshift.pperror = coregerror;
                [a1,~]=spm_jobman('run', matlabbatch);
                coregfile = cell2mat(a1{1}.D);               
                
                % Invert
                matlabbatch = [];
                matlabbatch{1}.spm.meeg.source.invertiter.D = {coregfile};
                matlabbatch{1}.spm.meeg.source.invertiter.val = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.whatconditions.all = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invfunc = 'Classic';
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.invtype = algorithm;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.woi = [0 300];                          % time window of interest; ms
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.foi = [0 80];                           % frequency window of interest; Hz
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.hanning = 1;
                if strcmp(algorithm,'EBB') || strcmp(algorithm, 'IID')
                    matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.npatches = 512;
                    matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.randpatch.niter = 1;
                elseif strcmp(algorithm, 'GS')
                    if anatModel == 1
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {['nestedpriors' num2str(ceil(run/10)) '.mat']};
                    elseif anatModel == 2
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedfile = {['nonnestedpriors' num2str(ceil(run/10)) '.mat']};
                        matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.isfixedpatch.fixedpatch.fixedrows = [];
                    end
                end
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.patchfwhm = 0.6;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.mselect = 0;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.nsmodes = 274;                          % number of spatial modes; all channels used
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.umodes = {'U.mat'};                     % spatial mode file = eye(nSpatialmodes)
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.ntmodes = 1;                            % number of temporal modes
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.priorsmask = {''};
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.priors.space = 1;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.locs = zeros(0, 3);
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.restrict.radius = 32;
                matlabbatch{1}.spm.meeg.source.invertiter.isstandard.custom.outinv = '';
                matlabbatch{1}.spm.meeg.source.invertiter.modality = {'All'};
                matlabbatch{1}.spm.meeg.source.invertiter.crossval = [0 1]; 
                [a2,~]=spm_jobman('run', matlabbatch);
                D = spm_eeg_load(cell2mat(a2{1}.D));                               
                
                
                results(1,run,anatModel) = D.inv{1}.inverse.F;         % Free energy
                results(2,run,anatModel) = D.inv{1}.inverse.R2;        % Variance explained
                results(3,run,anatModel) = find(D.inv{1}.inverse.qC==max(D.inv{1}.inverse.qC),1); % Peak in source distribution
                
                
                if run == 30
                    name = ['results_snr_' num2str(SNR) '_err_', num2str(coregerror)];
                    save(name,'results');
                end
                
            end
        end % run        
    end % coregerror
end % SNR
