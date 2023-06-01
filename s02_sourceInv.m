function s02_sourceInv(arrayInput)

%% Source inversion for single-trial EEG data

% This script uses the SPM12 toolbox for source inversion

% Here we take the single-trial EEG ERP data and perform source
% reconstruction

% Elise Rowe, Monash University, 2020
addpath(genpath(['../spm12']))

spm eeg %load the SPM EEG conditions
close all

%Set the participant names
load('PLIST.mat'); %list containing all participant names
thisP = PList(arrayInput,:);
sprintf(['This P selected:' num2str(thisP)])


%% INPUT: Settings and filenames for pre-processing
epochTimeWindow = [-100 500]; %time window of original ERP epoching
epochSourceTimeWindow = [0 500]; %time window for epoching
trialType = 'faces'; %define trial type to run: either 'faces'or 'random'

nTrials = 260; %set max number of trials (~260 for Faces and 400 for Random)
allPhases = {'Phase1'}; %set which phase to examine

% LOAD all required auxiliary files
load('chanFudicials_ShaftoFINAL.mat'); %load 'fudicials' for nas, lpa and rpa positions

%% Run source inversion process for each trial
% Steps include: head template, coregistrations, define and invert forward
% model, specify time/frequency window and whether to use taper, create
% images of the results

for jj = 1:length(allPhases)

    thisPhase = allPhases{jj};

    for ii = 1:nTrials

        spm('defaults', 'EEG')
        spm_jobman('initcfg');

        filename = ['bm_Sae_rm_final_' num2str(trialType) 'aEBn_noHP_' num2str(thisP) '_' ...
            num2str(thisPhase) '_'  num2str(epochTimeWindow(1)) ...
            'to' num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat'] % filename for FACE trials

        if ~exist(filename, 'file')
            continue
        end

        %Template, Coregister, Forward Model
        matlabbatch{1}.spm.meeg.source.headmodel.D = {[filename]}; %load this file
        matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.comment = 'Source'; %run source reconstruction
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1; %mesh template (1 = 'yes', 0 = 'no')
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2; % 1=coarse, 2 = normal, 3 = fine
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas'; %position of nasion
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fudicials(1,:);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa'; %position of left periauricular
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fudicials(2,:);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa'; %position of right periauricular
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fudicials(3,:);
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0; % 0 = 'no, 1 = 'yes'
        matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM'; %boundary element method forward model
        matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell'; % single shell forward model

        %Invert forward model and specify time (and frequency) window
        matlabbatch{2}.spm.meeg.source.invert.D = {[filename]};
        matlabbatch{2}.spm.meeg.source.invert.val = 1;
        matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1; %use all conditions
        matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS'; %Multiple sparse priors (greedy search)
        matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = epochSourceTimeWindow; %time window
        matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256]; %frequency window
        matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 0; %1 = Hanning taper at start and end of trial (0 = no taper, 1 = yes)
        matlabbatch{2}.spm.meeg.source.invert.modality = {'EEG'}; %modality = EEG

        %Show results within time (and frequency) window and create images
        matlabbatch{3}.spm.meeg.source.results.D = {[filename]};
        matlabbatch{3}.spm.meeg.source.results.val = 1;
        matlabbatch{3}.spm.meeg.source.results.woi = epochSourceTimeWindow; % time of interest
        matlabbatch{3}.spm.meeg.source.results.foi = [0 256]; % frequency window specify
        matlabbatch{3}.spm.meeg.source.results.ctype = 'trials'; % 'evoked' 'induced' or single 'trials'
        matlabbatch{3}.spm.meeg.source.results.space = 1; % 1=MNI or Native
        matlabbatch{3}.spm.meeg.source.results.format = 'image';
        matlabbatch{3}.spm.meeg.source.results.smoothing = 12; % mm %Smoothing mm^3

        spm_jobman('serial',matlabbatch);

        clear spm_jobman

    end

end

end
