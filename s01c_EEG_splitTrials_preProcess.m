%% Pre-Procesing of EEG data into individual single-trial datafiles

% Run this script after s01a and s01b for preliminary preprocessing and
% then removal of line nosie

% This script will extract each trial per condition (FACE or RANDOM) as a
% separate file per participant and experimental phase


% Elise Rowe, Monash University, 2020
clear all; close all
addpath(['../spm12'])
spm eeg
close all

%% INPUT: Settings and filenames for pre-processing
epochTimeWindow = [-100 500]; %epoch around this time window (ms)
samplingFreq = 500; %sampling frequency (Hz)
trialType = 'Random'; %either 'Faces' or 'Random'

Participant_Names = {'P0002'}

filepath = ['../preprocessingEEG/']

%% STEP 3: Extract each trial as a separate data-file
% Here, we first determine if any trials were marked as 'bad' (i.e. reject)
% using artefact detection (i.e. eyeblinks or 100uV noise). These trials
% are removed from further procesing and their labels are saved for future
% reference.

for bb = 1:length(Participant_Names)
    
    countA = 1; countB = 1; countC = 1;
    
    filename = [Participant_Names{bb} '_Move_Markers-edf.mat'];
    % We then extract each 'good' trial as a separate datafile (1 per trial)

    if strcmp(trialType, 'Faces')
        spm_jobman('initcfg');
        
        load([filepath 'ae_faces_aEBn_noHP_spmeeg_' filename]); %load dataset

        %Rename to Spartan location
        %If path needs to be re-written, uncomment here
        D.path = '../ShaftoDCM'; %Rename the file path in the EEG files to match the current location
        D.data.fname = ['../ShaftoDCM/ae_' num2str(trialType) '_aEBn_noHP_spmeeg_' Participant_Names{bb} '.dat'];
        save(filename, 'D')
        
        %Now lets look at splitting up the trials
        nTrials_Face = length(D.trials); %total number of trials
        
        for tt = 1:nTrials_Face  %Find bad trials (marked to reject)
            thisTrialReject_face(tt) = D.trials(tt).bad; % ('1' = bad, '0' = good)
        end
        
        %Save this 'bad' face trials for this participant (for use later)
        badTrials_Face = find(thisTrialReject_face == 1)
        count = 1; countA = 1; countB = 1; countC = 1;
        
        for ii = 1:nTrials_Face
            
            if any(badTrials_Face(:) == ii) %do not extract 'bad' trials
                continue
            end
            
            load([filepath 'ae_faces_aEBn_noHP_spmeeg_' filename]); %load dataset
            
            D_orig.trials = D.trials; %extract trial data under separate name
            D.trials = []; % delete the original trial structure (to be overwritten)
            
            for kk = ii+1:length(D_orig.trials)
                D_orig.trials(kk).label = 'Phase_OTHER_Faces'; % re-label all trials BEFORE the current trial
            end
            
            if ii > 1 % re-label all other trials AFTER the current trial
                fillUp = 1:ii-1;
                for jj = 1:max(fillUp)
                    D_orig.trials(jj).label = 'Phase_OTHER_Faces';
                end
            end
            
            D.trials = D_orig.trials; %re-define the trial types in this data file
            
            
            % DETERMINE WHICH TYPE OF TRIAL THIS ONE IS
            thisLabel = D.trials(ii).label; %select each label and compare
            
            if strcmp(thisLabel,'Phase1_Face')
                savename = ['Sae_faces_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase1_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countA) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countA = countA+1;
            elseif strcmp(thisLabel,'Phase2_Face')
                savename = ['Sae_faces_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase2_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countB) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countB = countB+1;
            elseif strcmp(thisLabel,'Phase3_Face')
                savename = ['Sae_faces_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase3_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countC) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countC = countC+1;
            end
            
            clear D D_orig
            
            
            %Baseline correct (sometimes unnecessary to do this step)
            matlabbatch{1}.spm.meeg.preproc.bc.D(1) = {[filepath savename]};
            matlabbatch{1}.spm.meeg.preproc.bc.timewin = [-100 0];
            matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'bm';
            
            spm_jobman('run',matlabbatch);
            
            clear spm_jobman
            clear matlabbatch
                        
            clear D_orig.trials; clear D; % clear for next loop
        end
        
    elseif strcmp(trialType, 'Random')
        % Do the same for the random trials
        load([filepath 'ae_random_aEBn_noHP_spmeeg_' filename])
        
        D.path = '../ShaftoDCM'; %Rename the file path in the EEG files to match the current location        
        D.data.fname = ['../ShaftoDCM/ae_' num2str(trialType) '_aEBn_noHP_spmeeg_' Participant_Names{bb} '.dat'];
        save(filename, 'D')
        
        %Now lets look at splitting up the trials
        nTrials_Random = length(D.trials); %nuber of trials for this trial type
        
        for tt = 1:nTrials_Random %Find bad trials (marked to reject)
            thisTrialReject_rand(tt) = D.trials(tt).bad;
        end
        
        %Save this 'bad' face trials for this participant (for later)
        badTrials_Random = find(thisTrialReject_rand == 1)% ('1' = bad, '0' = good)
        count = 1; countA = 1; countB = 1; countC = 1;
        
        for ii = 1:nTrials_Random
            
            if any(badTrials_Random(:) == ii) %do not extract 'bad' trials
                continue
            end
            
            load([filepath 'ae_random_aEBn_noHP_spmeeg_' filename]); %load preprocessed dataset
            
            D_orig.trials = D.trials; %extract trial data under separate name
            D.trials = []; % delete the original trial structure (to be overwritten)
            
            for kk = ii+1:length(D_orig.trials)
                D_orig.trials(kk).label = 'OTHER_Random'; % re-label all trials BEFORE the current trial
            end
            
            if ii > 1
                fillUp = 1:ii-1; %relabel all the trials AFTER the current trial
                for jj = 1:max(fillUp)
                    D_orig.trials(jj).label = 'OTHER_Random';
                end
            end
            
            D.trials = D_orig.trials; %re-define the trial types in this data file
            
            
            
            % DETERMINE WHICH TYPE OF TRIAL THIS ONE IS
            thisLabel = D.trials(ii).label; %select each label and compare
            
            if strcmp(thisLabel,'Phase1_Random')
                savename = ['Sae_random_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase1_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countA) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countA = countA+1;
            elseif strcmp(thisLabel,'Phase2_Random')
                savename = ['Sae_random_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase2_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countB) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countB = countB+1;
            elseif strcmp(thisLabel,'Phase3_Random')
                savename = ['Sae_random_aEBn_noHP_spmeeg_' Participant_Names{bb} '_Phase3_' num2str(epochTimeWindow(1)) 'to' ...
                    num2str(epochTimeWindow(2)) 'ms_trial_' num2str(countC) '.mat']
                save(savename, 'D'); % save this trial as a separate file
                countC = countC+1;
            end

            clear D D_orig

            %Baseline correct (sometimes unnecessary to do this step)
            matlabbatch{1}.spm.meeg.preproc.bc.D(1) = {[filepath savename]};
            matlabbatch{1}.spm.meeg.preproc.bc.timewin = [-100 0];
            matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'bm';
            
            spm_jobman('run',matlabbatch);
            
            clear spm_jobman
            clear matlabbatch

            clear D_orig.trials; clear D; %clear for next loop
        end
    end
end
