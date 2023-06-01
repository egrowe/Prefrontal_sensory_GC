% Remove the 60Hz line noise in the data (to clean data)

%%%%%%%% This script requires SPM and the Chronux add-on toolbox %%%%%%
addpath(genpath('../spm12'));

%%%%%%%% START LOADING THE PARTICIPANT DATA %%%%%%
thisP = 'P0912'; %this participant (change as necessary)

%load to determine number of trials
load(['ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
nTrials_face = size(D.trials,2);
clear D;

% RESAVE THE AE_X_AEBN files WITH RMLINE APPLIED
D = spm_eeg_load(['ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
data_new_faces = spm2fieldtrip(D)

%Remove SPM and add Chronux
rmpath(genpath('../spm12'));
addpath(genpath(['../chronux_2_11']));

%APPLY LINE NOISE REMOVAL
params.Fs = 500 %D.Fsample; %2.0345e+03; %sampling freq.. possibly change to sampResms?
params.tapers = [5 9];
params.pad=2;

%Remove from each trial (Ph3 only) and channel
for thisTrial = 1:nTrials_face

    thisTrialData = data_new_faces.trial{1,thisTrial};

    for ch = 1:size(thisTrialData,1)

        thisFChTrialData = squeeze(thisTrialData(ch,:));
        rmSignal_face(ch,:) = rmlinesc(thisFChTrialData,params,[],'n',60);
    end

    data_new_faces.trial{1,thisTrial} = rmSignal_face;
    clear rmSignal_face thisTrialData
end

%SAVE THESE RESULTS
addpath(genpath(['../PData_MASSIVE/spm12']));
spm_eeg_ft2spm(data_new_faces, ['rm_final_ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
clear D

load(['ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
D_orig = D; clear D;
load(['rm_final_ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])

D.trials = D_orig.trials;
D.channels = D_orig.channels;
D.sensors = D_orig.sensors;
D.fiducials = D_orig.fiducials;
D.transform = D_orig.transform;
D.condlist = D_orig.condlist;
D.montage = D_orig.montage;
D.history = D_orig.history;
D.other = D_orig.other;

save(['rm_final_ae_faces_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'],'D')
clear D nTrials_face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NOW FOR RANDOM 
%load to determine trial number
load(['ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
nTrials_random = size(D.trials,2);
clear D;

% RESAVE THE AE_X_AEBN files WITH RMLINE APPLIED
D = spm_eeg_load(['ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
data_new_random = spm2fieldtrip(D)

%Remove SPM and add Chronux
rmpath(genpath('../PData_MASSIVE/spm12'));
addpath(genpath(['../PData_MASSIVE/chronux_2_11']));

%APPLY LINE NOISE REMOVAL
params.Fs = 500 %D.Fsample; %2.0345e+03; %sampling freq.. possibly change to sampResms?
params.tapers = [5 9];
params.pad=2;

%Remove from each trial (Ph3 only) and channel
for thisTrial = 1:nTrials_random

    thisTrialData = data_new_random.trial{1,thisTrial};

    for ch = 1:size(thisTrialData,1)

        thisRChTrialData = squeeze(thisTrialData(ch,:));
        rmSignal_rand(ch,:) = rmlinesc(thisRChTrialData,params,[],'n',60);
    end

    data_new_random.trial{1,thisTrial} = rmSignal_rand;
    clear rmSignal_rand thisTrialData
end

%SAVE THESE RESULTS
addpath(genpath(['../PData_MASSIVE/spm12']));
spm_eeg_ft2spm(data_new_random, ['rm_final_ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
clear D

load(['ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])
D_orig = D; clear D;
load(['rm_final_ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'])

D.trials = D_orig.trials;
D.channels = D_orig.channels;
D.sensors = D_orig.sensors;
D.fiducials = D_orig.fiducials;
D.transform = D_orig.transform;
D.condlist = D_orig.condlist;
D.montage = D_orig.montage;
D.history = D_orig.history;
D.other = D_orig.other;

save(['rm_final_ae_random_aEBn_noHP_spmeeg_' num2str(thisP) '_Move_Markers-edf.mat'],'D')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

