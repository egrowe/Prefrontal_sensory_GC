%% Pre-Procesing of EEG data into individual single-trial datafiles

% This script uses the SPM12 toolbox for all preprocessing steps

% Batch script for all stages of EEG data pre-processing
% STEP 1: Convert data from EDF format to SPM readable dat/mat
% MANUAL STEP: Assigned channel types and locations, add event triggers
% STEP 2: Pre-processing beings: Highpass filter, mark eyeblink
%   artefacts, epoch depending on trial type, remove trials marked with
%   eyeblinks and with noise > 100 uV.
% STEP 3: Separate trials into individual files and save
% STEP 4: Robustly average data and lowpass filter

%  NOTE: Prior to running this script, file must be converted  to EDF format in EEGLAB
%  NOTE: After Step 1, the EEG channels need to be manually assigned and trial triggers added

% Elise Rowe, Monash University, 2020

clear all; close all

%% INPUT: Settings and filenames for pre-processing
epochTimeWindow = [-100 500]; %epoch around this time window (ms)
samplingFreq = 500; %sampling frequency (Hz)
trialType = 'Random'; %either 'Faces' or 'Random' (will extract one of the two trial types)

minStartEventTime = 2; %minimium trial triggers can begin (remove unnecessary event labels)
maxStartEventTime = 7; %likely max trigger time could begin

P_Names = {'7MJN'}; %this participant's ID (real)

%Rename every participant to de-identify the data and allow for public sharing
Rename_Names = {'P001'};

filepath = ['../EEG_preprocessing_Shafto/']
cd(filepath)

%Needed for eyeblink detection and source localisation
load(['.../chanFudicials_ShaftoFINAL.mat']); %load 'fudicials' for source re-coregistration nas, lpa and rpa positions

%% LOAD THE FILE AND RENAME IT TO REMOVE IDENTIFICATION %%
for gg = 1:length(P_Names)
    
    filename = [P_Names{gg} 'Move Markers-edf.edf'] %Data file
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: Convert the datafile to SPM readable format (from EDF to dat/mat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spm_jobman('initcfg');
    clear matlabbatch
    
    %
    % %%Convert (from EDF to dat/mat file)
    matlabbatch{1}.spm.meeg.convert.dataset = {[filepath filename]};
    matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
    matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.convert.outfile = '';
    matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
    matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
    matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
    matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
    matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
    
    spm_jobman('run',matlabbatch);
    
    %Rename the participant to de-identify the data (resave)
    load([filepath 'spmeeg_'  P_Names{gg} 'Move Markers-edf.mat'])
    D.fname = ['spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.mat']
    D.data.fname = [filepath 'spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.dat']
    save([filepath D.fname],'D')
    
    clear matlabbatch
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHANNELS: Manually setup the file by adding the channels using the 'Prepare' function in GUI
% (1) Manually setup the channels using the 'channel locations .xyz file'

% Need to go into this file and convert the X Y Z coordinates to two
% separate files ---- (1) for EEG channels and (2) the other for EOG channels

% MANUALLY enter into D.channels structure if not automatically in last
% step -- need to also select 'uV' as units
for gg = 1:length(Rename_Names)
    
    load([filepath 'spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.mat'])
    
    %Load the D.channels structure from known participant (apply to others)
    chanInfoFile = 'channel_info_Shafto.mat'
    load(chanInfoFile)
    
    D.channels = allChanInfo; D.fiducials = allFidicuials;
    D.transform = allTransform; D.sensors = allSensors;
    
    clear allChanInfo allFidicuials allTransform allSensors
    save([filepath D.fname],'D')

end
    
for gg = 1:length(Rename_Names)
    %% TRIGGERS: Manually need to add the triggers
    % For this you need to use an old SPM files D.trials.events structure and
    % enter your own labels that correspond to the information given in the
    % vmrk file! (this needs to be converted into Text File, then Excel then into
    % Matlab row/col format -- so from vmrk to '.txt' file

% (1) Open the vmrk file (e.g., 7MJN_Move Markers.vmrk) in TextEdit
% (2) Copy everything here and OPEN A NEW TextEdit file and set to Plain
        % Text and then SAVE this file (e.g., orig_7MJN_Move Markers.txt)
% (3) Import this text file into Excel: For this you will need to, 
        % Import text file, select 'Start Import at Row' 2 then click next
        % Select the delimiters as 'Tab' and 'Comma' (this should separate the data
        % into defined columns)
        % Click 'Next' and then 'Finish' to import to Excel in correct format
% (4) Save this file as a Text File (Tab delimited text) (e.g., 7MJN_Move_Markers.txt)
% (5) Reload the previous participant data file and add these into their D
        % structure in a format SPM can read

    load([filepath 'spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.mat'])
    
    load('emptyTrialStructure.mat');
    D.trials = trialStruct; clear trialStruct;

    D.trials.events = table2struct( readtable([filepath P_Names{gg} '_Move_Markers.txt']) )  
    
    %Rename the labels for the fields in this structuce to match SPM
    % For this, first copy Type with proper name 'type' and delete 'Type'!
    [D.trials.events.type] = D.trials.events.Type; 

    %Need to reorder the event fields to match SPM
    D.trials.events = orderfields(D.trials.events,[1:0,6,1:5]); 
    D.trials.events = rmfield(D.trials.events,'Type');    

    %Rename Description to 'value'
    [D.trials.events.value] = D.trials.events.Description; 
    D.trials.events = orderfields(D.trials.events,[1:1,6,2:5]); 
    D.trials.events = rmfield(D.trials.events,'Description');
    
    %Rename Position to 'time'
    [D.trials.events.time] = D.trials.events.Position; 
    D.trials.events = orderfields(D.trials.events,[1:2,6,3:5]); 
    D.trials.events = rmfield(D.trials.events,'Position');

    %Rename Length to 'duration' -- and made all empty
    [D.trials.events.duration] = D.trials.events.Length; 
    D.trials.events = orderfields(D.trials.events,[1:3,6,4:5]); 
    D.trials.events = rmfield(D.trials.events,'Length');
    
    for jj = 1:length(D.trials.events)
        D.trials.events(jj).duration = [];
    end

   %Rename Channel to 'offset' and set to equal 0
   [D.trials.events.offset] = D.trials.events.Channel; 
   D.trials.events = orderfields(D.trials.events,[1:4,6,5:5]); 
   D.trials.events = rmfield(D.trials.events,'Channel');

    for jj = 1:length(D.trials.events)
        D.trials.events(jj).offset = 0;
    end
   
   %Determine when the first trigger 'event' occurs
    for hh = minStartEventTime:maxStartEventTime
        trigPres = D.trials.events(hh).value;
        if isempty(trigPres) == 0
            break
        end
    end
    
    % INPUT: Every trigger NAME into the 'value' column
    % % Run the script to replace every 'type' as "STATUS"
    for ii = hh:length(D.trials.events)
        D.trials.events(ii).type = 'STATUS';
    end
    
    D.trials.events(1).type = 'Epoch';
    D.trials.events = D.trials.events';

    %INPUT: Every trigger TIME into the 'time' column
    %Run the script to shift these raw time values into secs
    for ii = hh:length(D.trials.events)
        currTime = D.trials.events(ii).time;
        D.trials.events(ii).time = currTime/samplingFreq;
    end
    
    
    %Save these newly defined events/chans etc
    save([filepath 'spmeeg_'  P_Names{gg} 'Move Markers-edf.mat'], 'D')
    clear D
    
end


%% NOW RENAME AFTER TRIGGERS ADDED
% This is done so there are two copies of the data with the triggers added
% (one with the original name and another with the 'de-identified' name)
for gg = 1:length(P_Names)
    
    load([filepath 'spmeeg_'  P_Names{gg} 'Move Markers-edf.mat'])
    
    D.fname = ['spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.mat']
    %D.data.fname = [filepath 'spmeeg_'  Rename_Names{gg} '_Move_Markers-edf.dat']
    
    save([filepath D.fname],'D')
    clear D
end



%% STEP 2 BEGINS NOW!
for bb = 1:length(P_Names)
    
    filename = [Rename_Names{bb} '_Move_Markers-edf.mat'];
    %% STEP 2: Data pre-processing begins here: Filter and remove artefacts
    % (1) mark eyeblink artefacts, (2) epoch the data depending
    % on trial type (defned by trigger), (3) detect noisy  artefacts (> 100 Hz) and
    % remove trials marked with eyeblinks and noisy artefacts.
    
    spm_jobman('initcfg');
    clear matlabbatch
    
    %Setup HEAD MODEL for eyeblink detection (using source-level coregistration)
    matlabbatch{1}.spm.meeg.source.headmodel.D = {[filepath 'spmeeg_' filename]};
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = 'Source';
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2; % 1=coarse, 2 = normal, 3 = fine
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fudicials(1,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fudicials(2,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fudicials(3,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
    
    %Setup eyeblink detection (use VEOG channel to mark trials with eyeblinks)
    matlabbatch{2}.spm.meeg.preproc.artefact.D = {[filepath 'spmeeg_' filename]};
    matlabbatch{2}.spm.meeg.preproc.artefact.mode = 'mark';
    matlabbatch{2}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
    matlabbatch{2}.spm.meeg.preproc.artefact.append = true;
    matlabbatch{2}.spm.meeg.preproc.artefact.methods.channels{1}.chan = 'VEOG';
    matlabbatch{2}.spm.meeg.preproc.artefact.methods.fun.eyeblink.threshold = 4;
    matlabbatch{2}.spm.meeg.preproc.artefact.methods.fun.eyeblink.excwin = 0;
    matlabbatch{2}.spm.meeg.preproc.artefact.prefix = 'aEBn_noHP_';
    
    %Epoch by trialType (either Faces or Random)
    matlabbatch{3}.spm.meeg.preproc.epoch.D = {[filepath 'aEBn_noHP_spmeeg_' filename]};
    matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.timewin = epochTimeWindow;
    % Select trial definitions to use (depending on 'trialType' defined at start of script)
    if strcmp(trialType, 'Faces')
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'Phase1_Face';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S 11';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'Phase2_Face';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 'S 21';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).conditionlabel = 'Phase3_Face';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventvalue = 'S 31';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.prefix = 'e_faces_';
    elseif strcmp(trialType, 'Random')
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'Phase1_Random';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 'S 12';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'Phase2_Random';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 'S 22';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).conditionlabel = 'Phase3_Random';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventtype = 'STATUS';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventvalue = 'S 32';
        matlabbatch{3}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).trlshift = 0;
        matlabbatch{3}.spm.meeg.preproc.epoch.prefix = 'e_random_';
    end
    matlabbatch{3}.spm.meeg.preproc.epoch.bc = 1;
    matlabbatch{3}.spm.meeg.preproc.epoch.eventpadding = 0;
    
    
    %Detect artefacts at 100uV and REMOVE trials marked with eyeblinks
    if strcmp(trialType, 'Faces')
        matlabbatch{4}.spm.meeg.preproc.artefact.D = {[filepath 'e_faces_aEBn_noHP_spmeeg_' filename]};
    elseif strcmp(trialType, 'Random')
        matlabbatch{4}.spm.meeg.preproc.artefact.D = {[filepath 'e_random_aEBn_noHP_spmeeg_' filename]};
    end
    
    matlabbatch{4}.spm.meeg.preproc.artefact.mode = 'reject'; %Reject channels according to methods below
    matlabbatch{4}.spm.meeg.preproc.artefact.badchanthresh = 0.2; %bad channel threshold (default = 0.2; usually used)
    matlabbatch{4}.spm.meeg.preproc.artefact.append = true;
    matlabbatch{4}.spm.meeg.preproc.artefact.methods(1).channels{1}.chan = 'VEOG'; %find VEOG channels w eyeblinks
    matlabbatch{4}.spm.meeg.preproc.artefact.methods(1).fun.events.whatevents.artefacts = 1; % 1 = "all"
    matlabbatch{4}.spm.meeg.preproc.artefact.methods(2).channels{1}.all = 'all'; %Also rejet channels w noise above 100 uV
    matlabbatch{4}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.threshold = 100; %100 uV
    matlabbatch{4}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.excwin = 1000; %default
    matlabbatch{4}.spm.meeg.preproc.artefact.prefix = 'a';
    
    spm_jobman('run',matlabbatch);
    
    clear spm_jobman
    clear matlabbatch
end