function s04_topCoord(arrayInput)

%% Determine which voxels have the highest t-stat value over time
%   (use this as the representative voxel for the ROI -- used for GC and
%   decodng)

% NOTE: Before executing this script, all 8196 voxel coordinates need to be
% assigned to an ROI. For this, use the batch coordinate input option of
% the AnatomyToolbox extension for SPM12. The output of this file will be
% the assigned cortical locations for each voxel.

% Next, manually sort these coordinates to ROIs. In this study, we selected
% Left and Right PFC (defined as posterior, superiorl, medial prefrontal
% without inferior or orbital prefrontal), FFA, ITG, parietal and
% occipital. Once manually sorted, the coordinates for each ROI need to be
% saved as separate .mat files that are loaded in this script to determine
% the 'top voxel' per ROI (defined as the voxel with the highest t-stat
% value when cmoparing face and random trials over the 0 to 500 ms time
% window -- top t-stat at any point within this time window).

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
listCoords = ['Coordinates_8196_voxels_byROI.mat'];
trialTypeComb = {'random'}%,'random'}; %if above = 1, these trial types

% MAIN PARAMETERS
%Set the participant name
load('PLIST.mat'); %participant fileame list
PName = PList(arrayInput,:);
sprintf(['This P selected:' num2str(PName)])

%Set the Phases to examine and the ROIs (the ROIs will have an associated
%list of coordinates)
allPhases = {'Phase3'};
allROIs = {'L_SP_PFC','L_FFA','R_FFA','L_Occ','R_Occ',...
    'R_SP_PFC','L_Inf_PFC','R_Inf_PFC'};

% Other parameters
zscore = 1; %1 = on, 0 = offclea
timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,251); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 1; %to show data visualisation =1, =0 for not
maxFreq = 55; %maximum frequency to examine

%Set reference Phase coordinates
allPhases = {'Phase3'}

%% Loop through each PHASE and COORDINATE and SAVE THE TOP VOXEL
for bb = 1:length(allPhases)

    Phase = allPhases{bb};

    for vv = 1:length(allROIs)

        coordList = load([listCoords], allROIs{vv})
        coordList = cell2mat(struct2cell(coordList))
        nameROI = allROIs{vv};

        %Set the savename
        savename = ['CoordSorted_rmline_' num2str(PName) '_Phase1_EBRem_' num2str(nameROI) '_by_TstatOverTime_' ...
            num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract waveforms for face trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(['Final_rmline_' (PName) '_Phase3_0to500ms_faces_SourceWaveform_ALL_Trials.mat'])

        for gg = 1:length(coordList)
            % Determine the INDEXES FOR THE COORDINATES
            coordDiffs = abs(listIdx_all-coordList(gg,:,:));
            [x,allCoord_Idxs(gg)] = min(abs(sum(coordDiffs,2)));
        end

        for voxel = 1:length(coordList)
            thisVoxel = allCoord_Idxs(voxel);
            GM_faces(voxel,:,:) = waveform_all(thisVoxel,:,:);% time x trial x feature
        end

        clear waveform_all; clear listIdx_all;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract waveforms for RANDOM trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(['Final_rmline_' (PName) '_Phase3_0to500ms_random_SourceWaveform_ALL_Trials.mat'])

        for voxel = 1:length(coordList)
            thisVoxel = allCoord_Idxs(voxel);
            GM_random(voxel,:,:) = waveform_all(thisVoxel,:,:);% time x trial x feature
        end

        clear waveform_all;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Gather all data within time window AND VISUALISE (if required)
        % GATHER DATA FOR DECODING
        GM_faces = GM_faces(:,:,useWind); %restrict to specified time window (if applicable)
        GM_random = GM_random(:,:,useWind); %restrict to specified time window (if applicable)

        %% Plot the timecourse
        for voxel = 1:length(coordList)
            %plot timecourse of tstat
            for ii = 1:size(GM_faces,3)

                thistimePoint_F = squeeze(GM_faces(voxel,:,ii));
                thisTimePoint_NF = squeeze(GM_random(voxel,:,ii));

                [h,p,i,stat] = ttest2(thistimePoint_F,thisTimePoint_NF);

                tstat(voxel,ii) = stat.tstat;
            end
        end

        % %Sort by maximum t-stat value over time
        maxEachRow = max(tstat')
        [valMax,y] = sort(maxEachRow, 'descend')
        sortedTstat = tstat(y,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SORT THE VOXEL COORDS BY HIGHEST T-STAT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Sort coordinates by tstat over time
        sortedCoords = coordList(y,:)
        topCoord = sortedCoords(1,:)

        %EXTRACT THE DATA FROM THIS COORDINATE
        extractThisIdx = ismember(listIdx_all,topCoord,'rows');

        %%%%%%%%%%%%% WE USE THESE VOXEL LOCATIONS IDENTIFIED IN PHASE 3 in
        %%%%%%%%%%%%% BOTH PHASE 2 and 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % NOW APPLY TO PHASE 1
        %EXTRACT THE DATA FROM THIS COORDINATE IN PHASE 1
        load(['Final_rmline_' (PName) '_Phase1_0to500ms_faces_SourceWaveform_ALL_Trials.mat']);
        faceData_thisCoord_Phase1 = squeeze(waveform_all(extractThisIdx,:,:));
        mean_faceData_thisCoord_Phase1 = mean(faceData_thisCoord_Phase1);
        clear waveform_all

        load(['Final_rmline_' (PName) '_Phase1_0to500ms_random_SourceWaveform_ALL_Trials.mat']);
        randData_thisCoord_Phase1 = squeeze(waveform_all(extractThisIdx,:,:));
        mean_randData_thisCoord_Phase1 = mean(randData_thisCoord_Phase1);
        clear waveform_all randData_thisCoord

        save(savename, 'sortedTstat','tstat','sortedCoords','topCoord',...
            'faceData_thisCoord_Phase1','randData_thisCoord_Phase1','mean_faceData_thisCoord_Phase1', ...
            'mean_randData_thisCoord_Phase1')

        % NOW APPLY TO PHASE 2
        %EXTRACT THE DATA FROM THIS COORDINATE IN PHASE 2
        load(['Final_rmline_' (PName) '_Phase2_0to500ms_faces_SourceWaveform_ALL_Trials.mat']);
        faceData_thisCoord_Phase2 = squeeze(waveform_all(extractThisIdx,:,:));
        mean_faceData_thisCoord_Phase2 = mean(faceData_thisCoord_Phase2);
        clear waveform_all

        load(['Final_rmline_' (PName) '_Phase2_0to500ms_random_SourceWaveform_ALL_Trials.mat']);
        randData_thisCoord_Phase2 = squeeze(waveform_all(extractThisIdx,:,:));
        mean_randData_thisCoord_Phase2 = mean(randData_thisCoord_Phase2);
        clear waveform_all

        save(savename, 'sortedTstat','tstat','sortedCoords','topCoord',...
            'faceData_thisCoord_Phase2','randData_thisCoord_Phase2','mean_faceData_thisCoord_Phase2', ...
            'mean_randData_thisCoord_Phase2')

        clear theseCoords nameROI savename coordList sortedCoords ...
            topCoord sortedTstat maxEachRow tstat faceData_thisCoord ...
            randData_thisCoord mean_randData_thisCoord_Phase2 mean_faceData_thisCoord_Phase2 ...
            mean_randData_thisCoord_Phase3 mean_faceData_thisCoord_Phase3

    end
end

end

