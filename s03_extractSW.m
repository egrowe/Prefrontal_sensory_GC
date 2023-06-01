function s03_extractSW(arrayInput)

%% Extract source reconstructed waveforms from SPM files (single-trial per file)
% Using specific input coordinates (coordList)
%addpath(['../spm12'])
spm eeg
close all

%Set the participant name
load('PLIST.mat'); %this list of participant names
thisP = PList(arrayInput,:);
sprintf(['This P selected:' num2str(thisP)])

%% Extract waveform from this file
%filepath = ['../ShaftoDCM/']; %set file path
con = 1; % Condition = 1 (always 1 for faces only or random only input files)
trials = 'faces'; %set the name of the experimental condition
allPhases = {'Phase1'}; %set the phase(s) to examine here


for ff = 1:length(allPhases)
    phase = allPhases{ff}; %this Phase

    nTrials = 250; %max number of trials (end when final trial reached)
    rangeCoords = 1:8196; %range of coordinates to extract (always 8196 using Multipe Sparse Priors in SPM)

    saveCoordName = ['CoordIdxs_for_' thisP '_all8196_' phase]; %how to save these coordinates to later extract

    % % Obtain voxel coordinates for ALL 8196!!
    load(['bm_Sae_rm_final_' trials 'aEBn_noHP_' thisP '_' phase '_-100to500ms_trial_26.mat'])
    model = D.other.inv{1,1};
    vert   = model.mesh.tess_mni.vert; %find voxel coordinates
    coordList = vert;
    save(saveCoordName,'coordList')

    %Setup trial type and counters
    count = 1;


    for jj = 1:length(rangeCoords)

        countB = 1;
        useCoord = rangeCoords(jj); %which coordinate index to select
        PST = coordList(useCoord,:); %the coordinate itself

        %Loop throough eact trial and extract waveform
        for thisTrial = 1:nTrials

            filename = ['bm_Sae_rm_final_' trials 'aEBn_noHP_' thisP ...
                '_' phase '_-100to500ms_trial_' num2str(thisTrial) '.mat']; % set filename to load

            if exist(filename, 'file')
                load([filename]);

                % D - SPM data structure (extract inverted data)
                %==========================================================================
                model = D.other.inv{1,1};

                % get solution and spatiotemporal basis
                %--------------------------------------------------------------------------
                J      = model.inverse.J; %J{con} = 8196 x 5 (voxels x temporal modes) CONDITIONAL EXPECTATION
                T      = model.inverse.T; %51 x 5 (time x  temporal modes ) TEMPORAL PROJECTOR
                Is     = model.inverse.Is; % 1 x 8196 (voxels) INDICES OF ACTIVE DIPOLES/VOXELS
                pst    = model.inverse.pst; % 1 x 51 (time) (PERISTIMULUS TIME)
                R2     = model.inverse.R2; % single value (VARIANCE IN SUBSPACES ACCOUNTED FOR BY MODEL (%))
                F      = model.inverse.F; % single value (LOG EVIDENCE)
                Nd     = model.inverse.Nd; % single value (total number dipoles/voxels)
                VE     = model.inverse.VE; % single value (VARIANCE EXPLAINED IN SPATIAL/TEMPORAL SUBSPACE (%))

                % - project J onto pst
                %--------------------------------------------------------------------------
                J      = J{con}*T'; % Gives J = 8196 x 5

                % - determine the coordinates within the cortical mesh
                vert   = model.mesh.tess_mni.vert; %find voxel coordinates

                % Find response at XYZ
                %--------------------------------------------------------------------------
                [i,js] = min(sum([vert(Is,1)-PST(1), vert(Is,2)-PST(2), vert(Is,3)-PST(3)].^2,2));   %js = Find real voxel coords (closest to input X,Y,Z coords)
                [i,jt] = max(abs(J(js,:)));     %jt = maximum resp TIMEPOINT (regardless of +ve or -ve) at XYZ

                % % gather response over time
                Jt    = J(js,:);                     % extract waveform over time at this voxel
                Js    = J(:,jt);                     % extract max RESPONSE from ALL voxels at THIS TIME POINT (deemed as max for this voxel)
                XYZ   = round(vert(Is(js),:));       % true XYZ coods for voxel (given the index)
                Jmax  = abs(sparse(Is,1,Js,Nd,1));   % maximum response per voxel (collapsed over time)

                % % gather confidence intervals
                % %----------------------------------------------------------------------
                qC  = model.inverse.qC(js).*diag(model.inverse.qV)';
                ci  = 1.64*sqrt(abs(qC));

                % SAVE THESE RESULTS
                % %----------------------------------------------------------------------
                waveform(jj,countB,:) = Jt; %this waveform
                clear D;
                countB = countB + 1;
            end
        end

        listIndexes(count,:) = XYZ; %save the 'aligned coordinates' that SPM pulled waveform from
        count = count + 1;
    end

    resultsName = ['Final_' thisP '_' phase '_rmline_0to500ms_' num2str(trials) '_Source_Waveform_at_for_ALL_trials_Coords_from_' ...
        num2str(rangeCoords(1)) '_to_' num2str(rangeCoords(end)) '.mat']; %filename

    save([resultsName], 'listIndexes', 'waveform'); %save these variables

    clear waveform resultsName

end
end
