function s05_applyGC_setupDecode(arrayInput)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMBINED GC ESTIMATION AND SVM CLASSIFICATION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script requires the fieldtrip add-on, the bsmart toolbox and
%decode_libSVM

restoredefaultpath;
addpath(genpath(['../fieldtrip'])); ft_defaults;
addpath(genpath(['../bsmart']))
addpath(genpath(['../decode_libSVM']))

%Define parameters
load('PLIST.mat'); %list of participant filenames
PName = PList(arrayInput,:);
sprintf(['This P selected:' num2str(PName)])
Phases = 'Phase1';

timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,251); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
maxFreq = 55; %maximum freq bin to examine (from 0 to maxFreq)

%Which trials are we examining?
possTrials = {'faces'; 'random'}; %'random' or 'faces'
%Training split proportion (here 70%)
trainProp = 0.7; % percentage to split as training set (to be used in s06 script)
% Set nReps if not specified
nRep = 10;
%Set z score to = 1 if not specified
zVal = 1;

%Save configs together
configs.timeWindow = timeWindow; configs.maxFreq = maxFreq;
configs.trainProp = trainProp;
configs.zscoreSVM = zVal; configs.nReps = nRep;

%LOAD FACE AND RANDOM DATA ALL
allPFCROIs = {'L_SP_PFC','R_SP_PFC','L_Inf_PFC','R_Inf_PFC'};
allSensoryROIs = {'L_FFA','R_FFA','L_Occ','R_Occ'};

%DETERMINE WHICH ROI LOCATIONS TO COMPARE
for mm = 1:length(allPFCROIs)

    for ss = 1:length(allSensoryROIs)

        PFC_ROI = allPFCROIs{mm};
        Sensory_ROI = allSensoryROIs{ss};

        %Set the loadnames
        loadPFCData = ['CoordSorted_rmline_' (PName) '_' (Phases) '_EBRem_' num2str(PFC_ROI) '_by_TstatOverTime_' ...
            num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];
        loadSensoryData = ['CoordSorted_rmline_' (PName) '_' (Phases) '_EBRem_' num2str(Sensory_ROI) '_by_TstatOverTime_' ...
            num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];

        %Load data at PFC ("Node Y") location
        load(loadPFCData)

        clear faceData_thisCoord_Phase3 randData_thisCoord_Phase3

        allFACE_data_NodePFC = faceData_thisCoord_Phase1;
        allRANDOM_data_NodePFC = randData_thisCoord_Phase1;

        clear faceData_thisCoord_Phase1 randData_thisCoord_Phase1

        %Load data at sensory ("Node X") location
        load(loadSensoryData)

        clear faceData_thisCoord_Phase3 randData_thisCoord_Phase3

        allFACE_data_NodeSensory = faceData_thisCoord_Phase1;
        allRANDOM_data_NodeSensory = randData_thisCoord_Phase1;

        clear faceData_thisCoord_Phase1 randData_thisCoord_Phase1

        for rr = 1:nRep

            %Set the savename
            savename =  ['GC_' (PName) '_' (Phases) '_thisRep' num2str(rr) ...
                '_zscore=_' num2str(zVal) '_bw_' num2str(PFC_ROI) ...
                '_and_' num2str(Sensory_ROI) '_usingTrainProp_' num2str(trainProp) '_and_' ...
                num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];

            sprintf(['Savename has been set for loop' num2str(rr)])

            %Randomise the random number generator each loop
            rng('shuffle')

            %Randomly pick nTrials for RANDOM that equal nTrials for FACE
            if size(allFACE_data_NodePFC,1) ~= size(allRANDOM_data_NodePFC,1)
                allRANDOM_data_NodePFC_thisRep = allRANDOM_data_NodePFC(randperm(size(allRANDOM_data_NodePFC,1),size(allFACE_data_NodePFC,1)),:);
                allRANDOM_data_NodeSensory_thisRep = allRANDOM_data_NodeSensory(randperm(size(allRANDOM_data_NodeSensory,1),size(allFACE_data_NodeSensory,1)),:);
            end

            %Set trials in test/train
            trainSplit = ceil(size(allFACE_data_NodePFC,1)*trainProp); %number of trials for 70%

            %Partition the FACE data at NODE X and Y randomly (randomise trial order)
            randomiseTrials_Face = randperm(size(allFACE_data_NodeSensory,1));  %select trial order randomly
            Face_Sensory_train = allFACE_data_NodeSensory(randomiseTrials_Face(1:trainSplit),:); % select 70% train trials: Obs x Features
            Face_Sensory_test = allFACE_data_NodeSensory(randomiseTrials_Face(trainSplit+1:end),:); % select 30% test trials: Obs x Features

            %ENSURE BELOW WORKS
            Face_PFC_train = allFACE_data_NodePFC(randomiseTrials_Face(1:trainSplit),:); % select 70% train trials: Obs x Features
            Face_PFC_test = allFACE_data_NodePFC(randomiseTrials_Face(trainSplit+1:end),:); % select 30% test trials: Obs x Features


            %Partition the RANDOM data NODE X and Y randomly (randmomise trial order)
            randomiseTrials_Rand = randperm(size(allRANDOM_data_NodeSensory_thisRep,1));  %select trial order randomly
            Random_Sensory_train = allRANDOM_data_NodeSensory_thisRep(randomiseTrials_Rand(1:trainSplit),:); % select 70% train trials: Obs x Features
            Random_Sensory_test = allRANDOM_data_NodeSensory_thisRep(randomiseTrials_Rand(trainSplit+1:end),:); % select 30% test trials: Obs x Features

            % ENUSRE THE BELOW WORKS
            Random_PFC_train = allRANDOM_data_NodePFC_thisRep(randomiseTrials_Rand(1:trainSplit),:); % select 70% train trials: Obs x Features
            Random_PFC_test = allRANDOM_data_NodePFC_thisRep(randomiseTrials_Rand(trainSplit+1:end),:); % select 30% test trials: Obs x Features

            % NORMALISE THE TRAINING AND TEST DATA SEPARATELY (demean/divide stdev)
            % APPLY DEMEAN PROCESS (NEW VERSION)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %PFC TRAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_ERP_PFC_train = [Face_PFC_train; Random_PFC_train];
            all_ERP_PFC_train = all_ERP_PFC_train - mean(all_ERP_PFC_train); %demean
            ERPforGC_PFC_train = all_ERP_PFC_train./ std(all_ERP_PFC_train); % std across trials

            Face_PFC_train = ERPforGC_PFC_train(1:size(Face_PFC_train,1),:);
            Random_PFC_train = ERPforGC_PFC_train(size(Face_PFC_train,1)+1:end,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %PFC TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_ERP_PFC_test = [Face_PFC_test; Random_PFC_test];
            all_ERP_PFC_test = all_ERP_PFC_test - mean(all_ERP_PFC_test); %demean
            ERPforGC_PFC_test = all_ERP_PFC_test./ std(all_ERP_PFC_test); % std across trials

            Face_PFC_test = ERPforGC_PFC_test(1:size(Face_PFC_test,1),:);
            Random_PFC_test = ERPforGC_PFC_test(size(Face_PFC_test,1)+1:end,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SENSORY TRAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_ERP_Sensory_train = [Face_Sensory_train; Random_Sensory_train];
            all_ERP_Sensory_train = all_ERP_Sensory_train - mean(all_ERP_Sensory_train);%demean
            ERPforGC_Sensory_train = all_ERP_Sensory_train./ std(all_ERP_Sensory_train); %std across trials

            Face_Sensory_train = ERPforGC_Sensory_train(1:size(Face_Sensory_train,1),:);
            Random_Sensory_train = ERPforGC_Sensory_train(size(Random_Sensory_train,1)+1:end,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % OTHER TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_ERP_Sensory_test = [Face_Sensory_test; Random_Sensory_test];
            all_ERP_Sensory_test = all_ERP_Sensory_test - mean(all_ERP_Sensory_test);%demean
            ERPforGC_Sensory_test = all_ERP_Sensory_test./ std(all_ERP_Sensory_test); %std across trials

            Face_Sensory_test = ERPforGC_Sensory_test(1:size(Face_Sensory_test,1),:);
            Random_Sensory_test = ERPforGC_Sensory_test(size(Random_Sensory_test,1)+1:end,:);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% APPLY GC ESTIMATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %LOAD UP SOME PARAMETERS FOR GC ESTIMATION
            nNodes = 2; %number of nodes in the GC network
            coordList = []; %empty field for coordinate list
            nTimePts = 251; %number of time points
            timeIdx = linspace(timeWindow(1),timeWindow(2),nTimePts); % find index for time points in data
            useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
            sampFreq = 500; %sampling frequency

            for nCond = 1:length(possTrials) %perform GC for these R then F trials

                theseTrials = possTrials{nCond};

                %% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%  Set nTrials: Depending on the type of trials we are analysing
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(theseTrials, 'faces')
                    nTrials = size(Face_Sensory_train,1);
                elseif strcmp(theseTrials, 'random')
                    nTrials = size(Random_Sensory_train,1);
                end

                count = 1; %setup the counter for the next loop
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Run Non-Parametric Granger Causality Analysis on a trial-by-trial basis
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for ii = 1:nTrials

                    %% Put data into fieldtrip format
                    data.time{1} = linspace(0,timeWindow(2),nTimePts)/1000; %time point indices
                    data.fsample = sampFreq; % sampling freq

                    if strcmp(theseTrials, 'random') %if 'random' selected, run on random trials
                        data.trial{1,1}(1,:) = Random_Sensory_train(ii,:); %data node 1 (i.e. lower in hierarchy)
                        data.trial{1,1}(2,:) = Random_PFC_train(ii,:); % data node 2 (i.e. higher in hierarchy)
                        data.label{1,:} = 'Random_Sensory'; %label for node 1
                        data.label{2,:} = 'Random_PFC'; %label for node 2
                    elseif strcmp(theseTrials, 'faces') %if 'faces' selected, run on face trials
                        data.trial{1,1}(1,:) = Face_Sensory_train(ii,:)
                        data.trial{1,1}(2,:) = Face_PFC_train(ii,:)
                        data.label{1,:} = 'Face_Sensory'
                        data.label{2,:} = 'Face_PFC'
                    end

                    %Skip this trial if it is invalid (i.e. all zeros)
                    if sum(data.trial{1,1}(1,:)) == 0 || sum(data.trial{1,1}(2,:)) == 0 %is this a valid trial?
                        continue
                    end
                    %% Computation of the multivariate autoregressive model (parametric)

                    cfg         = [];
                    cfg.fsample     = data.fsample; %sampling frequency
                    cfg.ntrials     = 1; %number of trials
                    cfg.triallength = 0.5; %length of trial
                    cfg.nsignal     = 2; %number of signals
                    cfg.method      = 'ar'; %method = auto-regressive
                    cfg.method = 'bsmart'; %use toolbox
                    mdata       = ft_mvaranalysis(cfg, data); %apply multivariate analysis

                    cfg.method = 'mvar';
                    mfreq      = ft_freqanalysis(cfg, mdata); %obtain spectral estimates

                    %% NON-PARAMETRIC POWER CALCULATION
                    cfg.method    = 'mtmfft'; %multitaper fourier transform
                    cfg.taper     = 'dpss'; %taper type = multi-taper DPSS
                    cfg.output    = 'fourier'; %output in fourier transform
                    cfg.tapsmofrq = 5; %number of tapers (i.e. smoothing window, here, +/-3 Hz = 6 Hz)

                    freq_NP          = ft_freqanalysis(cfg, data);
                    fd_NP            = ft_freqdescriptives(cfg, freq_NP); %gather power descriptives


                    %% Computation and inspection of the connectivity measures
                    cfg.method    = 'coh';
                    coherence_NP     = ft_connectivityanalysis(cfg, freq_NP); %calculate coherence measures

                    %% Computation of GRANGER
                    cfg           = [];
                    cfg.method    = 'granger'; %apply Granger causality measures
                    granger       = ft_connectivityanalysis(cfg, freq_NP);

                    %% GATHER ALL DATA
                    allGC_to_PFC_TRAIN(nCond,count,:) = squeeze(granger.grangerspctrm(2,1,:));
                    allGC_from_PFC_TRAIN(nCond,count,:) = squeeze(granger.grangerspctrm(1,2,:));

                    count = count + 1;

                    close all
                    clear fd_NP coherence_NP granger
                end

                clear data cfg
                %% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%  Set nTrials: Depending on the type of trials we are analysing
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(theseTrials, 'faces')
                    nTrials = size(Face_Sensory_test,1);
                elseif strcmp(theseTrials, 'random')
                    nTrials = size(Random_Sensory_test,1);
                end

                count = 1; %setup the counter for the next loop
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Run Non-Parametric Granger Causality Analysis on a trial-by-trial basis
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for ii = 1:nTrials

                    %% Put data into fieldtrip format
                    data.time{1} = linspace(0,timeWindow(2),nTimePts)/1000; %time point indices
                    data.fsample = sampFreq; % sampling freq

                    if strcmp(theseTrials, 'random') %if 'random' selected, run on random trials
                        data.trial{1,1}(1,:) = Random_Sensory_test(ii,:); %data node 1 (i.e. lower in hierarchy)
                        data.trial{1,1}(2,:) = Random_PFC_test(ii,:); % data node 2 (i.e. higher in hierarchy)
                        data.label{1,:} = 'Random_Sensory'; %label for node 1
                        data.label{2,:} = 'Random_PFC'; %label for node 2
                    elseif strcmp(theseTrials, 'faces') %if 'faces' selected, run on face trials
                        data.trial{1,1}(1,:) = Face_Sensory_test(ii,:);
                        data.trial{1,1}(2,:) = Face_PFC_test(ii,:);
                        data.label{1,:} = 'Face_Sensory'
                        data.label{2,:} = 'Face_PFC'
                    end

                    %Skip this trial if it is invalid (i.e. all zeros)
                    if sum(data.trial{1,1}(1,:)) == 0 || sum(data.trial{1,1}(2,:)) == 0 %is this a valid trial?
                        continue
                    end
                    %% Computation of the multivariate autoregressive model (parametric)

                    cfg         = [];
                    cfg.fsample     = data.fsample; %sampling frequency
                    cfg.ntrials     = 1; %number of trials
                    cfg.triallength = 0.5; %length of trial
                    cfg.nsignal     = 2; %number of signals
                    cfg.method      = 'ar'; %method = auto-regressive
                    cfg.method = 'bsmart'; %use toolbox
                    mdata       = ft_mvaranalysis(cfg, data); %apply multivariate analysis

                    cfg.method = 'mvar';
                    mfreq      = ft_freqanalysis(cfg, mdata); %obtain spectral estimates

                    %% NON-PARAMETRIC POWER CALCULATION
                    cfg.method    = 'mtmfft'; %multitaper fourier transform
                    cfg.taper     = 'dpss'; %taper type = multi-taper DPSS
                    cfg.output    = 'fourier'; %output in fourier transform
                    cfg.tapsmofrq = 5; %number of tapers (i.e. smoothing window, here, +/-3 Hz = 6 Hz)

                    freq_NP          = ft_freqanalysis(cfg, data);
                    fd_NP            = ft_freqdescriptives(cfg, freq_NP); %gather power descriptives


                    %% Computation and inspection of the connectivity measures
                    cfg.method    = 'coh';
                    coherence_NP     = ft_connectivityanalysis(cfg, freq_NP); %calculate coherence measures

                    %% Computation of GRANGER
                    cfg           = [];
                    cfg.method    = 'granger'; %apply Granger causality measures
                    granger       = ft_connectivityanalysis(cfg, freq_NP);

                    %% GATHER ALL DATA
                    allGC_to_PFC_TEST(nCond,count,:) = squeeze(granger.grangerspctrm(2,1,:));
                    allGC_from_PFC_TEST(nCond,count,:) = squeeze(granger.grangerspctrm(1,2,:));

                    if ii == 1
                        %First look at coherence per trial
                        freq_roi =  0 <= granger(1).freq & granger(1).freq < maxFreq;
                        freqWindow = granger(1).freq(freq_roi);
                    end

                    count = count + 1;

                    close all
                    clear fd_NP coherence_NP granger

                end

            end

            %Save all GC estimate data in structure
            allData{rr}.GCtoPFC_TEST = allGC_to_PFC_TEST;
            allData{rr}.PFCtoGC_TEST= allGC_from_PFC_TEST;
            allData{rr}.GCtoPFC_TRAIN= allGC_to_PFC_TRAIN;
            allData{rr}.PFCtoGC_TRAIN = allGC_from_PFC_TRAIN;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% BEGIN SVM CLASSIFICATION PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Restrict to the frequency window of interest
            allGC_to_PFC_TRAIN = allGC_to_PFC_TRAIN(:,:,1:length(freqWindow));
            allGC_to_PFC_TEST = allGC_to_PFC_TEST(:,:,1:length(freqWindow));

            %Gather together to normalise
            allTrain_Sensory_to_PFC = [squeeze(allGC_to_PFC_TRAIN(1,:,:));squeeze(allGC_to_PFC_TRAIN(2,:,:))];
            allTest_Sensory_to_PFC = [squeeze(allGC_to_PFC_TEST(1,:,:));squeeze(allGC_to_PFC_TEST(2,:,:))];

            allTrainTest.trainToPFC = allTrain_Sensory_to_PFC;
            allTrainTest.testToPFC = allTest_Sensory_to_PFC;

            save(savename)

            sprintf('Saved, now loop starting again')


            clearvars -except loadPFCData loadSensoryData ss mm rr nRep PName Phases zVal PFC_ROI ...
                Sensory_ROI timeWindow allRANDOM_data_NodePFC allRANDOM_data_NodeSensory ...
                allFACE_data_NodePFC allFACE_data_NodeSensory configs possTrials trainProp ...
                maxFreq timeIdx useWind allPFCROIs allSensoryROIs

        end

        clear savename allFACE_data_NodePFC allRANDOM_data_NodePFC ...
            allFACE_data_NodeSensory allRANDOM_data_NodeSensory
    end

end

end
