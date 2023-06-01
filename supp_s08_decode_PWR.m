%% Visualise and decode NP coherence, power and Granger Causality results
function supp_s08_decode_PWR(arrayInput)


% This script takes the POWER at the TOP VOXEL within each ROI and runs
% classificaton analysis using 70/30 (train/test) split and 10-fold cross
% validation (CV).

% The output decoding accuracy represents the ability of the classifying to
% distinguish between a FACE and RANDOM trial

% The script runs for each participant and experimental phase

% This script requires libSVM toolbox to run and the fieldtrip and bsmart
% add-ons

% Add fieldtrip to path
restoredefaultpath;
addpath(genpath(['../fieldtrip'])); ft_defaults;
addpath(genpath(['../bsmart']))

%Add the libSVM decoding library to the path
cd(['../decode_libSVM/libsvm/matlab'])
make

addpath(genpath(['../decode_libSVM/libsvm/matlab']))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countA = 1; countB = 1; countC = 1;

sprintf(['This is the current array input number:' num2str(arrayInput)])

% SELECT THE TWO NODES TO DECODE
possPFCCoords = {'L_SP_PFC','R_SP_PFC','L_Inf_PFC','R_Inf_PFC'}% all PFC ROIs
% (switch the above to all sensory for sensory ROI analysis)
possSensoryCoords = {'L_FFA'}; %select a random sensory ROI (as this isn't examined)

load('PLIST.mat') %participant list
PName = PList(arrayInput,:);

sprintf(['This P selected:' num2str(PName)])

thisPhase = 'Phase1'; % this Phase
trainProp = '0.7'; %this proportion of training trials
timeWindow = [0,500]; %this time window

useThisName = ['_ALL_PFC_']; %use savename relevant to above coords
maxFreq = 55; %maximum freq bin to examine (from 0 to maxFreq)
nCVReps = 10; %number of CV repetitions
zVal = 1; %is zscoring turned on

%Setup some counters (so we only take power values ONCE)
count = 1; cA = 1; cB = 1; cA_PFC = 1; cB_PFC = 1;
cC = 1; cD = 1; cC_PFC = 1; cD_PFC = 1;

for thisCVRep = 1:nCVReps

    ROIcount = 0;

    for PFCIdx = 1:length(possPFCCoords)

        usePFCCoord = possPFCCoords{PFCIdx};

        for sensIdx = 1:length(possSensoryCoords)

            ROIcount = ROIcount + 1;
            useSensCoord = possSensoryCoords{sensIdx};

            loadname =  ['GC_PWR_rmline_' (PName) '_' (thisPhase) '_CVmethod_thisRep' num2str(thisCVRep) ...
                '_NoFilt_wBC_zscore=_' num2str(zVal) '_bw_' num2str(usePFCCoord) ...
                '_and_' num2str(useSensCoord) '_by_TstatOverTime_usingTrainProp_' num2str(trainProp) '_and_' ...
                num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];

            savename_Results = ['results_PWR_rmline_' (PName) '_' (thisPhase) ...
                'CVmethod_10REPS_' num2str(useThisName) '_0to500ms_thisRep' num2str(thisCVRep) '.mat'] %savename

            % ASSIGN DECODING PARAMETERS (libSVM)
            nReps = 1; % number of repetitions of decoding steps
            costRange = [-10:5:10]; %cost range = 1 value (%-20:5:15)
            zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')

            %% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load([loadname]); %load face data
            freqBand = find(freq_NP.freq<maxFreq); % freq band (dep. on maxFreq)
            freqBandEnd = freqBand(end); %end freq band (dep. on maxFreq)

            % %Gather GRANGER for BOTH hemispheres
            all_Train_GC_at_PFC_faces(:,:,ROIcount) = squeeze(allGC_at_PFC_TRAIN(1,:,1:freqBandEnd));
            all_Test_GC_at_PFC_faces(:,:,ROIcount) = squeeze(allGC_at_PFC_TEST(1,:,1:freqBandEnd));

            %% RANDOM: SQUEEZE DATA INTO REQUIRED FORMAT
            all_Train_GC_at_PFC_rand(:,:,ROIcount) = squeeze(allGC_at_PFC_TRAIN(2,:,1:freqBandEnd));
            all_Test_GC_at_PFC_rand(:,:,ROIcount) = squeeze(allGC_at_PFC_TEST(2,:,1:freqBandEnd));


        end
    end
    %% ORGANISE THE DATA ONCE ALL GATHERED
    all_Train_Face = reshape(all_Train_GC_at_PFC_faces,[size(all_Train_GC_at_PFC_faces,1),...
        size(all_Train_GC_at_PFC_faces,2)*size(all_Train_GC_at_PFC_faces,3)]);
    all_Test_Face = reshape(all_Test_GC_at_PFC_faces,[size(all_Test_GC_at_PFC_faces,1),...
        size(all_Test_GC_at_PFC_faces,2)*size(all_Test_GC_at_PFC_faces,3)]);

    all_Train_Rand = reshape(all_Train_GC_at_PFC_rand,[size(all_Train_GC_at_PFC_rand,1),...
        size(all_Train_GC_at_PFC_rand,2)*size(all_Train_GC_at_PFC_rand,3)]);
    all_Test_Rand = reshape(all_Test_GC_at_PFC_rand,[size(all_Test_GC_at_PFC_rand,1),...
        size(all_Test_GC_at_PFC_rand,2)*size(all_Test_GC_at_PFC_rand,3)]);

    allTrain_Sensory_to_PFC = [all_Train_Face; all_Train_Rand];
    allTest_Sensory_to_PFC = [all_Test_Face; all_Test_Rand];


    %% NORMALISE ONCE ALL DATA FOR THIS REP IS EXTRACTED
    % Gather all LABELS for TRAIN and TEST sets
    classLabel_toPFC_train = [repmat(1, size(allTrain_Sensory_to_PFC,1)/2,1); repmat(2, size(allTrain_Sensory_to_PFC,1)/2,1)];
    classLabel_toPFC_test = [repmat(1, size(allTest_Sensory_to_PFC,1)/2,1); repmat(2, size(allTest_Sensory_to_PFC,1)/2,1)];

    if zVal == 1
        % Normalise training data
        means = mean(allTrain_Sensory_to_PFC, 1);
        stds = std(allTrain_Sensory_to_PFC, [], 1);
        means_mat = repmat(means, [size(allTrain_Sensory_to_PFC, 1), 1]);
        stds_mat = repmat(stds, [size(allTrain_Sensory_to_PFC, 1), 1]);
        allTrain_Sensory_to_PFC = (allTrain_Sensory_to_PFC - means_mat) ./ stds_mat;

        % Normalise validation data (using same parameters as for
        % training data)
        means_mat = repmat(means, [size(allTest_Sensory_to_PFC, 1), 1]);
        stds_mat = repmat(stds, [size(allTest_Sensory_to_PFC, 1), 1]);
        allTest_Sensory_to_PFC = (allTest_Sensory_to_PFC - means_mat) ./ stds_mat;
    end


    %% Optimise the cost parameter through grid search
    % bestcv = [];
    thisChanCV = 0;
    for costCounter = 1:length(costRange)
        log2c = 2^costRange(costCounter);

        %Setup the model parameters and define data + labels
        cmdCV = ['-s 0 -t 0 -c ' num2str(log2c) ' -v 10']; %
        labels = classLabel_toPFC_train;
        data = sparse(allTrain_Sensory_to_PFC);

        % 10-fold CROSS VALIDATION TO TUNE COST PARAMETER
        cv = svmtrain(labels, data , cmdCV); %CV accuracy
        cv_acc(costCounter) = cv; %save the CV accuracy

        % Determine the best cost parameter (i.e. highest CV accuracy)
        if (cv > thisChanCV)
            thisChanCV = cv; bestcv = costRange(costCounter);
        end
    end
    allBestCostParams = bestcv; %Save each of the best cost parameters (per repetition of outer loop)
    useCostFinal = 2^allBestCostParams; %choose best C parameter

    %% TRAIN THE CLASSIFIER
    cmdModel = ['-s 0 -t 0 -c ' num2str(useCostFinal) ' -q']; %

    %TRAIN THE CLASSIFIER
    labels = classLabel_toPFC_train;
    data = sparse(allTrain_Sensory_to_PFC);
    model = svmtrain(labels, data , cmdModel);

    %% TEST THE CLASSIFIER
    [predicted_label, accuracy, decision_values] = svmpredict(classLabel_toPFC_test,sparse(allTest_Sensory_to_PFC),model);
    decodability(thisCVRep) = accuracy(1,:)

end

%MEAN ACCURACY
mean_accuracy =  mean(decodability);
std_accuracy = std(decodability);

%% SAVE DATA AND RESULTS
save(savename_Results)

end
